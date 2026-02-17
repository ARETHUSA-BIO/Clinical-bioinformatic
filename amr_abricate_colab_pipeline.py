"""
AMR deep analysis pipeline for Google Colab using ABRicate.

How to run in Colab:
1. Upload this script and run: !python amr_abricate_colab_pipeline.py
2. Follow prompts to upload or paste a bacterial genome FASTA.

The pipeline will:
- Install all Python dependencies with pip.
- Install ABRicate and required tools with micromamba/conda.
- Perform DNA FASTA QC + IUPAC checks (protein sequences are rejected).
- Compute metrics (GC%, Shannon entropy, sequence length, complexity).
- Run ABRicate across ALL resistance databases.
- Preserve EVERY ABRicate output column:
    #FILE  SEQUENCE  START  END  STRAND  GENE  COVERAGE  COVERAGE_MAP
    GAPS   %COVERAGE  %IDENTITY  DATABASE  ACCESSION  PRODUCT  RESISTANCE
- Filter top hits (per gene, ranked by %IDENTITY then %COVERAGE).
- Report all hits + dedicated top-hits section with inferred mechanisms.
"""

from __future__ import annotations

import csv
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Dependency bootstrap — pip installs happen before any import
# ---------------------------------------------------------------------------
REQUIRED_PY_MODULES = [
    ("pandas", "pandas"),
    ("tabulate", "tabulate"),
]


def pip_install_if_needed() -> None:
    """Install Python dependencies with pip inside Colab/runtime."""
    for import_name, pip_name in REQUIRED_PY_MODULES:
        try:
            __import__(import_name)
        except ImportError:
            print(f"[INFO] Installing missing Python module via pip: {pip_name}")
            subprocess.run(
                [sys.executable, "-m", "pip", "install", "-q", pip_name],
                check=True,
            )


pip_install_if_needed()

import pandas as pd  # noqa: E402
from tabulate import tabulate  # noqa: E402


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# IUPAC DNA alphabet (allow N + ambiguity codes)
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
STANDARD_DNA = set("ACGTN")

# ABRicate canonical column order (matches --noheader + default header line)
ABRICATE_COLUMNS: List[str] = [
    "#FILE",
    "SEQUENCE",
    "START",
    "END",
    "STRAND",
    "GENE",
    "COVERAGE",
    "COVERAGE_MAP",
    "GAPS",
    "%COVERAGE",
    "%IDENTITY",
    "DATABASE",
    "ACCESSION",
    "PRODUCT",
    "RESISTANCE",
]

# Top-hit filtering thresholds
TOP_HIT_MIN_COVERAGE: float = 80.0   # %COVERAGE >= this value
TOP_HIT_MIN_IDENTITY: float = 80.0   # %IDENTITY >= this value


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SequenceQC:
    seq_id: str
    length: int
    gc_percent: float
    shannon_entropy: float
    complexity_k4: float
    invalid_chars: List[str]
    non_standard_iupac: List[str]


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd: List[str], env: Optional[Dict[str, str]] = None) -> None:
    """Run a shell command, streaming stdout/stderr."""
    print(f"[CMD] {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(str(c) for c in cmd)}")


def detect_colab() -> bool:
    """Return True when executing inside Google Colab."""
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Conda / micromamba + ABRicate installation
# ---------------------------------------------------------------------------

def install_micromamba_and_abricate(base_dir: Path) -> Tuple[Path, Dict[str, str]]:
    """
    Install micromamba into *base_dir/micromamba*, then create a conda env
    with ABRicate, BLAST, any2fasta (required by ABRicate), and seqkit.

    Returns (env_prefix, env_dict) where env_dict has the updated PATH so
    every subsequent subprocess call finds the correct binaries.
    """
    mm_dir = base_dir / "micromamba"
    mm_bin = mm_dir / "bin" / "micromamba"

    # --- Install micromamba if missing ---
    if not mm_bin.exists():
        mm_dir.mkdir(parents=True, exist_ok=True)
        print("[INFO] Downloading micromamba …")
        # The official endpoint returns a .tar.bz2 archive; extract only bin/micromamba
        dl_cmd = (
            "curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest "
            f"| tar -xvj -C '{mm_dir}' bin/micromamba"
        )
        subprocess.run(["bash", "-c", dl_cmd], check=True)

    if not mm_bin.exists():
        raise RuntimeError(
            "micromamba binary not found after download attempt. "
            "Check network access in your Colab session."
        )

    env_prefix = base_dir / "envs" / "amr-abricate"

    # --- Create conda environment with ABRicate + dependencies ---
    if not env_prefix.exists():
        print("[INFO] Creating conda environment with ABRicate (this may take ~5 min) …")
        run_cmd([
            str(mm_bin), "create", "-y",
            "-p", str(env_prefix),
            "-c", "conda-forge",
            "-c", "bioconda",
            # ABRicate is a Perl script — include perl explicitly for robustness
            "perl",
            "perl-bioperl",      # BioPerl required by ABRicate internals
            "abricate",
            "blast",
            "any2fasta",         # ABRicate hard-requires this helper
            "seqkit",
        ])

    env = os.environ.copy()
    env["PATH"] = f"{env_prefix / 'bin'}:{env.get('PATH', '')}"
    # Set PERL5LIB so ABRicate finds BioPerl inside the env
    env["PERL5LIB"] = str(env_prefix / "lib" / "perl5" / "site_perl")

    # --- Populate ABRicate databases ---
    print("[INFO] Setting up ABRicate databases …")
    run_cmd(["abricate", "--setupdb"], env=env)

    return env_prefix, env


# ---------------------------------------------------------------------------
# FASTA parsing and QC
# ---------------------------------------------------------------------------

def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Simple FASTA parser — returns list of (header, sequence_uppercase)."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    seq_chunks: List[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_chunks).upper()))
            header = line[1:].strip() or "unnamed"
            seq_chunks = []
        else:
            seq_chunks.append(line.replace(" ", "").upper())

    if header is not None:
        records.append((header, "".join(seq_chunks).upper()))

    return records


def is_probable_protein(seq: str) -> bool:
    """
    Heuristic: reject sequences that look like amino-acid FASTA.
    Flags any character that only exists in the protein alphabet
    (E, F, I, L, P, Q, Z, X, J, O, U) OR a high AA-fraction with
    low DNA-fraction.
    """
    if not seq:
        return False
    aa_exclusive = set("EFILPQZXJOU")
    if any(ch in aa_exclusive for ch in seq):
        return True

    dna_count = sum(ch in IUPAC_DNA for ch in seq)
    aa_count = sum(ch in set("ACDEFGHIKLMNPQRSTVWY") for ch in seq)
    dna_frac = dna_count / len(seq)
    aa_frac = aa_count / len(seq)
    return dna_frac < 0.85 and aa_frac > 0.90


def shannon_entropy(seq: str) -> float:
    counts = Counter(seq)
    total = len(seq)
    if total == 0:
        return 0.0
    return -sum((n / total) * math.log2(n / total) for n in counts.values())


def kmer_complexity(seq: str, k: int = 4) -> float:
    if len(seq) < k:
        return 0.0
    kmers = {seq[i: i + k] for i in range(len(seq) - k + 1)}
    max_possible = min(4 ** k, len(seq) - k + 1)
    return len(kmers) / max_possible if max_possible else 0.0


def qc_record(seq_id: str, seq: str) -> SequenceQC:
    invalid = sorted({ch for ch in seq if not ch.isalpha() or ch not in IUPAC_DNA})
    non_standard = sorted({ch for ch in seq if ch in IUPAC_DNA and ch not in STANDARD_DNA})
    gc = 100.0 * sum(ch in "GC" for ch in seq) / len(seq) if seq else 0.0
    return SequenceQC(
        seq_id=seq_id,
        length=len(seq),
        gc_percent=gc,
        shannon_entropy=shannon_entropy(seq),
        complexity_k4=kmer_complexity(seq, k=4),
        invalid_chars=invalid,
        non_standard_iupac=non_standard,
    )


# ---------------------------------------------------------------------------
# User input (Colab / local)
# ---------------------------------------------------------------------------

def get_fasta_from_user(workdir: Path) -> Path:
    """
    In Colab: offer file-upload widget or paste.
    Outside Colab: prompt for local path or paste.
    """
    print("\nChoose FASTA input method:")
    print("  1) Upload genome FASTA file")
    print("  2) Paste FASTA content directly")
    choice = input("Enter 1 or 2 [default 1]: ").strip() or "1"

    fasta_path = workdir / "input_genome.fasta"

    if choice == "2":
        print("Paste your FASTA content below. When finished, enter a line with only: END")
        lines: List[str] = []
        while True:
            line = input()
            if line.strip() == "END":
                break
            lines.append(line)
        fasta_text = "\n".join(lines).strip()
        if not fasta_text.startswith(">"):
            raise ValueError("Pasted content is not valid FASTA (must start with '>').")
        fasta_path.write_text(fasta_text + "\n")
        return fasta_path

    if detect_colab():
        from google.colab import files  # type: ignore

        print("Please upload a genome FASTA file (.fa / .fasta / .fna)")
        uploaded = files.upload()
        if not uploaded:
            raise RuntimeError("No file was uploaded.")
        first_name = next(iter(uploaded))
        fasta_path.write_bytes(uploaded[first_name])
        return fasta_path

    # Non-Colab fallback
    local_path = input("Enter local path to FASTA file: ").strip()
    src = Path(local_path)
    if not src.exists():
        raise FileNotFoundError(f"File not found: {src}")
    shutil.copy(src, fasta_path)
    return fasta_path


# ---------------------------------------------------------------------------
# ABRicate runner
# ---------------------------------------------------------------------------

def run_abricate_all_dbs(
    fasta_path: Path,
    out_dir: Path,
    env: Dict[str, str],
) -> Path:
    """
    Run ABRicate against every available database and merge results into a
    single TSV that preserves ALL canonical ABRicate output columns.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Discover available databases
    db_proc = subprocess.run(
        ["abricate", "--list"],
        env=env, check=True, text=True, capture_output=True,
    )
    dbs: List[str] = []
    for line in db_proc.stdout.splitlines()[1:]:  # skip header row
        parts = re.split(r"\s+", line.strip())
        if parts and parts[0]:
            dbs.append(parts[0])

    if not dbs:
        raise RuntimeError("ABRicate returned no databases. Run 'abricate --setupdb' manually.")

    print(f"[INFO] Databases found: {', '.join(dbs)}")

    combined_tsv = out_dir / "abricate_combined.tsv"
    header_written = False

    with combined_tsv.open("w") as out_fh:
        for db in dbs:
            print(f"[INFO] Screening against database: {db}")
            run = subprocess.run(
                ["abricate", "--db", db, str(fasta_path)],
                env=env, check=True, text=True, capture_output=True,
            )
            lines = run.stdout.strip().splitlines()
            if not lines:
                continue

            # lines[0] is the header; lines[1:] are data rows
            if not header_written:
                out_fh.write(lines[0] + "\n")
                header_written = True

            for row in lines[1:]:
                if row.strip():
                    out_fh.write(row + "\n")

    return combined_tsv


# ---------------------------------------------------------------------------
# Mechanism inference
# ---------------------------------------------------------------------------

def infer_mechanism(product: str, resistance: str) -> str:
    """
    Rule-based lookup: map gene/product keywords to a resistance mechanism
    description.  Returns a generic note when no keyword matches.
    """
    text = f"{product} {resistance}".lower()
    rules = [
        ("beta-lactamase",   "Enzymatic hydrolysis of the beta-lactam ring (beta-lactamase)."),
        ("bla",              "Beta-lactam hydrolysis by beta-lactamase enzyme."),
        ("aminoglycoside",   "Enzymatic modification (acetylation/phosphorylation/adenylation) or 16S rRNA methylation reducing aminoglycoside binding."),
        ("efflux",           "Active efflux pump expels antibiotic, lowering intracellular concentration."),
        ("tet",              "Ribosomal protection protein or efflux-mediated tetracycline resistance."),
        ("qnr",              "Pentapeptide repeat protein protects DNA gyrase/topoisomerase IV from quinolone binding."),
        ("mec",              "Altered penicillin-binding protein (PBP2a) with reduced beta-lactam affinity."),
        ("van",              "Cell-wall precursor remodelling (D-Ala→D-Lac or D-Ser) reduces glycopeptide binding."),
        ("sul",              "Alternative dihydropteroate synthase confers sulfonamide resistance."),
        ("dfr",              "Alternative dihydrofolate reductase confers trimethoprim resistance."),
        ("erm",              "rRNA N6-adenosine methylase reduces macrolide/lincosamide/streptogramin B binding."),
        ("mph",              "Macrolide phosphotransferase inactivates macrolide antibiotics."),
        ("cat",              "Chloramphenicol acetyltransferase inactivates chloramphenicol."),
        ("int",              "Integron-associated gene cassette — context-dependent resistance."),
        ("oxa",              "OXA-type beta-lactamase with variable spectrum (including carbapenemase activity)."),
        ("kpc",              "Klebsiella pneumoniae carbapenemase — hydrolysis of carbapenems and other beta-lactams."),
        ("ndm",              "New Delhi metallo-beta-lactamase — broad-spectrum carbapenemase."),
        ("vim",              "VIM-type metallo-beta-lactamase — carbapenem hydrolysis."),
        ("imp",              "IMP-type metallo-beta-lactamase — carbapenem hydrolysis."),
        ("mcr",              "Phosphoethanolamine transferase modifies lipid A, reducing colistin binding."),
    ]
    for keyword, annotation in rules:
        if keyword in text:
            return annotation
    return "Putative AMR determinant detected; mechanism requires manual curation."


# ---------------------------------------------------------------------------
# Top-hits selection
# ---------------------------------------------------------------------------

def extract_top_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    From the full hit table, select the BEST hit per unique gene:
      1. Filter rows where %COVERAGE >= TOP_HIT_MIN_COVERAGE
                        AND %IDENTITY >= TOP_HIT_MIN_IDENTITY.
      2. Within each GENE group, keep the single row with the highest
         %IDENTITY; break ties by %COVERAGE, then by DATABASE name.
      3. Sort the resulting table by %IDENTITY DESC, then %COVERAGE DESC.

    If no rows survive the thresholds, the function relaxes to the
    best-per-gene from the unfiltered table so the report is never empty.
    """
    pct_cov_col = "%COVERAGE"
    pct_id_col  = "%IDENTITY"
    gene_col    = "GENE"

    # Coerce numeric columns (ABRicate sometimes writes 'NA')
    for col in (pct_cov_col, pct_id_col):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    filtered = df[
        (df[pct_cov_col].fillna(0) >= TOP_HIT_MIN_COVERAGE) &
        (df[pct_id_col].fillna(0)  >= TOP_HIT_MIN_IDENTITY)
    ].copy() if (pct_cov_col in df.columns and pct_id_col in df.columns) else df.copy()

    if filtered.empty:
        print(
            f"[WARN] No hits pass coverage≥{TOP_HIT_MIN_COVERAGE}% / "
            f"identity≥{TOP_HIT_MIN_IDENTITY}% thresholds — "
            "reporting best-per-gene from full hit set instead."
        )
        filtered = df.copy()

    if gene_col not in filtered.columns or filtered.empty:
        return filtered

    # Best per GENE: sort so highest %IDENTITY (then %COVERAGE) comes first,
    # then keep first occurrence in each group.
    sort_cols  = [c for c in (pct_id_col, pct_cov_col) if c in filtered.columns]
    if sort_cols:
        filtered = filtered.sort_values(sort_cols, ascending=False)

    top = filtered.drop_duplicates(subset=[gene_col], keep="first")

    if sort_cols:
        top = top.sort_values(sort_cols, ascending=False)

    return top.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Results summariser
# ---------------------------------------------------------------------------

def summarize_results(tsv_path: Path, out_dir: Path) -> Path:
    """
    Read the combined ABRicate TSV, produce:
      - amr_all_hits_annotated.csv   — every hit with MECHANISM column
      - amr_top_hits.csv             — best-per-gene filtered hits
      - amr_report.txt               — human-readable summary
    """
    # ----- Edge cases -----
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        report = out_dir / "amr_report.txt"
        report.write_text("No ABRicate hits found in any configured database.\n")
        return report

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)

    if df.empty:
        report = out_dir / "amr_report.txt"
        report.write_text(
            "There is no gene responsible for antibiotic resistance "
            "based on ABRicate databases.\n"
        )
        return report

    # ----- Normalise column names -----
    # ABRicate header row starts with '#FILE'; strip leading '#' for safer access
    df.columns = [c.lstrip("#") for c in df.columns]
    # Re-add '#FILE' as an alias so we keep the original name in outputs
    if "FILE" in df.columns and "#FILE" not in df.columns:
        df.insert(0, "#FILE", df["FILE"])

    # Ensure all expected columns are present (fill missing with 'NA')
    for col in ABRICATE_COLUMNS:
        clean = col.lstrip("#")
        if clean not in df.columns and col not in df.columns:
            df[col] = "NA"

    # Guard optional columns
    for col in ("RESISTANCE", "PRODUCT", "ACCESSION", "COVERAGE_MAP", "GAPS"):
        if col not in df.columns:
            df[col] = "Unknown"

    # ----- Add MECHANISM column -----
    df["MECHANISM"] = [
        infer_mechanism(
            str(p) if pd.notna(p) else "",
            str(r) if pd.notna(r) else "",
        )
        for p, r in zip(df.get("PRODUCT", [""]*len(df)), df.get("RESISTANCE", [""]*len(df)))
    ]

    # ----- Build keep_cols preserving full ABRicate column order -----
    ordered_cols = [
        "#FILE", "FILE",    # one or the other
        "SEQUENCE", "START", "END", "STRAND",
        "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS",
        "%COVERAGE", "%IDENTITY",
        "DATABASE", "ACCESSION",
        "PRODUCT", "RESISTANCE",
        "MECHANISM",
    ]
    keep_cols = [c for c in ordered_cols if c in df.columns]
    # De-duplicate while preserving order
    seen: set = set()
    keep_cols = [c for c in keep_cols if not (c in seen or seen.add(c))]  # type: ignore[func-returns-value]

    all_df = df[keep_cols].copy()

    # ----- Coerce numeric for sorting -----
    for col in ("%COVERAGE", "%IDENTITY"):
        if col in all_df.columns:
            all_df[col] = pd.to_numeric(all_df[col], errors="coerce")

    all_df_sorted = all_df.sort_values(
        by=[c for c in ("GENE", "DATABASE", "%IDENTITY") if c in all_df.columns],
        ascending=[True, True, False],
        na_position="last",
    )

    # ----- Top hits -----
    top_df = extract_top_hits(all_df.copy())

    # ----- Save files -----
    all_csv  = out_dir / "amr_all_hits_annotated.csv"
    top_csv  = out_dir / "amr_top_hits.csv"
    all_df_sorted.to_csv(all_csv, index=False)
    top_df.to_csv(top_csv, index=False)

    # ----- Build report -----
    report_lines = [
        "=" * 70,
        "  AMR ABRicate Analysis Report",
        "=" * 70,
        f"  Total raw hits (all databases)  : {len(all_df_sorted)}",
        f"  Unique genes detected           : "
        f"{all_df_sorted['GENE'].nunique() if 'GENE' in all_df_sorted.columns else 'NA'}",
        f"  Top hits (post-filter, best/gene): {len(top_df)}",
        f"  Coverage threshold applied      : ≥{TOP_HIT_MIN_COVERAGE}%",
        f"  Identity threshold applied      : ≥{TOP_HIT_MIN_IDENTITY}%",
        "",
    ]

    if all_df_sorted.empty:
        report_lines.append(
            "  There is no gene responsible for antibiotic resistance "
            "based on ABRicate databases."
        )
    else:
        # ---- Section 1: ALL hits table ----
        report_lines += [
            "-" * 70,
            "  SECTION 1 — ALL ABRicate Hits (full column set)",
            "-" * 70,
        ]
        display_all_cols = [
            c for c in (
                "FILE", "#FILE", "SEQUENCE", "START", "END", "STRAND",
                "GENE", "COVERAGE", "COVERAGE_MAP", "GAPS",
                "%COVERAGE", "%IDENTITY",
                "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE",
            )
            if c in all_df_sorted.columns
        ]
        report_lines.append(
            tabulate(
                all_df_sorted[display_all_cols].fillna("NA"),
                headers="keys",
                tablefmt="github",
                showindex=False,
            )
        )
        report_lines.append("")

        # ---- Section 2: TOP hits table ----
        report_lines += [
            "-" * 70,
            f"  SECTION 2 — TOP HITS  "
            f"(best per gene | %COV≥{TOP_HIT_MIN_COVERAGE} | %ID≥{TOP_HIT_MIN_IDENTITY})",
            "-" * 70,
        ]
        if top_df.empty:
            report_lines.append("  No hits passed the top-hit thresholds.")
        else:
            display_top_cols = [
                c for c in (
                    "GENE", "DATABASE", "ACCESSION",
                    "%COVERAGE", "%IDENTITY",
                    "COVERAGE", "COVERAGE_MAP", "GAPS",
                    "PRODUCT", "RESISTANCE", "MECHANISM",
                )
                if c in top_df.columns
            ]
            report_lines.append(
                tabulate(
                    top_df[display_top_cols].fillna("NA"),
                    headers="keys",
                    tablefmt="github",
                    showindex=False,
                )
            )
            report_lines.append("")

        # ---- Section 3: Mechanism annotations ----
        report_lines += [
            "-" * 70,
            "  SECTION 3 — Inferred Resistance Mechanisms (top hits only)",
            "-" * 70,
        ]
        if top_df.empty:
            report_lines.append("  No top hits available for mechanism annotation.")
        else:
            for i, (_, row) in enumerate(top_df.iterrows(), start=1):
                gene       = row.get("GENE", "NA")
                resistance = row.get("RESISTANCE", "Unknown")
                product    = row.get("PRODUCT", "Unknown")
                mechanism  = row.get("MECHANISM", "Unknown")
                pct_cov    = row.get("%COVERAGE", "NA")
                pct_id     = row.get("%IDENTITY", "NA")
                db         = row.get("DATABASE", "NA")
                accession  = row.get("ACCESSION", "NA")
                report_lines.append(
                    f"  [{i:02d}] Gene      : {gene}\n"
                    f"       Database  : {db}  |  Accession : {accession}\n"
                    f"       %Coverage : {pct_cov}  |  %Identity : {pct_id}\n"
                    f"       Class     : {resistance}\n"
                    f"       Product   : {product}\n"
                    f"       Mechanism : {mechanism}\n"
                )

    report_lines += [
        "=" * 70,
        f"  Full annotated table : {all_csv}",
        f"  Top hits table       : {top_csv}",
        "=" * 70,
    ]

    report_path = out_dir / "amr_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    return report_path


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    base_dir = Path.cwd() / "amr_colab_workspace"
    base_dir.mkdir(exist_ok=True)

    print("\n" + "=" * 70)
    print("  AMR Deep Analysis Pipeline  —  ABRicate + Google Colab")
    print("=" * 70)

    # ---- Step 1: Collect FASTA input ----
    fasta_path = get_fasta_from_user(base_dir)
    fasta_text = fasta_path.read_text()
    records = parse_fasta(fasta_text)

    if not records:
        raise ValueError("No valid FASTA records detected in the supplied file.")

    # ---- Step 2: Protein rejection + per-sequence QC ----
    qc_rows: List[SequenceQC] = []
    for seq_id, seq in records:
        if is_probable_protein(seq):
            raise ValueError(
                f"Sequence '{seq_id}' appears to be protein-like. "
                "Please provide a nucleotide genome FASTA only."
            )
        row = qc_record(seq_id, seq)
        if row.invalid_chars:
            raise ValueError(
                f"Sequence '{seq_id}' contains invalid non-IUPAC DNA characters: "
                f"{', '.join(row.invalid_chars)}"
            )
        qc_rows.append(row)

    # ---- Step 3: Write per-contig QC table ----
    qc_tsv = base_dir / "qc_metrics.tsv"
    with qc_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "seq_id", "length_bp", "gc_percent",
            "shannon_entropy", "complexity_k4", "non_standard_iupac_codes",
        ])
        for row in qc_rows:
            writer.writerow([
                row.seq_id,
                row.length,
                round(row.gc_percent, 4),
                round(row.shannon_entropy, 4),
                round(row.complexity_k4, 4),
                ",".join(row.non_standard_iupac) if row.non_standard_iupac else "None",
            ])

    total_len      = sum(r.length for r in qc_rows)
    weighted_gc    = sum(r.length * r.gc_percent for r in qc_rows) / total_len
    mean_entropy   = statistics.mean(r.shannon_entropy for r in qc_rows)
    mean_complexity = statistics.mean(r.complexity_k4 for r in qc_rows)

    print("\n[QC SUMMARY]")
    print(f"  Contigs / sequences : {len(qc_rows)}")
    print(f"  Total length        : {total_len:,} bp")
    print(f"  Weighted GC%%       : {weighted_gc:.2f}")
    print(f"  Mean Shannon entropy: {mean_entropy:.4f}")
    print(f"  Mean k4 complexity  : {mean_complexity:.4f}")
    print(f"  Per-contig QC table : {qc_tsv}")

    # ---- Step 4: Install ABRicate (micromamba) ----
    print("\n[INSTALL] Checking / installing ABRicate environment …")
    _, env = install_micromamba_and_abricate(base_dir)

    # ---- Step 5: Screen against all ABRicate databases ----
    results_dir = base_dir / "results"
    abricate_tsv = run_abricate_all_dbs(fasta_path, results_dir, env)

    # ---- Step 6: Summarise & produce reports ----
    report_path = summarize_results(abricate_tsv, results_dir)

    print("\n" + "=" * 70)
    print("  [AMR REPORT]")
    print("=" * 70)
    print(report_path.read_text())

    print("\n[OUTPUT FILES]")
    print(f"  QC metrics          : {qc_tsv}")
    print(f"  Raw ABRicate TSV    : {abricate_tsv}")
    print(f"  All hits (CSV)      : {results_dir / 'amr_all_hits_annotated.csv'}")
    print(f"  Top hits (CSV)      : {results_dir / 'amr_top_hits.csv'}")
    print(f"  Text report         : {report_path}")
    print("\n✅  Pipeline completed successfully.")


if __name__ == "__main__":
    main()
