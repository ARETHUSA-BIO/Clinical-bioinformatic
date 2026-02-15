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
- Run ABRicate across resistance databases.
- Report resistance genes, putative antibiotic resistance classes,
  and inferred mechanism annotations.
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
from typing import Dict, List, Tuple


# --- Dependency bootstrap (pip for all Python modules) ---
REQUIRED_PY_MODULES = ["pandas"]


def pip_install_if_needed() -> None:
    """Install Python dependencies with pip inside Colab/runtime."""
    for module in REQUIRED_PY_MODULES:
        try:
            __import__(module)
        except ImportError:
            print(f"[INFO] Installing missing Python module via pip: {module}")
            subprocess.run([sys.executable, "-m", "pip", "install", "-q", module], check=True)


pip_install_if_needed()
import pandas as pd  # noqa: E402


# IUPAC DNA alphabet (allow N + ambiguity codes)
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
STANDARD_DNA = set("ACGTN")


@dataclass
class SequenceQC:
    seq_id: str
    length: int
    gc_percent: float
    shannon_entropy: float
    complexity_k4: float
    invalid_chars: List[str]
    non_standard_iupac: List[str]


def run_cmd(cmd: List[str], env: Dict[str, str] | None = None) -> None:
    """Run command and stream output."""
    print(f"[CMD] {' '.join(cmd)}")
    subprocess.run(cmd, check=True, env=env)


def detect_colab() -> bool:
    return "google.colab" in sys.modules


def install_micromamba_and_abricate(base_dir: Path) -> Tuple[Path, Dict[str, str]]:
    """Install micromamba, then ABRicate + dependencies in a conda env."""
    mm_dir = base_dir / "micromamba"
    mm_bin = mm_dir / "bin" / "micromamba"

    if not mm_bin.exists():
        mm_dir.mkdir(parents=True, exist_ok=True)
        run_cmd([
            "bash",
            "-lc",
            (
                "curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest "
                "| tar -xvj -C {dest} bin/micromamba"
            ).format(dest=str(mm_dir)),
        ])

    env_prefix = base_dir / "envs" / "amr-abricate"
    if not env_prefix.exists():
        run_cmd([
            str(mm_bin),
            "create",
            "-y",
            "-p",
            str(env_prefix),
            "-c",
            "conda-forge",
            "-c",
            "bioconda",
            "abricate",
            "blast",
            "seqkit",
        ])

    env = os.environ.copy()
    env["PATH"] = f"{env_prefix / 'bin'}:{env.get('PATH', '')}"

    # Ensure ABRicate databases are set up.
    run_cmd(["abricate", "--setupdb"], env=env)
    return env_prefix, env


def parse_fasta(text: str) -> List[Tuple[str, str]]:
    """Simple FASTA parser returning list of (header, sequence)."""
    records = []
    header = None
    seq_chunks: List[str] = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):  # new record
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
    """Reject clear protein sequences based on amino acid-only characters/frequency."""
    if not seq:
        return False
    aa_only = set("EFILPQZXJOU")
    if any(ch in aa_only for ch in seq):
        return True

    dna_count = sum(ch in IUPAC_DNA for ch in seq)
    aa_count = sum(ch in set("ACDEFGHIKLMNPQRSTVWY") for ch in seq)
    dna_fraction = dna_count / len(seq)
    aa_fraction = aa_count / len(seq)

    # if sequence looks much more like protein than DNA
    return dna_fraction < 0.85 and aa_fraction > 0.90


def shannon_entropy(seq: str) -> float:
    counts = Counter(seq)
    total = len(seq)
    if total == 0:
        return 0.0
    entropy = 0.0
    for n in counts.values():
        p = n / total
        entropy -= p * math.log2(p)
    return entropy


def kmer_complexity(seq: str, k: int = 4) -> float:
    if len(seq) < k:
        return 0.0
    kmers = {seq[i : i + k] for i in range(len(seq) - k + 1)}
    max_possible = min(4**k, len(seq) - k + 1)
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


def get_fasta_from_user(workdir: Path) -> Path:
    """Ask user to upload or paste FASTA in Colab; fallback to local path prompt."""
    print("\nChoose FASTA input method:")
    print("1) Upload genome FASTA file")
    print("2) Paste FASTA content")
    choice = input("Enter 1 or 2 [default 1]: ").strip() or "1"

    fasta_path = workdir / "input_genome.fasta"

    if choice == "2":
        print("Paste your FASTA content. End with a single line containing only: END")
        lines = []
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

        print("Please upload a genome FASTA file (.fa/.fasta/.fna)")
        uploaded = files.upload()
        if not uploaded:
            raise RuntimeError("No file uploaded.")

        first_name = next(iter(uploaded))
        content = uploaded[first_name]
        fasta_path.write_bytes(content)
        return fasta_path

    # Non-colab fallback
    local_path = input("Enter local path to FASTA file: ").strip()
    src = Path(local_path)
    if not src.exists():
        raise FileNotFoundError(f"File not found: {src}")
    shutil.copy(src, fasta_path)
    return fasta_path


def run_abricate_all_dbs(fasta_path: Path, out_dir: Path, env: Dict[str, str]) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)

    db_list_cmd = subprocess.run(
        ["abricate", "--list"],
        env=env,
        check=True,
        text=True,
        capture_output=True,
    )

    dbs = []
    for line in db_list_cmd.stdout.splitlines()[1:]:
        parts = re.split(r"\s+", line.strip())
        if parts and parts[0]:
            dbs.append(parts[0])

    combined_tsv = out_dir / "abricate_combined.tsv"
    with combined_tsv.open("w") as out_fh:
        header_written = False
        for db in dbs:
            print(f"[INFO] Running ABRicate DB: {db}")
            run = subprocess.run(
                ["abricate", "--db", db, str(fasta_path)],
                env=env,
                check=True,
                text=True,
                capture_output=True,
            )
            lines = run.stdout.strip().splitlines()
            if not lines:
                continue
            if not header_written:
                out_fh.write(lines[0] + "\n")
                header_written = True
            for row in lines[1:]:
                if row.strip():
                    out_fh.write(row + "\n")

    return combined_tsv


def infer_mechanism(product: str, resistance: str) -> str:
    text = f"{product} {resistance}".lower()
    rules = [
        ("beta-lactamase", "Enzymatic inactivation of beta-lactam antibiotics by hydrolysis."),
        ("aminoglycoside", "Enzymatic modification or target alteration reducing aminoglycoside binding."),
        ("efflux", "Active efflux pump lowers intracellular antibiotic concentration."),
        ("tet", "Ribosomal protection protein or efflux-mediated tetracycline resistance."),
        ("qnr", "Protection of DNA gyrase/topoisomerase from quinolone inhibition."),
        ("mec", "Altered penicillin-binding protein with reduced beta-lactam affinity."),
        ("van", "Cell wall precursor remodeling reduces glycopeptide binding."),
        ("sul", "Alternative dihydropteroate synthase confers sulfonamide resistance."),
        ("dfr", "Alternative dihydrofolate reductase confers trimethoprim resistance."),
        ("erm", "rRNA methylation reduces macrolide/lincosamide/streptogramin binding."),
        ("bla", "Beta-lactam hydrolysis by beta-lactamase enzyme."),
    ]
    for key, annotation in rules:
        if key in text:
            return annotation
    return "Putative AMR determinant detected; specific mechanism requires manual curation."


def summarize_results(tsv_path: Path, out_dir: Path) -> Path:
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        report = out_dir / "amr_report.txt"
        report.write_text("No ABRicate hits found in any configured database.\n")
        return report

    df = pd.read_csv(tsv_path, sep="\t")
    if df.empty:
        report = out_dir / "amr_report.txt"
        report.write_text("There is no gene responsible for antibiotic resistance based on ABRicate databases.\n")
        return report

    if "RESISTANCE" not in df.columns:
        df["RESISTANCE"] = "Unknown"
    if "PRODUCT" not in df.columns:
        df["PRODUCT"] = "Unknown"

    df["MECHANISM"] = [infer_mechanism(p, r) for p, r in zip(df["PRODUCT"], df["RESISTANCE"])]

    keep_cols = [c for c in ["#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE", "PRODUCT", "RESISTANCE", "%COVERAGE", "%IDENTITY", "DATABASE", "MECHANISM"] if c in df.columns]
    final_df = df[keep_cols].sort_values(by=["GENE", "DATABASE"], na_position="last")

    csv_path = out_dir / "amr_hits_annotated.csv"
    final_df.to_csv(csv_path, index=False)

    report_lines = [
        "AMR ABRicate Analysis Report",
        "=" * 30,
        f"Total AMR hits: {len(final_df)}",
        f"Unique genes: {final_df['GENE'].nunique() if 'GENE' in final_df.columns else 'NA'}",
        "",
    ]

    if len(final_df) == 0:
        report_lines.append("There is no gene responsible for antibiotic resistance based on ABRicate databases.")
    else:
        report_lines.append("Detected resistance genes and annotations:")
        for _, row in final_df.iterrows():
            report_lines.append(
                f"- Gene: {row.get('GENE', 'NA')} | Antibiotic/Class: {row.get('RESISTANCE', 'Unknown')} "
                f"| Product: {row.get('PRODUCT', 'Unknown')} | Mechanism: {row.get('MECHANISM', 'Unknown')}"
            )

    report_path = out_dir / "amr_report.txt"
    report_path.write_text("\n".join(report_lines) + "\n")
    return report_path


def main() -> None:
    base_dir = Path.cwd() / "amr_colab_workspace"
    base_dir.mkdir(exist_ok=True)

    print("\n=== AMR Deep Analysis Pipeline (ABRicate + Colab) ===")
    fasta_path = get_fasta_from_user(base_dir)
    fasta_text = fasta_path.read_text()
    records = parse_fasta(fasta_text)

    if not records:
        raise ValueError("No valid FASTA records detected.")

    # Protein sequence rejection and QC checks
    qc_rows: List[SequenceQC] = []
    for seq_id, seq in records:
        if is_probable_protein(seq):
            raise ValueError(
                f"Sequence '{seq_id}' appears to be protein-like. Please provide nucleotide genome FASTA only."
            )
        row = qc_record(seq_id, seq)
        if row.invalid_chars:
            raise ValueError(
                f"Sequence '{seq_id}' contains invalid non-IUPAC DNA characters: {','.join(row.invalid_chars)}"
            )
        qc_rows.append(row)

    qc_tsv = base_dir / "qc_metrics.tsv"
    with qc_tsv.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "seq_id",
            "length",
            "gc_percent",
            "shannon_entropy",
            "complexity_k4",
            "non_standard_iupac_codes",
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

    total_len = sum(r.length for r in qc_rows)
    weighted_gc = sum(r.length * r.gc_percent for r in qc_rows) / total_len
    mean_entropy = statistics.mean(r.shannon_entropy for r in qc_rows)
    mean_complexity = statistics.mean(r.complexity_k4 for r in qc_rows)

    print("\n[QC SUMMARY]")
    print(f"Contigs/sequences: {len(qc_rows)}")
    print(f"Total length: {total_len:,} bp")
    print(f"Weighted GC%: {weighted_gc:.2f}")
    print(f"Mean Shannon entropy: {mean_entropy:.3f}")
    print(f"Mean complexity (k=4): {mean_complexity:.3f}")
    print(f"Detailed QC table: {qc_tsv}")

    # Install ABRicate and run analysis
    _, env = install_micromamba_and_abricate(base_dir)
    results_dir = base_dir / "results"
    abricate_tsv = run_abricate_all_dbs(fasta_path, results_dir, env)
    report_path = summarize_results(abricate_tsv, results_dir)

    print("\n[AMR OUTPUT]")
    print(report_path.read_text())
    print(f"Annotated table (if hits): {results_dir / 'amr_hits_annotated.csv'}")
    print("\nPipeline completed successfully.")


if __name__ == "__main__":
    main()
