# Clinical Bioinformatics AMR Pipeline (ABRicate)

Welcome! This project is a **student-friendly clinical bioinformatics pipeline** for detecting antimicrobial resistance (AMR) genes from **bacterial whole-genome FASTA** files using **ABRicate**.

## Why this project matters in clinical bioinformatics

In clinical bioinformatics, we analyze microbial genomes to support:
- faster understanding of infection risks,
- detection of resistance genes,
- better antibiotic decision support in research and healthcare.

AMR is a major global health challenge. A genome-based AMR workflow helps students and researchers learn how computational analysis can identify genes linked to antibiotic resistance.

## What is ABRicate?

**ABRicate** is a genome screening tool that compares assembled contigs/genomes against curated resistance gene databases (for example CARD, ResFinder, NCBI AMR, and others depending on setup).

This pipeline automates ABRicate usage and adds quality-control and interpretation steps.

## What this pipeline does

The script: `amr_abricate_colab_pipeline.py`

### 1) Input from user
- Asks user to either:
  - upload a genome FASTA file, or
  - paste FASTA text directly.
- Accepts bacterial **nucleotide** genome sequences only.
- Rejects probable **protein** sequences.

### 2) QC + IUPAC validation
Before AMR screening, the pipeline checks sequence quality:
- FASTA format parsing,
- IUPAC DNA validity,
- rejection of invalid characters.

### 3) Genome metrics
For each sequence/contig, it calculates:
- sequence length,
- GC content (%),
- Shannon entropy,
- k-mer complexity (k=4).

Outputs QC table:
- `amr_colab_workspace/qc_metrics.tsv`

### 4) ABRicate AMR screening
- Installs ABRicate environment (micromamba/conda + bioconda packages).
- Runs `abricate --setupdb` to configure databases.
- Scans the genome across available ABRicate databases.

### 5) AMR interpretation report
If resistance genes are detected, it reports:
- gene name,
- resistance class / antibiotic association,
- product annotation,
- inferred mechanism of resistance.

If no AMR gene is detected, it returns a clear no-hit message.

Main outputs:
- `amr_colab_workspace/results/abricate_combined.tsv`
- `amr_colab_workspace/results/amr_hits_annotated.csv`
- `amr_colab_workspace/results/amr_report.txt`

## Quick start (Google Colab)

1. Upload `amr_abricate_colab_pipeline.py` to Colab.
2. Run:

```bash
!python amr_abricate_colab_pipeline.py
```

3. Choose upload/paste input and follow prompts.
4. Review QC and AMR outputs in `amr_colab_workspace/`.

## Notes for students

- Use **assembled bacterial genome FASTA** (not proteins).
- A detected gene suggests potential resistance, but final clinical interpretation requires expert review and phenotype correlation.
- Database coverage and update status affect results.

## Author

Prepared for students by **Eng. Taha Bilel Chalbi**.
