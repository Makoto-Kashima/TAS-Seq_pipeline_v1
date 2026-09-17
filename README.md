# tas-seq-pipeline

A single-file, self-contained processing pipeline for TAS-Seq /
BD Rhapsody Enhanced-bead scRNA-seq data: cell-barcode & UMI
extraction, adapter/quality trimming, alignment, UMI deduplication,
per-cell demultiplexing, and transcript quantification.

> This README describes the pipeline as implemented in
> `tas-seq-pipeline` (the Python script in this repository). Adjust
> paths/examples below to your own environment before running.

## Overview

The pipeline replaces a multi-script workflow (separate barcode
extraction, SAM demultiplexing, SAM→BAM conversion, and salmon
quantification scripts) with one Python file that runs the full flow
end to end, with resumable per-step and per-file checkpoints:

```
extract (CB/UMI extraction + correction)
  -> fastp (poly-X / adapter trim)
  -> bwa mem (alignment) -> sorted, indexed BAM
  -> directional UMI dedup (self-implemented, parallel)
  -> demultiplex SAM by cell barcode
  -> SAM -> BAM (parallel)
  -> salmon quant (alignment-based, parallel, per cell)
  -> per-sample summary TSV
  -> cleanup of intermediates (keeps salmon/, fastp report, read counts)
```

## Requirements

External tools (must be on `PATH`):

- [`fastp`](https://github.com/OpenGene/fastp)
- [`bwa`](https://github.com/lh3/bwa) (`bwa mem`)
- [`samtools`](https://github.com/samtools/samtools)
- [`salmon`](https://github.com/COMBINE-lab/salmon)

Python: 3.8+, with [`pysam`](https://github.com/pysam-developers/pysam)
installed (only required for the UMI-dedup step).

Optional (used automatically if present, for faster gzip I/O):
`pigz` / `unpigz` / `igzip`, or the `isal` Python package.

## Installation

```bash
# Place the script anywhere on your PATH, e.g.:
cp tas-seq-pipeline /usr/local/bin/    # or any directory already on PATH
chmod +x /usr/local/bin/tas-seq-pipeline

pip install pysam
```

## Input files

- Paired-end FASTQ files, gzip-compressed, with `R1`/`R2` in the
  filenames (R1 = barcode/UMI read, R2 = cDNA read). Only the R1 path
  is passed on the command line; R2 is inferred by substituting
  `R1` -> `R2` in the filename.
- A reference FASTA (`-f/--ref`) used both as the `bwa mem` alignment
  target and as the `salmon quant -t` target.
- Three cell-label reference lists (`--cls1`, `--cls2`, `--cls3`), one
  per 9-nt barcode segment, one sequence per line (BD Rhapsody
  Enhanced-bead cell-label lists).
- (Optional) a custom adapter FASTA for fastp (`--adapter`).

## Usage

```bash
tas-seq-pipeline \
  -i sample_R1.fastq.gz \
  -f reference.fa \
  --cls1 BD_CLS1_enh.txt --cls2 BD_CLS2_enh.txt --cls3 BD_CLS3_enh.txt \
  --adapter TAS-adapter.fa \
  --threads 32
```

Output is written to a new `TAS-Seq_<sample>/` directory in the
current working directory.

### Options

| Flag | Default | Description |
|---|---|---|
| `-i, --input` | *(required)* | R1 FASTQ(.gz); R2 is inferred by replacing `R1`→`R2` in the filename |
| `-f, --ref` | *(required)* | Reference FASTA (bwa index prefix & salmon target) |
| `--cls1` / `--cls2` / `--cls3` | *(paths must be set for your environment)* | Cell-label reference lists for barcode segments 1/2/3 |
| `--adapter` | *(path must be set for your environment)* | fastp adapter FASTA |
| `--mismatch` | `1` | Max mismatches allowed when correcting each 9-nt barcode segment against its reference list |
| `--threads` | `nproc` | CPU threads (auto-detected) |
| `--min-len` | `31` | fastp minimum read length after trimming |
| `--chunk-size` | `100000` | Reads per chunk for the parallel extraction step |
| `--linker-mm` | off | Allow 1 mismatch in each linker sequence (GTGA/GACA) during extraction |
| `--keep-intermediates` | off | Keep intermediate BAM/SAM/gz files instead of deleting them at the end |
| `--resume` | off | Resume an interrupted run: skip steps/files already completed |

> The script ships with default paths for `--cls1/2/3` and
> `--adapter` pointing at a specific host's file layout — pass your
> own paths explicitly, or edit the defaults in the script for your
> deployment.

## Barcode structure (BD Rhapsody Enhanced beads)

```
R1: [0-3 nt spacer][9 nt CLS1][4 nt linker "GTGA"][9 nt CLS2][4 nt linker "GACA"][9 nt CLS3][8 nt UMI]
```

The pipeline locates the two linkers (exact match, or ≤1 mismatch
each with `--linker-mm`) to determine the spacer length, extracts the
27-nt cell barcode (CLS1+CLS2+CLS3) and 8-nt UMI, then independently
corrects each 9-nt segment against its reference list (exact match,
falling back to a single allowed mismatch; ambiguous 1-mismatch calls
are discarded). Reads whose linkers can't be located, or whose
segments can't be resolved, are dropped. The corrected barcode + UMI
are encoded into the R2 read name for all downstream steps.

## UMI deduplication

Duplicates are removed with a self-implemented **directional
deduplication** algorithm (Smith, Heger & Sudbery, *Genome Res.* 2017),
applied per `(reference contig, cell barcode)` group and parallelized
across groups of contigs:

- Build an edge UMI A → UMI B when `hamming(A, B) == 1` and
  `count(A) >= 2 * count(B) - 1`.
- Connected components (started from the highest-count "hub" UMI) are
  collapsed into one molecule.
- One representative read per molecule is kept, chosen deterministically
  by `(leftmost alignment position, read name)`.

This is algorithmically equivalent to
`umi_tools dedup --per-contig --per-gene --per-cell`, reimplemented in
Python/`pysam` for parallel execution and to avoid the `umi_tools`
dependency.

## Output

```
TAS-Seq_<sample>/
├── fastp.html, fastp.json         # trimming QC report
├── R2-extracted-read.txt          # read count per cell barcode (post-trim)
├── summary_<sample>.tsv           # per-sample QC summary (see below)
└── salmon/
    └── <cell_barcode>/
        └── quant.sf               # per-cell transcript quantification
```

Intermediate files (extracted/trimmed FASTQ, alignment BAM, dedup SAM,
per-cell SAM/BAM before quantification) are deleted after a successful
run unless `--keep-intermediates` is passed.

`summary_<sample>.tsv` reports: total reads, cell-barcode-assigned
reads (count & %), barcode-corrected reads, no-linker reads, discarded
reads, mapped reads entering dedup, unique molecules after dedup, UMI
conversion rate (%), duplication rate (%), number of cells, and mean
molecules per cell.

## Resuming an interrupted run

Re-run the exact same command with `--resume` added. Step-level
progress is tracked via marker files (`.step_<name>.done`) in the
sample's working directory; the SAM→BAM and salmon-quant steps
additionally resume at the level of individual files/cells.

## Citation

If you use this pipeline, please cite:

- Smith T, Heger A, Sudbery I. UMI-tools: modeling sequencing errors
  in Unique Molecular Identifiers to improve quantification accuracy.
  *Genome Res.* 2017;27(3):491-499. (directional deduplication method)
- Shichino S, Ueha S, Hashimoto S, et al. TAS-Seq is a robust and sensitive amplification method for bead-based scRNA-seq. Commun Biol. 2022;5:602. (terminator-assisted solid-phase cDNA amplification and sequencing method)
- Hasegawa M, Oshita M, Naruse K, et al. Oxytocin regulates TN-GnRH3 circuit maturation and mate preference through C1q-dependent synaptic mechanisms. bioRxiv [Preprint]. 2026. doi:10.64898/2026.05.17.725056.
Also cite the underlying tools: `fastp`, `bwa`, `samtools`, `salmon`.

## License

[Add a LICENSE file and state the license here, e.g. MIT, before
making the repository public.]
