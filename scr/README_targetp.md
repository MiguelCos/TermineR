# TargetP2 Wrapper Script

This repository provides `run_targetp2.sh`, a thin wrapper around the TargetP 2.0 binary deployed at `/home/shared/proteomics/targetp2/targetp-2.0/bin/targetp`. The script captures common defaults for the rat reference proteome FASTA while letting you override any TargetP flags per run.

## Prerequisites

- Bash 4+
- Access to the Linux environment where TargetP 2.0 is installed
- Permission to read the desired FASTA file(s) and write output files

## Installation

Copy the script to the Linux server (if you edited it locally) and make it executable:

```bash
chmod +x run_targetp2.sh
```

Optionally export defaults for the executable, FASTA, or prefix:

```bash
export TARGETP_BIN=/home/shared/proteomics/targetp2/targetp-2.0/bin/targetp
export TARGETP_FASTA=/path/to/your.fasta
export TARGETP_PREFIX=my_run
export TARGETP_OUT_DIR=/path/to/results
```

## Basic Usage

Invoke the script with any combination of the supported options; omit flags to keep the defaults shown below:

```bash
./run_targetp2.sh \
  --exec /home/shared/proteomics/targetp2/targetp-2.0/bin/targetp \
  --fasta /home/mcosenza/media/huesgen_share/0_people/mcosenza/fastas_and_organims_annotations/fastas/ebi_reference_proteomes_uniprot_format/UP000002494_10116_rat.fasta/UP000002494_10116_rat.fasta \
  --out-dir /home/mcosenza/targetp_results \
  --prefix rat_reference_run \
  --format long \
  --plot png \
  --gff3 --mature
```

Key flags map directly to the TargetP CLI:

| Wrapper flag | TargetP flag | Description |
| --- | --- | --- |
| `--fasta PATH` | `-fasta PATH` | Input FASTA file |
| `--batch INT` | `-batch INT` | Number of sequences processed simultaneously |
| `--format {short,long}` | `-format` | Output verbosity |
| `--org {non-pl,pl}` | `-org` | Organism type |
| `--plot {png,eps,none}` | `-plot` | Plot output format (when `long`) |
| `--out-dir DIR` | n/a (wrapper helper) | Directory where output files are written |
| `--prefix NAME` | `-prefix` | Output file prefix |
| `--tmp DIR` | `-tmp` | Temporary directory |
| `--gff3` / `--no-gff3` | `-gff3` | Emit processed GFF3 |
| `--mature` / `--no-mature` | `-mature` | Emit mature sequence FASTA |
| `--stdout` / `--no-stdout` | `-stdout` | Stream summary to STDOUT |
| `--quiet` | `-verbose=false` | Silence verbose logs |

Add `--dry-run` to print the fully-expanded TargetP command without executing it.

## Example Output Location

TargetP writes files into the working directory by default. Use `--out-dir` to choose a destination directory and `--prefix` to distinguish multiple runs. For example, `--out-dir /data/runs/2026-05-06 --prefix rat_reference_long` produces files under `/data/runs/2026-05-06`.

If you need to place temporary files elsewhere, pass `--tmp /scratch/$USER/targetp`.

## Troubleshooting

- "TargetP executable not found": Verify the path passed via `--exec` or `TARGETP_BIN`, and confirm the binary is executable.
- "FASTA file not found": Double-check the `--fasta` path and that the user account has read permissions.
- Need to see the command before running: use `--dry-run` to preview.

## Example execution

```
chmod +x run_targetp2.sh
./run_targetp2.sh \
  --exec /home/shared/proteomics/targetp2/targetp-2.0/bin/targetp \
  --tmp /home/mcosenza/media/huesgen_share/0_people/mcosenza/fastas_and_organims_annotations/fastas/ebi_reference_proteomes_uniprot_format/UP000002494_10116_rat.fasta/temp \
  --fasta /home/mcosenza/media/huesgen_share/0_people/mcosenza/fastas_and_organims_annotations/fastas/ebi_reference_proteomes_uniprot_format/UP000002494_10116_rat.fasta/UP000002494_10116_rat.fasta \
  --prefix rat_reference \
  --format short
```