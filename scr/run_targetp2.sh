#!/usr/bin/env bash
# Wrapper script to run TargetP 2.0 with convenient defaults and overrides.

set -euo pipefail

DEFAULT_TARGETP_BIN="/home/shared/proteomics/targetp2/targetp-2.0/bin/targetp"
DEFAULT_FASTA="/home/mcosenza/media/huesgen_share/0_people/mcosenza/fastas_and_organims_annotations/fastas/uniprot/UP000002494_10116_rat.fasta/uniprotkb_yeast_AND_model_organism_5592_2026_05_06.fasta"
DEFAULT_PREFIX="targetp_run"
DEFAULT_OUT_DIR=""

usage() {
  cat <<EOF
Usage: run_targetp2.sh [options]

Options:
  --exec PATH            Path to the TargetP binary (default: ${DEFAULT_TARGETP_BIN} or \$TARGETP_BIN env).
  --fasta PATH           Input FASTA file (default: rat proteome reference FASTA).
  --batch INT            Number of sequences processed simultaneously (default: 100).
  --format FORMAT        Output format: short|long (default: short).
  --org ORG              Organism flag: non-pl|pl (default: non-pl).
  --plot TYPE            Plot output when using long format: png|eps|none (default: png).
  --prefix NAME          Prefix for output files (default: targetp_run).
  --out-dir DIR          Output directory for generated files (default: current working directory).
  --tmp DIR              Directory for temporary files.
  --gff3 / --no-gff3     Enable/disable GFF3 output (default: disabled).
  --mature / --no-mature Enable/disable mature sequence FASTA output (default: disabled).
  --stdout / --no-stdout Print prediction summary to STDOUT (default: disabled).
  --verbose / --quiet    Toggle verbose logging (default: verbose).
  --dry-run              Print the command without executing it.
  -h, --help             Show this help message and exit.
  -- <args>              Everything after "--" is passed directly to TargetP.

Environment overrides:
  TARGETP_BIN      Override default TargetP executable path.
  TARGETP_FASTA    Default FASTA input file.
  TARGETP_PREFIX   Default prefix for output files.
  TARGETP_OUT_DIR  Default output directory for generated files.
EOF
}

# Defaults (can be overridden by env vars)
targetp_bin="${TARGETP_BIN:-$DEFAULT_TARGETP_BIN}"
fasta="${TARGETP_FASTA:-$DEFAULT_FASTA}"
prefix="${TARGETP_PREFIX:-$DEFAULT_PREFIX}"
out_dir="${TARGETP_OUT_DIR:-$DEFAULT_OUT_DIR}"
batch=100
format="short"
org="non-pl"
plot="png"
tmp_dir=""
gff3=false
mature=false
stdout_flag=false
verbose=true
dry_run=false
extra_args=()

if [[ $# -eq 0 ]]; then
  usage
  exit 0
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exec)
      targetp_bin="$2"; shift 2;;
    --fasta|-f)
      fasta="$2"; shift 2;;
    --batch)
      batch="$2"; shift 2;;
    --format)
      format="$2"; shift 2;;
    --org)
      org="$2"; shift 2;;
    --plot)
      plot="$2"; shift 2;;
    --prefix)
      prefix="$2"; shift 2;;
    --out-dir)
      out_dir="$2"; shift 2;;
    --tmp)
      tmp_dir="$2"; shift 2;;
    --gff3)
      gff3=true; shift;;
    --no-gff3)
      gff3=false; shift;;
    --mature)
      mature=true; shift;;
    --no-mature)
      mature=false; shift;;
    --stdout)
      stdout_flag=true; shift;;
    --no-stdout)
      stdout_flag=false; shift;;
    --verbose)
      verbose=true; shift;;
    --quiet|--no-verbose)
      verbose=false; shift;;
    --dry-run)
      dry_run=true; shift;;
    -h|--help)
      usage; exit 0;;
    --)
      shift
      extra_args+=("$@")
      break;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ ! -x "$targetp_bin" ]]; then
  echo "[ERROR] TargetP executable not found or not executable: $targetp_bin" >&2
  exit 1
fi

if [[ ! -f "$fasta" ]]; then
  echo "[ERROR] FASTA file not found: $fasta" >&2
  exit 1
fi

if ! [[ "$batch" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --batch must be a positive integer: $batch" >&2
  exit 1
fi

if [[ "$batch" -lt 1 ]]; then
  echo "[ERROR] --batch must be a positive integer: $batch" >&2
  exit 1
fi

effective_prefix="$prefix"
if [[ -n "$out_dir" ]]; then
  mkdir -p "$out_dir"
  effective_prefix="$out_dir/$prefix"
fi

cmd=("$targetp_bin" "-fasta" "$fasta" "-batch" "$batch" "-format" "$format" "-org" "$org" "-plot" "$plot")

if [[ -n "$effective_prefix" ]]; then
  cmd+=("-prefix" "$effective_prefix")
fi

if [[ -n "$tmp_dir" ]]; then
  cmd+=("-tmp" "$tmp_dir")
fi

if $gff3; then
  cmd+=("-gff3")
fi

if $mature; then
  cmd+=("-mature")
fi

if $stdout_flag; then
  cmd+=("-stdout")
fi

if ! $verbose; then
  cmd+=("-verbose=false")
fi

if (( ${#extra_args[@]} )); then
  cmd+=("${extra_args[@]}")
fi

echo "[INFO] Running TargetP from: $targetp_bin"
echo "[INFO] FASTA: $fasta"
if [[ -n "$out_dir" ]]; then
  echo "[INFO] Output directory: $out_dir"
else
  echo "[INFO] Output directory: current working directory"
fi
echo "[INFO] Output prefix: $effective_prefix"

if $dry_run; then
  printf '[DRY-RUN] %q ' "${cmd[@]}"
  printf '\n'
  exit 0
fi

"${cmd[@]}"
