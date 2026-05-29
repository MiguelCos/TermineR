#!/usr/bin/env bash
# Download UniProt FASTA files and run TargetP 2.0 locally on Linux.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUNNER="${SCRIPT_DIR}/run_targetp2.sh"
OUT_ROOT="${REPO_DIR}/data-raw/targetp"
FORMAT="short"
BATCH="100"
TMP_DIR=""
DRY_RUN=false
FORCE_DOWNLOAD=false
RUN_PREPARE=false
ORGANISMS=()

# ORGANISM REGISTRY -----
# Columns: organism | taxonomy_id | targetp_org | output_prefix
TARGETP_REGISTRY=$(cat <<'EOF'
medicago_trucantula|3880|pl|medicago_trucantula_reference
rhizobium_melitoli|266834|non-pl|rhizobium_melitoli_reference
pig|9823|non-pl|pig_reference
human_iso|9606|non-pl|human_iso_reference
ecoli|83333|non-pl|ecoli_reference
c_elegans|6239|non-pl|c_elegans_reference
synechocystis|1111708|non-pl|synechocystis_reference
EOF
)

usage() {
  cat <<EOF
Usage: scr/run_targetp_all_linux.sh [options]

Downloads UniProt FASTA files into data-raw/targetp/<organism>/ and runs
scr/run_targetp2.sh for each selected organism.

Options:
  --organism NAME       Run only one organism. Can be supplied multiple times.
  --out-root DIR        Output root directory (default: data-raw/targetp).
  --format FORMAT       TargetP output format: short|long (default: short).
  --batch INT           Number of sequences processed simultaneously (default: 100).
  --tmp DIR             TargetP temporary directory.
  --targetp-bin PATH    TargetP executable path, exported as TARGETP_BIN.
  --force-download      Re-download FASTA even if it already exists.
  --prepare             Run scr/prepare_targetp_processing_annotations.R after TargetP.
  --dry-run             Print commands without executing them.
  -h, --help            Show this help message.

Default organisms:
  medicago_trucantula rhizobium_melitoli pig human_iso ecoli c_elegans synechocystis

The TargetP executable is normally resolved by scr/run_targetp2.sh. Override it
with --targetp-bin or by setting TARGETP_BIN before running this script.
EOF
}

registry_row() {
  local organism_name="$1"
  local row

  row="$(printf '%s\n' "$TARGETP_REGISTRY" | awk -F'|' -v organism="$organism_name" '$1 == organism {print; exit}')"

  if [[ -z "$row" ]]; then
    echo "[ERROR] Organism not found in TargetP Linux registry: $organism_name" >&2
    exit 1
  fi

  printf '%s\n' "$row"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --organism)
      ORGANISMS+=("$2"); shift 2;;
    --out-root)
      OUT_ROOT="$2"; shift 2;;
    --format)
      FORMAT="$2"; shift 2;;
    --batch)
      BATCH="$2"; shift 2;;
    --tmp)
      TMP_DIR="$2"; shift 2;;
    --targetp-bin)
      export TARGETP_BIN="$2"; shift 2;;
    --force-download)
      FORCE_DOWNLOAD=true; shift;;
    --prepare)
      RUN_PREPARE=true; shift;;
    --dry-run)
      DRY_RUN=true; shift;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "[ERROR] Unknown option: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ ${#ORGANISMS[@]} -eq 0 ]]; then
  ORGANISMS=(
    "medicago_trucantula"
    "rhizobium_melitoli"
    "pig"
    "human_iso"
    "ecoli"
    "c_elegans"
    "synechocystis"
  )
fi

run_cmd() {
  if $DRY_RUN; then
    printf '[DRY-RUN]'
    printf ' %q' "$@"
    printf '\n'
  else
    "$@"
  fi
}

download_fasta() {
  local url="$1"
  local fasta="$2"

  if [[ -f "$fasta" && "$FORCE_DOWNLOAD" == false ]]; then
    echo "[INFO] FASTA exists, skipping download: $fasta"
    return 0
  fi

  run_cmd mkdir -p "$(dirname "$fasta")"

  if command -v curl >/dev/null 2>&1; then
    run_cmd curl -L --fail --output "$fasta" "$url"
  elif command -v wget >/dev/null 2>&1; then
    run_cmd wget -O "$fasta" "$url"
  else
    echo "[ERROR] Neither curl nor wget is available for FASTA download." >&2
    exit 1
  fi
}

if [[ ! -f "$RUNNER" ]]; then
  echo "[ERROR] TargetP runner not found: $RUNNER" >&2
  exit 1
fi

if ! [[ "$BATCH" =~ ^[0-9]+$ ]]; then
  echo "[ERROR] --batch must be a positive integer: $BATCH" >&2
  exit 1
fi

if [[ "$BATCH" -lt 1 ]]; then
  echo "[ERROR] --batch must be a positive integer: $BATCH" >&2
  exit 1
fi

run_cmd chmod +x "$RUNNER"
run_cmd mkdir -p "$OUT_ROOT"

for organism in "${ORGANISMS[@]}"; do
  registry_info="$(registry_row "$organism")"
  IFS='|' read -r organism_name taxonomy_id targetp_org prefix <<< "$registry_info"
  organism_dir="${OUT_ROOT}/${organism}"
  fasta="${organism_dir}/${organism}_uniprotkb_taxonomy_${taxonomy_id}.fasta"
  url="https://rest.uniprot.org/uniprotkb/stream?download=true&format=fasta&query=%28%28taxonomy_id%3A${taxonomy_id}%29%29"

  echo "[INFO] Processing ${organism}"
  echo "[INFO] TargetP org: ${targetp_org}"
  echo "[INFO] UniProt URL: ${url}"

  download_fasta "$url" "$fasta"
  run_cmd mkdir -p "$organism_dir"

  targetp_args=(
    "--fasta" "$fasta"
    "--out-dir" "$organism_dir"
    "--prefix" "$prefix"
    "--format" "$FORMAT"
    "--batch" "$BATCH"
    "--org" "$targetp_org"
    "--quiet"
  )

  if [[ -n "$TMP_DIR" ]]; then
    targetp_args+=("--tmp" "$TMP_DIR")
  fi

  run_cmd "$RUNNER" "${targetp_args[@]}"
done

if $RUN_PREPARE; then
  if command -v Rscript >/dev/null 2>&1; then
    run_cmd Rscript "${REPO_DIR}/scr/prepare_targetp_processing_annotations.R"
  else
    echo "[ERROR] Rscript is not available; cannot run the TargetP prep script." >&2
    exit 1
  fi
fi

echo "[INFO] TargetP batch finished."
