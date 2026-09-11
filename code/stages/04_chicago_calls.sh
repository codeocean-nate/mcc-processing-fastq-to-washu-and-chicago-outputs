#!/usr/bin/env bash
# Stage 4: chinput generation (bam2chicago.sh) and CHiCAGO interaction calling
# (runChicago.R), exporting WashU long-range text + tabix track, interBed/ibed,
# seqMonk and the .Rds database. Optionally writes a filtered ibed.
set -euo pipefail

chicago_bam=""; baitmap=""; rmap=""; design_dir=""; cutoff="5"
export_formats="interBed,washU_text,seqMonk,washU_track"
sample=""; workdir="/scratch/mcc"; results="/results"; filter="true"
settings_file=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chicago-bam=*) chicago_bam="${1#*=}" ;;
    --baitmap=*) baitmap="${1#*=}" ;;
    --rmap=*) rmap="${1#*=}" ;;
    --design-dir=*) design_dir="${1#*=}" ;;
    --cutoff=*) cutoff="${1#*=}" ;;
    --export-formats=*) export_formats="${1#*=}" ;;
    --sample=*) sample="${1#*=}" ;;
    --workdir=*) workdir="${1#*=}" ;;
    --results=*) results="${1#*=}" ;;
    --filter=*) filter="${1#*=}" ;;
    --settings-file=*) settings_file="${1#*=}" ;;
    -h|--help) sed -n '2,5p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
  shift
done
for v in chicago_bam baitmap rmap design_dir sample; do
  [[ -n "${!v}" ]] || { echo "04_chicago_calls.sh: --${v//_/-} is required" >&2; exit 2; }
done
[[ -s "$chicago_bam" ]] || { echo "04_chicago_calls.sh: missing bam $chicago_bam" >&2; exit 2; }

tools="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/vendor/chicagoTools"

# CHiCAGO (bioconductor-chicago) lives in its own conda env so its R stack does not
# clash with the base python environment.
RSCRIPT="${CHICAGO_RSCRIPT:-/opt/conda/envs/chicago-r/bin/Rscript}"
[[ -x "$RSCRIPT" ]] || RSCRIPT="$(command -v Rscript || true)"
[[ -n "$RSCRIPT" ]] || { echo "04_chicago_calls.sh: no Rscript with the Chicago package found" >&2; exit 2; }
mkdir -p "$workdir/chinput" "$results"

chin_prefix="$workdir/chinput/${sample}_chinput"
rm -rf "$chin_prefix"
bash "$tools/bam2chicago.sh" "$chicago_bam" "$baitmap" "$rmap" "$chin_prefix"

chinput=$(find "$chin_prefix" -maxdepth 1 -type f -name '*.chinput' | sort | head -1)
[[ -n "$chinput" ]] || { echo "bam2chicago.sh produced no chinput file" >&2; exit 3; }
cp -f "$chinput" "$results/${sample}.chinput"

out_name="${sample}_chicago_calls"
out_prefix="$results/$out_name"
rm -rf "$out_prefix"
# CHiCAGO writes design-file caches next to the design dir, so use a writable copy.
design_work="$workdir/design/$(basename "$design_dir")"
mkdir -p "$design_work"
cp -f "$design_dir"/* "$design_work"/ 2>/dev/null || true

settings_args=()
if [[ -n "$settings_file" ]]; then
  [[ -s "$settings_file" ]] || { echo "04_chicago_calls.sh: settings file $settings_file not found" >&2; exit 2; }
  settings_args=(--settings-file "$settings_file")
fi

# runChicago.R requires a bare output prefix plus an --output-dir.
"$RSCRIPT" "$tools/runChicago.R" --design-dir "$design_work" --output-dir "$out_prefix" \
  --cutoff "$cutoff" --export-format "$export_formats" \
  "${settings_args[@]+"${settings_args[@]}"}" \
  "$chinput" "$out_name"

ibed=$(find "$out_prefix" -type f -name '*.ibed' | sort | head -1 || true)
if [[ "$filter" == "true" && -n "$ibed" ]]; then
  awk 'NR>1 && $1 == $5 && \
    (($6 > $3 && ($6 - $3) < 2000000) || ($6 < $3 && ($2 - $7) < 2000000)) && \
    (($6 > $3 && ($6 - $3) >   10000) || ($6 < $3 && ($2 - $7) >   10000)) && \
    $9 >= 5 { print }' "$ibed" > "$results/${sample}.filtered_interactions.ibed"
  echo "wrote $results/${sample}.filtered_interactions.ibed"
fi

echo "CHiCAGO outputs:"
find "$out_prefix" -type f | sort
