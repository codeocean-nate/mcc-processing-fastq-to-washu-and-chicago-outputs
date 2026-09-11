#!/usr/bin/env bash
# Stage 2: library QC from the pairtools stats file, plus target-enrichment QC
# (on-target read pairs and probe coverage) when a probe bed file is available.
set -euo pipefail

stats=""; bam=""; probes_bed=""; threads="4"; out_qc=""; out_dir=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --stats=*) stats="${1#*=}" ;;
    --bam=*) bam="${1#*=}" ;;
    --probes-bed=*) probes_bed="${1#*=}" ;;
    --threads=*) threads="${1#*=}" ;;
    --out-qc=*) out_qc="${1#*=}" ;;
    --out-dir=*) out_dir="${1#*=}" ;;
    -h|--help) sed -n '2,4p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
  shift
done
for v in stats bam out_qc out_dir; do
  [[ -n "${!v}" ]] || { echo "02_library_qc.sh: --${v//_/-} is required" >&2; exit 2; }
done
[[ -s "$stats" ]] || { echo "02_library_qc.sh: missing stats file $stats" >&2; exit 2; }
mkdir -p "$out_dir" "$(dirname "$out_qc")"

qc_script="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/vendor/get_qc.py"
python3 "$qc_script" -p "$stats" > "$out_qc"

if [[ -n "$probes_bed" && -s "$probes_bed" ]]; then
  echo "" >> "$out_qc"
  echo "Target enrichment (probes: $(basename "$probes_bed"))" >> "$out_qc"
  on_target=$(samtools view "$bam" -L "$probes_bed" -@ "$threads" \
    | awk -F '\t' '{print "@"$1}' | sort -u | wc -l)
  nodup=$(awk '/^No-Dup Read Pairs/ {gsub(",","",$4); print $4; exit}' "$out_qc")
  printf 'On-Target Read Pairs %s\n' "$on_target" >> "$out_qc"
  if [[ -n "${nodup:-}" && "$nodup" -gt 0 ]]; then
    awk -v a="$on_target" -v b="$nodup" 'BEGIN{printf "On-Target Rate %.2f%%\n", (a/b)*100}' >> "$out_qc"
  fi
  if command -v mosdepth >/dev/null 2>&1; then
    mosdepth -t "$threads" -b "$probes_bed" -x "$out_dir/probe_coverage" -n "$bam"
    echo "" >> "$out_qc"
    echo "Coverage depth (mosdepth summary)" >> "$out_qc"
    head -1 "$out_dir/probe_coverage.mosdepth.summary.txt" >> "$out_qc"
    tail -2 "$out_dir/probe_coverage.mosdepth.summary.txt" >> "$out_qc"
  fi
else
  echo "" >> "$out_qc"
  echo "Target enrichment QC skipped: no probe bed file supplied." >> "$out_qc"
fi
echo "wrote $out_qc"
