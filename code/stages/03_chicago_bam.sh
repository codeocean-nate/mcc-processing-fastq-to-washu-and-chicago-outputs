#!/usr/bin/env bash
# Stage 3: CHiCAGO-compatible bam. CHiCAGO expects HiCUP-style input, so
# supplementary alignments (-F 2048) are dropped and sorting switches from
# coordinate to read-name order.
set -euo pipefail

in_bam=""; out_bam=""; threads="4"; tmpdir="/scratch/mcc/tmp"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-bam=*) in_bam="${1#*=}" ;;
    --out-bam=*) out_bam="${1#*=}" ;;
    --threads=*) threads="${1#*=}" ;;
    --tmpdir=*) tmpdir="${1#*=}" ;;
    -h|--help) sed -n '2,5p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
  shift
done
[[ -s "$in_bam" ]] || { echo "03_chicago_bam.sh: missing input bam $in_bam" >&2; exit 2; }
[[ -n "$out_bam" ]] || { echo "03_chicago_bam.sh: --out-bam is required" >&2; exit 2; }
mkdir -p "$tmpdir" "$(dirname "$out_bam")"

set -o pipefail
samtools view -@ "$threads" -Shu -F 2048 "$in_bam" \
  | samtools sort -n -T "$tmpdir/chicago_sort" --threads "$threads" -o "$out_bam" -
echo "wrote $out_bam"
