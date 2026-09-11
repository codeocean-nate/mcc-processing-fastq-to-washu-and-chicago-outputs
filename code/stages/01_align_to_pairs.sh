#!/usr/bin/env bash
# Stage 1: fastq -> deduplicated valid-pairs bam (.PT.bam), .pairs and dup stats.
# Piped exactly as in the Dovetail capture protocol:
#   bwa mem -5SP -T0 | pairtools parse | sort | dedup | split | samtools sort
set -euo pipefail

ref=""; genome=""; r1=""; r2=""; threads="4"; nproc_io="2"; min_mapq="40"
tmpdir="/scratch/mcc/tmp"; out_bam=""; out_pairs=""; out_stats=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref=*) ref="${1#*=}" ;;
    --genome=*) genome="${1#*=}" ;;
    --r1=*) r1="${1#*=}" ;;
    --r2=*) r2="${1#*=}" ;;
    --threads=*) threads="${1#*=}" ;;
    --nproc-io=*) nproc_io="${1#*=}" ;;
    --min-mapq=*) min_mapq="${1#*=}" ;;
    --tmpdir=*) tmpdir="${1#*=}" ;;
    --out-bam=*) out_bam="${1#*=}" ;;
    --out-pairs=*) out_pairs="${1#*=}" ;;
    --out-stats=*) out_stats="${1#*=}" ;;
    -h|--help) sed -n '2,5p' "$0"; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
  shift
done

for v in ref genome r1 r2 out_bam out_pairs out_stats; do
  [[ -n "${!v}" ]] || { echo "01_align_to_pairs.sh: --${v//_/-} is required" >&2; exit 2; }
done
for f in "$ref" "$genome" "$r1" "$r2"; do
  [[ -s "$f" ]] || { echo "01_align_to_pairs.sh: missing or empty input $f" >&2; exit 2; }
done

mkdir -p "$tmpdir" "$(dirname "$out_bam")"
tmpdir=$(cd "$tmpdir" && pwd)   # pairtools sort requires an absolute path
rm -f "$out_stats"              # pairtools dedup appends to an existing stats file

set -o pipefail
bwa mem -5SP -T0 -t"$threads" "$ref" "$r1" "$r2" \
  | pairtools parse --min-mapq "$min_mapq" --walks-policy 5unique \
      --max-inter-align-gap 30 --nproc-in "$nproc_io" --nproc-out "$nproc_io" \
      --chroms-path "$genome" \
  | pairtools sort --tmpdir="$tmpdir" --nproc "$threads" \
  | pairtools dedup --nproc-in "$nproc_io" --nproc-out "$nproc_io" --mark-dups \
      --output-stats "$out_stats" \
  | pairtools split --nproc-in "$nproc_io" --nproc-out "$nproc_io" \
      --output-pairs "$out_pairs" --output-sam - \
  | samtools view -bS -@"$threads" \
  | samtools sort -@"$threads" -T "$tmpdir/sort" -o "$out_bam"
samtools index -@"$threads" "$out_bam"

echo "wrote $out_bam, $out_pairs, $out_stats"
