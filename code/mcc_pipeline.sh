#!/usr/bin/env bash
# MCC / Dovetail capture processing: fastq -> valid-pairs bam -> CHiCAGO-compatible
# bam -> chinput -> CHiCAGO interaction calls (WashU long-range text + tabix track,
# interBed/ibed, seqMonk, Rds). Follows the Dovetail capture analysis protocol.
set -euo pipefail

DATA_DIR="${DATA_DIR:-/data}"
RESULTS_DIR="${RESULTS_DIR:-/results}"
WORK_DIR="${WORK_DIR:-/scratch/mcc}"
CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

sample_name=""; species="human"; design_resolution="10kb"
chicago_cutoff="5"; min_mapq="40"; threads="$(nproc)"
filter_interactions="true"
export_formats="interBed,washU_text,seqMonk,washU_track"
r1_fastq=""; r2_fastq=""; reference_fasta=""; design_dir=""; probes_bed=""
chicago_settings_file=""

usage() {
  cat <<'USAGE'
MCC processing: fastq -> WashU + CHiCAGO outputs.
Inputs are discovered under /data when the matching flag is not given.
  --sample_name=NAME               output prefix (default: derived from the R1 file)
  --r1_fastq=PATH --r2_fastq=PATH  read 1 / read 2 fastq(.gz)
  --reference_fasta=PATH           reference fasta (bwa/faidx indexes built if absent)
  --design_dir=PATH                CHiCAGO design dir holding .rmap and .baitmap
  --probes_bed=PATH                padded probe bed for on-target / coverage QC
  --species=human|mouse            reported in the run summary (default human)
  --design_resolution=5kb|10kb|20kb  picks a design dir by name (default 10kb)
  --chicago_cutoff=N               CHiCAGO score cutoff (default 5)
  --min_mapq=N                     pairtools parse --min-mapq (default 40)
  --threads=N                      default: all available cores
  --export_formats=LIST            CHiCAGO --export-format list
  --filter_interactions=true|false  also write cis 10kb-2Mb, >=5 read ibed
  --chicago_settings_file=PATH      optional CHiCAGO settings file (e.g. a lower
                                   minNPerBait for shallow pilot libraries)
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample_name=*) sample_name="${1#*=}" ;;
    --r1_fastq=*) r1_fastq="${1#*=}" ;;
    --r2_fastq=*) r2_fastq="${1#*=}" ;;
    --reference_fasta=*) reference_fasta="${1#*=}" ;;
    --design_dir=*) design_dir="${1#*=}" ;;
    --probes_bed=*) probes_bed="${1#*=}" ;;
    --species=*) species="${1#*=}" ;;
    --design_resolution=*) design_resolution="${1#*=}" ;;
    --chicago_cutoff=*) chicago_cutoff="${1#*=}" ;;
    --min_mapq=*) min_mapq="${1#*=}" ;;
    --threads=*) threads="${1#*=}" ;;
    --export_formats=*) export_formats="${1#*=}" ;;
    --filter_interactions=*) filter_interactions="${1#*=}" ;;
    --chicago_settings_file=*) chicago_settings_file="${1#*=}" ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
  shift
done

log() { printf '[mcc %s] %s\n' "$(date -u +%H:%M:%S)" "$*" >&2; }

find_one() {
  local dir="$1"; shift; local pat found
  for pat in "$@"; do
    found=$(find -L "$dir" -type f -iname "$pat" 2>/dev/null | sort | head -1)
    [[ -n "$found" ]] && { printf '%s\n' "$found"; return 0; }
  done
  return 1
}

[[ -d "$DATA_DIR" ]] || { echo "Missing input directory $DATA_DIR" >&2; exit 2; }
mkdir -p "$RESULTS_DIR" "$WORK_DIR/tmp" "$WORK_DIR/ref"

# --- inputs -----------------------------------------------------------------
if [[ -z "$r1_fastq" ]]; then
  r1_fastq=$(find_one "$DATA_DIR" '*_R1*.fastq.gz' '*_R1*.fq.gz' '*R1*.fastq.gz' \
    '*_1.fastq.gz' '*_R1*.fastq' '*R1*.fastq' '*_1.fastq') || {
    echo "No read 1 fastq found under $DATA_DIR" >&2; exit 2; }
fi
if [[ -z "$r2_fastq" ]]; then
  for a in R1 r1 _1.; do
    b=${a/1/2}
    cand="${r1_fastq//$a/$b}"
    [[ "$cand" != "$r1_fastq" && -f "$cand" ]] && { r2_fastq="$cand"; break; }
  done
fi
[[ -n "$r2_fastq" && -f "$r2_fastq" ]] || { echo "No read 2 fastq matching $r1_fastq" >&2; exit 2; }

if [[ -z "$reference_fasta" ]]; then
  reference_fasta=$(find_one "$DATA_DIR" '*.fa' '*.fasta' '*.fna' '*.fa.gz' '*.fasta.gz') || {
    echo "No reference fasta found under $DATA_DIR" >&2; exit 2; }
fi
if [[ "$reference_fasta" == *.gz ]]; then
  plain="$WORK_DIR/ref/$(basename "${reference_fasta%.gz}")"
  [[ -s "$plain" ]] || { log "decompressing reference"; gunzip -c "$reference_fasta" > "$plain"; }
  reference_fasta="$plain"
fi

if [[ -z "$design_dir" ]]; then
  baitmap=$(find -L "$DATA_DIR" -type f -iname '*.baitmap' | grep -i "$design_resolution" | sort | head -1 || true)
  [[ -n "$baitmap" ]] || baitmap=$(find -L "$DATA_DIR" -type f -iname '*.baitmap' | sort | head -1 || true)
  [[ -n "$baitmap" ]] || { echo "No CHiCAGO .baitmap found under $DATA_DIR" >&2; exit 2; }
  design_dir=$(dirname "$baitmap")
fi
baitmap=$(find -L "$design_dir" -maxdepth 1 -type f -iname '*.baitmap' | sort | head -1)
rmap=$(find -L "$design_dir" -maxdepth 1 -type f -iname '*.rmap' | sort | head -1)
[[ -n "$baitmap" && -n "$rmap" ]] || { echo "design dir $design_dir needs a .rmap and a .baitmap" >&2; exit 2; }
[[ -n "$probes_bed" ]] || probes_bed=$(find_one "$DATA_DIR" '*probes*200bp*.bed' '*probes*.bed' || true)

if [[ -z "$sample_name" ]]; then
  sample_name=$(basename "$r1_fastq"); sample_name=${sample_name%%.*}
  sample_name=$(sed -E 's/([._-]R?1)$//; s/([._-]R1[._-].*)$//' <<<"$sample_name")
fi

nproc_io=$(( threads / 2 )); (( nproc_io < 1 )) && nproc_io=1

log "sample=$sample_name species=$species threads=$threads"
log "R1=$r1_fastq"; log "R2=$r2_fastq"; log "reference=$reference_fasta"
log "design_dir=$design_dir (rmap=$(basename "$rmap") baitmap=$(basename "$baitmap"))"

# --- reference indexes ------------------------------------------------------
ref="$reference_fasta"
if [[ ! -f "$ref.bwt" || ! -f "$ref.fai" ]]; then
  work_ref="$WORK_DIR/ref/$(basename "$ref")"
  [[ -e "$work_ref" ]] || ln -sfn "$ref" "$work_ref"
  for ext in amb ann bwt pac sa fai; do
    [[ -f "$ref.$ext" && ! -e "$work_ref.$ext" ]] && ln -sfn "$ref.$ext" "$work_ref.$ext"
  done
  [[ -f "$work_ref.fai" ]] || { log "indexing reference (samtools faidx)"; samtools faidx "$work_ref"; }
  [[ -f "$work_ref.bwt" ]] || { log "indexing reference (bwa index; this is slow but done once)"; bwa index "$work_ref"; }
  ref="$work_ref"
fi
genome_file="$WORK_DIR/ref/$(basename "$ref").genome"
cut -f1,2 "$ref.fai" > "$genome_file"

# --- stages -----------------------------------------------------------------
pt_bam="$RESULTS_DIR/${sample_name}.PT.bam"
pairs="$RESULTS_DIR/${sample_name}.mapped.pairs.gz"
stats="$RESULTS_DIR/${sample_name}.pairtools.stats.txt"
chicago_bam="$RESULTS_DIR/${sample_name}.chicago.bam"

log "stage 1/4 aligning reads and calling valid pairs"
bash "$CODE_DIR/stages/01_align_to_pairs.sh" --ref="$ref" --genome="$genome_file" \
  --r1="$r1_fastq" --r2="$r2_fastq" --threads="$threads" --nproc-io="$nproc_io" \
  --min-mapq="$min_mapq" --tmpdir="$WORK_DIR/tmp" --out-bam="$pt_bam" \
  --out-pairs="$pairs" --out-stats="$stats"

log "stage 2/4 library and target-enrichment QC"
bash "$CODE_DIR/stages/02_library_qc.sh" --stats="$stats" --bam="$pt_bam" \
  --probes-bed="$probes_bed" --threads="$threads" \
  --out-qc="$RESULTS_DIR/${sample_name}.qc_report.txt" --out-dir="$RESULTS_DIR/qc"

log "stage 3/4 building CHiCAGO-compatible bam"
bash "$CODE_DIR/stages/03_chicago_bam.sh" --in-bam="$pt_bam" --out-bam="$chicago_bam" \
  --threads="$threads" --tmpdir="$WORK_DIR/tmp"

log "stage 4/4 chinput + CHiCAGO interaction calling (WashU exports)"
bash "$CODE_DIR/stages/04_chicago_calls.sh" --chicago-bam="$chicago_bam" \
  --baitmap="$baitmap" --rmap="$rmap" --design-dir="$design_dir" \
  --cutoff="$chicago_cutoff" --export-formats="$export_formats" \
  --sample="$sample_name" --workdir="$WORK_DIR" --results="$RESULTS_DIR" \
  --filter="$filter_interactions" --settings-file="$chicago_settings_file"

{
  echo "sample: $sample_name"
  echo "species: $species"
  echo "read 1: $r1_fastq"
  echo "read 2: $r2_fastq"
  echo "reference: $reference_fasta"
  echo "design dir: $design_dir"
  echo "rmap: $(basename "$rmap")"
  echo "baitmap: $(basename "$baitmap")"
  echo "min mapq: $min_mapq"
  echo "chicago cutoff: $chicago_cutoff"
  echo "export formats: $export_formats"
  echo "threads: $threads"
} > "$RESULTS_DIR/${sample_name}.run_parameters.txt"

log "done; outputs in $RESULTS_DIR"
