# MCC processing: FASTQ to WashU and CHiCAGO outputs

Implements the Dovetail capture (Micro-Capture-C / capture Micro-C) analysis
protocol end to end: paired FASTQ reads in, promoter-interaction calls out,
including the WashU browser tracks and the CHiCAGO files.

## Stages

1. `stages/01_align_to_pairs.sh` — `bwa mem -5SP -T0` piped through
   `pairtools parse` (`--min-mapq 40 --walks-policy 5unique
   --max-inter-align-gap 30`) → `sort` → `dedup --mark-dups` → `split`,
   then coordinate-sorted and indexed. Outputs `*.PT.bam`, `*.mapped.pairs.gz`,
   `*.pairtools.stats.txt`.
2. `stages/02_library_qc.sh` — proximity-ligation QC summary from the pairtools
   stats file; on-target read-pair rate and probe coverage depth (mosdepth) when
   a probe BED file is supplied.
3. `stages/03_chicago_bam.sh` — CHiCAGO-compatible BAM: supplementary
   alignments removed (`-F 2048`) and read-name sorted, as CHiCAGO expects
   HiCUP-style input.
4. `stages/04_chicago_calls.sh` — `bam2chicago.sh` builds the `.chinput` from the
   design `.rmap`/`.baitmap`, then `runChicago.R` calls significant interactions
   and exports `interBed`, `washU_text`, `washU_track` (bgzipped + `.tbi`) and
   `seqMonk`, plus the `.Rds` database, diagnostic plots and example plots.

## Inputs (attach as data assets under `/data`)

| Input | Discovery pattern | Notes |
| --- | --- | --- |
| Paired FASTQ | `*_R1*.fastq.gz` / `*_R2*.fastq.gz` | gzipped or plain |
| Reference FASTA | `*.fa`, `*.fasta`, `*.fna` (`.gz` accepted) | `bwa index` / `samtools faidx` built into `/scratch` if absent |
| CHiCAGO design dir | any dir holding `*.rmap` + `*.baitmap` | Dovetail 5 kb / 10 kb / 20 kb design files |
| Probe BED (optional) | `*probes*.bed` | enables target-enrichment QC |
| CHiCAGO settings file (optional) | pass `--chicago_settings_file=` | overrides CHiCAGO defaults, e.g. a lower `minNPerBait` for shallow pilot libraries |

Nothing is hardcoded: every input is located by pattern at run time, and any
input can be overridden with the matching flag.

## Parameters

`--sample_name --r1_fastq --r2_fastq --reference_fasta --design_dir --probes_bed
--species --design_resolution --chicago_cutoff --min_mapq --threads
--export_formats --filter_interactions --chicago_settings_file` (see `./run --help`).

## Outputs (`/results`)

- `<sample>.PT.bam` (+ `.bai`), `<sample>.mapped.pairs.gz`
- `<sample>.pairtools.stats.txt`, `<sample>.qc_report.txt`, `qc/`
- `<sample>.chicago.bam`, `<sample>.chinput`
- `<sample>_chicago_calls/data/` — `*_washU_text.txt`,
  `*_washU_track.txt.gz` + `.tbi`, `*.ibed`, `*_seqMonk.txt`, `*.Rds`
- `<sample>_chicago_calls/diag_plots/`, `<sample>_chicago_calls/examples/`
- `<sample>.filtered_interactions.ibed` — cis, 10 kb–2 Mb, ≥5 reads
- `<sample>.run_parameters.txt`

Load a WashU track in the browser as a **longrange** local track, supplying both
the `*_washU_track.txt.gz` and its `.tbi`; the `*_washU_text.txt` file loads as a
**long-range text** local text track.

## Depth expectations

CHiCAGO calls interactions at the bait level and by default discards baits with
fewer than 250 reads (`minNPerBait`), so a full library (the protocol recommends
at least 150M read pairs) is needed for meaningful calls at the 10 kb design.
Shallower pilot libraries run through unchanged up to the calling step; supply
`--chicago_settings_file` with a lower `minNPerBait` (one `name<TAB>value` line
per setting) if you want calls from a pilot run, and treat them as diagnostic.

## Validation

Built and exercised end to end on the published Dovetail NSC capture test
library (1M read pairs) against a single-chromosome (chr8) hg38 reference and the
human 10 kb design files: alignment and pair calling, QC report, on-target rate
and probe coverage, CHiCAGO-compatible BAM, chinput (14,871 bait/other-end
records), and CHiCAGO export of `washU_text`, `washU_track` (bgzipped + tabix
indexed and queryable), `interBed`/`.ibed`, `seqMonk` and `.Rds`.
