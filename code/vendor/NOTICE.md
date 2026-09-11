# Third-party scripts bundled here

- `chicagoTools/` — from the Dovetail Genomics fork of CHiCAGO
  (`dovetail-genomics/chicago`, which the Dovetail capture protocol recommends
  because it carries minor bug fixes over the upstream chicagoTeam release).
  Used here: `bam2chicago.sh` (chinput generation) and `runChicago.R`
  (interaction calling and export of interBed / washU_text / washU_track /
  seqMonk). Artistic License 2.0, as in the CHiCAGO project.
- `get_qc.py` — from `dovetail-genomics/capture`; summarises the pairtools
  dedup stats file into the proximity-ligation QC report.

The CHiCAGO R package itself is installed from bioconda
(`bioconductor-chicago`, with `bioconductor-rsamtools` for the bgzipped/tabixed
WashU track) into a separate conda environment at `/opt/conda/envs/chicago-r`,
built by `environment/postInstall`. Override the interpreter with the
`CHICAGO_RSCRIPT` environment variable if needed.

Local patch to `chicagoTools/runChicago.R`: the two `packageVersion("argparser")`
comparisons against bare numeric literals (`< 0.3`, `< 0.4`) are quoted, because
R >= 4.4 refuses to coerce a double to a version. No other change.
