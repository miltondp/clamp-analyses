# Fresh-checkout validation — 2026-10-19

Changes were prepared on `clamp-analyses-testing` in a temporary worktree.
The original checkout, scientific inputs, outputs, and installed environments
were not modified. No scientific fits, timing campaigns, GPU clustering jobs,
large data downloads, or cluster jobs were executed.

## Passed

- Installed `envs/snakemake.lock` into an isolated temporary Conda prefix:
  Python 3.11, Snakemake 8.30.0, SLURM executor plugin 2.8.0. Lock pins match YAML.
- Shell syntax checks, Python compilation, local documentation links, and
  `git diff --check`.
- Fresh checkout dry-runs fail clearly on genuine absent manual inputs.
- In a separate disposable copy, supplied **46 empty external-input placeholders**
  and no inherited scientific outputs/Snakemake metadata. Full DAGs construct:
  `benchmark` (415 jobs), `archs4_biology` (2394 jobs), `all` (2776 jobs).
  These checks validate dependency wiring, **not scientific input validity**.
  No missing workflow-produced outputs were replaced by placeholders.
- Local, portable SLURM, and Pico profile parsing and target dry-runs.
- Pico read-only checks: `defq`, 900000 MB RAM, 168 CPUs, unlimited partition
  time limit, one NVIDIA L40S (~46 GB VRAM). No job submission was tested.
- Installed the required CLAMP commit into a temporary R library and verified
  its `RemoteSha`. The pre-existing local installation lacks verified revision
  metadata; it was left unchanged.
- `scripts/check_setup.py` passed using the isolated runner, existing analysis/GPU
  environments, and temporary pinned CLAMP library: method revisions, Python/R
  imports, notebook tools, workflow parsing and Inkscape availability.
- Synthetic 80-cell, 60-gene, eight-donor dataset: cohort selection, pseudobulk
  construction, exact normalized aggregation, cell-type fractions, CLAMP
  preprocessing/rank estimation, and an R notebook loading pinned methods.
- Synthetic ORA aggregation regression: canonical and coverage directory layouts
  produce identical expected long, cross-dataset, and panel summaries.

## Fixes found by the checks

- Declared GTEx global-alignment figure inputs as outputs of their producing rule.
- Added explicit GTEx notebook model, tissue-reference and xCell dependencies.
- Restricted canonical aggregation wildcards and supported the canonical ORA
  directory layout in the aggregator.
- Added GPU scheduling to GTEx clustering and single-GPU scheduling pools.
- Profiles now select the Snakefile, activate environments, tolerate shared-disk
  latency, and serialize timing fits through a shared timing slot.
- Setup now locks Snakemake and installs pinned non-Conda comparison dependencies.
- Historical GTEx CoGAPS evidence is tracked rather than assumed under `output/`.

## Limits and remaining access requirements

A full reinstall of the large analysis/GPU environments and all source dependency
builds was not repeated. Their installed imports/revisions were checked, and the
CLAMP source installation and runner installation were tested independently.
Pip dependencies have version pins, not a complete transitive hash lock.

The PhenoPLIER source returned HTTP 404 anonymously. Its existing local environment
is a different release, so v0.5.2 installation/GLS was not validated. Setup now uses
Git authentication to obtain the installer from the requested release. Reviewers
need source access and the reference bundle from the authors.

Manual cohort/reference files with missing provenance remain author-supplied
requirements; the input checklist states this explicitly. Memory guidance records
configured allocations and Pico capacity, not a newly measured minimum or peak RSS.
Complete scientific reproduction and successful cluster submission remain untested.
