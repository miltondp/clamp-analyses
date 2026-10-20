# PhenoPLIER, drug–disease and CRISPR inputs

**Many inputs are public, but the CLI's `full_data` download is not a complete
replacement for the processed reference bundle used by these analyses.** The
pinned sources/checksums are in
[phenoplier_inputs.json](../workflow/config/phenoplier_inputs.json).

## What can be downloaded?

| Input | Source and current status |
| --- | --- |
| LINCS drug signatures, drug–disease gold standard, trait/EFO and ontology mappings | Download automatically from [PhenoPLIER commit 5894a6b](https://github.com/greenelab/phenoplier/tree/5894a6bd9b091f5cf1f4806f1f3e20fb6caa91cc). LINCS/gold-standard downloads use Git LFS media URLs, not pointer files. All five match the current analysis inputs. |
| S-MultiXcan z-scores | Available from the same repository's pinned Git LFS object (~738 MB). Its SHA-256 matches the current analysis input. The historical CRISPR-specific path is linked to this shared file. |
| Gene-name/Ensembl mapping | Download automatically from [Zenodo 10944491](https://doi.org/10.5281/zenodo.10944491), the source used by phenoplier-cli v0.5.2. |
| Processed phenotype metadata | The exact current `phenotypes_info.tsv` is bundled. It is distinct from the raw phenotype metadata that the CLI downloads. |
| CRISPR `lipid_DEG.csv` | The original small file is bundled. It matches the checksum in PhenoPLIER's download script; its historical Box URL returns HTTP 404. |
| Processed S-MultiXcan P values and fastENLOC colocalization probabilities | Current files match PhenoPLIER's recorded checksums, but the recorded Box URLs return HTTP 404. Supply matching files from an existing PhenoPLIER workspace/extracted bundle. |
| 49 processed S-PrediXcan tissue pickles | Supply the existing processed `gene_assoc/spredixcan/pkl` directory. Historical raw/HDF5 archives are different inputs and are not silently substituted for these pickles. |
| GLS reference correlations, GWAS/TWAS and prediction-model resources | Use the matching PhenoPLIER workspace/reference bundle. Some underlying resources have CLI/Zenodo downloads; the CLI's registered `full_data` resources do not guarantee every processed GLS/reference file required by this workflow. |

PhenoPLIER paper and CRISPR screen:
[DOI:10.1038/s41467-023-41057-4](https://doi.org/10.1038/s41467-023-41057-4).
The original project's archived code/data are described at
[DOI:10.5281/zenodo.8071382](https://doi.org/10.5281/zenodo.8071382).
The archive was not assumed to be a drop-in replacement for the CLI workspace.

## Stage the inputs

Set the root of a workspace or an extracted bundle containing `data/phenomexcan`:

```bash
conda activate snakemake
export PHENOPLIER_DATA_ROOT=/shared/phenoplier_full_data
snakemake --profile workflow/profiles/local --dry-run phenoplier_inputs
snakemake --profile workflow/profiles/local phenoplier_inputs
```

This downloads live public files and links matching processed files from that
root. It does not run GLS or model fitting. Existing files are checked and never
replaced. Large inputs remain outside Git. A missing/unavailable processed input
fails with staging instructions rather than being replaced by another dataset.
For cluster jobs, keep the source root on shared storage. A configuration override
`phenoplier_input_root` can supply the same path without environment variables.

To stage only a subset:

```bash
snakemake --profile workflow/profiles/local drug_disease_inputs
snakemake --profile workflow/profiles/local crispr_inputs
```

The drug–disease target also requests the existing pinned DrugBank download and
the 49 processed tissue pickles. The CRISPR target includes its PhenomeXcan
metadata/association inputs. Large files may take substantial time to download.

To verify the manifest's file inputs after staging:

```bash
python scripts/inputs/reference_inputs.py --manifest workflow/config/phenoplier_inputs.json --check
```

This verifies file SHA-256 values; the S-PrediXcan staging helper separately
checks all 49 tissue-pickle filenames and byte sizes against the recorded
inventory; it does not claim content-hash verification for these 21 GB of files. It
links individual files into a new owned directory, so scheduler metadata never
needs to be written inside the original reference directory.

## Initialize GLS separately

Keep phenoplier-cli **v0.5.2**, as required by the analysis configuration.
See [installation/workspace setup](../scripts/phenoplier/README.md#3-phenoplier-cli-environment-and-workspace).
For a new workspace, the CLI supports:

```bash
conda activate phenoplier-cli-neo
export PHENOPLIER_HOME=/shared/my-phenoplier-workspace
phenoplier workspace init
phenoplier workspace download full_data --list
```

`--list` only inspects available downloads. Removing it downloads the CLI's
registered resources, but does not establish that the full processed bundle is
complete. For the existing bundle, follow the `workspace link` instructions in
the PhenoPLIER guide and use a bundle copy you own. `PHENOPLIER_HOME` selects the
GLS workspace; `PHENOPLIER_DATA_ROOT` above selects input files for this repository.

## Timing results are recorded evidence

The author-generated timing CSVs are included unchanged under
`references/runtime/{pseudobulk,gtex}/`, with source paths and hashes in
`references/runtime/recorded_timings.json`. The default `benchmark` and `all`
workflows restore those tables; they do **not** rerun timing fits or memory
measurements. Supplementary Figure 4 consumes the restored results.

Model benchmarking elsewhere in phase 1 still fits the scientific models when
needed; the saved timing results do not replace those analyses. Existing explicit
timing-fit rules remain available for deliberate new measurements. Recomputing
both timing reports requires an explicit `--config timing_mode=recompute` override
and the full scientific input prerequisites; reviewers do not need this.
