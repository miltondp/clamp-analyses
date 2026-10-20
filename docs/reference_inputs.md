# Small reference inputs

The small pathway, cell/tissue and gene-length snapshots are included in Git.
GTEx annotations have download rules. Run from the repository root:

```bash
conda activate snakemake
snakemake --profile workflow/profiles/local reference_inputs
python scripts/inputs/reference_inputs.py --check
```

The first command restores missing bundled files or downloads missing GTEx
annotations; the second verifies every installed file against its SHA-256.
Downstream analyses also request missing references automatically. Existing
files are checked and never overwritten. A changed upstream response fails the
checksum check rather than silently changing the analysis.

## Sources and citations

The complete filename, version, source URL, DOI list and SHA-256 for each file
are in [the machine-readable manifest](../workflow/config/reference_inputs.json).
Verification dates below are **not historical download dates**.

| Reference | Included/downloaded from | Citation and verification |
| --- | --- | --- |
| MSigDB `c2.all`, `c2.cp`, `c2.cp.reactome`, `c5.go.cc`, `c8.all`, `h.all` | Included unchanged; [official 2026.1.Hs release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/) | All six files byte-identical to upstream on 2026-10-20. [GSEA, DOI:10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102), [MSigDB, DOI:10.1093/bioinformatics/btr260](https://doi.org/10.1093/bioinformatics/btr260), [hallmarks, DOI:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004). |
| CellMarker human workbook | Included unchanged; [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/), [documented download endpoint](https://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx) | [CellMarker 2.0, DOI:10.1093/nar/gkac947](https://doi.org/10.1093/nar/gkac947). The endpoint timed out on 2026-10-20; original acquisition date/URL was not recorded. Use the bundled workbook, sheet `human`, columns `cell_name` and `Symbol`. |
| Allen Brain Atlas 10x scRNA 2021 | Included unchanged; [Enrichr library download](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Allen_Brain_Atlas_10x_scRNA_2021) | Byte-identical on 2026-10-20; [Allen Brain Atlas](https://portal.brain-map.org/), [Enrichr, DOI:10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377). |
| Azimuth 2023 / ARCHS4 Cell-lines | Included snapshots; Enrichr endpoints in the manifest | [Azimuth](https://azimuth.hubmapconsortium.org/), [ARCHS4](https://maayanlab.cloud/archs4/), [Enrichr, DOI:10.1002/cpz1.90](https://doi.org/10.1002/cpz1.90). Preserve the pinned snapshots if the current libraries differ. |
| `GTEx_Tissues_pathMat.rds` | Included derived snapshot; [Enrichr GTEx_Tissues_V8_2023](https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GTEx_Tissues_V8_2023) | 8,076 genes × 511 terms; all memberships match the upstream GMT after adding the `GTEx_Tissues_` prefix. Verified 2026-10-20. Historical serialization script is unavailable. [GTEx, DOI:10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776), [Enrichr, DOI:10.1093/nar/gkw377](https://doi.org/10.1093/nar/gkw377). |
| `gene_lengths.rds` | Included existing analysis snapshot: 44,737 named gene lengths, in base pairs | The [historical notebook](https://github.com/pivlab/clamp-analyses/blob/6935e139dc9ef6c4d17ea4ff17aa3d7f603051fe/nbs/01_model_building/04_archs4/01_archs4_preprocess.ipynb) and its config query Ensembl 107 BioMart, then take the maximum transcript length per `(hgnc_symbol, ensembl_gene_id)` pair. [Ensembl release 107](https://ftp.ensembl.org/pub/release-107/) and [Ensembl 2022, DOI:10.1093/nar/gkab1049](https://doi.org/10.1093/nar/gkab1049) describe that recipe. The acquisition date/serialization environment is unknown; the bundled vector has not been replaced by a new query. |
| GTEx sample attributes / subject phenotypes | Downloaded from [official GTEx v8 metadata](https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=annotations/v8/metadata-files/) using manifest URLs; sample attributes also remain bundled for compatibility | Both files byte-identical on 2026-10-20. [GTEx, DOI:10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776). Subject phenotypes are a companion reference; current biology rules do not require them. |
| GTEx seven-cell-type xCell scores | [Official GTEx v8 download](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/cell-type-iqtl/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt.gz), decompressed automatically | Installed bytes match the existing file exactly. [GTEx, DOI:10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776), [xCell, DOI:10.1186/s13059-017-1349-1](https://doi.org/10.1186/s13059-017-1349-1). |

MSigDB files retain their release and original description columns. Third-party
terms and notices, including KEGG/BioCarta conditions, are retained in
[data/pathways/NOTICE.md](../data/pathways/NOTICE.md). Cite underlying gene-set
sources as directed by MSigDB. No scientific inputs were upgraded in this change.

## Restoring or checking selected files

```bash
# Offline restoration of bundled references from the current Git commit
python scripts/inputs/reference_inputs.py --offline --only archs4_gene_lengths cellmarker

# Download only the GTEx companion annotations and xCell scores
python scripts/inputs/reference_inputs.py --only gtex_sample_attributes gtex_subject_phenotypes gtex_xcell
```

To test an upstream URL against the pinned checksum without touching analysis
inputs, use `--upstream --output-root /tmp/reference-download-check --only <key>`.
Use a new empty directory for each test; existing files are always checked and
reused. Bundled derived RDS inputs have no original download URL and are restored
from Git. Git snapshots, hashes and provenance records make those historical
inputs portable without claiming an unrecorded reconstruction method.

## Historical ARCHS4 gene-length recipe

To reproduce the recovered calculation into a separate file (requires the raw
ARCHS4 v2.5 HDF5 and access to Ensembl's version-107 BioMart archive):

```bash
conda run -n clamp-analyses Rscript scripts/inputs/build_archs4_gene_lengths.R \
  --h5 data/archs4/human_gene_v2.5.h5 --out /tmp/archs4-gene-lengths-rebuilt.rds
```

This is optional; fresh clones already contain the exact reference vector used
by the analyses. The recipe does not alter the current vector. Re-querying the
archive was not needed for the input automation checks; compare names, values
and order before adopting a newly serialized vector.
