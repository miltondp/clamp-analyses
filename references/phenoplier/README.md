# Processed S-PrediXcan references

`spredixcan_download.json` pins the public archive, tissue filenames and exact
numerical fingerprints used by automatic conversion. The 49 converted matrices
match the author inputs in labels, ordering and values; pickle bytes can differ
between pandas versions. `spredixcan_inventory.json` is the historical filename
and size inventory. `input_equivalence_audit.json` records previous comparisons.

Remaining CLI workspace and downstream input issues are described in the
[PhenoPLIER handoff](../../scripts/phenoplier/README.md).
