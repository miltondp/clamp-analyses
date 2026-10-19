# Final Models GLS

The active workflow is `workflow/rules/phenoplier.smk`. Run from the repository
root after setting up phenoplier-cli v0.5.2 and its reference workspace:

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 12 phenoplier_canonical
```

Coverage and saturation use CLAMPbase only. Final CLAMPfull models use the
canonical prior. The historical BP-only launch scripts have been removed;
`scripts/phenoplier/RUN_SUMMARY.md` retains the historical run record.
