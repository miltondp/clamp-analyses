rule quickstart_example:
    input:
        script="scripts/quickstart/run.py",
        builder="scripts/pseudobulk/build_pseudobulk.py",
        common="scripts/pseudobulk/common.py",
    output:
        report="output/00_quickstart/report.json",
        expression="output/00_quickstart/expression.csv",
        fractions="output/00_quickstart/fractions.csv",
        notebook="output/00_quickstart/quickstart.ipynb",
        raw="output/00_quickstart/synthetic.h5ad",
    log:
        "output/00_quickstart/run.log",
    threads: 1
    resources:
        mem_mb=2000,
        runtime=5,
    conda: "clamp-analyses"
    shell:
        "python {input.script:q} --output-dir output/00_quickstart > {log:q} 2>&1"

rule quickstart:
    input:
        rules.quickstart_example.output,
