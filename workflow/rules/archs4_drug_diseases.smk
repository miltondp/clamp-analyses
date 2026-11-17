import shlex

A4_DD_CFG = A4_CFG["drug_diseases"]
A4_DD_DATA = A4_DD_CFG["data_dir"]
A4_DD_SPREDIXCAN = A4_DD_CFG["spredixcan_dir"]
A4_DD_COMPUTE = A4_DD_CFG["compute_root"]
A4_DD_REPORT = A4_DD_CFG["report_root"]
A4_DD_NLVS = ",".join(map(str, A4_DD_CFG["n_top_lvs"]))
A4_DD_NGENES = ",".join(map(str, A4_DD_CFG["n_top_genes"]))

A4_DD_METHODS = A4_DD_CFG["methods"]
A4_DD_COMPARISON = list(A4_DD_CFG["comparison_methods"])
A4_DD_METHOD_PATTERN = "|".join(A4_DD_METHODS)
A4_DD_COMPARISON_ROOT = f"{A4_DD_REPORT}/method_comparison"
A4_DD_COMPARISON_TABLES = ["per_tissue_metrics", "max_aggregate_reference", "max_aggregate_ci",
                           "ordering_stability", "paired_tissue_tests", "prediction_scores_max"]
A4_DD_L2 = "$(python -c \"import json; print(repr(json.load(open('{input.manifest}'))['L2']))\")"


def dd_compendium(method):
    return A4_DD_METHODS[method]["compendium"]


def dd_model_dir(method):
    m = A4_DD_METHODS[method]
    root = CON_ROOT if m["model"] == "full" else CON_BASE_ROOT
    return f"{root}/{m['compendium']}"


def dd_model_inputs(method):
    return {"z": f"{dd_model_dir(method)}/Z.csv", "manifest": f"{dd_model_dir(method)}/manifest.json"}


def dd_prediction_dir(method):
    return f"{A4_DD_COMPUTE}/predictions/module_based_{method}"


def dd_scores(path):
    return f"{path}/lincs/predictions/dotprod_neg"


def dd_comparison_args(methods):
    specs = [f"{A4_DD_METHODS[m]['label']}={dd_scores(dd_prediction_dir(m))}=lv" for m in methods]
    specs.append(f"Gene-based={dd_scores(f'{A4_DD_COMPUTE}/predictions/gene_based')}=gene")
    return " ".join(f"--method {shlex.quote(spec)}" for spec in specs)


# ============================================================
# Step 1: project traits and drug signatures
# ============================================================

rule project_spredixcan_drug_diseases:
    input:
        unpack(lambda wc: dd_model_inputs(wc.method)),
        spredixcan=A4_DD_SPREDIXCAN,
        script="scripts/archs4/drug_diseases/project_spredixcan.py"
    output:
        directory(f"{A4_DD_COMPUTE}/projections/{{method}}/spredixcan")
    wildcard_constraints:
        method=A4_DD_METHOD_PATTERN
    resources:
        mem_mb=32000, runtime=720
    params:
        compendium=lambda wc: dd_compendium(wc.method)
    conda: "clamp-analyses"
    shell:
        "python {input.script} --z {input.z} --l2 " + A4_DD_L2 + " --spredixcan-dir {input.spredixcan} "
        "--compendium {params.compendium} --output-dir {output}"

rule project_lincs_drug_diseases:
    input:
        unpack(lambda wc: dd_model_inputs(wc.method)),
        lincs=f"{A4_DD_DATA}/lincs-data.pkl",
        script="scripts/archs4/drug_diseases/project_lincs.py"
    output:
        f"{A4_DD_COMPUTE}/lincs/{{method}}.pkl"
    wildcard_constraints:
        method=A4_DD_METHOD_PATTERN
    resources:
        mem_mb=16000, runtime=360
    params:
        compendium=lambda wc: dd_compendium(wc.method)
    conda: "clamp-analyses"
    shell:
        "python {input.script} --z {input.z} --l2 " + A4_DD_L2 + " --lincs-file {input.lincs} "
        "--compendium {params.compendium} --output {output}"

rule fetch_drugbank_drug_diseases:
    output:
        A4_DD_CFG["drugbank_tsv"]
    log:
        A4_DD_CFG["drugbank_tsv"] + ".log"
    shell:
        "mkdir -p $(dirname {output:q}) && "
        "curl --fail --location --retry 3 '{A4_DD_CFG[drugbank_url]}' --output {output}.download 2> {log} && "
        "echo '{A4_DD_CFG[drugbank_sha256]}  {output}.download' | sha256sum --check --status && "
        "mv {output}.download {output}"


# ============================================================
# Step 2: predict drug-disease associations
# ============================================================

rule predict_module_drug_diseases:
    input:
        lincs=rules.project_lincs_drug_diseases.output,
        spredixcan=rules.project_spredixcan_drug_diseases.output,
        gold=f"{A4_DD_DATA}/gold_standard.pkl", ukb=f"{A4_DD_DATA}/phenomexcan_traits_fullcode_to_efo.tsv",
        efo=f"{A4_DD_DATA}/term_id_xrefs.tsv.gz", do=f"{A4_DD_DATA}/xrefs-prop-slim.tsv",
        script="scripts/archs4/drug_diseases/predict_module_based.py"
    output:
        directory(f"{A4_DD_COMPUTE}/predictions/module_based_{{method}}")
    wildcard_constraints:
        method=A4_DD_METHOD_PATTERN
    resources:
        mem_mb=16000, runtime=720
    params:
        compendium=lambda wc: dd_compendium(wc.method)
    conda: "clamp-analyses"
    shell:
        "python {input.script} --compendium {params.compendium} --lincs-projection {input.lincs} --spredixcan-proj-dir {input.spredixcan}/proj "
        "--gold-standard {input.gold} --ukb-efo {input.ukb} --efo-xrefs {input.efo} --do-xrefs {input.do} --n-top-lvs {A4_DD_NLVS} --output-dir {output}"

rule predict_gene_drug_diseases:
    input:
        spredixcan_raw=f"{A4_DD_COMPUTE}/projections/archs4_full/spredixcan",
        lincs=f"{A4_DD_DATA}/lincs-data.pkl",
        gold=f"{A4_DD_DATA}/gold_standard.pkl", ukb=f"{A4_DD_DATA}/phenomexcan_traits_fullcode_to_efo.tsv",
        efo=f"{A4_DD_DATA}/term_id_xrefs.tsv.gz", do=f"{A4_DD_DATA}/xrefs-prop-slim.tsv",
        script="scripts/archs4/drug_diseases/predict_single_gene.py"
    output:
        directory(f"{A4_DD_COMPUTE}/predictions/gene_based")
    resources:
        mem_mb=16000, runtime=360
    conda: "clamp-analyses"
    shell:
        "python {input.script} --spredixcan-raw-dir {input.spredixcan_raw}/raw --lincs-file {input.lincs} "
        "--gold-standard {input.gold} --ukb-efo {input.ukb} --efo-xrefs {input.efo} --do-xrefs {input.do} --n-top-genes {A4_DD_NGENES} --output-dir {output}"


# ============================================================
# Step 3: compare prediction methods
# ============================================================

rule aggregate_drug_diseases:
    input:
        modules=[dd_prediction_dir(m) for m in A4_DD_COMPARISON],
        gene=f"{A4_DD_COMPUTE}/predictions/gene_based",
        gold=f"{A4_DD_DATA}/gold_standard.pkl",
        script="scripts/archs4/drug_diseases/aggregate_predictions.py"
    output:
        directory(f"{A4_DD_REPORT}/aggregate/max")
    params:
        methods=dd_comparison_args(A4_DD_COMPARISON),
    resources:
        mem_mb=32000, runtime=720
    conda: "clamp-analyses"
    shell:
        "python {input.script} {params.methods} --gold-standard {input.gold} "
        "--n-top-lvs {A4_DD_NLVS} --n-top-genes {A4_DD_NGENES} --output-dir {output}"

rule drug_disease_method_comparison:
    input:
        modules=[dd_prediction_dir(m) for m in A4_DD_COMPARISON],
        gene=f"{A4_DD_COMPUTE}/predictions/gene_based",
        gold=f"{A4_DD_DATA}/gold_standard.pkl",
        script="scripts/archs4/drug_diseases/prepare_method_comparison.py"
    output:
        stats=f"{A4_DD_COMPARISON_ROOT}/data/figure_statistics.txt",
        tables=expand(f"{A4_DD_COMPARISON_ROOT}/data/{{name}}.csv", name=A4_DD_COMPARISON_TABLES),
    params:
        methods=dd_comparison_args(A4_DD_COMPARISON),
        out=lambda wc, output: os.path.dirname(output.stats),
    resources:
        mem_mb=16000, runtime=360
    conda: "clamp-analyses"
    shell:
        "python {input.script} {params.methods} --gold-standard {input.gold} "
        "--n-top-lvs {A4_DD_NLVS} --n-top-genes {A4_DD_NGENES} --n-boot {A4_DD_CFG[bootstrap][n_boot]} "
        "--seed {A4_DD_CFG[bootstrap][seed]} --output-dir {params.out}"

rule drug_disease_topk_comparison:
    input:
        modules=[dd_prediction_dir(m) for m in A4_DD_COMPARISON],
        gene=f"{A4_DD_COMPUTE}/predictions/gene_based",
        gold=f"{A4_DD_DATA}/gold_standard.pkl",
        script="scripts/archs4/drug_diseases/prepare_method_comparison.py"
    output:
        stats=f"{A4_DD_COMPARISON_ROOT}/top{{k}}/figure_statistics.txt",
        tables=expand(f"{A4_DD_COMPARISON_ROOT}/top{{{{k}}}}/{{name}}.csv", name=A4_DD_COMPARISON_TABLES),
    params:
        methods=dd_comparison_args(A4_DD_COMPARISON),
        out=lambda wc, output: os.path.dirname(output.stats),
    wildcard_constraints:
        k="|".join(map(str, A4_DD_CFG["n_top_lvs"]))
    resources:
        mem_mb=16000, runtime=360
    conda: "clamp-analyses"
    shell:
        "python {input.script} {params.methods} --gold-standard {input.gold} "
        "--n-top-lvs {wildcards.k} --n-top-genes {A4_DD_NGENES} --n-boot {A4_DD_CFG[bootstrap][n_boot]} "
        "--seed {A4_DD_CFG[bootstrap][seed]} --output-dir {params.out}"

A4_DD_TOPK = expand(f"{A4_DD_COMPARISON_ROOT}/top{{k}}/max_aggregate_reference.csv", k=A4_DD_CFG["n_top_lvs"])


# ============================================================
# Step 4: render reports
# ============================================================

rule drug_disease_method_comparison_report:
    input:
        stats=rules.drug_disease_method_comparison.output.stats,
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/02_method_comparison.ipynb"
    output:
        complete=touch(f"{A4_DD_COMPARISON_ROOT}/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/02_method_comparison.executed.ipynb"
    resources:
        mem_mb=4000, runtime=120
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/02_method_comparison.ipynb"

rule drug_diseases_prediction_performance_report:
    input:
        aggregate=f"{A4_DD_REPORT}/aggregate/max",
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.ipynb"
    output:
        complete=touch(f"{A4_DD_REPORT}/00_prediction_performance/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.executed.ipynb"
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/00_prediction_performance.ipynb"

rule drug_diseases_prediction_plots_report:
    input:
        aggregate=f"{A4_DD_REPORT}/aggregate/max",
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.ipynb"
    output:
        complete=touch(f"{A4_DD_REPORT}/01_prediction_performance_plots/notebook.complete")
    log:
        notebook=f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.executed.ipynb"
    conda: "clamp-analyses"
    notebook:
        f"{A4_BIO_NB}/02_drug_diseases_associations/01_prediction_performance_plots.ipynb"

rule archs4_drug_diseases:
    input:
        rules.aggregate_drug_diseases.output,
        rules.drug_disease_method_comparison.output,
        A4_DD_TOPK,
        rules.drug_disease_method_comparison_report.output,
        rules.drug_diseases_prediction_performance_report.output,
        rules.drug_diseases_prediction_plots_report.output,
