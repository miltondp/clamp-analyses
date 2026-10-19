# Define cell-type markers and test module enrichment.
CELLMARKER_PATH <- snakemake@input[["cell_marker_file"]]
ALLEN_PATH      <- snakemake@input[["allen_brain_gmt_file"]]
AZIMUTH_PATH    <- snakemake@input[["azimuth_file"]]

load_cellmarker_sets <- function(path) {
    df <- readxl::read_excel(path, sheet = "human")
    df <- df[df$species == "Human" & !is.na(df$Symbol) & nzchar(df$Symbol), ]
    lapply(split(df$Symbol, df$cell_name), unique)
}

read_gmt_sets <- function(path) {
    sp <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
    sp <- sp[lengths(sp) >= 3]
    setNames(lapply(sp, function(x) unique(x[-c(1, 2)])),
             vapply(sp, `[`, character(1), 1))
}

cellmarker_t2g <- load_cellmarker_sets(CELLMARKER_PATH)
allen_raw <- read_gmt_sets(ALLEN_PATH)
allen_t2g <- allen_raw[grepl("^Human ", names(allen_raw)) &
                        grepl(" up$", names(allen_raw))]
BRAIN_DS <- unlist(snakemake@config$references$brain_datasets)

normalize_marker_text <- function(x) {
    stringr::str_squish(gsub("[^a-z0-9]+", " ", tolower(x)))
}

CELL_TYPE_ALIASES <- list(
    Ast = c("astrocyte"), Exc = c("excitatory neuron", "glutamatergic neuron"),
    Inh = c("inhibitory neuron", "gabaergic neuron"),
    Mic = c("microglia", "microglial"), Oli = c("oligodendrocyte"),
    Opc = c("oligodendrocyte precursor", "oligodendrocyte progenitor", "opc"),
    B_cell = c("b cell", "b lymphocyte"), NK = c("natural killer", "nk cell"),
    CD4_T = c("cd4 t cell", "cd4 positive t cell"),
    CD8_T = c("cd8 t cell", "cd8 positive t cell"),
    CD14_Mono = c("cd14 monocyte", "classical monocyte", "monocyte"),
    CD16_Mono = c("cd16 monocyte", "non classical monocyte", "monocyte"),
    DC = c("dendritic cell"), gd_T = c("gamma delta t cell", "gd t cell"),
    Plasma_B = c("plasma cell", "plasmablast"),
    Myeloid = c("myeloid", "monocyte", "macrophage"),
    T_cell = c("t cell", "t lymphocyte"),
    Macrophage = c("macrophage", "myeloid"),
    Endothelial = c("endothelial cell"), Pericyte = c("pericyte"),
    VentricularCM = c("ventricular cardiomyocyte", "cardiomyocyte"),
    AtrialCM = c("atrial cardiomyocyte", "cardiomyocyte"),
    Fibroblast = c("fibroblast"), Adipocyte = c("adipocyte"),
    Endocardial = c("endocardial cell", "endothelial cell"),
    Epicardial = c("epicardial cell", "mesothelial cell"),
    LymphaticEndothelial = c("lymphatic endothelial"),
    Lymphocyte = c("lymphocyte", "t cell", "b cell", "natural killer"),
    Mast = c("mast cell"), Neuronal = c("neuron", "neuronal"),
    VSMC = c("vascular smooth muscle", "smooth muscle cell"),
    `Blood vessels` = c("endothelial cell", "vascular endothelial"),
    Lymphoid = c("lymphocyte", "t cell", "b cell", "natural killer"),
    `Alveolar epithelium` = c("alveolar epithelial", "type i pneumocyte", "type ii pneumocyte"),
    `Airway epithelium` = c("ciliated cell", "club cell", "goblet cell", "basal cell"),
    `Fibroblast lineage` = c("fibroblast", "myofibroblast"),
    `Lymphatic EC` = c("lymphatic endothelial"),
    `Submucosal Gland` = c("submucosal gland", "serous cell", "mucous cell")
)

cellmarker_terms_for <- function(cell_type) {
    aliases <- CELL_TYPE_ALIASES[[cell_type]]
    if (is.null(aliases)) aliases <- ct_label(cell_type)
    aliases <- normalize_marker_text(aliases)
    terms_norm <- normalize_marker_text(names(cellmarker_t2g))
    names(cellmarker_t2g)[vapply(terms_norm, function(term)
        any(vapply(aliases, function(alias)
            grepl(paste0(" ", alias, " "), paste0(" ", term, " "),
                  fixed = TRUE), logical(1))), logical(1))]
}

ALLEN_ALIAS_MAP <- list(
    Ast = "Astro", Exc = "Exc", Inh = "Inh", Mic = "Micro",
    Oli = "Oligo", Opc = "OPC"
)

allen_terms_for <- function(dataset, cell_type) {
    if (!dataset %in% BRAIN_DS) return(character(0))
    prefixes <- ALLEN_ALIAS_MAP[[cell_type]]
    if (is.null(prefixes)) return(character(0))
    pattern <- paste0("^Human (", paste(prefixes, collapse = "|"), ") ")
    names(allen_t2g)[grepl(pattern, names(allen_t2g))]
}

read_azimuth_sets <- function(path) {
    sp <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
    sp <- sp[lengths(sp) >= 3]
    setNames(lapply(sp, function(x) unique(x[-c(1, 2)])),
             vapply(sp, `[`, character(1), 1))
}
azimuth_t2g <- read_azimuth_sets(AZIMUTH_PATH)

AZIMUTH_TERM_MAP <- list(
    Brain_Mathys2023 = list(
        Oli = "Motor Cortex-subclass-oligodendrocyte",
        Exc = "Motor Cortex-class-Glutamatergic Neuron",
        Ast = "Motor Cortex-subclass-Astrocyte",
        Inh = "Motor Cortex-class-GABAergic Neuron",
        Mic = "Motor Cortex-subclass-microglia / Perivascular Macrophage",
        Opc = "Motor Cortex-subclass-oligodendrocyte Precursor Cell"
    ),
    Brain_Xiong2023 = list(
        Oli = "Motor Cortex-subclass-oligodendrocyte",
        Exc = "Motor Cortex-class-Glutamatergic Neuron",
        Ast = "Motor Cortex-subclass-Astrocyte",
        Inh = "Motor Cortex-class-GABAergic Neuron",
        Mic = "Motor Cortex-subclass-microglia / Perivascular Macrophage",
        Opc = "Motor Cortex-subclass-oligodendrocyte Precursor Cell"
    ),
    Heart_Datar2026 = list(
        Macrophage = "Heart-L2-Macrophage",
        Endothelial = "Heart-L2-Endothelial",
        Pericyte = "Heart-L2-Pericyte",
        VentricularCM = "Heart-L2-Ventricular Cardiomyocyte",
        Fibroblast = "Heart-L2-Fibroblast",
        Adipocyte = "Heart-L2-Adipocyte",
        AtrialCM = "Heart-L2-Atrial Cardiomyocyte",
        Endocardial = "Heart-L2-Endocardial",
        Epicardial = "Heart-L2-Mesothelial",
        LymphaticEndothelial = "Heart-L2-Lymphatic Endothelial",
        Lymphocyte = c("Heart-L2-B", "Heart-L2-ILC", "Heart-L2-NK", "Heart-L2-T"),
        Mast = "Heart-L2-Mast",
        Neuronal = "Heart-L2-Neuronal",
        VSMC = "Heart-L2-Smooth Muscle"
    ),
    PBMC_1k1k = list(
        B_cell = "PBMC-L1-B Cell",
        NK = "PBMC-L1-natural Killer Cell",
        CD4_T = "PBMC-L1-CD4+ T Cell",
        CD8_T = "PBMC-L1-CD8+ T Cell",
        CD14_Mono = "PBMC-L2-CD14+ Monocyte",
        CD16_Mono = "PBMC-L2-CD16+ Monocyte",
        DC = "PBMC-L1-dendritic Cell",
        gd_T = "PBMC-L2-gamma-delta T",
        Plasma_B = "PBMC-L2-Plasmablast"
    ),
    PBMC_Perez2022 = list(
        Myeloid = c("PBMC-L1-Monocyte", "PBMC-L1-dendritic Cell"),
        T_cell = c("PBMC-L1-CD4+ T Cell", "PBMC-L1-CD8+ T Cell",
                   "PBMC-L1-other T Cell"),
        B_cell = "PBMC-L1-B Cell",
        NK = "PBMC-L1-natural Killer Cell"
    ),
    Lung_Sikkema2023 = list(
        "Blood vessels" = "Lung V2 (HLCA)-ann Level 2-Blood Vessels",
        Lymphoid = "Lung V2 (HLCA)-ann Level 2-Lymphoid",
        "Alveolar epithelium" = "Lung V2 (HLCA)-ann Level 2-Alveolar Epithelium",
        "Airway epithelium" = "Lung V2 (HLCA)-ann Level 2-Airway Epithelium",
        Myeloid = "Lung V2 (HLCA)-ann Level 2-Myeloid",
        "Fibroblast lineage" = "Lung V2 (HLCA)-ann Level 2-Fibroblast Lineage",
        "Lymphatic EC" = "Lung V2 (HLCA)-ann Level 2-Lymphatic EC",
        "Submucosal Gland" = "Lung V2 (HLCA)-ann Level 2-Submucosal Gland"
    )
)

mapped_terms <- unique(unlist(AZIMUTH_TERM_MAP, recursive = TRUE, use.names = FALSE))
missing_terms <- setdiff(mapped_terms, names(azimuth_t2g))
if (length(missing_terms))
    stop("Azimuth mapping contains missing terms: ", paste(missing_terms, collapse = "; "))

azimuth_terms_for <- function(dataset, cell_type) {
    ds_map <- AZIMUTH_TERM_MAP[[dataset]]
    if (is.null(ds_map) || is.null(ds_map[[cell_type]]))
        stop("No Azimuth marker mapping for ", dataset, " / ", cell_type)
    ds_map[[cell_type]]
}

azimuth_dataset_terms <- function(dataset, cell_types = names(AZIMUTH_TERM_MAP[[dataset]])) {
    unique(unlist(AZIMUTH_TERM_MAP[[dataset]][cell_types], use.names = FALSE))
}

prefix_marker_sets <- function(sets, source) {
    stats::setNames(sets, paste0(source, "::", names(sets)))
}

prefix_marker_terms <- function(terms, source) {
    if (length(terms) == 0) return(character(0))
    paste0(source, "::", terms)
}

reference_sets_for <- function(dataset, cell_types = names(AZIMUTH_TERM_MAP[[dataset]])) {
    az_terms <- azimuth_dataset_terms(dataset, cell_types)
    out <- c(prefix_marker_sets(cellmarker_t2g, "CellMarker"),
             prefix_marker_sets(azimuth_t2g[az_terms], "Azimuth"))
    if (dataset %in% BRAIN_DS)
        out <- c(out, prefix_marker_sets(allen_t2g, "Allen"))
    out
}

target_reference_terms <- function(dataset, cell_type) {
    unique(c(prefix_marker_terms(cellmarker_terms_for(cell_type), "CellMarker"),
             prefix_marker_terms(azimuth_terms_for(dataset, cell_type), "Azimuth"),
             prefix_marker_terms(allen_terms_for(dataset, cell_type), "Allen")))
}

build_marker_sets <- function(dataset, cell_types, gene_universe, min_size = 5,
                              promiscuity_max_frac = 0.5) {
    references <- reference_sets_for(dataset, cell_types)
    sets <- lapply(cell_types, function(ct) {
        terms <- intersect(target_reference_terms(dataset, ct), names(references))
        intersect(unique(unlist(references[terms], use.names = FALSE)), gene_universe)
    })
    names(sets) <- cell_types
    too_small <- names(sets)[lengths(sets) < min_size]
    if (length(too_small))
        stop("Combined marker sets below min_size for ", dataset, ": ",
             paste(too_small, collapse = ", "))

    if (length(sets) > 1) {
        gene_counts <- table(unlist(sets))
        promiscuous <- names(gene_counts)[gene_counts > promiscuity_max_frac * length(sets)]
        sets <- lapply(sets, function(g) setdiff(g, promiscuous))
        too_small <- names(sets)[lengths(sets) < min_size]
        if (length(too_small))
            stop("Combined marker sets below min_size after promiscuity filter for ",
                 dataset, ": ", paste(too_small, collapse = ", "))
    }
    sets
}


get_top_genes <- function(Z, lv, top_pct, positive_only = FALSE) {
    if (!lv %in% colnames(Z)) return(character(0))
    loadings <- Z[[lv]]; names(loadings) <- rownames(Z)
    if (positive_only) return(names(loadings)[is.finite(loadings) & loadings > 0])
    n_genes  <- max(1L, ceiling(length(loadings) * top_pct))
    ord      <- order(loadings, decreasing = TRUE)
    names(loadings)[ord][seq_len(n_genes)]
}


build_term2gene <- function(marker_sets) {
    dplyr::bind_rows(lapply(names(marker_sets), function(ct) {
        data.frame(term = ct, gene = marker_sets[[ct]], stringsAsFactors = FALSE)
    }))
}

run_ora_enrichment <- function(selected_lvs_df, Z, term2gene, universe, top_pct,
                               positive_only = FALSE) {
    do.call(rbind, lapply(seq_len(nrow(selected_lvs_df)), function(i) {
        row       <- selected_lvs_df[i, ]
        top_genes <- get_top_genes(Z, row$LV, top_pct, positive_only = positive_only)
        res <- tryCatch(
            clusterProfiler::enricher(gene = top_genes, universe = universe,
                                      TERM2GENE = term2gene, pvalueCutoff = 1,
                                      qvalueCutoff = 1, minGSSize = 1, maxGSSize = Inf),
            error = function(e) NULL)
        if (is.null(res) || nrow(res@result) == 0) return(NULL)
        rdf <- res@result
        data.frame(lv_cell_type = row$cell_type, LV = row$LV,
                  marker_cell_type = rdf$ID, Count = rdf$Count,
                  pvalue = rdf$pvalue, FDR = rdf$p.adjust, stringsAsFactors = FALSE)
    }))
}

