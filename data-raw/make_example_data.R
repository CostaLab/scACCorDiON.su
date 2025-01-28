##
dat <- data.table::fread(file.path("data-raw", "PDAC_2910-cci_matrix.csv"))

# clean data for r processing
pdac_cci_data <-
        dat |>
        tibble::as_tibble() |>
        dplyr::mutate(LR = V1) |>
        dplyr::select(-V1) |>
        dplyr::relocate(LR) |>
        (\(x){
                gsub("@", "_", colnames(x)) |>
                        make.names() |>
                        setNames(object = x)
        })()


# query_TCGA = TCGAbiolinks::GDCquery(
#         project = "TCGA-PAAD",
#         data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
#         data.type = "Gene Expression Quantification",
#         experimental.strategy = "RNA-Seq",
#         workflow.type = "STAR - Counts",
#         sample.type = c("Primary Tumor", "Solid Tissue Normal")
# )
#
# TCGAbiolinks::getResults(query_TCGA)
# TCGAbiolinks::GDCdownload(query = query_TCGA)
#
# tcga_data <- TCGAbiolinks::GDCprepare(query_TCGA)
tcga_data <- qs::qread(file.path("data-raw", "PAAD_tcga_data.qs"))

# normalize and filter
# TODO: check effect
tcga_postfilt <-
        tcga_data |>
        TCGAbiolinks::TCGAanalyze_Preprocessing() |>
        TCGAbiolinks::TCGAanalyze_Normalization(geneInfo = TCGAbiolinks::geneInfoHT, method = "gcContent") |>
        TCGAbiolinks::TCGAanalyze_Filtering(method = "quantile")

# add gene symbol info
paad_tcga_expression_data <-
        tcga_postfilt |>
        as.data.frame() |>
        tibble::rownames_to_column("ENSG") |>
        tibble::as_tibble() |>
        dplyr::inner_join(
                annotables::grch38 |>
                        dplyr::select(ensgene, symbol),
                by = c("ENSG" = "ensgene")
        ) |>
        dplyr::relocate(symbol) |>
        dplyr::select(-ENSG)


select_clinical_cols <- c(
        "patient", "vital_status", "days_to_death", "days_to_last_follow_up", "gender", "age_at_index",
        "ajcc_pathologic_stage", "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m"
)

paad_tcga_clinical_data <-
        tcga_data@colData |>
        as.data.frame() |>
        dplyr::filter(definition == "Primary solid Tumor") |>
        dplyr::select(dplyr::all_of(select_clinical_cols)) |>
        dplyr::mutate(deceased = vital_status == "Dead") |>
        dplyr::mutate(overall_survival = ifelse(deceased, days_to_death, days_to_last_follow_up)) |>
        # if OS is zero and status is alive, it prob means theres some issue with the clinical data
        dplyr::filter(!(overall_survival == 0 & !deceased)) |>
        dplyr::mutate(time = overall_survival, event = deceased, age = age_at_index) |>
        dplyr::mutate(stage = gsub("A$|B$", "", ajcc_pathologic_stage)) |>
        dplyr::mutate(stage = gsub(" ", "_", stage)) |>
        dplyr::mutate(staging_tumor = gsub("a$|b$", "", ajcc_pathologic_t)) |>
        dplyr::mutate(staging_nodes = gsub("a$|b$", "", ajcc_pathologic_n)) |>
        dplyr::mutate(staging_metastasis = gsub("a$|b$", "", ajcc_pathologic_m)) |>
        dplyr::select(-ajcc_pathologic_stage, -ajcc_pathologic_t, -ajcc_pathologic_n, -ajcc_pathologic_m) |>
        tibble::rownames_to_column("sample_id") |>
        tibble::as_tibble()


comparison_of_interest <- "Ductal.cell.type.2_Ductal.cell.type.1"




data.table::fwrite(pdac_cci_data, file = "data-raw/pdac_cci_data.csv")
data.table::fwrite(paad_tcga_clinical_data, file = "data-raw/paad_tcga_clinical_data.csv")
data.table::fwrite(paad_tcga_expression_data, file = "data-raw/paad_tcga_expression_data.csv")

usethis::use_data(pdac_cci_data, overwrite = TRUE, compress = "xz")
usethis::use_data(paad_tcga_clinical_data, overwrite = TRUE, compress = "xz")
usethis::use_data(paad_tcga_expression_data, overwrite = TRUE, compress = "xz")
usethis::use_data(comparison_of_interest, overwrite = TRUE, compress = "xz")
