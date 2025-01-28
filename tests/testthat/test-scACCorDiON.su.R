test_that("Runs with default example.", {
        library(scACCorDiON.su)

        data(pdac_cci_data)
        data(paad_tcga_clinical_data)
        data(paad_tcga_expression_data)
        data(comparison_of_interest)

        results <-
                run_survival(
                        comparison_oi = comparison_of_interest,
                        lr_data = pdac_cci_data,
                        tcga_exp_data = paad_tcga_expression_data,
                        tcga_clinical_data = paad_tcga_clinical_data,
                        selection_method = c("limma"),
                        custom_selection = NULL,
                        n_lr_selected = 5,
                        which_get = c("ligand", "receptor"),
                        is_signif = TRUE,
                        clinical_vars_in_model = c("stage")
                )

        testthat::expect_type(results, "list")
        testthat::expect_length(results, 2)
        testthat::expect_s3_class(results$model_data, "data.frame")
        testthat::expect_type(results$cox_models, "list")
        testthat::expect_length(results$cox_models, 4)
        testthat::expect_type(results$cox_models$coxph_lr_pairs, "list")
        testthat::expect_type(results$cox_models$coxph_gene_expression, "list")
        testthat::expect_type(results$cox_models$coxph_clinical.lr_pairs, "list")
        testthat::expect_type(results$cox_models$coxph_clinical.gene_expression, "list")
})
