#' Survival analysis for scACCorDiON results
#'
#' @param  comparison_oi TODO
#'
#' @param  lr_data TODO
#'
#' @param  tcga_exp_data TODO
#'
#' @param  tcga_clinical_data TODO
#'
#' @param  selection_method TODO
#'
#' @param  custom_selection TODO
#'
#' @param  n_lr_selected TODO
#'
#' @param  which_get TODO
#'
#' @param  is_signif TODO
#'
#' @param  clinical_vars_in_model TODO
#'
#' @return TODO
#'
#' @export
make_comparison <- function(
  comparison_oi,
  lr_data,
  tcga_exp_data,
  tcga_clinical_data,
  selection_method = c("limma", "prauc"),
  custom_selection = NULL,
  n_lr_selected = 10,
  which_get = c("ligand", "receptor"),
  is_signif = TRUE,
  clinical_vars_in_model = NULL
){

  symbol <- NULL

  assertthat::assert_that(comparison_oi %in% colnames(lr_data))
  assertthat::assert_that("LR" %in% colnames(lr_data))
  assertthat::assert_that(length(selection_method) == 1)
  assertthat::assert_that(selection_method %in% c("limma", "prauc"))
  assertthat::assert_that(all(which_get %in% c("ligand", "receptor")))
  assertthat::assert_that(is.logical(is_signif))
  assertthat::assert_that("symbol" %in% colnames(tcga_exp_data))
  assertthat::assert_that(all(c("sample_id","patient") %in% colnames(tcga_clinical_data)))
  assertthat::assert_that(all(c("time","event") %in% colnames(tcga_clinical_data)))
  assertthat::assert_that(
    is.null(clinical_vars_in_model) | all(clinical_vars_in_model %in% colnames(tcga_clinical_data))
  )

  genes_oi <- NULL

  if(!is.null(custom_selection)){

    assertthat::assert_that(all(colnames(custom_selection) == c("ligand","receptor")))
    genes_oi <- custom_selection |> dplyr::select(dplyr::all_of(which_get))

  }else{

    lr_oi <- get_lr_oi(mat = lr_data |> tibble::column_to_rownames("LR"), col_oi = comparison_oi)
    genes_oi <- get_lr_genes(
      lr_oi_obj = lr_oi,
      set_name = selection_method,
      n_lr = n_lr_selected,
      which_get = which_get,
      is_ofinterest = is_signif
    )

  }

  # assert that our genes of interest are present in the data
  # otherwise these should be removed apriori from pool of genes of interest
  assertthat::assert_that(all(unlist(genes_oi) %in% tcga_exp_data$symbol))
  #assertthat::assert_that(all(unlist(genes_oi) %in% tcga_exp_data$symbol))

  tcga_survready <-
    tcga_exp_data |>
    dplyr::filter(symbol %in% (genes_oi |> unlist() |> unique())) |>
    tibble::column_to_rownames("symbol") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("sample_id") |>
    tibble::as_tibble()

  #########
  # log2 normalization
  #########
  tcga_survready <-
    tcga_survready |>
    dplyr::mutate(
      dplyr::across(tidyselect::where(is.numeric), (\(x) log2(x + 1)))
    )

  # joining clinical and expression data together
  survgene_data <-
    tcga_clinical_data |> dplyr::left_join(tcga_survready, by = "sample_id")


  ### this shouldnt be necessary anymore as we assert it at an earlier state
  # excluding interaction because gene not in data
  # genes_oi <- genes_oi |> dplyr::filter(!ligand %in% "MIF")

  # if we are only getting ligand OR receptor, then we dont build lr scores
  if (!identical(which_get, c("ligand", "receptor"))) {
    res_cox <- compute_cox_ph_models(
      data = survgene_data,
      clinical_vars = clinical_vars_in_model,
      genes_oi = genes_oi,
      which_get = which_get
    )

    return(list(model_data = survgene_data, cox_models = res_cox))
  }

  survgenelr_data <-
    build_lr_score(
      survgene_dat = survgene_data,
      lr_genes_df = genes_oi
    )

  res_cox <- compute_cox_ph_models(
    data = survgenelr_data,
    clinical_vars = clinical_vars_in_model,
    genes_oi = genes_oi,
    which_get = which_get
  )

  return(list(model_data = survgenelr_data, cox_models = res_cox))
}


#' @rdname make_comparison
#' @export
run_survival <- make_comparison

