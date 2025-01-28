

get_lr_oi <- function(mat, col_oi){

  adj.P.Val <- truth <- prauc <- NULL

  # get ligand-receptor of interest w/ limma and auprc
  cols <- colnames(mat)
  target_vec <- ifelse(col_oi == cols, "X1", "X0")

  colnames(mat) <- target_vec

  dsg_mat <- model.matrix(~0 + comp, data = data.frame(comp = target_vec))
  colnames(dsg_mat) <- colnames(dsg_mat) |> gsub(pattern = "comp", replacement = "")

  fit <- limma::lmFit(mat, dsg_mat)
  contrast_matrix <- limma::makeContrasts(
    contrasts = "X1-X0",
    levels = c("X0", "X1")
  )

  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  res <- limma::topTable(fit2, number = Inf, sort.by = "none")

  res <- res |> dplyr::arrange(adj.P.Val) |>
    tibble::rownames_to_column("LR") |>
    tibble::as_tibble()

  tmp_dat <- 
    t(mat) |> as.data.frame() |> tibble::rownames_to_column("truth") |> 
    dplyr::mutate(truth = factor(truth, levels = c("X1", "X0")))

  res_pr <-
    tmp_dat |> 
    dplyr::summarize(dplyr::across(dplyr::where(is.numeric), (\(x){yardstick::pr_auc_vec(truth = tmp_dat$truth, estimate = x)}))) |>
    t() |> as.data.frame() |> tibble::rownames_to_column("LR") |> tibble::as_tibble() |> setNames(c("LR","prauc")) |>
    dplyr::arrange(-prauc)

  return(list(limma = res, prauc = res_pr))
}

get_lr_genes <- function(lr_oi_obj, set_name, n_lr = 1,  which_get = c("ligand","receptor"), is_ofinterest = TRUE){

  adj.P.Val <- logFC <- LR <- prauc <- NULL

  # get ligand-receptor gene-pairs
  assertthat::assert_that(set_name %in% names(lr_oi_obj))
  assertthat::assert_that("LR" %in% colnames(lr_oi_obj[[set_name]]))

  df <- lr_oi_obj[[set_name]]
  if(set_name == "limma" & is_ofinterest){
    df <- df |> dplyr::filter(adj.P.Val < 0.05) |> dplyr::arrange(logFC)
  }
  if(set_name == "prauc" & is_ofinterest){
    df <- df |> dplyr::filter(prauc > 0.9)
  }

  lr_genes_df <-
    df[seq_len(min(n_lr, nrow(df))), "LR"] |>
    tidyr::separate_wider_delim(LR, delim = ":", names = c("ligand", "receptor")) |>
    dplyr::select(dplyr::all_of(which_get))
  return(lr_genes_df)
}

build_lr_score <- function(survgene_dat, lr_genes_df){

  # Building  LR score with geometric mean
  survgene_dat <- as.data.frame(survgene_dat)

  # testing if all LRs elements are present in data
  assertthat::assert_that(all(lr_genes_df$ligand %in% colnames(survgene_dat)))
  assertthat::assert_that(all(lr_genes_df$receptor %in% colnames(survgene_dat)))

  res_list <- purrr::map(seq_len(nrow(lr_genes_df)), dat = survgene_dat, function(i, dat) {
    lig <- lr_genes_df[i, "ligand"] |> unlist()
    rec <- lr_genes_df[i, "receptor"] |> unlist()
    lr_df <- data.frame(sample_id = dat[, "sample_id"], lig = dat[, lig], rec = dat[, rec])
    lr_name <- paste0(lig, "_", rec)
    res <-
      lr_df |>
      # dplyr::mutate(lig = log2(1 + lig)) |>
      # dplyr::mutate(rec = log2(1 + rec)) |>
      dplyr::rowwise() |>
      dplyr::mutate({{lr_name}} := DescTools::Gmean(c(lig, rec))) |>
      dplyr::select(-lig, -rec) |>
      dplyr::ungroup()
    return(res)
  })

  lr_tbl <- purrr::reduce(res_list, dplyr::left_join, by = "sample_id")
  f_res <- survgene_dat |> dplyr::left_join(lr_tbl, by = "sample_id")
  return(f_res)
}


compute_cox_ph_models <- function(data, genes_oi, clinical_vars = NULL, which_get = c("ligand","receptor")){

  ligand <- receptor <- lr <- NULL

  assertthat::assert_that(all(c("time", "event") %in% colnames(data)))
  assertthat::assert_that(all(which_get %in% colnames(genes_oi)))
  # assertthat::assert_that(all(genes_oi$ligand %in% colnames(data)))
  # assertthat::assert_that(all(genes_oi$receptor %in% colnames(data)))

  # There's 4 possibilies if using LRs:
  # RUN COX W/ CLINICAL DATA AND INDIVIDUAL LR GENES
  # RUN COX W/ CLINICAL DATA AND LR PAIRS
  # RUN COX W/ INDIVIDUAL LR GENES
  # RUN COX W/ LR PAIRS
  BASE_FORMULA_STRING <- "survival::Surv(time=time, event=event) ~ "

  res_list <- list()

  ind_genes <- genes_oi |> unlist() |> unique()

  genes_formula_str_token <- paste(ind_genes, collapse = "+")

  # formula for individual LR genes
  coxph_expr_formula <- as.formula(paste(BASE_FORMULA_STRING, genes_formula_str_token))

  # creating cox proportional hazards models
  # ind genes
  coxph_expr <- survival::coxph(coxph_expr_formula, data = data)
  coxph_expr$call$formula <- coxph_expr_formula

  res_list <- append(x = res_list, values = list("coxph_gene_expression"=coxph_expr))

  if (length(which_get) > 1) {
    assertthat::assert_that(all(genes_oi$ligand %in% colnames(data)))
    assertthat::assert_that(all(genes_oi$receptor %in% colnames(data)))

    lr_pairs <- genes_oi |> dplyr::mutate(lr = paste0(ligand, "_", receptor)) |> dplyr::pull(lr)

    assertthat::assert_that(all(lr_pairs %in% colnames(data)))
    lr_formula_str_token <- paste(lr_pairs, collapse = "+")
    # formula for LR pairs
    coxph_lr_formula <- as.formula(paste(BASE_FORMULA_STRING, lr_formula_str_token))
    # creating cox proportional hazards models
    # lr paurs
    coxph_lr <- survival::coxph(coxph_lr_formula, data = data)
    coxph_lr$call$formula <- coxph_lr_formula
    res_list <- append(x = res_list, values = list("coxph_lr_pairs"=coxph_lr))
  }


  if (!is.null(clinical_vars)) {

    clinical_formula_str_token <- paste0(clinical_vars, collapse = " + ")

    # formula for clinical data + individual LR genes
    coxph_clinexpr_formula <- as.formula(paste(BASE_FORMULA_STRING, clinical_formula_str_token, " + ", genes_formula_str_token))

    # creating cox proportional hazards models
    # clinical + ind genes
    coxph_clinexpr <- survival::coxph(coxph_clinexpr_formula, data = data)
    coxph_clinexpr$call$formula <- coxph_clinexpr_formula
    res_list <- append(x = res_list, values = list("coxph_clinical.gene_expression" = coxph_clinexpr))

    if (length(which_get) > 1) {
      # formula for clinical data + LR pairs
      coxph_clinlr_formula <- as.formula(paste(BASE_FORMULA_STRING, clinical_formula_str_token, " + ", lr_formula_str_token))

      # creating cox proportional hazards models
      # clinical + lr paurs
      coxph_clinlr <- survival::coxph(coxph_clinlr_formula, data = data)
      coxph_clinlr$call$formula <- coxph_clinlr_formula

      res_list <- append(x = res_list, values = list("coxph_clinical.lr_pairs" = coxph_clinlr))
    }

    return(res_list)
  }

  return(res_list)
}

