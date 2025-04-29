
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scACCorDiON.su
![image](https://github.com/user-attachments/assets/e0015cf5-9eda-4ec8-a59e-38e0430d8d7a)

<!-- badges: start -->

<!-- badges: end -->

The goal of scACCorDiON.su is to run survival analysis on the results of
scACCorDiON.

## Installation

You can install the development version of scACCorDiON.su like so:

``` r
pak::pak("CostaLab/scACCorDiON.su")
```

## Example

This is a basic example of how you could run this pipeline:

``` r
library(scACCorDiON.su)
#> ---------------------------------
#> scACCorDiON.su version 0.0.0.9000
#> ---------------------------------

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
#> Warning: Zero sample variances detected, have been offset away from zero


results$cox_models$coxph_clinical.lr_pairs |>
  broom::tidy() |>
  print()
#> # A tibble: 8 × 5
#>   term              estimate std.error statistic p.value
#>   <chr>                <dbl>     <dbl>     <dbl>   <dbl>
#> 1 stageStage_II       0.244      0.430     0.567   0.571
#> 2 stageStage_III     -0.497      1.09     -0.458   0.647
#> 3 stageStage_IV       0.194      0.832     0.234   0.815
#> 4 PCSK9_VLDLR         0.139      0.102     1.36    0.172
#> 5 TNFSF11_TNFRSF11A   0.0823     0.123     0.668   0.504
#> 6 WNT7B_FZD1          0.154      0.129     1.19    0.233
#> 7 DKK1_KREMEN1        0.191      0.118     1.62    0.105
#> 8 L1CAM_EPHB2         0.0170     0.104     0.163   0.871
```
