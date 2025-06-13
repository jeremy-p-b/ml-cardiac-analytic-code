# Combine and add FDR --------------------------------------

for (outcome in OUTCOME_VARS) {
  # Read in data
  maternal_rx_pval_results_top500_lgr <- read_csv(glue("output/maternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
  maternal_dx_pval_results_top500_lgr <- read_csv(glue("output/maternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
  paternal_rx_pval_results_top500_lgr <- read_csv(glue("output/paternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
  paternal_dx_pval_results_top500_lgr <- read_csv(glue("output/paternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
  
  # Combine data
  pval_results_top500_lgr <- bind_rows(
    maternal_rx_pval_results_top500_lgr %>% mutate(code=as.character(code), parent="Mother", type="rx"),
    maternal_dx_pval_results_top500_lgr %>% mutate(parent="Mother", type="dx"),
    paternal_rx_pval_results_top500_lgr %>% mutate(code=as.character(code), parent="Father", type="rx"),
    paternal_dx_pval_results_top500_lgr %>% mutate(parent="Father", type="dx")
  ) 
  
  # Add indicator for significant at different FDRs
  pval_results_top500_lgr <- pval_results_top500_lgr %>% 
    mutate(fdr1 = find_significant_benjamini_hochberg(p.value, 0.01),
           fdr5 = find_significant_benjamini_hochberg(p.value, 0.05),
           fdr10 = find_significant_benjamini_hochberg(p.value, 0.1))
  
  # Add bootstrap corrected odds ratios
  pval_results_top500_lgr <- pval_results_top500_lgr %>% 
    mutate(rsid=row_number())
  
  br_ests <- pval_results_top500_lgr %>% mutate(beta=ln_aor, se=stderr_ln_aor) %>% 
    select(rsid, beta, se) %>% se_adjust(method="BR_ss") 
  
  pval_results_top500_lgr <- br_ests %>% 
    rename(ln_aor_corr=beta_BR_ss, stderr_ln_aor_corr=adj_se) %>% 
    select(rsid, ln_aor_corr, stderr_ln_aor_corr) %>% 
    inner_join(pval_results_top500_lgr, by="rsid") %>% 
    select(!rsid) %>% relocate(ln_aor_corr, .after=ln_aor) %>% 
    relocate(stderr_ln_aor_corr, .after=stderr_ln_aor) 
  
  # Order by p-value
  pval_results_top500_lgr <- pval_results_top500_lgr %>%
    arrange(p.value)
  
  # Save outputs
  write_csv(pval_results_top500_lgr, glue("output/pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
  file.remove(c(glue("output/maternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"),
            glue("output/maternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"),
            glue("output/paternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"),
            glue("output/paternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv")))
}
