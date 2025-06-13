# Calculate p-values for drug associations using logistic regression ---------------------------------------------------

# Specify list of codes to screen
maternal_rx_codes <- maternal_medications %>% filter(!is.na(THERDTL)) %>% {.[["THERDTL"]]} %>% unique()
paternal_rx_codes <- paternal_medications %>% filter(!is.na(THERDTL)) %>% {.[["THERDTL"]]} %>% unique()

top_maternal_rx_codes <- maternal_medications %>%
  group_by(pregid, THERDTL) %>% 
  summarise(code_count = n()) %>%
  mutate(above_threshold = as.integer(code_count >= threshold_no)) %>%
  group_by(THERDTL) %>%
  summarise(no_exposed = n(), no_above_threshold=sum(above_threshold)) %>% 
  arrange(desc(no_above_threshold)) %>% 
    filter(row_number() <= 500) %>%
  filter(no_above_threshold >= 10) %>%
    {.[["THERDTL"]]} %>% unique() 
  
top_paternal_rx_codes <- paternal_medications %>%
  group_by(pregid, THERDTL) %>% 
  summarise(code_count = n()) %>%
  mutate(above_threshold = as.integer(code_count >= threshold_no)) %>%
  group_by(THERDTL) %>%
  summarise(no_exposed = n(), no_above_threshold=sum(above_threshold)) %>%
  arrange(desc(no_above_threshold)) %>% 
  filter(row_number() <= 500) %>%
  filter(no_above_threshold >= 10) %>%
  {.[["THERDTL"]]} %>% unique()

calculate_and_write_maternal_rx_pvals <- function(outcome) {
  maternal_rx_pval_results_top_lgr <- calculate_pvals_lgr(maternal_wide_medications, THERDTL, THRDTDS, top_maternal_rx_codes, rx_detail_look_up, outcome)
  write_csv(maternal_rx_pval_results_top_lgr, glue("output/maternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
}
calculate_and_write_paternal_rx_pvals <- function(outcome) {
  paternal_rx_pval_results_top_lgr <- calculate_pvals_lgr(paternal_wide_medications, THERDTL, THRDTDS, top_paternal_rx_codes, rx_detail_look_up, outcome)
  write_csv(paternal_rx_pval_results_top_lgr, glue("output/paternal_rx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
}

# Calculate p-values Rx
maternal_medications <- maternal_medications %>% filter(THERDTL %in% top_maternal_rx_codes)
gc()
maternal_wide_medications <- make_wide(maternal_medications, THERDTL, pregnancies, threshold_no)
datadist(maternal_wide_medications)
mclapply(OUTCOME_VARS, calculate_and_write_maternal_rx_pvals, mc.cores=3)
rm(maternal_wide_medications)
gc()

paternal_medications <- paternal_medications %>% filter(THERDTL %in% top_paternal_rx_codes)
gc()
paternal_wide_medications <- make_wide(paternal_medications, THERDTL, pregnancies %>% filter(!is.na(enrolid_dad)), threshold_no)
datadist(paternal_wide_medications)
mclapply(OUTCOME_VARS, calculate_and_write_paternal_rx_pvals, mc.cores=4)
rm(paternal_wide_medications)
gc()