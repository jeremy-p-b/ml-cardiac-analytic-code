# Calculate p-values for medical associations using logistic regression ------------------------------------------------

if (ICD_VERSION == "both") {
  maternal_dx <- maternal_icd_combined %>% mutate(code=section_code) 
  paternal_dx <- paternal_icd_combined %>% mutate(code=section_code) 
  lookup <- icd9_3digit_look_up 
  pregnancies_subset <- pregnancies
} else if (ICD_VERSION == "icd9") {
  maternal_dx <- maternal_medical %>% filter(SVCDATE < ymd("2015/10/01"), DXVER == 9) %>% mutate(code = substr(codedp, 1, NO_DIGITS_ICD))
  paternal_dx <- paternal_medical %>% filter(SVCDATE < ymd("2015/10/01"), DXVER == 9) %>% mutate(code = substr(codedp, 1, NO_DIGITS_ICD))
  if (NO_DIGITS_ICD == 3) {
    lookup <- icd9_3digit_look_up
  } else {
    lookup <- icd9_4digit_look_up
  }
  pregnancies_subset <- pregnancies %>% filter(lmp + 90 < ymd("2015/10/01"))
} else if (ICD_VERSION == "icd10") {
  maternal_dx <- maternal_medical %>% filter(SVCDATE >= ymd("2015/10/01"), DXVER == 0) %>% mutate(code = substr(codedp, 1, NO_DIGITS_ICD))
  paternal_dx <- paternal_medical %>% filter(SVCDATE >= ymd("2015/10/01"), DXVER == 0) %>% mutate(code = substr(codedp, 1, NO_DIGITS_ICD))
  if (NO_DIGITS_ICD == 3) {
    lookup <- icd10_3digit_look_up
  } else {
    lookup <- icd10_4digit_look_up
  }
  pregnancies_subset <- pregnancies %>% filter(lmp - 180 >= ymd("2015/10/01"))
} 

rm(maternal_medical, paternal_medical, maternal_icd_combined, paternal_icd_combined)
gc()

# Specify list of codes to screen
maternal_dx_codes <- maternal_dx$code %>% unique()
paternal_dx_codes <- paternal_dx$code %>% unique()

top_maternal_dx_codes <- maternal_dx %>% 
  fgroup_by(pregid, code) %>% 
  fsummarise(code_count = fnobs(pregid)) %>%
  fmutate(above_threshold = as.integer(code_count >= threshold_no)) %>%
  fgroup_by(code) %>%
  fsummarise(no_exposed = fnobs(pregid), no_above_threshold=fsum(above_threshold)) %>% 
  arrange(desc(no_above_threshold)) %>% 
  filter(row_number() <= 500) %>%
  filter(no_above_threshold > 10) %>%
  {.[["code"]]} %>% unique() 
  
top_paternal_dx_codes <- pregnancies_subset %>% 
  filter(!is.na(enrolid_dad)) %>% 
  select(pregid) %>%
  inner_join(paternal_dx, by="pregid") %>% 
  fgroup_by(pregid, code) %>% 
  fsummarise(code_count = fnobs(pregid)) %>%
  mutate(above_threshold = as.integer(code_count >= threshold_no)) %>%
  fgroup_by(code) %>%
  fsummarise(no_exposed = fnobs(pregid), no_above_threshold=fsum(above_threshold)) %>%
  arrange(desc(no_above_threshold)) %>% 
  filter(row_number() <= 500) %>%
  filter(no_above_threshold > 10) %>%
  {.[["code"]]} %>% unique() 

calculate_and_write_maternal_dx_pvals <- function(outcome) {
  maternal_dx_pval_results_top_lgr <- calculate_pvals_lgr(maternal_dx_wide, code, description, top_maternal_dx_codes, lookup, outcome)
  write_csv(maternal_dx_pval_results_top_lgr, glue("output/maternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
}
calculate_and_write_paternal_dx_pvals <- function(outcome) {
  paternal_dx_pval_results_top_lgr <- calculate_pvals_lgr(paternal_dx_wide, code, description, top_paternal_dx_codes, lookup, outcome)
  write_csv(paternal_dx_pval_results_top_lgr, glue("output/paternal_dx_pval_results_top500_{outcome}_lgr{singleton_stub}{window_stub}{icd_stub}{train_stub}{diabetes_stub}{threshold_stub}{prior_malformation_stub}{chrom_stub}{paternal_stub}.csv"))
}

# Calculate p-values Dx
maternal_dx <- maternal_dx %>% filter(code %in% top_maternal_dx_codes)
gc()
maternal_dx_wide <- make_wide(maternal_dx, code, pregnancies_subset, threshold_no)
datadist(maternal_dx_wide)
for (outcome in OUTCOME_VARS) {
  calculate_and_write_maternal_dx_pvals(outcome)
}
mclapply(OUTCOME_VARS, calculate_and_write_maternal_dx_pvals, mc.cores=8)
rm(maternal_dx_wide)
gc()

paternal_dx <- paternal_dx %>% filter(code %in% top_paternal_dx_codes)
gc()
paternal_dx_wide <- make_wide(paternal_dx, code, pregnancies_subset %>% filter(!is.na(enrolid_dad)), threshold_no)
datadist(paternal_dx_wide)
mclapply(OUTCOME_VARS, calculate_and_write_paternal_dx_pvals, mc.cores=4)
rm(paternal_dx_wide)
gc()