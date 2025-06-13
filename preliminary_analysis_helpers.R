# Create query matrix
create_query_matrix <- function(all_cols, cols_of_interest) {
  query_indices <- (1:length(all_cols))[all_cols %in% cols_of_interest]
  query_matrix <- matrix(0, ncol=length(all_cols), nrow=length(cols_of_interest)) # rows represent queries
  for (i in 1:length(query_indices)) {
    query_matrix[i,query_indices[i]] <- 1
  } 
  return(query_matrix)
}

# Create term-frequency tibble
count_code_frequency <- function(medications, clinical) {
  medications <- medications %>% 
    mutate(code=paste0("rx", THERDTL)) %>% select(pregid, code) 
  
  clinical <- clinical %>% mutate(code=paste0("dx", section_code)) %>% 
    select(pregid, code)

  freqs <- bind_rows(medications, clinical) %>%
    group_by(pregid, code) %>% summarise(occurrence=n()) %>% ungroup() 
  
  # split data to handle integer overflow
  ids <- unique(freqs$pregid)
  first_half_ids <- sample(ids, size=round(length(ids)/2), replace=FALSE)
  second_half_ids <- ids[!(ids %in% first_half_ids)]
  first_half_freqs <- freqs %>% filter(pregid %in% first_half_ids)
  second_half_freqs <- freqs %>% filter(pregid %in% second_half_ids)
  first_half_freqs <- first_half_freqs %>%
    arrange(code) %>%
    pivot_wider(names_from=code, values_from=occurrence, values_fill=0) 
  second_half_freqs <- second_half_freqs %>% 
    arrange(code) %>% 
    pivot_wider(names_from=code, values_from=occurrence, values_fill=0)
  freqs <- bind_rows(first_half_freqs, second_half_freqs) %>%
    right_join(pregnancies %>% select(pregid), by="pregid") %>%
    mutate(across(c(starts_with("rx"), starts_with("dx")), ~ if_else(is.na(.x), 0L, .x))) 
  return(freqs)
}

# Tabulate exposure-outcome 2x2 table
tabulate_exposure_malformations <- function(coded_data, malformations, code_var, code, malformation) {
  tabulated <- coded_data %>% filter({{code_var}} == code) %>% 
    group_by(pregid) %>% filter(row_number() == 1) %>% ungroup() %>%
    mutate(exposed = 1) %>% right_join(malformations, by="pregid") %>% mutate(exposed = replace_na(exposed,0)) %>% 
    tabyl(exposed, .data[[malformation]]) 
  return(tabulated)
}

# Calculate deviance
calculate_deviance <- function(preds, labels) {
  deviance <- -2*(1/length(preds))*sum(labels*log(preds) + (1-labels)*log(1-preds))
  return(deviance)
}

# Calculate Brier score
calculate_brier <- function(preds, labels) {
  brier <- (1/length(preds))*sum((preds - labels)^2)
  return(brier)
}

# Calculate MSE
calculate_mse <- function(preds, labels) {
  mse <- mean((preds - labels)^2)
  return(mse)
}

# Calculate MAE 
calculate_mae <- function(preds, labels) {
  mae <- mean(abs(preds-labels))
  return(mae)
}

# Calculate performance 
calculate_performance <- function(preds, labels) {
  perf <- list()
  perf$deviance <- calculate_deviance(preds, labels)
  perf$brier <- calculate_brier(preds, labels)
  perf$auc <- Metrics::auc(labels, preds)
  perf$mse <- calculate_mse(preds, labels)
  perf$mae <- calculate_mae(preds, labels)
  return(perf)
}


# Calculate p-values for given codes
calculate_pvals <- function(coded_data, malformations, code_var, description_var, codes, look_up, malformation, exact=FALSE, midp=FALSE, boschloo_test=FALSE) {
  # specify p-value method
  if (exact) {
    if (midp) {
      pvalue_test <- function(x) {x[,2:3] %>% as.matrix() %>% exact2x2::fisher.exact(midp=TRUE)}
    } else {
      pvalue_test <- function(x) {x[,2:3] %>% as.matrix() %>% exact2x2::fisher.exact(tsmethod="central")}  
    }
  } else if (boschloo_test) {
      pvalue_test <- function(x) {boschloo(x[2,3], x[2,3] + x[2,2], x[1,3], x[1,3] + x[1,2], midp=TRUE)}
  } else {
    pvalue_test <- janitor::chisq.test
  }
  # calculate p-values
  pval_results <- NULL
  for (i in 1:length(codes)) {
    code <- codes[i]
    code_description <- {{look_up}} %>% filter({{code_var}} == code) %>% pull({{description_var}})
    if (length(code_description) < 1) {
      code_description <- "Missing description"
    }
    tabulated <- tabulate_exposure_malformations(coded_data, malformations, {{code_var}}, code, malformation) 

    pval <- tabulated %>% pvalue_test() %>% {.[["p.value"]]}
    new_row <- tibble_row(code = code, description = code_description, x1y1=tabulated[2,3],x1=(tabulated[2,3] + tabulated[2,2]), perc_x1=100*tabulated[2,3]/(tabulated[2,3] + tabulated[2,2]),
                          perc_x0=100*tabulated[1,3]/(tabulated[1,3] + tabulated[1,2]),
                          abs_ln_rr=abs(log(perc_x1/perc_x0)),
                          p.value=pval)
    pval_results <- bind_rows(pval_results, new_row) 
  }
  return(pval_results)
}

# Make data wide
make_wide <- function(coded_data, code_var, malformations, threshold=1L) {
  coded_data <- coded_data %>% 
    mutate(coded_var = str_trim(format({{code_var}}, scientific=FALSE))) %>% 
    select(coded_var, pregid) %>%
    group_by(coded_var, pregid) %>% 
    summarise(occurrence=n()) %>%
    ungroup() %>%
    mutate(occurrence = as.integer(occurrence >= threshold)) %>%
    pivot_wider(names_from=coded_var, values_from=occurrence, names_prefix="code_") %>%
    right_join(malformations,by="pregid") %>% mutate(across(starts_with("code_"), ~ if_else(is.na(.x), 0, .x)))
  return(coded_data)
}

# Calculate p-values for given codes with age
calculate_pvals_lgr <- function(coded_data, code_var, description_var, codes, look_up, malformation) {
  # calculate p-values
  pval_results <- NULL
  
  reduced_model <- speedglm(as.formula(glue("{malformation} ~ rcs(AGE, 3)")), data=coded_data %>% select(AGE, all_of(malformation)), family=binomial(link="logit"))
  
  estimate_pvals <- function(codename) {
  
    code_description <- {{look_up}} %>% filter({{code_var}} == codename) %>% pull({{description_var}})
    if (length(code_description) < 1) {
      code_description <- "Missing description"
    }
    full_model <-speedglm(as.formula(glue("{malformation} ~ rcs(AGE, 3) + code_{format(codename, scientific=FALSE)}")), data=coded_data %>% select(AGE, all_of(c(malformation, glue("code_{format(codename, scientific=FALSE)}")))), 
                           family=binomial(link="logit"))
    crude_model <- speedglm(as.formula(glue("{malformation} ~ code_{format(codename, scientific=FALSE)}")), data=coded_data %>% select(AGE, all_of(c(malformation, glue("code_{format(codename, scientific=FALSE)}")))), 
                       family=binomial(link="logit"))
    
    if ((full_model$convergence) & (reduced_model$convergence)) {
      lr_stat <- -2 * (reduced_model$logLik - full_model$logLik)
      df_diff <- length(full_model$coefficients) - length(reduced_model$coefficients)
      pval <- pchisq(lr_stat, df=df_diff, lower.tail=FALSE)
    } else {
      pval <- NA_real_
    }
    
    x1y1 <- coded_data[coded_data[malformation] == 1 & coded_data[glue("code_{format(codename, scientific=FALSE)}")] == 1,] %>% nrow()
    x0y1 <- coded_data[coded_data[malformation] == 1 & coded_data[glue("code_{format(codename, scientific=FALSE)}")] == 0,] %>% nrow()
    x1 <- coded_data[coded_data[glue("code_{format(codename, scientific=FALSE)}")] == 1,] %>% nrow() 
    x0 <- coded_data[coded_data[glue("code_{format(codename, scientific=FALSE)}")] == 0,] %>% nrow()
    
    ln_or <- crude_model %>% broom::tidy() %>% filter(str_detect(term, glue("code_{format(codename, scientific=FALSE)}"))) %>% {.[["estimate"]]}
    ln_aor <- full_model %>% broom::tidy() %>% filter(str_detect(term, glue("code_{format(codename, scientific=FALSE)}"))) %>% {.[["estimate"]]}
    stderr_ln_or <- crude_model %>% broom::tidy() %>% filter(str_detect(term, glue("code_{format(codename, scientific=FALSE)}"))) %>% {.[["std.error"]]}
    stderr_ln_aor <- full_model %>% broom::tidy() %>% filter(str_detect(term, glue("code_{format(codename, scientific=FALSE)}"))) %>% {.[["std.error"]]}

    new_row <- tibble_row(code = codename, description = code_description, x1y1=x1y1, x1=x1, perc_x1=100*x1y1/x1,
                          perc_x0=100*x0y1/x0,
                          ln_or=ln_or, ln_aor=ln_aor,
                          stderr_ln_or=stderr_ln_or, stderr_ln_aor=stderr_ln_aor,
                          p.value=pval)
    return(new_row)
  }
  pval_results <- mclapply(codes, estimate_pvals, mc.cores=6) %>% bind_rows()
  return(pval_results)
}


# Create data table of p-values
create_pval_table <- function(pval_results) {
  pval_table <- pval_results %>% arrange(p.value) %>% select(code, description, x1, x1y1, perc_x1, perc_x0, p.value, p.adjusted) %>%
    datatable(class = 'cell-border stripe', rownames=FALSE, colnames=c("Code", "Description", "No. exposed", "No. exposed with outcome", "Percentage of exposed", "Percentage of unexposed", "p.value", "p.adjusted")) %>%
    formatRound(columns=c('perc_x1', 'perc_x0'), digits=2) %>% formatRound(columns=c("p.value", "p.adjusted"), digit=7)
  return(pval_table)
}

# Create static Manhattan plot
create_static_manhattan <- function(labelled_results, group_var, group_lab, fdrs, scale_y = 20) {
  labelled_results <- labelled_results 
  
  static_manhattan <- ggplot(data = labelled_results, aes(x=fct_rev({{group_var}}), y=-log10(p.value), stroke=1.2)) + 
    ylab("-log10(p-value)") + 
    geom_jitter(alpha=0.75, aes(colour={{group_var}}), show.legend=TRUE) +
    guides(colour="none")
  
  if (fdrs[1] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[1]), linetype="1%"), color="#b51963")
  }
  if (fdrs[2] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[2]), linetype="5%"), color="#0073e6")
  }
  if (fdrs[3] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[3]), linetype="10%"), color="#89ce00")
  }
  
  if (sum(fdrs) > 0) {
    static_manhattan <- static_manhattan + scale_linetype_manual(name = "FDR", values = c("1%" = 2,"5%" = 2, "10%" = 2)[fdrs > 0], 
                                                                     guide = guide_legend(override.aes = list(color = c("#b51963", "#0073e6", "#89ce00")[fdrs > 0])))
  }
    
  static_manhattan <- static_manhattan  +
    theme_minimal() + 
    #theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=-0.05)) +
    xlab(group_lab) 

  static_manhattan <- static_manhattan + 
      ylim(NA, scale_y)
  
  static_manhattan <- static_manhattan + scale_x_discrete(drop=FALSE) + coord_flip()
  return(static_manhattan)
}

# Create sized static Manhattan plot
create_sized_static_manhattan <- function(labelled_results, group_var, group_lab, padjusted=FALSE) {
  if (padjusted) {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.adjusted))) + 
      ylab("-log10(p-value)") + 
      geom_hline(aes(yintercept = -log10(0.05)), color = "red")
  } else {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.value))) + 
      ylab("-log10(p-value)")
  }
  static_manhattan <- static_manhattan + 
    geom_jitter(alpha=0.75, aes(colour={{group_var}}, size=abs_ln_rr), show.legend=FALSE) +
    theme_minimal() + 
    theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=-0.05)) +
    xlab(group_lab)
  return(static_manhattan)
}

# Create sized static Manhattan plot with threshold
create_sized_threshold_static_manhattan <- function(labelled_results, group_var, group_lab, padjusted=FALSE) {
  labelled_results <- labelled_results %>% mutate(rr = (perc_x1/perc_x0)) %>% 
    mutate(rrabove2 = if_else(rr > 2, "Yes", "No")) %>% 
    mutate(rrabove2 = factor(rrabove2, levels=c("Yes", "No")))
  
  if (padjusted) {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.adjusted))) + 
      ylab("-log10(p-value)") + 
      geom_hline(aes(yintercept = -log10(0.05)), color = "red")
  } else {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.value))) + 
      ylab("-log10(p-value)")
  }
  static_manhattan <- static_manhattan + 
    geom_jitter(alpha=0.75, stroke=1.5, aes(colour={{group_var}}, size=abs_ln_rr, shape=rrabove2), show.legend=TRUE) +
    scale_shape_manual(values=c(19, 21), name="Risk ratio above 2") + 
    theme_minimal() + 
    theme(legend.position="top", axis.text.x=element_text(angle = -90, hjust = 0, vjust=-0.05)) +
    xlab(group_lab) + 
    guides(colour="none", size="none")
  return(static_manhattan)
}

# Create size threshold static Manhattan plot
create_size_threshold_static_manhattan <- function(labelled_results, group_var, group_lab, padjusted=FALSE) {
  labelled_results <- labelled_results %>% mutate(rr = (perc_x1/perc_x0)) %>% 
    mutate(rrabove2 = if_else(rr > 2, "Yes", "No")) %>% 
    mutate(rrabove2 = factor(rrabove2, levels=c("Yes", "No")))
    
  if (padjusted) {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.adjusted))) + 
      ylab("-log10(p-value)") + 
      geom_hline(aes(yintercept = -log10(0.05)), color = "red")
  } else {
    static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=-log10(p.value))) + 
      ylab("-log10(p-value)")
  }
  static_manhattan <- static_manhattan + 
    geom_jitter(alpha=0.75, stroke=1.5, size=1, aes(colour={{group_var}}, shape=rrabove2), show.legend=TRUE) +
    scale_shape_manual(values=c(19, 21), name="Risk ratio above 2") + 
    theme_minimal() + 
    theme(legend.position="top", axis.text.x=element_text(angle = -90, hjust = 0, vjust=-0.05)) +
    xlab(group_lab) + 
    guides(colour="none")
  return(static_manhattan)
}

# Create static Manhattan subway plot
create_static_manhattan_subway <- function(labelled_results, group_var, group_lab, fdrs, scale_y=20) {
  
  labelled_results <- labelled_results %>% 
    mutate(minuslog10pval = -log10(p.value)) %>%
    mutate(large = minuslog10pval > scale_y) %>%
    mutate(minuslog10pval = sign(ln_aor)*if_else(minuslog10pval > scale_y, scale_y, minuslog10pval))

  static_manhattan <- ggplot(data = labelled_results, aes(x=fct_rev({{group_var}}), y=minuslog10pval)) + 
    geom_jitter(alpha=0.75, aes(colour={{group_var}}, shape=large, size=large), show.legend=FALSE) + 
    theme_minimal() + 
    theme(axis.text.y = element_text(size = 12), axis.title = element_text(size=12), legend.text = element_text(size =12)) + 
    scale_shape_manual(values=c(19, 8)) + 
    scale_size_manual(values=c(1.5,3)) +
    #theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
    labs(x=group_lab, y="sign(ln(aOR))*-log10(p-value)")
  
  if (fdrs[1] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[1]), linetype="1%"), color="#b51963")
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=log10(fdrs[1]), linetype="1%"), color="#b51963")
  }
  if (fdrs[2] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[2]), linetype="5%"), color="#0073e6")
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=log10(fdrs[2]), linetype="5%"), color="#0073e6")
  }
  if (fdrs[3] > 0) {
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=-log10(fdrs[3]), linetype="10%"), color="#89ce00")
    static_manhattan <- static_manhattan + geom_hline(aes(yintercept=log10(fdrs[3]), linetype="10%"), color="#89ce00")
  }
  
  if (sum(fdrs) > 0) {
    static_manhattan <- static_manhattan + scale_linetype_manual(name = "FDR", values = c("1%" = 2,"5%" = 2, "10%" = 2)[fdrs > 0], labels=c("1%", "5%", "10%")[fdrs > 0],
                                                                 guide = guide_legend(override.aes = list(color = c("#b51963", "#0073e6", "#89ce00")[fdrs > 0])))
  }
  

    static_manhattan <- static_manhattan + 
      ylim(-10, scale_y + 1)

  static_manhattan <- static_manhattan + scale_x_discrete(drop=FALSE) + coord_flip()
  return(static_manhattan)
}

warn_if <- function(condition, msg = "Condition failed") {
  if (!condition) {
    warning(msg)
    return(FALSE)
  }
  return(TRUE)
}


# Plot p-value against ratio
create_pval_ratio_plot <- function(pval_data, fdrs, logratio, ratio_description="risk", short_ratio_description="RR", max_ratio=6, min_ratio=0.3, max_log10_pval=20) {
  pval_data <- pval_data %>% 
    mutate(logratio = {{logratio}}) %>%
    mutate(abs_logratio = abs(logratio),
           minuslog10_pvalue = -log10(p.value)) %>%
    mutate(large = minuslog10_pvalue > max_log10_pval,
           logratio = if_else(is.infinite(abs_logratio), NA_real_, logratio),
           abs_logratio = if_else(is.infinite(abs_logratio), NA_real_, abs_logratio),
           minuslog10_pvalue = if_else(large, max_log10_pval, minuslog10_pvalue),
           ratio = exp(logratio)) %>%
    mutate(description = if_else(-log10(p.value) > -log10(fdrs[1]), str_wrap(description, 30), NA_character_))
  

    pval_ratio_plot <- ggplot(data=pval_data, aes(x=ratio, y=minuslog10_pvalue)) + 
      geom_point(aes(color=codetype, shape=large)) +
      theme_minimal() + 
      labs(y="-log10(p-value)",x=glue("{str_to_sentence(ratio_description)} ratio")) +
      guides(colour=guide_legend(position="inside"), linetype=guide_legend(position="inside", override.aes = list(color = c("#b51963","#0073e6", "#5ba300")[fdrs > 0])))

  
  if (fdrs[1] > 0) {
    pval_ratio_plot <- pval_ratio_plot + geom_hline(aes(yintercept=-log10(fdrs[1]), linetype="1%"), color="#b51963")
  }
  if (fdrs[2] > 0) {
    pval_ratio_plot <- pval_ratio_plot + geom_hline(aes(yintercept=-log10(fdrs[2]), linetype="5%"), color="#0073e6")
  }
  if (fdrs[3] > 0) {
    pval_ratio_plot <- pval_ratio_plot + geom_hline(aes(yintercept=-log10(fdrs[3]), linetype="10%"), color="#5ba300")
  }
  
  if (sum(fdrs) > 0) {
    pval_ratio_plot <- pval_ratio_plot + scale_linetype_manual(name = "FDR", values = c("1%" = 2,"5%" = 2, "10%" = 2)[fdrs > 0], labels=c("1%", "5%", "10%")[fdrs > 0]
                                                         )
  }
    
  pval_ratio_plot <- pval_ratio_plot + 
    ylim(0, max_log10_pval) + 
    scale_x_continuous(trans='log',breaks=c(1:max_ratio), limits=c(min_ratio, max_ratio)) +
    geom_text_repel(aes(label=description),  size=8/.pt, max.overlaps=20) + 
    scale_shape_manual(values=c(19, 8))  + 
    guides(shape = "none") +
    labs(color="Category") + 
    theme(legend.position.inside = c(0.83, 0.83),
          legend.box.background = element_rect(color = "black", linewidth = 1),
          legend.text = element_text(size = 8),    # Change the size of the legend items
          legend.title = element_text(size = 10),
          legend.box.margin = margin(1, 1, 1, 1))
  
  warn_if(max(pval_data$logratio, na.rm=TRUE) <= log(max_ratio),"max ratio outside range") 
  warn_if(min(pval_data$logratio, na.rm=TRUE) >= log(min_ratio), "min ratio outside range")
  return(pval_ratio_plot)
}

# Create sized Manhattan subway plot
create_sized_static_manhattan_subway <- function(labelled_results, group_var, group_lab) {
  labelled_results <- labelled_results %>% mutate(log10p.value = -sign(ln_aor)* log10(p.value))
  static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=log10p.value)) + geom_jitter(alpha=0.75, aes(colour={{group_var}}, size=abs_ln_rr), show.legend=FALSE) + 
    theme_minimal() + 
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
    labs(x=group_lab, y="-sign(ln(aOR))*log10(p-value)")
  return(static_manhattan)
}

# Create size threshold Manhattan subway plot
create_size_threshold_static_manhattan_subway <- function(labelled_results, group_var, group_lab) {
  labelled_results <- labelled_results %>% mutate(log10p.value = -sign(ln_aor)*log10(p.value)) %>% mutate(rrabove2 = (perc_x1/perc_x0) > 2)
  static_manhattan <- ggplot(data = labelled_results, aes(x={{group_var}}, y=log10p.value)) + geom_jitter(alpha=0.75, aes(colour={{group_var}}, shape=rrabove2), show.legend=FALSE) + 
    scale_shape_manual(values=c(16,17), name="Significant")+
    theme_minimal() + 
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) + 
    labs(x=group_lab, y="-sign(ln(aOR))*log10(p-value)")
  return(static_manhattan)
}

# Create dynamic Manhattan plot Rx
create_dynamic_manhattan_rx <- function(labelled_results) {
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~description, jitter=1, x=~THRGRDS, y=~-log10(p.value), color=~THRGRDS, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="Therapeutic group"), yaxis=list(title="-log10(p-value)")) 
  return(dynamic_manhattan)
}

# Create dynamic Manhattan plot ICD-9
create_dynamic_manhattan_icd9 <- function(labelled_results) {
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~description, jitter=1, x=~chapter, y=~-log10(p.value), color=~chapter, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="ICD-9 Chapter"), yaxis=list(title="-log10(p-value)")) 
  return(dynamic_manhattan)
}

# Create dynamic Manhattan plot ICD-10
create_dynamic_manhattan_icd10 <- function(labelled_results) {
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~category, jitter=1, x=~chapter, y=~-log10(p.value), color=~chapter, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="ICD-10 Chapter"), yaxis=list(title="-log10(p-value)")) 
  return(dynamic_manhattan)
}

# Create dynamic Manhattan subway plot Rx
create_dynamic_manhattan_subway_rx <- function(labelled_results) {
  labelled_results <- labelled_results %>% mutate(log10p.value = -sign(ln_aor)*log10(p.value))
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~description, jitter=1, x=~THRGRDS, y=~log10p.value, color=~THRGRDS, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="Therapeutic group"), yaxis=list(title="-sign(ln(aOR))*log10(p-value)"))
  return(dynamic_manhattan)
}

# Create dynamic Manhattan subway plot ICD-9
create_dynamic_manhattan_subway_icd9 <- function(labelled_results) {
  labelled_results <- labelled_results %>% mutate(log10p.value = if_else(perc_x1 > perc_x0, -log10(p.value), log10(p.value)))
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~description, jitter=1, x=~chapter, y=~log10p.value, color=~chapter, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="ICD-9 chapter"), yaxis=list(title="-sign(ln(aOR))*log10(p-value)"))
  return(dynamic_manhattan)
}

# Create dynamic Manhattan subway plot ICD-10
create_dynamic_manhattan_subway_icd10<- function(labelled_results) {
  labelled_results <- labelled_results %>% mutate(log10p.value = -sign(ln_aor)*log10(p.value))
  dynamic_manhattan <- plot_ly(data=labelled_results, type="box", boxpoints="all", text=~description, jitter=1, x=~chapter, y=~log10p.value, color=~chapter, line = list(color = 'rgba(255,255,255,0)'), fillcolor = list(color = 'rgba(255,255,255,0)')) %>%
    layout(showlegend=FALSE, xaxis=list(title="ICD-10 chapter"), yaxis=list(title="-sign(ln(aOR))*log10(p-value)"))
  return(dynamic_manhattan)
}

# create ICD9 section from ICD9 code 
categorise_icd9_section <- function(icd9_codes) {
  icd9_section <- if_else(str_detect(icd9_codes, "^E"), substr(icd9_codes, 1, 4), substr(icd9_codes, 1, 3))
  return(icd9_section)
}

# convenience function for formatting numbers
format_number <- function(a_number) {
  formatted_number <- scales::number(a_number, big.mark=",")
  return(formatted_number)
}

# convenience function for formatting percentages
format_percent <- function(a_number) {
  formatted_number <- scales::number(a_number, accuracy=0.1)
  return(formatted_number)
}

# function to read in ICD10 mapping file
read_in_icd10_mapping <- function(mapping_file) {
  mapping_table <- read.table(mapping_file, sep="" ,stringsAsFactors = FALSE, colClasses=c(rep("character",3))) %>% 
    as_tibble() %>% rename(icd10=V1, icd9=V2, flags=V3) %>%
    mutate(approximate=as.integer(substr(flags,1,1)),
           nomap=as.integer(substr(flags,2,2)),
           combination=as.integer(substr(flags,3,3)),
           scenario=as.integer(substr(flags,4,4)),
           choicelist=as.integer(substr(flags,5,5)))
  return(mapping_table)
}

# function to create code map to ICD-9 section from ICD-10 code
create_map_to_icd9_section <- function(icd10_to_icd9, medical_diagnoses_icd10) {
  # select ICD-10 codes in data
  icd10_codes_in_data <- medical_diagnoses_icd10 %>%
    group_by(codedp) %>%
    filter(row_number() == 1) %>%
    ungroup() %>% 
    select(codedp)
  
  # keep ICD-10 codes that can map to one or more ICD-9 sections representing distinct concepts
  code_map <- icd10_codes_in_data %>% 
    inner_join(icd10_to_icd9, by=c("codedp" = "icd10")) %>%
    filter(nomap == 0) %>% 
    mutate(section_code=categorise_icd9_section(icd9)) %>%
    group_by(codedp, section_code) %>% 
    filter(row_number() == 1) %>% ungroup() %>% 
    group_by(codedp) %>% 
    mutate(no_scenarios = length(unique(scenario))) %>% 
    filter(no_scenarios == 1 & (combination == 1 | n() ==1)) %>%
    select(codedp, section_code)
  return(code_map)
}

find_max_benjamini_hochberg <- function(pvalues, alpha) {
  
  # Sort p-values
  sorted_pvalues <- sort(pvalues)
  
  # Get length of vector 
  no_pvalues <- length(pvalues)
  
  # Create index vector
  index_vector <- 1:no_pvalues
  
  # Create adjusted p-value
  adjusted_pvalues <- no_pvalues*sorted_pvalues/index_vector
  
  # Identify p-values below threshold
  sig_pvalues <- sorted_pvalues[adjusted_pvalues < alpha]
  
  # Find max p-value
  if (length(sig_pvalues) < 1) {
    max_pvalue <- 0
  } else {
    max_pvalue <- max(sig_pvalues)
  }
  
  return(max_pvalue)
}

find_significant_benjamini_hochberg <- function(pvalues, alpha) {
  
  # Sort p-values
  sorted_pvalues <- sort(pvalues)
  
  # Get length of vector 
  no_pvalues <- length(pvalues)
  
  # Create index vector
  index_vector <- 1:no_pvalues
  
  # Create adjusted p-value
  adjusted_pvalues <- no_pvalues*sorted_pvalues/index_vector
  
  # Identify p-values below threshold
  sig_pvalues <- sorted_pvalues[adjusted_pvalues < alpha]
  
  # Find max p-value
  if (length(sig_pvalues) < 1) {
    max_pvalue <- 0
  } else {
    max_pvalue <- max(sig_pvalues)
  }
  
  # Return boolean vector
  return(pvalues <= max_pvalue)
}

# max keep na
max_keep_na <- function(a_vec) {
  if (all(is.na(a_vec))) {
    return(NA_real_)
  } else {
    return(max(a_vec, na.rm=TRUE))
  }
}

# cosine distance between two vectors
cosine_distance <- function(vec1, vec2) {
  cs_distance <- vec1 %*% vec2 /(sqrt(vec1 %*% vec1) * sqrt(vec2 %*% vec2))
  return(cs_distance[1,1])
}

binom_lower <- function(x, num) {
  return(mapply(function(x, num) binom.wilson(x, num)["lower"][1,1], x, num))
}
binom_upper <- function(x, num) {
  return(mapply(function(x, num) binom.wilson(x, num)["upper"][1,1], x, num))
}