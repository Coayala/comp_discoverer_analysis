#
#
# Christian Ayala
# Functions for 3_Differential_analysis.Rmd
#
#
# -------------------------------------------------------------------------
get_samples <- function(metadata.df, Treatment, value){
  # Get value to filter samples
  selector <- syms({{Treatment}})

  samples <- metadata.df %>% 
    filter((!!! selector) == value)
  
  # Get only sample names
  samples <- samples$SampleID
  
  return(samples)
}

# -------------------------------------------------------------------------
get_vectors <- function(df, filter_by, value, get_col){
  # Column where the value will be filtered
  filter_col <- syms({{filter_by}})
  
  # Column that will be retrieved
  get <- syms({{get_col}})
  
  vector <- df %>% 
    filter((!!! filter_col) == value) %>% 
    pull((!!! get))
  
  return(vector)
}

# -------------------------------------------------------------------------
get_diff_table <- function(auc_matrix, control.sample_list, treatment.sample_list){
  
  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>% 
    select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control)
  
  temp.df_treatment <- auc_matrix %>% 
    select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment)
  
  diff_table <- as.data.frame(cbind(control_means, treatment_means))
  diff_table <- diff_table %>% 
    mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
    mutate(log2FC = log2(ratio)) # calculate log2FC
  
  # Initialize pvalues matrix
  pvalues <- data.frame(row.names = rownames(auc_matrix), pval = rep(0, length(rownames(auc_matrix))))

  #Calculate pvalue per each of the features
  for(i in 1:nrow(pvalues)){
    t.test <- t.test(as.numeric(temp.df_control[i,]), as.numeric(temp.df_treatment[i,]), paired = FALSE)
    pvalues$pval[i] <- t.test$p.value
  }

  diff_table <- merge(diff_table, pvalues, by = 'row.names')
  diff_table$pval.adj <- p.adjust(diff_table$pval, method = 'fdr')
  diff_table <- diff_table %>%
    rename(FeatureID = Row.names)

  return(diff_table)
  
}

# -------------------------------------------------------------------------
get_diff_table_no_pval <- function(auc_matrix, control.sample_list, treatment.sample_list){
  
  # Get the AUC values per each sample and calculate the means per feature
  temp.df_control <- auc_matrix %>% 
    select(all_of(control.sample_list))
  control_means <- rowMeans(temp.df_control)
  
  temp.df_treatment <- auc_matrix %>% 
    select(all_of(treatment.sample_list))
  treatment_means <- rowMeans(temp.df_treatment)
  
  diff_table <- as.data.frame(cbind(control_means, treatment_means))
  diff_table <- diff_table %>% 
    mutate(ratio = treatment_means/control_means) %>% # get the control/treatment ratio
    mutate(log2FC = log2(ratio)) # calculate log2FC
  
  return(diff_table)
  
}

# -------------------------------------------------------------------------
plot_col <- function(df, my_x, my_y, color_by1, color_by2, dodge = FALSE){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_col(aes(fill = {{color_by1}},
                 color = {{color_by2}}),
             size = 2,
             width = 0.75,
             position = ifelse(dodge == TRUE, 'dodge', 'stack')) +
    scale_fill_jco() +
    scale_color_jama() +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}

# -------------------------------------------------------------------------
plot_volcano <- function(df, log2FC, pval, log2FC.threshold, pval.threshold){
  
  #Generate label for the plot
  
  ggplot(df,
         aes(x = {{log2FC}},
             y = -log10({{pval}}))) +
    geom_point(color = ifelse(abs(df$log2FC) > {{log2FC.threshold}} & 
                                -log10(df$pval.adj) > -log10({{pval.threshold}}), "#FF0000", "#000000")) +
    geom_vline(xintercept = c(-{{log2FC.threshold}}, {{log2FC.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    geom_hline(yintercept = -log10({{pval.threshold}}),
               linetype = 'dotted',
               size = 1,
               color = 'blue') +
    geom_text_repel(data = . %>%
                      mutate(label = ifelse(abs(df$log2FC) > log2FC.threshold &
                                              -log10(df$pval.adj) > -log10(pval.threshold), name4plot, NA)),
                    aes(label = label),
                    size = 3,
                    max.overlaps = 20,
                    force = 2) +
    theme_bw() +
    labs(title = 'Volcano plot', 
         x = expression("Log"[2]*" Fold Change"), 
         y = expression("-Log"[10]*" pvalue")) +
    theme(plot.title = element_text(hjust = 0.5, 
                                    face = 'bold', 
                                    size = 18),
          plot.subtitle = element_text(hjust = 0.5, 
                                       face = 'bold', 
                                       size = 15))
}

# -------------------------------------------------------------------------
plot_venn <- function(my_list, my_colors){
  venn(my_list,
       zcolor = {{my_colors}},
       ilcs = 1,
       sncs = 1)
}

# -------------------------------------------------------------------------
plot_vank <- function(df, color_by, facet_by = NULL, facet_by2 = NULL){
  ggplot(df,
         aes(x = O_to_C,
             y = H_to_C,
             color = {{color_by}})) +
    geom_point(size = 2) +
    scale_color_igv() +
    theme_bw() +
    labs(title = 'Van Krevelen Diagram',
         x = 'O/C',
         y = 'H/C') +
    theme(plot.title = element_text(face = 'bold',
                                    hjust = 0.5)) +
    facet_grid(rows = vars({{facet_by}}),
               cols = vars({{facet_by2}}))
  
}

# -------------------------------------------------------------------------
plot_boxplot <- function(df, my_x, my_y, color_by, my_colors){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}},
             fill = {{color_by}})) +
    geom_boxplot() +
    scale_fill_manual(values = {{my_colors}}) +
    theme_bw() +
    theme(plot.title = element_text(face = 'bold', hjust = 0.5))
}