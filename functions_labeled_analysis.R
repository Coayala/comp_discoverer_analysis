#
#
# Christian Ayala
# Functions for Labeled_analysis.Rmd
#
#

# -------------------------------------------------------------------------
get_elements <- function(df){
  ## Obtain element list in the formula
  element_list <- gsub('[[:digit:]]+', '', df$Formula)
  element_list <- paste(element_list, collapse = ' ')
  element_list <- str_split(element_list, ' ')[[1]]
  element_list <- unique(element_list)
  
  return(element_list)
}

# -------------------------------------------------------------------------
### OLD NOT AUTOMATED VERSION OF THE separate_formula FUNCTION, IT WORKS ONLY FOR SOME ELEMENTS
# separate_formula <- function(df){
#   # This function will split the Formula column
#   # into columns with the number of each element
#   
#   ## Get formula in a new df
#   
#   new_df <- select(df, Formula)
#   
#   ## Separate Formula into the elements
#   new_df <- separate(new_df, 
#                  Formula,
#                  c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7'),
#                  sep = ' ',
#                  remove = FALSE)
#   
#   ## Get the indices of each element
#   new_df <- new_df %>% 
#     mutate(C = ifelse(grepl('C', S1), ifelse(grepl('C\\d\\d|C\\d', S1), str_replace(S1, 'C', ''), 1), 0)) %>% 
#     mutate(H = ifelse(grepl('H', S2), ifelse(grepl('H\\d\\d|H\\d', S2), str_replace(S2, 'H', ''), 1), 0)) %>% 
#     mutate(O = ifelse(grepl('O', S1), ifelse(grepl('O\\d\\d|O\\d', S1), str_replace(S1, 'O', ''), 1),
#                       ifelse(grepl('O', S2), ifelse(grepl('O\\d\\d|O\\d', S2), str_replace(S2, 'O', ''), 1),
#                              ifelse(grepl('O', S3), ifelse(grepl('O\\d\\d|O\\d', S3), str_replace(S3, 'O', ''), 1),
#                                     ifelse(grepl('O', S4), ifelse(grepl('O\\d\\d|O\\d', S4), str_replace(S4, 'O', ''), 1),
#                                            ifelse(grepl('O', S5), ifelse(grepl('O\\d\\d|O\\d', S5), str_replace(S5, 'O', ''), 1),0)))))) %>% 
#     mutate(N = ifelse(grepl('N', S1), ifelse(grepl('N\\d\\d|N\\d', S1), str_replace(S1, 'N', ''), 1),
#                       ifelse(grepl('N', S2), ifelse(grepl('N\\d\\d|N\\d', S2), str_replace(S2, 'N', ''), 1),
#                              ifelse(grepl('N', S3), ifelse(grepl('N\\d\\d|N\\d', S3), str_replace(S3, 'N', ''), 1),
#                                     ifelse(grepl('N', S4), ifelse(grepl('N\\d\\d|N\\d', S4), str_replace(S4, 'N', ''), 1),
#                                            ifelse(grepl('N', S5), ifelse(grepl('N\\d\\d|N\\d', S5), str_replace(S5, 'N', ''), 1),0)))))) %>% 
#     mutate(P = ifelse(grepl('P', S1), ifelse(grepl('P\\d\\d|P\\d', S1), str_replace(S1, 'P', ''), 1),
#                       ifelse(grepl('P', S2), ifelse(grepl('P\\d\\d|P\\d', S2), str_replace(S2, 'P', ''), 1),
#                              ifelse(grepl('P', S3), ifelse(grepl('P\\d\\d|P\\d', S3), str_replace(S3, 'P', ''), 1),
#                                     ifelse(grepl('P', S4), ifelse(grepl('P\\d\\d|P\\d', S4), str_replace(S4, 'P', ''), 1),
#                                            ifelse(grepl('P', S5), ifelse(grepl('P\\d\\d|P\\d', S5), str_replace(S5, 'P', ''), 1),0)))))) %>% 
#     mutate(S = ifelse(grepl('S', S1), ifelse(grepl('S\\d\\d|S\\d', S1), str_replace(S1, 'S', ''), 1),
#                       ifelse(grepl('S', S2), ifelse(grepl('S\\d\\d|S\\d', S2), str_replace(S2, 'S', ''), 1),
#                              ifelse(grepl('S', S3), ifelse(grepl('S\\d\\d|S\\d', S3), str_replace(S3, 'S', ''), 1),
#                                     ifelse(grepl('S', S4), ifelse(grepl('S\\d\\d|S\\d', S4), str_replace(S4, 'S', ''), 1),
#                                            ifelse(grepl('S', S5), ifelse(grepl('S\\d\\d|S\\d', S5), str_replace(S5, 'S', ''), 1),0)))))) %>% 
#     mutate(Cl = ifelse(grepl('Cl', S1), ifelse(grepl('Cl\\d\\d|Cl\\d', S1), str_replace(S1, 'Cl', ''), 1),
#                        ifelse(grepl('Cl', S2), ifelse(grepl('Cl\\d\\d|Cl\\d', S2), str_replace(S2, 'Cl', ''), 1),
#                               ifelse(grepl('Cl', S3), ifelse(grepl('Cl\\d\\d|Cl\\d', S3), str_replace(S3, 'Cl', ''), 1),
#                                      ifelse(grepl('Cl', S4), ifelse(grepl('Cl\\d\\d|Cl\\d', S4), str_replace(S4, 'Cl', ''), 1),
#                                             ifelse(grepl('Cl', S5), ifelse(grepl('Cl\\d\\d|Cl\\d', S5), str_replace(S5, 'Cl', ''), 1),0))))))
#   
#   ## Eliminate extra columns
#   new_df <- select(new_df, Formula, C, H, O, N, P, S, Cl)
#   
#   ## Bind new_df to original df
#   
#   df <- left_join(df, new_df, by = "Formula")
#   
#   return(df)
# }

# -------------------------------------------------------------------------
separate_formula <- function(df){
  # This function will split the Formula column
  # into columns with the number of each element
  
  ## Get the elements that make each of the compounds
  element_list <- get_elements(joint_table)
  
  ## Get formula in a df were results will be stored
  result_df <- select(joint_table, Formula)
  
  ## Separate Formula into the elements
  new_df <- separate(result_df, 
                     Formula,
                     into = {{element_list}},
                     sep = ' ')
  
  ## Initialize a temporal df and an accumulator for the for loops
  temp_df <- tibble(.rows = nrow(new_df))
  j <- 1
  
  ## For loop to obtain the coefficients of each element
  for(el in element_list){
    my_exp <- paste0(el, "\\d\\d|", el, "\\d")
    for(i in 1:length(element_list)){
      temp_df[,i] <- ifelse(grepl(el, deframe(new_df[,i])), # because new_df is a tibble, deframe allows it to be sliced as a vector
                            ifelse(grepl(my_exp, deframe(new_df[,i])), str_replace(deframe(new_df[,i]), paste0(el), ''), 1), 0)
      temp_df[,i] <- as.numeric(unlist(temp_df[,i])) # unlist produces atomic components that can be changed into numeric values
    }
    temp_df$sum <- rowSums(temp_df)
    result_df[,j + 1] <- temp_df$sum
    j <- j + 1
    temp_df <- tibble(.rows = nrow(new_df))
  }
  
  ## Put back the names of the elements in each column
  colnames(result_df) <- c('Formula', element_list)
  
  ## Merge with original matrix
  df <- left_join(df, result_df, by = "Formula")
  
  
}

# -------------------------------------------------------------------------
calc_ratios_n_idxs <- function(df){
  # This function will calculate H/c and O/C ratios
  # as well as other thermodynamics index
  
  df$C <- as.numeric(df$C)
  df$H <- as.numeric(df$H)
  df$O <- as.numeric(df$O)
  df$N <- as.numeric(df$N)
  df$P <- as.numeric(df$P)
  df$S <- as.numeric(df$S)
  
  ## Get ratios
  df <- df %>% 
    mutate(H_to_C = H / C) %>% 
    mutate(O_to_C = O / C)
  
  ## Calculate thermodynamic indices
  df <- df %>% 
    mutate(NOSC = -((4*C + H - 3*N - 2* O + 5*P - 2*S) / C) + 4) %>% 
    mutate(GFE = 60.3 - 28.5*NOSC) %>% 
    mutate(DBE = 1 + 0.5 * (2*C - H + N + P)) %>% 
    mutate(DBE_O = DBE - O) %>% 
    mutate(AI = (1 + C - O - S - ((H + P + N) * 0.5)) / (C - O - S - N - P)) %>% 
    mutate(AI_mod = (1 + C - (O * 0.5) - S - ((H + P + N) * 0.5)) / (C - (O * 0.5) - S - N - P)) %>% 
    mutate(DBE_AI = 1 + C -O -S - (0.5 * (H + N + P)))
}

# -------------------------------------------------------------------------
calc_classes <- function(df){
  df <- df %>% 
    mutate(Class = ifelse(between(O_to_C, 0, 0.3)&between(H_to_C, 1.5, 2.5), 'Lipid',
                          ifelse(between(O_to_C, 0, 0.125)&between(H_to_C, 0.8, 1.5), 'Unsaturated HC',
                                 ifelse(between(O_to_C, 0, 0.95)&between(H_to_C, 0.2, 0.8), 'Condensed HC',
                                        ifelse(between(O_to_C, 0.3, 0.55)&between(H_to_C, 1.5, 2.3), 'Protein',
                                               ifelse(between(O_to_C, 0.55, 0.7)&between(H_to_C, 1.5, 2.2), 'Amino Sugar',
                                                      ifelse(between(O_to_C, 0.7, 1.5)&between(H_to_C, 1.5, 2.5), 'Carbohydrate',
                                                             ifelse(between(O_to_C, 0.125, 0.65)&between(H_to_C, 0.8, 1.5), 'Lignin',
                                                                    ifelse(between(O_to_C, 0.65, 1.1)&between(H_to_C, 0.8, 1.5), 'Tannin', 'Other')))))))))
  return(df)
}

# -------------------------------------------------------------------------
plot_vank <- function(df, color_by, facet_by, facet_by2 = NULL){
  ggplot(df,
         aes(x = O_to_C,
             y = H_to_C,
             color = {{color_by}})) +
    geom_point(size = 2) +
    scale_color_jco() +
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
plot_col <- function(df, my_x, my_y, color_by1, color_by2){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}})) +
    geom_col(aes(fill = {{color_by1}},
                 color = {{color_by2}}),
             size = 2,
             width = 0.75) +
    scale_fill_jco() +
    scale_color_jama() +
    theme_bw()
}

# -------------------------------------------------------------------------
plot_boxplot <- function(df, my_x, my_y, my_comparisons = NULL){
  ggplot(df,
         aes(x = {{my_x}},
             y = {{my_y}},
             fill = {{my_x}})) +
    geom_boxplot() +
    scale_fill_jama() +
    stat_compare_means(comparisons = {{my_comparisons}},
                       method = 't.test',
                       label = 'p.signif') +
    theme_bw()
}

# -------------------------------------------------------------------------
plot_venn <- function(my_list, my_colors){
  venn(my_list,
       zcolor = {{my_colors}},
       ilcs = 1,
       sncs = 1)
}
