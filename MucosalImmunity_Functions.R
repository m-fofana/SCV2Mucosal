# Functions to generate and plot bootstrapped correlogram 

# mycorrboot ------------------------------------------------------------------

mycorboot <- function(
    data,
    method = "pearson",
    cor_matrix = NULL,
    p_matrix = NULL,
    diag = T,
    numrep = 999,
    boot = T,
    varorder = NULL,
    cimethod = "basic"){
  # Function to run and compile correlation coefficients
  # And prepare dataframe suitable for plotting function  
  
  # -- required packages -------------------------------------------------------
  
  require(reshape2, quietly = TRUE)
  require(tibble, quietly = TRUE)
  require(rstatix, quietly = TRUE)
  require(TOSTER, quietly = TRUE)
  
  # -- check data columns ------------------------------------------------------
  
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      data = as.data.frame(data)
    }
    x = which(!sapply(data, is.numeric))
    if (length(x) > 0) {
      warning(paste("data in column(s)",
                    paste0(paste0("'", names(data)[x], "'"), collapse = ", "),
                    "are not numeric and were ignored"))
      data = data[, -x ]
    }
  }
  
  if (is.null(varorder)) {
    varorder = colnames(data)} else if (!is.null(varorder)) {
      data = select(data, varorder, everything())
    }
  
  # -- correlation matrix and p-values-----------------------------------------
  
  if (is.null(cor_matrix)) {
    cor_matrix = cor_mat(data, method = method) 
  }
  
  mc = cor_matrix %>% column_to_rownames()
  mp = cor_matrix %>% cor_get_pval() %>% column_to_rownames()
  
  # -- correlation data.frame --------------------------------------------------
  
  mc$row_names = colnames(mc)
  mc = reshape2::melt(mc, id.vars = "row_names") %>%
    rename(coeff = value)
  
  mp$row_names = colnames(mp)
  mp = reshape2::melt(mp, id.vars = "row_names") %>%
    rename(pval = value)
  
  mt = merge(mc, mp, by = c("row_names", "variable")) %>%
    mutate(row_names = as.character(row_names),
           variable = as.character(variable)) %>%
    mutate(row_names = factor(row_names, levels = varorder),
           variable = factor(variable, levels = varorder)) %>%
    mutate(x1 = as.numeric(row_names),
           y = as.numeric(variable)) 
  
  # Bootstrap computations -----------------------------------------
  resdat <- expand.grid(row_names = colnames(data), variable = colnames(data),
                        stringsAsFactors = F, KEEP.OUT.ATTRS = F) %>%
    mutate(coeff = NA, pval = NA, lci = NA, uci = NA)
  
  for(i in 1:nrow(resdat)){
    dattmp = data.frame(x = c(select(data, resdat$row_names[i])), 
                        y = c(select(data, resdat$variable[i]))) %>%
      filter(complete.cases(.) == T) 
    
    colnames(dattmp) = c("x", "y")
    
    # Run bootstrapped correlation
    restmp = boot_cor_test(x = dattmp$x, y = dattmp$y, boot_ci = cimethod,
                           method = "pearson", alternative = "two.sided",
                           R = numrep, null = 0)
    
    resdat$coeff[i] = restmp$estimate
    resdat$pval[i] = restmp$p.value
    resdat$lci[i] = restmp$conf.int[1]
    resdat$uci[i] = restmp$conf.int[2]
  }
  
  # Combine together -------------------------------------------------
  datfinal = merge(mt, resdat, by = c("row_names", "variable"),
                   suffixes = c("", "_boot")) %>%
    # Keep lower triangle
    subset(y <= x1) %>%
    # Remove nonsignificant and missing values
    mutate(across(.cols = c(pval, coeff),
                  .fns = ~replace(., coeff == 0 | is.na(coeff) | pval >= 0.05, 
                                  NA))) %>%
    mutate(across(.cols = c(pval_boot, coeff_boot),
                  .fns = ~replace(., coeff_boot == 0 | is.na(coeff_boot) | pval_boot >= 0.05, 
                                  NA)))
  
  return(datfinal)
}




# mycorrplot ------------------------------------------------------------------

mycorrplot <- function(
    mt,
    nbreaks = NULL,
    digits = 2,
    name = "", # Name of color legend
    high = "#67001F",
    mid = "#EEEEEE",
    low = "#053061",
    midpoint = 0,
    size_legend = F, # "legend" to include size legend
    size_col = "gray75", # color of size legend symbol
    glines = T,
    min_size = 1,
    max_size = 7,
    varorder = NULL,
    colseq = NULL, # color of text labels
    xlab = T,
    ylab = T,
    diaglab = F,
    labsize = 12, # axis label size
    diagsize = 12, # diagonal label size
    gridseq = c(6, 12), # where to draw thick lines
    flip_x = T,
    legend.position = "right",
    legend.size = 9, 
    siggbox = T, ...){
  
  # Customized function to generate correlogram with grid of circles
  # Pairwise complete, Pearson correlation
  # correlation coefficient represented by shade
  # p-value represented by size of circles (proportional to log(1/p))
  # Modified from package ggcor (https://github.com/briatte/ggcorr/blob/master/ggcorr.R)
  
  # -- required packages -------------------------------------------------------
  
  require(ggplot2, quietly = TRUE)
  
  # -- check data columns ------------------------------------------------------
  
  if (!is.data.frame(mt)) {data = as.data.frame(mt)}
  
  if(is.null(varorder)){varorder = unique(mt$row_names)}
  
  if(is.factor(mt$row_names)){mt$row_names = as.character(mt$row_names)} 
  if(is.factor(mt$variable)){mt$variable = as.character(mt$variable)} 
  
  # -- edit correlation matrix -----------------------------------------
  mt = mt %>%
    mutate(row_names = factor(row_names, levels = varorder),
           variable = factor(variable, levels = varorder)) 
  
  # Flip x axis
  mt$x2 = max(mt$x) - mt$x + 1
  mt  = mt %>%
    mutate(boxcol = "No")
  
  if(siggbox == T){
    mt = mt %>%
      mutate(boxcol = replace(boxcol,
                              grepl("Stotal.IgG", row_names) &
                                grepl("Stotal.IgG", variable) &
                                variable != row_names,
                              "Yes"))
  }
  
  mt = mt %>% mutate(boxcol = factor(boxcol))
  
  if(flip_x == T) {mt$x = mt$x2} else {mt$x = mt$x1}
  
  # Set point size
  # proportional to log10(1-p), with upper limit at p<1e-10, scale to max = 1
  # minimum size is log10(1/0.05) / 10 = 0.13 (only significant results shown)
  ptransf = function(x){min(log10(1/x), 10) / 10}
  mt$size = sapply(mt$pval, ptransf) 
  
  mt = mt %>%
    mutate(var_lab = str_replace_all(as.character(variable), 
                                     c("[.]" = " ", "Stotal" = "S total")),
           row_lab = str_replace_all(as.character(row_names), 
                                     c("[.]" = " ", "Stotal" = "S total")))
  
  # -- plot structure ----------------------------------------------------------
  
  p = ggplot(na.omit(mt, coeff), aes(x, y)) +
    # Background grid
    geom_tile(data = mt, value = 1, col = "gray75", fill = NA, linewidth = 1) +
    geom_tile(data = mt, aes(linetype = boxcol), col = "black",
              value = 1, fill = NA, linewidth = 1) +
    scale_linetype_manual("", breaks = c("No", "Yes"), values = c(0, 1), 
                          guide = "none") +
    # Draw points
    geom_point(aes(size = size, color = coeff)) +
    # Size scale   
    scale_size_continuous("p-value",
                          range = c(min_size, max_size),
                          limits=c(0, 1),
                          labels = c(0.05, 0.01, 0.001, "0.0001"),
                          breaks = sapply(c(0.05, 0.01, 0.001, 0.0001), ptransf)) +
    
    guides(size = guide_legend(override.aes = list(col = size_col))) +
    guides(size = size_legend) +
    # Color scale
    scale_color_gradient2(name, low = low, mid = mid, high = high,
                          midpoint = midpoint, limits = c(-1, 1))
  
  # -- Square grid and add labels ------------------------------------
  
  # Axis labels
  if(xlab == T) {
    axis.text.x = element_text(color = colseq[as.numeric(factor(mt$row_names))],
                               angle = 90, vjust = 0.5, hjust = 1, size = labsize)
  } else {axis.text.x = element_blank()}
  
  if(ylab == T) {
    axis.text.y = element_text(color = colseq[as.numeric(factor(mt$variable))],
                               size = labsize)
  } else {axis.text.y = element_blank()}
  
  if(is.null(colseq)){
    colseq = rep("gray25", ncol(data))
  }
  
  dnudge = diaglab * 5
  
  p = p  +
    coord_equal() +
    # Allow for some extra space for diagonal labels
    scale_x_continuous("", labels = mt$row_lab, 
                       breaks = mt$x,
                       lim = c((0 - (1-flip_x)*dnudge), 1 + max(mt$x) + flip_x*dnudge)) +
    scale_y_continuous("", labels = mt$var_lab, 
                       breaks = mt$y,
                       lim = c(0, 1 + max(mt$y) + dnudge)) +
    theme(
      axis.text.x = axis.text.x,
      axis.text.y = axis.text.y,
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      legend.key = element_blank(),
      legend.position = legend.position,
      legend.title = element_text(size = legend.size),
      legend.text = element_text(size = legend.size)
    )
  
  # Diagonal labels
  mtdiag <- subset(mt, row_names == variable) %>%
    mutate(xgrid = x + 0.5, ygrid = y + 0.5) 
  
  hjust = ifelse(flip_x == T, 0, 1)
  angle = ifelse(flip_x == T, 45, -45)
  
  if(diaglab == T) {
    p = p +
      new_scale_color() +
      geom_text(data = mtdiag, 
                aes(label = var_lab, x = x, y = ygrid, col = variable), 
                hjust = hjust, vjust = 0, angle = angle, nudge_y = 0.2, 
                size = diagsize/.pt) +
      scale_color_manual("", breaks = mtdiag$variable, 
                         values = colseq[as.numeric(factor(mtdiag$variable))],
                         guide = NULL)
  } 
  
  # Grid lines
  if(glines == T) {
    tmpgrid = data.frame(tmp = gridseq, 
                         orig1 = 0.5, 
                         end1 = NA,
                         orig2 = gridseq + 0.5,
                         end2 = gridseq + 0.5)
    
    for(i in 1:length(gridseq)){
      tmpgrid$end1[i] <- mtdiag$ygrid[mtdiag$x == tmpgrid$tmp[i]]
    }
    
    if (flip_x == T) {
      p = p +
        geom_segment(data = tmpgrid, aes(x = orig1, xend = end1,
                                         y = orig2, yend = end2)) +
        geom_segment(data = tmpgrid, aes(x = orig2, xend = end2,
                                         y = orig1, yend = end1))
    } else if (flip_x == F) {
      p = p +
        geom_segment(data = tmpgrid, aes(xend = nrow(mtdiag)+0.5, x = end1,
                                         y = orig2, yend = end2)) +
        geom_segment(data = tmpgrid, aes(x = orig2, xend = end2,
                                         y = orig1, yend = end1))
    }
  }
  
  return(p)
}