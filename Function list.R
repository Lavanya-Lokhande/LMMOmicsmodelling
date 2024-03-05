GGscatterPlot <- function(data, mapping, ..., 
                          method = "spearman") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method)
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D , 
                                       df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    ggplot2::scale_color_viridis_c(direction = colDirection) +
    #                scale_color_gradient(low = "#0091ff", high = "#f0650e") +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel, 
                             y = ylabel, 
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_minimal()
  
  return(pp)
}

univ_cox_ph <- function(xx){
  aa <- summary(xx)
  tmp1 <- aa$conf.int
  tmp2 <- aa$coefficients
  ntot <- aa$n
  hr = tmp1[1,1]
  hr_low <- tmp1[1,3]
  hr_high <- tmp1[1,4]
  coeff <- tmp2[1,1]
  pval <- tmp2[1,5] ## for hazard ratio
  cc <- cox.zph(xx)
  tmp <- cc$table
  varnames <- rownames(tmp);
  p_check <- tmp[1,3]
  vname <- varnames[1]
  obj <- data.frame(cbind(vname,ntot,coeff, hr, hr_low, hr_high, pval, p_check))
  return(obj)
}

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

volcano_plots_binary <- function(data, subset = FALSE,
                                 ROI_type,
                                 formula,
                                 ctrl_biom=FALSE,
                                 biom_trans=NULL,
                                 p_sig, coeff_thresh, reflevel, contrastlevel, subset_cat, binary_group, name){
  data <- data
  
  if (subset == TRUE) {data <- dplyr::filter(data, data[[subset_cat]] == ROI_type)}
  
  data <- data[!data[[binary_group]]=="",]
  data[[binary_group]] <- relevel(factor(data[[binary_group]]), reflevel)
  
  biomarkers <- names(data[,s:e])
  
  ## removes the biomarkers 
  if (ctrl_biom==FALSE){biomarkers <- biomarkers[!biomarkers %in% control]}
  
  
  ## log transformation
  if (!is.null(biom_trans)){data[,names(data) %in% biomarkers] <- biom_trans(data[,names(data) %in% biomarkers])}
  
  data <- data[complete.cases(data[,binary_group]),]
  
  p_marker <- data.frame(biomarker=biomarkers)
  p_marker$coeff = p_marker$pvalue <- NA
  
  ### CHANGING ACORDING TO BEST MODEL
  for (biomarker in biomarkers){
    form <- as.formula(paste(biomarker,formula))
    one_var_regress <- lmerTest::lmer(form, data=data)
    
    cf <- coefficients(summary(one_var_regress))
    p_marker[p_marker$biomarker==biomarker,"coeff"] = cf[2,1]
    p_marker[p_marker$biomarker==biomarker,"pvalue"] = cf[2,5]
  }
  
  p_marker$FDR <- p.adjust(p_marker$pvalue,"fdr")
  p_marker$significance <- -log10(p_marker$FDR)
  p_marker <- p_marker[!is.na(p_marker$coeff),]
  p_marker <-  p_marker %>% dplyr::mutate(DE = cut(FDR, breaks=c(0, p_sig, Inf), labels=c("DE","NA"))) 
  p_marker$threshold <- ifelse(p_marker$coeff > coeff_thresh & p_marker$FDR < p_sig, paste0("Enriched in ", contrastlevel), ifelse(p_marker$coeff < -(coeff_thresh) & p_marker$FDR < p_sig, paste0("Enriched in ", reflevel), ifelse(p_marker$FDR < p_sig,  "DE", "ns")))
  p_marker <-p_marker[order(p_marker$FDR),]
  
  ## YOU CAN FIND THE ENTIRE LIST WITH THE COEFF, PVAL, FDR etc IN YOUR GLOBAL ENVIRONMENT
  assign(paste0("DR_proteins_for_", name), p_marker, envir = .GlobalEnv)
  
  plot <- ggplot(p_marker, aes(x=coeff, y=significance, label=biomarker, color = threshold, size = significance)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values=cols) +
    geom_hline(yintercept = -log10(p_sig), linetype="longdash", col="black") +
    geom_vline(xintercept = c(-coeff_thresh, coeff_thresh), linetype="dashed", col="black") +
    geom_text_repel(data = subset(p_marker, FDR<p_sig), size = 3, show.legend = F,
                    point.padding = .3, box.padding = .3,
                    min.segment.length = .2, color="black",
                    max.overlaps = Inf) +
    labs(title=paste0(ROI_type," ROI -", " Significant change in ", contrastlevel , " with ", reflevel," as reference and ",'-log[10] ', p_sig," as significance line for ", name),  x="Coefficient", y=expression('-log'[10]*' (FDR adjusted P-value)')) +
    theme(legend.position = "top",legend.text = element_text(size = 9), panel.border = element_rect(color = "black",fill = NA, size = 0.2), plot.title = element_text(size=10), panel.background = element_rect(fill = "white"))
  
  print_df <- p_marker[,1:5] %>% dplyr::filter(FDR<=p_sig) 
  
  colnames(print_df)[1] <- c(paste0(contrastlevel , " against ", reflevel))
  
  print(plot)
  #print(print_df%>% as.data.frame()%>% kbl(digits = 2, format = "html",  align=c(rep('c',times=8))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir"))
  
  exp_data <- data[s:e] %>% t() %>% as.data.frame()
  
  p_marker_subset <- p_marker %>% dplyr::filter(FDR<=p_sig) %>% as.data.frame()
  # probelist <- dplyr::select(p_marker_subset, "biomarker") %>% unlist()
  # assign("test1", probelist, envir = .GlobalEnv)
  rownames(p_marker_subset) <- p_marker_subset$biomarker
  # 
  # 
  if (nrow(p_marker_subset)>1) {
    heatmap_input <- exp_data %>% rownames_to_column(., var = "biomarker")
    heatmap_input <- left_join(p_marker_subset, heatmap_input, by = "biomarker")
    heatmap_input <- heatmap_input %>% column_to_rownames(var = "biomarker")
   
    
    #heatmap_input2 <-  data[unlist(probelist[1])]
    #assign("test2", heatmap_input2, envir = .GlobalEnv)
    
    n1 <- nrow(data)
    e1 <- ncol(heatmap_input)
    s1 <- e1-n1+1
    
    
    
    factor_rows <-  data.frame("Direction" = heatmap_input$threshold, "FDR" = heatmap_input$FDR)
    rownames(factor_rows) <- rownames(heatmap_input) # name matching
    factor_cols <- data[annot_file]
   #  #heatmap_input <- apply(heatmap_input[,s1:e1], 1, scale)
   #  
   # # print(as.data.frame(heatmap_input[,s1:e1]))
   #  
   #  mat_scaled = t(scale(t(heatmap_input[,s1:e1])))
   #  
   #  #heatmap_input2 <- cbind(heatmap_input[1:s1-1], apply(heatmap_input[s1:e1],1, cal_z_score))
   #  #  
   #  # #  
   #  # #  plot2 <-pheatmap(heatmap_input[,s1:e1], annotation_col = factor_cols, annotation_row = factor_rows, cluster_cols = T, show_colnames = F, fontsize = 5, cluster_rows = T, main = paste0(ROI_type," ROI -", "Significant change in ", contrastlevel , " with ", reflevel," as reference and ",'-log[10] ', p_sig," as significance line"),  x="coefficient", y=expression('-log'[10]*' (FDR adjusted P-value)'))
   #  # # # print(plot2)
   #  
   #  #pdf(qq(paste0("Heatmap_", ROI_type," ROI -", "Significant change in ", contrastlevel , " with ", reflevel," as reference and ",'-log[10] ', p_sig," as significance line for ", name)), width = 8, height = 8)
   #  ha <- HeatmapAnnotation(df = factor_cols, show_legend=FALSE, annotation_height = unit(5, "mm"))
   #  pa <- rowAnnotation(df = factor_rows, show_legend=FALSE, annotation_width = unit(0.7, "cm"))
   #  
   #  plot2 <- Heatmap(as.matrix(mat_scaled), row_names_gp = gpar(fontsize = 6), cluster_column_slices = TRUE, row_gap = unit(1, "mm"), border = TRUE, row_title_rot = 0, column_names_gp = gpar(fontsize = 0), column_split = data[binary_group], col = pal2(20), top_annotation = ha,  column_title = paste0(ROI_type," ROI -", "Significant change in ", contrastlevel , " with ", reflevel," as reference and ",'-log[10] ', p_sig," as significance line"),left_annotation = pa, show_heatmap_legend = FALSE, cluster_rows = TRUE)
   #  #draw(plot2, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
   #  #dev.off()
   #  #print(plot2)
   #  
   #  
   #  #plot3_input <- cbind(data[binary_group],t(heatmap_input[, 7:ncol(heatmap_input)]) ) %>% raster::as.data.frame() %>% melt(id.vars = binary_group)
   # 
   #  plot3 <- ggplot(plot3_input, aes(x=plot3_input[,1], y=value, fill=plot3_input[,1])) +
   #    geom_boxplot(outlier.colour="red",outlier.size = 0.3) +
   #    geom_jitter(color="black", size=0.3, alpha=0.5) +
   #    ylab("") +
   #    xlab("") +
   #    #scale_fill_manual(values=c("skyblue", "red")) +
   #    facet_wrap( ~ variable, scales="free") +
   #    labs(title = paste0(ROI_type," ROI -", "Significant change in ", contrastlevel , " with ", reflevel," as reference and ",'-log[10] with ', p_sig," as significance line for ", name))+
   #    theme(panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA),axis.title.x=element_blank(),
   #          axis.text.x=element_text(size = 8),
   #          axis.ticks.x=element_blank(), axis.text.y=element_text(size = 8), plot.title = element_text(size=10), legend.title = element_blank())
   #  #print(plot3)


    exp_data <- t(data[,s:e])

    gcor <- cor(t(heatmap_input[,s1:e1]), method="spearman")
    plot4 <- corrplot.mixed(gcor,outline = F, order="hclust",addrect = 2, rect.col = "black", rect.lwd = 1, tl.col = "black", tl.cex = 0.5, cl.cex = 0.5, sig.level = 0.05, tl.pos = 'lt', lower = 'color',upper = 'number', number.cex = 0.3, upper.col = "black")
    print(plot4)
  }


}


Compare_models <- function(baseModel, coxphModel,coxmeModel) {
  if (!class(coxphModel) == "coxph" | !class(coxmeModel) == "coxme") {
    stop("Wrong models")
  }
  
  anova_compare <- anova(baseModel, coxphModel, coxmeModel)
  print(anova_compare)
  
  ## Degrees of freedom
  coxphDf <- sum(!is.na(coef(coxphModel)))
  coxmeDf <- coxmeModel$df
  names(coxmeDf) <- c("Integrated","Penalized")
  
  ## DF differnces
  dfDiff <- coxmeDf - coxphDf
  ## Log likelihoods
  coxphLogLik <- coxphModel$loglik[2]  ## Log likelihood
  coxmeLogLik <- coxmeModel$loglik + c(0, 0, coxmeModel$penalty)
  ## Chi-squared value comparing the penalized model to the null model
  
  Loglikcompareph <-  -2 * (baseModel$loglik - logLik(coxphModel))
  Loglikcompareme <-  -2 * (coxmeLogLik["NULL"] - coxmeLogLik["Penalized"])
  
  outdf0 <- data.frame(c(Loglikcompareph,Loglikcompareme))
  rownames(outdf0) <- c("coxph to NULL model", "coxme to NULL model")
  colnames(outdf0) <- "Chi square value difference" 
  print(outdf0)
  
  ##______TESTING FOR RANDOM EFFECTS_____##
  ## -2 logLik difference
  logLikNeg2Diff <- c(-2 * (coxphLogLik - coxmeLogLik["Integrated"]),
                      -2 * (coxphLogLik - coxmeLogLik["Penalized"]))
  
  ## p-values
  pVals <- pchisq(q = logLikNeg2Diff, df = dfDiff, lower.tail = FALSE)
  
  ## Combine
  outDf <- data.frame(dfDiff, logLikNeg2Diff, pVals)
  colnames(outDf) <- c("Degrees of freedom difference", "-2logLik diff", "p-values")
  print(outDf)
  #formattable::formattable(anova_compare, align =c("l","c","c","c", "c","c","c", "c"))
}


univ_cox_me <- function(xx){
  beta <- xx$coefficients #$fixed is not needed
  hr <- exp(beta)
  z1 <- qnorm((1 + 0.95)/2, 0, 1)
  nvar <- length(beta)
  nfrail <- nrow(xx$var) - nvar
  se <- sqrt(diag(xx$var)[nfrail + 1:nvar])
  z <- round(beta/se, 2) 
  hr_low <- exp(beta - z1 * se)
  hr_high <- exp(beta + z1 * se)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  cc <- cox.zph(xx)
  tmp <- cc$table
  check_pval <- tmp[1,3]
  obj <- data.frame(cbind(nvar,nfrail,beta,hr,hr_low,hr_high,se,z,p,check_pval))
  return(obj)
}

univ_cox_ph <- function(xx){
  aa <- summary(xx)
  tmp1 <- aa$conf.int
  tmp2 <- aa$coefficients
  ntot <- aa$n
  hr = tmp1[1,1]
  hr_low <- tmp1[1,3]
  hr_high <- tmp1[1,4]
  coeff <- tmp2[1,1]
  pval <- tmp2[1,5] ## for hazard ratio
  cc <- cox.zph(xx)
  tmp <- cc$table
  varnames <- rownames(tmp);
  p_check <- tmp[1,3]
  vname <- varnames[1]
  obj <- data.frame(cbind(vname,ntot,coeff, hr, hr_low, hr_high, pval, p_check))
  return(obj)
}

distribution_plots <- function(dataset, x_title, type){
  plot1 <- dataset %>% 
    select(combined_unique_coreID, Segment_name_label, Within_the_tumor, cell_number, Nuclei_count) %>%  
    group_by(combined_unique_coreID, Segment_name_label, Within_the_tumor) %>% 
    summarise(mean_FCI_NC = mean(cell_number, na.rm=TRUE), mean_O_NC = mean(Nuclei_count, na.rm=TRUE)) %>%  na.omit() %>% 
    
    ggplot(aes(x =mean_FCI_NC, y=mean_O_NC, fill = Segment_name_label, colour = Segment_name_label)) + 
    facet_wrap(~Within_the_tumor, scales = "free")+
    geom_point( alpha = 0.5, shape=1) +
    theme(axis.text.x = element_text(size = 8, angle=0, vjust=0.5), axis.ticks.length.x = unit(0.1, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7)) +
    xlab(paste0(x_title)) +
    ylab("Original GeoMX DSP nuclei count") +
    geom_rug(aes(color =Segment_name_label)) +
    geom_smooth(aes(color = Segment_name_label, fill = Segment_name_label), method = "lm", fullrange = FALSE, alpha = 0.4, formula = y~x) +
    stat_poly_eq(formula = y~x, aes(label = paste(..eq.label..)), parse = TRUE)  +
    scale_color_manual(values=brewer.pal(n=4, name="Set2")) +
    scale_fill_manual(values=brewer.pal(n=4, name="Set2")) 
  print(plot1)
  
  plot2 <- dataset %>% select(combined_unique_coreID, Segment_name_label, Within_the_tumor, cell_number, Nuclei_count) %>%  
    group_by(combined_unique_coreID, Segment_name_label, Within_the_tumor) %>% 
    summarise(typer = mean(cell_number, na.rm=TRUE), "Geomx-DSP" = mean(Nuclei_count, na.rm=TRUE)) %>%  rename(c(typer = type)) %>% na.omit() %>% data.table::melt() %>%
    ggplot(aes(x = reorder(combined_unique_coreID, value), y = value, fill = Segment_name_label, colour = Segment_name_label)) + 
    facet_grid(variable~Within_the_tumor, scales = "free")+
    geom_bar(aes(), stat="identity", position="fill", colour = "black", size = 0) +
    # scale_fill_manual(values = c("darkgoldenrod1", "darkorchid", "red", "green", "darkblue")) +
    scale_fill_manual(values = colours)+
    guides(fill=guide_legend(nrow=1)) +
    xlab("Core based ID") +
    ylab("")+
    theme(axis.text.x = element_text(size = 0), axis.ticks.length.x = unit(0.0, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7)) +
    scale_y_continuous(labels = scales::percent)
 print(plot2)
}

cell_size_distribution <- function(dataset, cell_size, type, return_paired){
  
  dataset <- dataset %>% mutate(cell_size = cell_size)
  
  dataset_paired <- inner_join(dataset %>% dplyr::filter(Within_the_tumor =="Tumor-deficient zone") %>% dplyr::select(PATID, CoreID, Within_the_tumor, Segment_name_label, cell_size), dataset %>% dplyr::filter(Within_the_tumor =="Tumor-rich zone")  %>% dplyr::select(PATID, CoreID, Within_the_tumor, Segment_name_label, cell_size), by = c("PATID", "Segment_name_label"))%>% dplyr::rename(c("Tumor-deficient zone"="cell_size.x", "Tumor-rich zone"=cell_size.y)) %>% dplyr::select(-c(CoreID.x, CoreID.y, Within_the_tumor.x, Within_the_tumor.y))%>% as.data.table
  
  assign(paste0(return_paired), dataset_paired, envir = .GlobalEnv) 
  
  X <- aggregate(dataset_paired$`Tumor-deficient zone`,by=list(name=dataset_paired$PATID,dataset_paired$Segment_name_label),data=dataset_paired,FUN=mean) 
  Y <- aggregate(dataset_paired$`Tumor-rich zone`,by=list(name=dataset_paired$PATID,dataset_paired$Segment_name_label),data=dataset_paired,FUN=mean) 
  
  X <- X %>% as.data.frame %>% dplyr::rename(c("Tumor-deficient zone"=x, "Segment_name_label"=Group.2,  "PATID"=name))
  Y <- Y %>% as.data.frame %>% dplyr::select(x) %>% dplyr::rename(c("Tumor-rich zone"=x)) 
  final <- cbind(X,Y)
  
  cbind(final %>% dplyr::select(Segment_name_label) %>% group_by(., Segment_name_label) %>% count() %>% `colnames<-`(c("Cell Types", "Count")), final %>% dplyr::select(Segment_name_label,"Tumor-deficient zone") %>%  as.data.frame() %>% group_by(Segment_name_label)  %>% drop_na(`Tumor-deficient zone`)%>% dplyr::summarise(Mean_cell_size = mean(`Tumor-deficient zone`), sd_cell_size = sd(`Tumor-deficient zone`)) %>% mutate(`Avg cell size in TD` = paste0(round(Mean_cell_size,2),"(",round(sd_cell_size,2), ")")) %>% dplyr::select(`Avg cell size in TD`), final %>% dplyr::select(Segment_name_label,"Tumor-rich zone") %>%  as.data.frame() %>% group_by(Segment_name_label) %>% drop_na(`Tumor-rich zone`)%>% dplyr::summarise(Mean_cell_size = mean(`Tumor-rich zone`),sd_cell_size = sd(`Tumor-rich zone`)) %>% mutate(`Avg cell size in TR` = paste0(round(Mean_cell_size,2),"(",round(sd_cell_size,2), ")")) %>% dplyr::select(`Avg cell size in TR`)) %>% kbl(digits = 2, format = "html", row.names = TRUE, align=c(rep('c',times=4))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir") %>% print()
  
  test <- cbind(final %>% dplyr::select(Segment_name_label,"Tumor-deficient zone") %>%  as.data.frame() %>% group_by(Segment_name_label) %>% drop_na(`Tumor-deficient zone`) %>% dplyr::summarise(Mean_cell_size = mean(`Tumor-deficient zone`), sd_cell_size = sd(`Tumor-deficient zone`)))
  
  cbind(final %>% dplyr::select(Segment_name_label) %>% group_by(., Segment_name_label) %>% count() %>% `colnames<-`(c("Cell Types", "Count")), final %>% dplyr::select(Segment_name_label,"Tumor-deficient zone") %>%  as.data.frame() %>% group_by(Segment_name_label)%>% drop_na(`Tumor-deficient zone`) %>% dplyr::summarise(Median_cell_size = median(`Tumor-deficient zone`), sd_cell_size = sd(`Tumor-deficient zone`)) %>% mutate(`Median cell size in TD` = paste0(round(Median_cell_size,2),"(",round(sd_cell_size,2), ")")) %>% dplyr::select(`Median cell size in TD`), final %>% dplyr::select(Segment_name_label,"Tumor-rich zone") %>%  as.data.frame() %>% group_by(Segment_name_label)%>% drop_na(`Tumor-rich zone`) %>% dplyr::summarise(Median_cell_size = median(`Tumor-rich zone`),sd_cell_size = sd(`Tumor-rich zone`)) %>% mutate(`Median cell size in TR` = paste0(round(Median_cell_size,2),"(",round(sd_cell_size,2), ")")) %>% dplyr::select(`Median cell size in TR`)) %>% kbl(digits = 2, format = "html", row.names = TRUE, align=c(rep('c',times=4))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir")%>% print()
  
  dataset_paired2 <- final %>% reshape::melt(id.vars = c("PATID", "Segment_name_label"))
  
  plot <- ggpaired(dataset_paired2, x = "variable", y = "value",
           color = "variable", line.color = "gray", line.size = 0.3, palette = wes_palette("Darjeeling1", n = 2))+
    ggtitle(paste0(type)) +
    stat_compare_means(paired = TRUE, aes(group = variable), method = "t.test", size = 3, label.x = 1.3, hide.ns = F, label.y.npc = 0.8 )+
    stat_compare_means(paired = TRUE, aes(group = variable), method = "wilcox.test", size = 3, label.x = 1.3, hide.ns = F, label.y.npc = 0.9 )+facet_wrap(~Segment_name_label,scale="free") +theme(axis.text.x = element_text(size = 8, angle=0, vjust=0.5), axis.ticks.length.x = unit(0.1, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "none") +ylab(expression("Cell size (\U003BCm"^2*")")) + xlab("")
  
  print(plot)
  }

cell_density_plots <- function(dataset, cell_density, type){
  dataset<- dataset %>%mutate(cell_density = cell_density)
  final <- dataset %>% as.data.frame(.)
  
  cbind(final %>% dplyr::select(Segment_name_label,Within_the_tumor, cell_density) %>% dplyr::filter(Within_the_tumor=="Tumor-deficient zone")%>%  as.data.frame() %>% group_by(Segment_name_label) %>% drop_na(cell_density) %>% dplyr::summarise(Mean_cell_density = mean(`cell_density`), sd_cell_density = sd(`cell_density`)) %>% mutate(`Avg cell density in TD (%)` = paste0(round(Mean_cell_density,2),"(",round(sd_cell_density,2), ")")) %>% dplyr::select(Segment_name_label,`Avg cell density in TD (%)`), final %>% dplyr::select(Segment_name_label,Within_the_tumor, cell_density) %>% dplyr::filter(Within_the_tumor=="Tumor-rich zone") %>% as.data.frame() %>% group_by(Segment_name_label)%>% drop_na(cell_density) %>% dplyr::summarise(Mean_cell_density = mean(`cell_density`),sd_cell_density = sd(`cell_density`)) %>% mutate(`Avg cell density in TR (%)` = paste0(round(Mean_cell_density,2),"(",round(sd_cell_density,2), ")")) %>% dplyr::select(`Avg cell density in TR (%)`)) %>% dplyr::rename(c("Cell Types"=Segment_name_label)) %>% kbl(digits = 2, format = "html", row.names = TRUE, align=c(rep('c',times=4))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir") %>% print()
  
  cbind(final %>% dplyr::select(Segment_name_label,Within_the_tumor, cell_density) %>% dplyr::filter(Within_the_tumor=="Tumor-deficient zone")%>%  as.data.frame() %>% group_by(Segment_name_label)%>% drop_na(cell_density) %>% dplyr::summarise(Median_cell_density = median(`cell_density`), sd_cell_density = sd(`cell_density`)) %>% mutate(`Median cell density in TD (%)` = paste0(round(Median_cell_density,2),"(",round(sd_cell_density,2), ")")) %>% dplyr::select(Segment_name_label,`Median cell density in TD (%)`), final %>% dplyr::select(Segment_name_label,Within_the_tumor, cell_density) %>% dplyr::filter(Within_the_tumor=="Tumor-rich zone") %>% as.data.frame() %>% group_by(Segment_name_label)%>% drop_na(cell_density) %>% dplyr::summarise(Median_cell_density = median(`cell_density`),sd_cell_density = sd(`cell_density`)) %>% mutate(`Median cell density in TR (%)` = paste0(round(Median_cell_density,2),"(",round(sd_cell_density,2), ")")) %>% dplyr::select(`Median cell density in TR (%)`)) %>% dplyr::rename(c("Cell Types"=Segment_name_label))%>% kbl(digits = 2, format = "html", row.names = TRUE, align=c(rep('c',times=4))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir") %>% print()
  
  ggplot(dataset, aes(x=Within_the_tumor, y=cell_density, fill = Within_the_tumor)) + 
    facet_wrap(~Segment_name_label, scales = "free") +
    geom_violin(alpha=0.7, size=0.5,color="black") +
    #geom_boxplot(outlier.size = -1, color="black",lwd=0.5, alpha = 0.7) +
    geom_jitter(alpha=0.5, outlier.size = 0.5, outlier.color="red",width = 0.2, height = 0.1) + 
    stat_compare_means(paired = FALSE, aes(group = Within_the_tumor), method = "t.test", size = 3, label.x = 1.4, hide.ns = F, label.y.npc = 0.8 )+
    stat_compare_means(paired = FALSE, aes(group = Within_the_tumor), method = "wilcox.test", size = 3, label.x = 1.35, hide.ns = F, label.y.npc = 0.9) +theme(axis.text.x = element_text(size = 8, angle=0, vjust=0.5), axis.ticks.length.x = unit(0.1, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "none") +ylab("Cell density (%)") + xlab("")+scale_fill_manual(values=wes_palette(n=3, name="Chevalier1")) + ggtitle(paste0(type))
}

coverage_plots <- function(dataset, type, cell_density){
  
  dataset_wider <- dataset %>% dplyr::select(MERGING_IDalt, Segment_name_label, Within_the_tumor, cell_density) 
  dataset_wider <- aggregate(dataset_wider[[cell_density]], by = list(name=dataset_wider$MERGING_IDalt,dataset_wider$Segment_name_label, dataset_wider$Within_the_tumor), data=dataset_wider,FUN=mean) %>% spread(key = Group.2, value = x) 
  
  dataset_wider <- cbind(dataset_wider, dataset_wider %>% dplyr::select(`Early Cytotoxic`,`Early Non-Cytotoxic`, `Late Cytotoxic`,`Late Non-Cytotoxic`) %>% mutate(Total_area_coverage = rowSums(., na.rm = TRUE))%>% mutate(Total_Tumor_coverage = 100-Total_area_coverage) %>% dplyr::select(Total_area_coverage, Total_Tumor_coverage))
  
  dataset_wider2 <- dataset %>% dplyr::select(PATID,Segment_name_label, Within_the_tumor, cell_density)  
  dataset_wider2 <- aggregate(dataset_wider2[[cell_density]],by=list(name=dataset_wider2$PATID,dataset_wider2$Segment_name_label, dataset_wider2$Within_the_tumor),data=dataset_wider2,FUN=mean) %>% `colnames<-`(c("PATID", "Cell Type", "Tumor Zone", "Cell Density")) %>%  spread(key = "Tumor Zone", value = "Cell Density") %>% replace(is.na(.), 0) %>%mutate(sum = `Tumor-deficient zone`+`Tumor-rich zone`) %>% na_if(0)
  
  dataset_wider3 <-dataset_wider2 %>% dplyr::filter(`Cell Type`== "Early Cytotoxic" | `Cell Type`== "Early Non-Cytotoxic") 
  
  plot <- ggplot(dataset_wider3, aes(x=`Tumor-deficient zone`, y=`Tumor-rich zone`, fill = `Cell Type`)) + 
    theme(axis.text.x = element_text(size = 8, angle=0, vjust=0.5), axis.ticks.length.x = unit(0.1, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "bottom") +ylab("Cell density (%) in TR zone") + xlab("Cell density (%) in TD zone") +
    geom_rug(aes(color =`Cell Type`)) +
    geom_smooth(aes(color = `Cell Type`, fill = `Cell Type`), method = "lm", fullrange = FALSE, alpha = 0.4) +
    stat_poly_eq(formula = y~x, aes(label = paste(..eq.label..)), parse = TRUE)  +
    geom_point(aes(size=sum), alpha = 0.8, stroke = 0.9, shape = 21, show.legend = FALSE) +
    scale_color_manual(values=brewer.pal(n=4, name="Set2")) +
    scale_fill_manual(values=brewer.pal(n=4, name="Set2")) +
    scale_size(range = c(.1, 10)) +
    ggtitle(paste0(type))
  print(plot)
  
  dataset_wider4 <- dataset_wider %>% dplyr::filter(Group.3 == "Tumor-rich zone") %>% dplyr::select(-Total_area_coverage) %>%  dplyr::arrange(Total_Tumor_coverage) 
  #measurevars <- names(dataset_wider4)[grepl("[1-8]$",names(dataset_wider4))]
  dataset_wider4$name <- factor(dataset_wider4$name, levels = unique(dataset_wider4$name))
  # measure_list <- split(measurevars, split_on)
  # measudplyr::renames <- unique(dataset_wider4$Total_Tumor_coverage)
  dataset_wider4 <- reshape::melt(setDT(dataset_wider4))
  
  colours = c( "#8856a7", "#a6bddb", "#fa9fb5","#fec44f","#1c9099") 
  
  plot <- ggplot(dataset_wider4, aes(x=reorder(name,value), y=value, fill=variable)) +
    geom_bar(aes(), stat="identity", position="fill", colour = "black", size = 0.1) +
    #scale_fill_manual(values = c("darkgoldenrod1", "darkorchid", "red", "green", "darkblue")) +
    scale_fill_manual(values = colours)+
    guides(fill=guide_legend(nrow=1)) +
    xlab("") +
    ylab("% Coverage")+
    ggtitle("Tumor-rich zone") +
    theme(axis.text.x = element_text(size = 0), axis.ticks.length.x = unit(0.0, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7)) +
    scale_y_continuous(labels = scales::percent)+
    ggtitle(paste0(type))
  print(plot)
  
  dataset_wider4 <- dataset_wider %>% dplyr::filter(Group.3 == "Tumor-deficient zone") %>% replace(is.na(.), 0) %>% dplyr::select(-Total_area_coverage) %>% dplyr::arrange(Total_Tumor_coverage) %>% dplyr::rename(c("Unidentified"="Total_Tumor_coverage")) 
  #measurevars <- names(dataset_wider4)[grepl("[1-8]$",names(dataset_wider4))]
  dataset_wider4$name <- factor(dataset_wider4$name, levels = unique(dataset_wider4$name))
  # measure_list <- split(measurevars, split_on)
  # measudplyr::renames <- unique(dataset_wider4$Total_Tumor_coverage)
  dataset_wider4 <- reshape::melt(setDT(dataset_wider4))
  
  colours = c( "#8856a7", "#a6bddb", "#fa9fb5","#fec44f","#1c9099") 
  
  plot <- ggplot(dataset_wider4, aes(x=reorder(name,value), y=value, fill=variable)) +
    geom_bar(aes(), stat="identity", position="fill", colour = "black", size = 0.3) +
    #scale_fill_manual(values = c("darkgoldenrod1", "darkorchid", "red", "green", "darkblue")) +
    scale_fill_manual(values = colours)+
    guides(fill=guide_legend(nrow=1)) +
    xlab("") +
    ylab("% Coverage")+
    ggtitle("Tumor-deficient zone") +
    theme(axis.text.x = element_text(size = 0), axis.ticks.length.x = unit(0.0, "cm"),axis.ticks.x = element_line(size = 0.1),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white"),panel.border = element_rect(fill = NA), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 7)) +
    scale_y_continuous(labels = scales::percent)+
    ggtitle(paste0(type))
  print(plot)
}


coxme_model_immunesubtype <- function(dataset, type){
  
  dataset_wider <- dataset %>% dplyr::select(combined_unique_coreID, PATID, MERGING_IDalt, Segment_name_label, Within_the_tumor, cell_density) %>% left_join(my_data_aoi %>% select(`OS_(0,1)`, `OS_(years)`, "PATID", "MERGING_IDalt", "Within_the_tumor", "combined_unique_coreID"), by = c("combined_unique_coreID","PATID", "Within_the_tumor", "MERGING_IDalt")) %>% distinct()%>% spread(key = Segment_name_label, value = cell_density) 
  
  names(dataset_wider)<-str_replace_all(names(dataset_wider), c(" " = "_", "-"="_"))
  
  dataset_wider_subset <- cbind(dataset_wider, dataset_wider %>% dplyr::select(Early_Non_Cytotoxic, Early_Cytotoxic, Late_Non_Cytotoxic,Late_Cytotoxic) %>% mutate(Total_area_coverage = rowSums(., na.rm = TRUE)) %>% mutate(Total_Tumor_coverage = 100-Total_area_coverage) %>% dplyr::select(Total_area_coverage, Total_Tumor_coverage)) %>% dplyr::filter(Within_the_tumor =="Tumor-rich zone")
  
  vl <- names(dataset_wider_subset[,7:ncol(dataset_wider_subset)])
  len <- length(vl)
  final_uni_table <- data.frame()
  for (j in 1:len) {
    dataset <-  dataset_wider_subset %>% select(PATID,`OS_(years)`, `OS_(0,1)`, paste(vl[j]) ) %>% drop_na()
    f1 = as.formula(paste("Surv(as.numeric(as.character(`OS_(years)`)),as.numeric(as.character(`OS_(0,1)`))) ~ ", paste0(vl[j]) ,"+ (1|PATID)"))
    md <- coxme(f1, data=dataset_wider_subset)
    obj <- univ_cox_me(md)
    final_uni_table <- final_uni_table %>% rbind(obj)
  }
  
  final_uni_table  %>% kbl(digits = 2, format = "html", caption = paste0(type),row.names = TRUE, align=c(rep('c',times=4))) %>%kable_styling(c("striped", "hover","condensed"),full_width = F,font_size = 12,position = "left", html_font = "Avenir") 
}

