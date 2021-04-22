library(Repitools)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)
library(tidyverse)
library(plyr)
library(data.table)
library(survival)
library(survminer)

Pediatric_DEGs <- read.csv("Pediatric_DEGs.csv", header = FALSE)
Pediatric_ProjCodes <- read.csv("Pediatric_ProjCodes.csv", header = TRUE)

for(i in 1:8) {
  proj <- Pediatric_ProjCodes[i, "proj"]
  sample_type <- Pediatric_ProjCodes[i, "sample_type"]
  name <- Pediatric_ProjCodes[i, "name"]
  
  query <- GDCquery(project = proj,
                    data.category = "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts")
  
  if(proj == "TARGET-NBL") {
    query.2 = query
    temp = query.2$results[[1]]
    temp = temp[which(!duplicated(temp$cases)),]
    query.2$results[[1]] = temp
    GDCdownload(query.2)
    TARGET_RNA_Seq_Counts <- GDCprepare(query.2)
    TARGET_RNA_Seq_Counts_df <- GDCprepare(query.2, summarizedExperiment = FALSE)
  } else {
    GDCdownload(query)
    TARGET_RNA_Seq_Counts <- GDCprepare(query)
    TARGET_RNA_Seq_Counts_df <- GDCprepare(query, summarizedExperiment = FALSE)
  }
  
  TARGET_RNA_Seq_Counts_subset_type <- TARGET_RNA_Seq_Counts[, TARGET_RNA_Seq_Counts$sample_type == sample_type]
  
  ## Subset for DEGs based on matched "original_ensembl_gene_id"
  RNA_Seq_Counts <- rowRanges(TARGET_RNA_Seq_Counts)
  RNA_Seq_Counts_subset <- RNA_Seq_Counts[(elementMetadata(RNA_Seq_Counts)[, "external_gene_name"] %in% Pediatric_DEGs$V1)]
  RNA_Seq_Counts_subset_df <- annoGR2DF(RNA_Seq_Counts_subset)
  
  TARGET_RNA_Seq_Counts_subset_df <- subset(TARGET_RNA_Seq_Counts_df, X1 %in% RNA_Seq_Counts_subset_df$original_ensembl_gene_id)
  
  ClinData <- colData(TARGET_RNA_Seq_Counts_subset_type)
  ClinData_df <- as(ClinData, "data.frame")
  
  colnames(TARGET_RNA_Seq_Counts_subset_df)[1] <- "original_ensembl_gene_id"
  TARGET_RNA_Seq_Counts_subset_gene_df <- merge(TARGET_RNA_Seq_Counts_subset_df, unique(RNA_Seq_Counts_subset_df)[, c("original_ensembl_gene_id", "external_gene_name")], by = "original_ensembl_gene_id")
  rownames(TARGET_RNA_Seq_Counts_subset_gene_df) <- TARGET_RNA_Seq_Counts_subset_gene_df$external_gene_name
  rm_cols <- c("original_ensembl_gene_id", "id", "external_gene_name")
  TARGET_RNA_Seq_Counts_subset_gene_df <- TARGET_RNA_Seq_Counts_subset_gene_df[ , !(names(TARGET_RNA_Seq_Counts_subset_gene_df) %in% rm_cols)]
  
  ## Subset for Primary Tumor samples based on previously subsetted ClinData df 
  TARGET_RNA_Seq_Counts_subset_gene_type_df <- TARGET_RNA_Seq_Counts_subset_df %>% select(c(ClinData_df$barcode))
  rownames(TARGET_RNA_Seq_Counts_subset_gene_type_df) <- RNA_Seq_Counts_subset_df$external_gene_name
  
  ## Begin iterating through Pediatric DEGs
  for(j in 1:89) {
    gene <- Pediatric_DEGs[j, "V1"]
    
    subset_gene <- TARGET_RNA_Seq_Counts_subset_gene_type_df[sprintf("%s", gene),]
    rownames(subset_gene)[1] <- gene
    t_subset_gene <- t(subset_gene)
    gene_df <- as.data.frame(t_subset_gene)
    
    # Calculate tertile ranges
    vTert = unique(quantile(t_subset_gene[,1], c(0:3/3)))
    
    if(length(vTert) == 4) {
      # Classify expression values
      gene_df$tert = with(gene_df,
                          cut(gene_df[,1],
                              vTert,
                              include.lowest = T,
                              labels = c("Low", "Medium", "High")))
      
      
      gene_df <- gene_df[gene_df$tert != "Medium", ]
      setDT(gene_df, keep.rownames = "barcode")
      gene_clin_df <- merge(gene_df, unique(ClinData_df)[, c("barcode", "age_at_diagnosis", "gender", "vital_status", "days_to_last_follow_up")], by = "barcode")
      
      # ALIVE = 1, DEAD = 2
      # MALE = 1, FEMALE = 2
      gene_clin_df$vital_status <- revalue(gene_clin_df$vital_status, c("Alive"=1))
      gene_clin_df$vital_status <- revalue(gene_clin_df$vital_status, c("Dead"=2))
      gene_clin_df$gender <- revalue(gene_clin_df$gender, c("male"=1))
      gene_clin_df$gender <- revalue(gene_clin_df$gender, c("female"=2))
      gene_clin_df$vital_status = as.numeric(as.character(gene_clin_df$vital_status))
      gene_clin_df$gender = as.numeric(as.character(gene_clin_df$gender))
      
      res.cox <- coxph(Surv(days_to_last_follow_up, vital_status) ~ tert, data = gene_clin_df)
      res.cox_summ <- summary(res.cox)
      res.cox_summ_coef <- res.cox_summ$coefficients
      res.cox_summ_coef <- as.data.frame(res.cox_summ_coef)
      res.cox_summ_coef <- res.cox_summ_coef[!(row.names(res.cox_summ_coef) == "tertMedium"), ]
      rownames(res.cox_summ_coef)[1] <- gene
      
      if(j==1) {
        all.res.cox_summ_coef <- res.cox_summ_coef
      } else {
        all.res.cox_summ_coef <- rbind(all.res.cox_summ_coef, res.cox_summ_coef)
      }            
      
      file_name = sprintf("Results/Pediatric/Survival Plots/%s/KM_plot_%s_%s.tiff", name, name, gene)
      fit <- survfit(Surv(days_to_last_follow_up, vital_status) ~ tert, data = gene_clin_df)
      
      ggsurv <- ggsurvplot(
        fit,
        data = gene_clin_df,
        title = sprintf("TARGET-AML (TBM)\n%s", gene),
        palette =
          c("#457B9D", "#E63946"),# custom color palettes
        conf.int = TRUE,          # Add confidence interval
        pval = TRUE,              # Add p-value
        legend.title = "",
        legend.labs =
          c("Low", "High"),    # Change legend labels
        xlab = "Overall survival time (days)",
        font.title = c(18, "bold"),
        font.x = c(18, "bold"),
        font.y = c(18, "bold"),
        font.tickslab = 16,
        font.legend = 16,
        pval.size = 7,
        pval.coord = c(500, 0.2),
        axes.offset = FALSE
      )
      
      print(ggsurv)
      ggsave(file_name, width = 10, height = 10)
      dev.off()
    }
    print(j)
  }
  
  res.cox.gender <- coxph(Surv(days_to_last_follow_up, vital_status) ~ gender, data = gene_clin_df)
  res.cox.gender_summ <- summary(res.cox.gender)
  res.cox.gender_summ_coef <- res.cox.gender_summ$coefficients
  res.cox.gender_summ_coef <- as.data.frame(res.cox.gender_summ_coef)
  rownames(res.cox.gender_summ_coef)[1] <- "Gender"
  
  res.cox.age_at_diagnosis <- coxph(Surv(days_to_last_follow_up, vital_status) ~ age_at_diagnosis, data = gene_clin_df)
  res.cox.age_at_diagnosis_summ <- summary(res.cox.age_at_diagnosis)
  res.cox.age_at_diagnosis_summ_coef <- res.cox.age_at_diagnosis_summ$coefficients
  res.cox.age_at_diagnosis_summ_coef <- as.data.frame(res.cox.age_at_diagnosis_summ_coef)
  rownames(res.cox.age_at_diagnosis_summ_coef)[1] <- "Age"
  
  all.res.cox_summ_coef <- rbind(all.res.cox_summ_coef, res.cox.gender_summ_coef)    
  all.res.cox_summ_coef <- rbind(all.res.cox_summ_coef, res.cox.age_at_diagnosis_summ_coef)    
  write.csv(AML.res.cox_summ_coef, sprintf("Results/Pediatric/Survival Plots/%s/Res.Cox_Summ_%s_%s.tiff", name, name, gene))
}
      