### THIS CODE RUNS HYPERGRAPH ANALAYSIS AND CALCULATES ROW SUM DISTRIBUTION ON PSEUDOBULK RNA SAMPLES - GROUPED PATHOLOGY DATA. THIS DATA CAN THEN BE TAKEN FORWARD FOR GSEA ANALYSIS
### TO IDENTIFY PATHWAYS THAT ARE LINKED TO GENES WHICH ARE HIGHLY CORRELATED IN THE HYPERGRAPH.


# Load in packages needed
#BiocManager::install("dendextend")
library(dendextend)
#BiocManager::install("gplots")
library(gplots)
#BiocManager::install("plyr")
library(plyr)
#BiocManager::install("dplyr")
#BiocManager::install("forcats")
library(forcats)
library(ggplot2)
library(ggrepel)

setwd("~/Wound_Healing/New_analysis/R_hypergraphs")

#Load expression data and spatial coordinate data
exp_data <- read.csv("~/Wound_Healing/New_analysis/R_hypergraphs/Filtered_wound_data_May_2024/X.csv",header=F, sep = ',')
spatial_obs<-read.csv("~/Wound_Healing/New_analysis/R_hypergraphs/Filtered_wound_data_May_2024/obs.csv",header=T)
rownames(exp_data)<-spatial_obs$X

# Load in the var data
var_data<- read.csv("~/Wound_Healing/New_analysis/R_hypergraphs/Filtered_wound_data_May_2024/var.csv",header=T)

## ROWS NEED TO BE SAMPLES, COLUMNS ARE GENES
# Refine dataset as needed. Assign gene names to rownames, remove any additional variables (eg sample metadata), ensure rows are samples, columns are genes.

colnames(exp_data)<-var_data$X
#sd scores for Clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0
#run hypergraphs

#### Acute samples ####
set.seed(3)

Acute_samples <- spatial_obs[spatial_obs$type == "Acute", ]

Acute_samples_X <- exp_data[na.omit(match(Acute_samples$X,rownames(exp_data))),]

sd.scores0<-apply(Acute_samples_X,2,sd) #2 = selects columns not rows
Acute_samples_X_sd<-Acute_samples_X[,which(sd.scores0>0)]

# Run the correlation
Acute_gene_correlation <- cor(Acute_samples_X_sd,Acute_samples_X_sd)

#Binarize correlation values using sd of correlation matrix
bin_acute <-abs(Acute_gene_correlation)
bin_acute[which(bin_acute>sd(Acute_gene_correlation))]<-1
bin_acute[which(bin_acute!=1)]<-0 # this is the hypergraph incidence matrix

#Matrix multiplication to generate hypergraph adjacency matrix
hyp_acute<-bin_acute %*% t(bin_acute) #adjacency matrix 


#ranking
Acute_rowsum_result<- rowSums(hyp_acute)
hist(Acute_rowsum_result)
Acute_row_sums <- as.data.frame(Acute_rowsum_result)

colnames(Acute_row_sums)[1] <- "rowsum"
Acute_row_sums$Pathology <- "Acute"  
Acute_row_sums$rank <- rank(Acute_row_sums$rowsum, ties.method = "average")
Acute_row_sums$norm_rank <- (Acute_row_sums$rank)/17770 ##### divide by total number of row sum genes

Acute_row_sums <- Acute_row_sums[order(Acute_row_sums$rank, decreasing = TRUE), ]

write.csv(Acute_row_sums, "~/Wound_Healing/New_analysis/R_hypergraphs/output_csvs/Full_samples_row_sum/Acute_row_sums_results_ranked_Oct_24.csv", row.names=TRUE)

#### Chronic samples ####
set.seed(3)
Chronic_samples <- spatial_obs[spatial_obs$type == "Chronic", ]
Chronic_samples_X <- exp_data[na.omit(match(Chronic_samples$X,rownames(exp_data))),]

sd.scores1<-apply(Chronic_samples_X,2,sd) #2 = selects columns not rows
Chronic_samples_X_sd<-Chronic_samples_X[,which(sd.scores1>0)]

# Run the correlation
Chronic_gene_correlation <- cor(Chronic_samples_X_sd, Chronic_samples_X_sd)

#Binarize correlation values using sd of correlation matrix
bin_chronic <-abs(Chronic_gene_correlation)
bin_chronic[which(bin_chronic>sd(Chronic_gene_correlation))]<-1
bin_chronic[which(bin_chronic!=1)]<-0 # this is the hypergraph incidence matrix

#Matrix multiplication to generate hypergraph adjacency matrix
hyp_chronic<-bin_chronic %*% t(bin_chronic) #adjacency matrix 


#ranking
Chronic_rowsum_result<- rowSums(hyp_chronic)
hist(Chronic_rowsum_result)
Chronic_row_sums <- as.data.frame(Chronic_rowsum_result)

colnames(Chronic_row_sums)[1] <- "rowsum"
Chronic_row_sums$Pathology <- "Chronic"  
Chronic_row_sums$rank <- rank(-Chronic_row_sums$rowsum, ties.method = "average")
Chronic_row_sums$norm_rank <- (Chronic_row_sums$rank)/17676 ##### divide by total number of row sum genes

Chronic_row_sums <- Chronic_row_sums[order(Chronic_row_sums$rank, decreasing = TRUE), ]

write.csv(Chronic_row_sums, "~/Wound_Healing/New_analysis/R_hypergraphs/output_csvs/Full_samples_row_sum/Chronic_row_sums_results_ranked_Oct_24.csv", row.names=TRUE)


#### Control samples ####
set.seed(3)
Control_samples <- spatial_obs[spatial_obs$type == "Control", ]
Control_samples_X <- exp_data[na.omit(match(Control_samples$X,rownames(exp_data))),]

sd.scores1<-apply(Control_samples_X,2,sd) #2 = selects columns not rows
Control_samples_X_sd<-Control_samples_X[,which(sd.scores1>0)]

# Run the correlation
Control_gene_correlation <- cor(Control_samples_X_sd,Control_samples_X_sd)

#Binarize correlation values using sd of correlation matrix
bin_control <-abs(Control_gene_correlation)
bin_control[which(bin_control>sd(Control_gene_correlation))]<-1
bin_control[which(bin_control!=1)]<-0 # this is the hypergraph incidence matrix

#Matrix multiplication to generate hypergraph adjacency matrix
hyp_control<-bin_control%*% t(bin_control) #adjacency matrix 


#ranking
Control_rowsum_result<- rowSums(hyp_control)
hist(Control_rowsum_result)
Control_row_sums <- as.data.frame(Control_rowsum_result)

colnames(Control_row_sums)[1] <- "rowsum"
Control_row_sums$Pathology <- "Control"  
Control_row_sums$rank <- rank(-Control_row_sums$rowsum, ties.method = "average")
Control_row_sums$norm_rank <- (Control_row_sums$rank)/17290 ##### divide by total number of row sum genes

Control_row_sums <- Control_row_sums[order(Control_row_sums$rank, decreasing = TRUE), ]

write.csv(Control_row_sums, "~/Wound_Healing/New_analysis/R_hypergraphs/output_csvs/Full_samples_row_sum/Control_row_sums_results_ranked.csv", row.names=TRUE)