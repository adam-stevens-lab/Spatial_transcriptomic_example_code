### This script runs a 10,0000 iteration hypergraph to calculate entropy per cluster, per pathology sample
### The data is filtered to only include samples that are Acute, or Chronic, or Control. A randomly selected 100 genes are correlated to a randomly selected 1000 genes, with the 
### entropy calculated each run of the for loop. The for loop is iterated 10,000 times and the data saved out.
### This can be run as an additional for loop per pathology

# load in packages
library(gplots)
library(plyr)
library(forcats)
library(entropy)
library(ggplot2)

setwd("")

#Load expression data and spatial coordinate data
exp_data <- read.csv("path_to_file/X.csv",header=F, sep = ',')
spatial_obs<-read.csv("~path_to_file/obs.csv",header=T)
rownames(exp_data)<-spatial_obs$X

spatial_obs$cell_type<-as.factor(spatial_obs$leiden_1.7) #need to change cell type col from a character to a factor
var_data<- read.csv("path_to_file/var.csv",header=T)


## ROWS NEED TO BE SAMPLES, COLUMNS ARE GENES

# Refine dataset as needed. Assign gene names to rownames, remove any additional variables (eg sample metadata), ensure rows are samples, columns are genes.

colnames(exp_data)<-var_data$X

#sd scores for Clusters # IMPORTANT - Correlation will not work on genes with no variance - remove genes with SD of 0

#####ACUTE SAMPLES#####
#- generate a acute only exp data
#then filter clusters

Acute_samples <- spatial_obs[spatial_obs$type == "Acute", ]

Acute_samples_X <- exp_data[na.omit(match(Acute_samples$X,rownames(exp_data))),]

#### Cluster_0 ####
set.seed(3)
Cluster_0 <- spatial_obs[spatial_obs$cell_type == "0", ]
Cluster_0_X <- Acute_samples_X[na.omit(match(Cluster_0$X,rownames(Acute_samples_X))),]


sd.scores0<-apply(Cluster_0_X,2,sd) #2 = selects columns not rows
Cluster_0_X_sd<-Cluster_0_X[,which(sd.scores0>0)]

Cluster_0_result_df <- data.frame()
cor_data_Cluster_0<-cor(Cluster_0_X_sd, Cluster_0_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_0_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_0_DEGs<- colnames(Cluster_0_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_0_2 <- cor_data_Cluster_0[na.omit(match(Cluster_0_DEGs,rownames(cor_data_Cluster_0))),-na.omit(match(Cluster_0_DEGs,colnames(cor_data_Cluster_0)))] #na omit the ones that arent in my selected 100
  Cluster_0_sub_sample <- colnames(cor_data_Cluster_0_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_0_2 <- cor_data_Cluster_0_2[,na.omit(match(Cluster_0_sub_sample,colnames(cor_data_Cluster_0_2)))] #take all the rows and sample the columns
  bin_Cluster_0<-abs(cor_data_Cluster_0_2)
  bin_Cluster_0[which(bin_Cluster_0>sd(cor_data_Cluster_0_2))]<-1
  bin_Cluster_0[which(bin_Cluster_0!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_0<-bin_Cluster_0 %*% t(bin_Cluster_0) #adjacency matrix - 
  Cluster_0_entropy_result<- entropy(hyp_Cluster_0)
  Cluster_0_result_df <- rbind(Cluster_0_result_df, Cluster_0_entropy_result)
}
write.csv(Cluster_0_result_df, "path_to_file/Acute_cluster_0_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_1 ####
set.seed(3)
Cluster_1 <- spatial_obs[spatial_obs$cell_type == "1", ]
Cluster_1_X <- Acute_samples_X[na.omit(match(Cluster_1$X,rownames(Acute_samples_X))),]

sd.scores1<-apply(Cluster_1_X,2,sd) #2 = selects columns not rows
Cluster_1_X_sd<-Cluster_1_X[,which(sd.scores1>0)]


Cluster_1_result_df <- data.frame()
cor_data_Cluster_1<-cor(Cluster_1_X_sd, Cluster_1_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_1_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_1_DEGs<- colnames(Cluster_1_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_1_2 <- cor_data_Cluster_1[na.omit(match(Cluster_1_DEGs,rownames(cor_data_Cluster_1))),-na.omit(match(Cluster_1_DEGs,colnames(cor_data_Cluster_1)))] #na omit the ones that arent in my selected 100
  Cluster_1_sub_sample <- colnames(cor_data_Cluster_1_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_1_2 <- cor_data_Cluster_1_2[,na.omit(match(Cluster_1_sub_sample,colnames(cor_data_Cluster_1_2)))] #take all the rows and sample the columns
  bin_Cluster_1<-abs(cor_data_Cluster_1_2)
  bin_Cluster_1[which(bin_Cluster_1>sd(cor_data_Cluster_1_2))]<-1
  bin_Cluster_1[which(bin_Cluster_1!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_1<-bin_Cluster_1 %*% t(bin_Cluster_1) #adjacency matrix - 
  Cluster_1_entropy_result<- entropy(hyp_Cluster_1)
  Cluster_1_result_df <- rbind(Cluster_1_result_df, Cluster_1_entropy_result)
}
write.csv(Cluster_1_result_df, "path_to_file/Acute_cluster_1_entropy_results_10000.csv", row.names=FALSE)


#### Cluster_2 ####
set.seed(3)
Cluster_2 <- spatial_obs[spatial_obs$cell_type == "2", ]
Cluster_2_X <- Acute_samples_X[na.omit(match(Cluster_2$X,rownames(Acute_samples_X))),]

sd.scores2<-apply(Cluster_2_X,2,sd) #2 = selects columns not rows
Cluster_2_X_sd<-Cluster_2_X[,which(sd.scores2>0)]



Cluster_2_result_df <- data.frame()
cor_data_Cluster_2<-cor(Cluster_2_X_sd, Cluster_2_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_2_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_2_DEGs<- colnames(Cluster_2_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_2_2 <- cor_data_Cluster_2[na.omit(match(Cluster_2_DEGs,rownames(cor_data_Cluster_2))),-na.omit(match(Cluster_2_DEGs,colnames(cor_data_Cluster_2)))] #na omit the ones that arent in my selected 100
  Cluster_2_sub_sample <- colnames(cor_data_Cluster_2_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_2_2 <- cor_data_Cluster_2_2[,na.omit(match(Cluster_2_sub_sample,colnames(cor_data_Cluster_2_2)))] #take all the rows and sample the columns
  bin_Cluster_2<-abs(cor_data_Cluster_2_2)
  bin_Cluster_2[which(bin_Cluster_2>sd(cor_data_Cluster_2_2))]<-1
  bin_Cluster_2[which(bin_Cluster_2!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_2<-bin_Cluster_2 %*% t(bin_Cluster_2) #adjacency matrix - 
  Cluster_2_entropy_result<- entropy(hyp_Cluster_2)
  Cluster_2_result_df <- rbind(Cluster_2_result_df, Cluster_2_entropy_result)
}
write.csv(Cluster_2_result_df, "path_to_file/Acute_cluster_2_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_3 ####
set.seed(3)
Cluster_3 <- spatial_obs[spatial_obs$cell_type == "3", ]
Cluster_3_X <- Acute_samples_X[na.omit(match(Cluster_3$X,rownames(Acute_samples_X))),]

sd.scores3<-apply(Cluster_3_X,2,sd) #2 = selects columns not rows
Cluster_3_X_sd<-Cluster_3_X[,which(sd.scores3>0)]


Cluster_3_result_df <- data.frame()
cor_data_Cluster_3<-cor(Cluster_3_X_sd, Cluster_3_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_3_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_3_DEGs<- colnames(Cluster_3_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_3_2 <- cor_data_Cluster_3[na.omit(match(Cluster_3_DEGs,rownames(cor_data_Cluster_3))),-na.omit(match(Cluster_3_DEGs,colnames(cor_data_Cluster_3)))] #na omit the ones that arent in my selected 100
  Cluster_3_sub_sample <- colnames(cor_data_Cluster_3_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_3_2 <- cor_data_Cluster_3_2[,na.omit(match(Cluster_3_sub_sample,colnames(cor_data_Cluster_3_2)))] #take all the rows and sample the columns
  bin_Cluster_3<-abs(cor_data_Cluster_3_2)
  bin_Cluster_3[which(bin_Cluster_3>sd(cor_data_Cluster_3_2))]<-1
  bin_Cluster_3[which(bin_Cluster_3!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_3<-bin_Cluster_3 %*% t(bin_Cluster_3) #adjacency matrix - 
  Cluster_3_entropy_result<- entropy(hyp_Cluster_3)
  Cluster_3_result_df <- rbind(Cluster_3_result_df, Cluster_3_entropy_result)
}
write.csv(Cluster_3_result_df, "path_to_file/Acute_cluster_3_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_4 ####
set.seed(3)
Cluster_4 <- spatial_obs[spatial_obs$cell_type == "4", ]
Cluster_4_X <- Acute_samples_X[na.omit(match(Cluster_4$X,rownames(Acute_samples_X))),]

sd.scores4<-apply(Cluster_4_X,2,sd) #2 = selects columns not rows
Cluster_4_X_sd<-Cluster_4_X[,which(sd.scores4>0)]


Cluster_4_result_df <- data.frame()
cor_data_Cluster_4<-cor(Cluster_4_X_sd, Cluster_4_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_4_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_4_DEGs<- colnames(Cluster_4_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_4_2 <- cor_data_Cluster_4[na.omit(match(Cluster_4_DEGs,rownames(cor_data_Cluster_4))),-na.omit(match(Cluster_4_DEGs,colnames(cor_data_Cluster_4)))] #na omit the ones that arent in my selected 100
  Cluster_4_sub_sample <- colnames(cor_data_Cluster_4_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_4_2 <- cor_data_Cluster_4_2[,na.omit(match(Cluster_4_sub_sample,colnames(cor_data_Cluster_4_2)))] #take all the rows and sample the columns
  bin_Cluster_4<-abs(cor_data_Cluster_4_2)
  bin_Cluster_4[which(bin_Cluster_4>sd(cor_data_Cluster_4_2))]<-1
  bin_Cluster_4[which(bin_Cluster_4!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_4<-bin_Cluster_4 %*% t(bin_Cluster_4) #adjacency matrix - 
  Cluster_4_entropy_result<- entropy(hyp_Cluster_4)
  Cluster_4_result_df <- rbind(Cluster_4_result_df, Cluster_4_entropy_result)
}
write.csv(Cluster_4_result_df, "path_to_file/Acute_cluster_4_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_5 ####
set.seed(3)
Cluster_5 <- spatial_obs[spatial_obs$cell_type == "5", ]
Cluster_5_X <- Acute_samples_X[na.omit(match(Cluster_5$X,rownames(Acute_samples_X))),]

sd.scores5<-apply(Cluster_5_X,2,sd) #2 = selects columns not rows
Cluster_5_X_sd<-Cluster_5_X[,which(sd.scores5>0)]


Cluster_5_result_df <- data.frame()
cor_data_Cluster_5<-cor(Cluster_5_X_sd, Cluster_5_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_5_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_5_DEGs<- colnames(Cluster_5_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_5_2 <- cor_data_Cluster_5[na.omit(match(Cluster_5_DEGs,rownames(cor_data_Cluster_5))),-na.omit(match(Cluster_5_DEGs,colnames(cor_data_Cluster_5)))] #na omit the ones that arent in my selected 100
  Cluster_5_sub_sample <- colnames(cor_data_Cluster_5_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_5_2 <- cor_data_Cluster_5_2[,na.omit(match(Cluster_5_sub_sample,colnames(cor_data_Cluster_5_2)))] #take all the rows and sample the columns
  bin_Cluster_5<-abs(cor_data_Cluster_5_2)
  bin_Cluster_5[which(bin_Cluster_5>sd(cor_data_Cluster_5_2))]<-1
  bin_Cluster_5[which(bin_Cluster_5!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_5<-bin_Cluster_5 %*% t(bin_Cluster_5) #adjacency matrix - 
  Cluster_5_entropy_result<- entropy(hyp_Cluster_5)
  Cluster_5_result_df <- rbind(Cluster_5_result_df, Cluster_5_entropy_result)
}
write.csv(Cluster_5_result_df, "path_to_file/Acute_cluster_5_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_6 ####
set.seed(3)
Cluster_6 <- spatial_obs[spatial_obs$cell_type == "6", ]
Cluster_6_X <- Acute_samples_X[na.omit(match(Cluster_6$X,rownames(Acute_samples_X))),]

sd.scores6<-apply(Cluster_6_X,2,sd) #2 = selects columns not rows
Cluster_6_X_sd<-Cluster_6_X[,which(sd.scores6>0)]


Cluster_6_result_df <- data.frame()
cor_data_Cluster_6<-cor(Cluster_6_X_sd, Cluster_6_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_6_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_6_DEGs<- colnames(Cluster_6_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_6_2 <- cor_data_Cluster_6[na.omit(match(Cluster_6_DEGs,rownames(cor_data_Cluster_6))),-na.omit(match(Cluster_6_DEGs,colnames(cor_data_Cluster_6)))] #na omit the ones that arent in my selected 100
  Cluster_6_sub_sample <- colnames(cor_data_Cluster_6_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_6_2 <- cor_data_Cluster_6_2[,na.omit(match(Cluster_6_sub_sample,colnames(cor_data_Cluster_6_2)))] #take all the rows and sample the columns
  bin_Cluster_6<-abs(cor_data_Cluster_6_2)
  bin_Cluster_6[which(bin_Cluster_6>sd(cor_data_Cluster_6_2))]<-1
  bin_Cluster_6[which(bin_Cluster_6!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_6<-bin_Cluster_6 %*% t(bin_Cluster_6) #adjacency matrix - 
  Cluster_6_entropy_result<- entropy(hyp_Cluster_6)
  Cluster_6_result_df <- rbind(Cluster_6_result_df, Cluster_6_entropy_result)
}
write.csv(Cluster_6_result_df, "path_to_file/Acute_cluster_6_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_7 ####
set.seed(3)
Cluster_7 <- spatial_obs[spatial_obs$cell_type == "7", ]
Cluster_7_X <- Acute_samples_X[na.omit(match(Cluster_7$X,rownames(Acute_samples_X))),]

sd.scores7<-apply(Cluster_7_X,2,sd) #2 = selects columns not rows
Cluster_7_X_sd<-Cluster_7_X[,which(sd.scores7>0)]


Cluster_7_result_df <- data.frame()
cor_data_Cluster_7<-cor(Cluster_7_X_sd, Cluster_7_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_7_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_7_DEGs<- colnames(Cluster_7_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_7_2 <- cor_data_Cluster_7[na.omit(match(Cluster_7_DEGs,rownames(cor_data_Cluster_7))),-na.omit(match(Cluster_7_DEGs,colnames(cor_data_Cluster_7)))] #na omit the ones that arent in my selected 100
  Cluster_7_sub_sample <- colnames(cor_data_Cluster_7_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_7_2 <- cor_data_Cluster_7_2[,na.omit(match(Cluster_7_sub_sample,colnames(cor_data_Cluster_7_2)))] #take all the rows and sample the columns
  bin_Cluster_7<-abs(cor_data_Cluster_7_2)
  bin_Cluster_7[which(bin_Cluster_7>sd(cor_data_Cluster_7_2))]<-1
  bin_Cluster_7[which(bin_Cluster_7!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_7<-bin_Cluster_7 %*% t(bin_Cluster_7) #adjacency matrix - 
  Cluster_7_entropy_result<- entropy(hyp_Cluster_7)
  Cluster_7_result_df <- rbind(Cluster_7_result_df, Cluster_7_entropy_result)
}
write.csv(Cluster_7_result_df, "path_to_file/Acute_cluster_7_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_8 ####
set.seed(3)
Cluster_8 <- spatial_obs[spatial_obs$cell_type == "8", ]
Cluster_8_X <- Acute_samples_X[na.omit(match(Cluster_8$X,rownames(Acute_samples_X))),]

sd.scores8<-apply(Cluster_8_X,2,sd) #2 = selects columns not rows
Cluster_8_X_sd<-Cluster_8_X[,which(sd.scores8>0)]


Cluster_8_result_df <- data.frame()
cor_data_Cluster_8<-cor(Cluster_8_X_sd, Cluster_8_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_8_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_8_DEGs<- colnames(Cluster_8_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_8_2 <- cor_data_Cluster_8[na.omit(match(Cluster_8_DEGs,rownames(cor_data_Cluster_8))),-na.omit(match(Cluster_8_DEGs,colnames(cor_data_Cluster_8)))] #na omit the ones that arent in my selected 100
  Cluster_8_sub_sample <- colnames(cor_data_Cluster_8_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_8_2 <- cor_data_Cluster_8_2[,na.omit(match(Cluster_8_sub_sample,colnames(cor_data_Cluster_8_2)))] #take all the rows and sample the columns
  bin_Cluster_8<-abs(cor_data_Cluster_8_2)
  bin_Cluster_8[which(bin_Cluster_8>sd(cor_data_Cluster_8_2))]<-1
  bin_Cluster_8[which(bin_Cluster_8!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_8<-bin_Cluster_8 %*% t(bin_Cluster_8) #adjacency matrix - 
  Cluster_8_entropy_result<- entropy(hyp_Cluster_8)
  Cluster_8_result_df <- rbind(Cluster_8_result_df, Cluster_8_entropy_result)
}
write.csv(Cluster_8_result_df, "path_to_file/Acute_cluster_8_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_9 ####
set.seed(3)
Cluster_9 <- spatial_obs[spatial_obs$cell_type == "9", ]
Cluster_9_X <- Acute_samples_X[na.omit(match(Cluster_9$X,rownames(Acute_samples_X))),]

sd.scores9<-apply(Cluster_9_X,2,sd) #2 = selects columns not rows
Cluster_9_X_sd<-Cluster_9_X[,which(sd.scores9>0)]


Cluster_9_result_df <- data.frame()
cor_data_Cluster_9<-cor(Cluster_9_X_sd, Cluster_9_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_9_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_9_DEGs<- colnames(Cluster_9_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_9_2 <- cor_data_Cluster_9[na.omit(match(Cluster_9_DEGs,rownames(cor_data_Cluster_9))),-na.omit(match(Cluster_9_DEGs,colnames(cor_data_Cluster_9)))] #na omit the ones that arent in my selected 100
  Cluster_9_sub_sample <- colnames(cor_data_Cluster_9_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_9_2 <- cor_data_Cluster_9_2[,na.omit(match(Cluster_9_sub_sample,colnames(cor_data_Cluster_9_2)))] #take all the rows and sample the columns
  bin_Cluster_9<-abs(cor_data_Cluster_9_2)
  bin_Cluster_9[which(bin_Cluster_9>sd(cor_data_Cluster_9_2))]<-1
  bin_Cluster_9[which(bin_Cluster_9!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_9<-bin_Cluster_9 %*% t(bin_Cluster_9) #adjacency matrix - 
  Cluster_9_entropy_result<- entropy(hyp_Cluster_9)
  Cluster_9_result_df <- rbind(Cluster_9_result_df, Cluster_9_entropy_result)
}
write.csv(Cluster_9_result_df, "path_to_file/Acute_cluster_9_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_10 ####
set.seed(3)

Cluster_10 <- spatial_obs[spatial_obs$cell_type == "10", ]
Cluster_10_X <- Acute_samples_X[na.omit(match(Cluster_10$X,rownames(Acute_samples_X))),]

sd.scores10<-apply(Cluster_10_X,2,sd) #2 = selects columns not rows
Cluster_10_X_sd<-Cluster_10_X[,which(sd.scores10>0)]


Cluster_10_result_df <- data.frame()
cor_data_Cluster_10<-cor(Cluster_10_X_sd, Cluster_10_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_10_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_10_DEGs<- colnames(Cluster_10_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_10_2 <- cor_data_Cluster_10[na.omit(match(Cluster_10_DEGs,rownames(cor_data_Cluster_10))),-na.omit(match(Cluster_10_DEGs,colnames(cor_data_Cluster_10)))] #na omit the ones that arent in my selected 100
  Cluster_10_sub_sample <- colnames(cor_data_Cluster_10_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_10_2 <- cor_data_Cluster_10_2[,na.omit(match(Cluster_10_sub_sample,colnames(cor_data_Cluster_10_2)))] #take all the rows and sample the columns
  bin_Cluster_10<-abs(cor_data_Cluster_10_2)
  bin_Cluster_10[which(bin_Cluster_10>sd(cor_data_Cluster_10_2))]<-1
  bin_Cluster_10[which(bin_Cluster_10!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_10<-bin_Cluster_10 %*% t(bin_Cluster_10) #adjacency matrix - 
  Cluster_10_entropy_result<- entropy(hyp_Cluster_10)
  Cluster_10_result_df <- rbind(Cluster_10_result_df, Cluster_10_entropy_result)
}
write.csv(Cluster_10_result_df, "path_to_file/Acute_cluster_10_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_11 <90 SPOTS SO EXCLUDE#### 

#### Cluster_12 ####
set.seed(3)
Cluster_12 <- spatial_obs[spatial_obs$cell_type == "12", ]
Cluster_12_X <- Acute_samples_X[na.omit(match(Cluster_12$X,rownames(Acute_samples_X))),]

sd.scores12<-apply(Cluster_12_X,2,sd) #2 = selects columns not rows
Cluster_12_X_sd<-Cluster_12_X[,which(sd.scores12>0)]


Cluster_12_result_df <- data.frame()
cor_data_Cluster_12<-cor(Cluster_12_X_sd, Cluster_12_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_12_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_12_DEGs<- colnames(Cluster_12_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_12_2 <- cor_data_Cluster_12[na.omit(match(Cluster_12_DEGs,rownames(cor_data_Cluster_12))),-na.omit(match(Cluster_12_DEGs,colnames(cor_data_Cluster_12)))] #na omit the ones that arent in my selected 100
  Cluster_12_sub_sample <- colnames(cor_data_Cluster_12_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_12_2 <- cor_data_Cluster_12_2[,na.omit(match(Cluster_12_sub_sample,colnames(cor_data_Cluster_12_2)))] #take all the rows and sample the columns
  bin_Cluster_12<-abs(cor_data_Cluster_12_2)
  bin_Cluster_12[which(bin_Cluster_12>sd(cor_data_Cluster_12_2))]<-1
  bin_Cluster_12[which(bin_Cluster_12!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_12<-bin_Cluster_12 %*% t(bin_Cluster_12) #adjacency matrix - 
  Cluster_12_entropy_result<- entropy(hyp_Cluster_12)
  Cluster_12_result_df <- rbind(Cluster_12_result_df, Cluster_12_entropy_result)
}
write.csv(Cluster_12_result_df, "path_to_file/Acute_cluster_12_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_13 ####
set.seed(3)

Cluster_13 <- spatial_obs[spatial_obs$cell_type == "13", ]
Cluster_13_X <- Acute_samples_X[na.omit(match(Cluster_13$X,rownames(Acute_samples_X))),]

sd.scores13<-apply(Cluster_13_X,2,sd) #2 = selects columns not rows
Cluster_13_X_sd<-Cluster_13_X[,which(sd.scores13>0)]

Cluster_13_result_df <- data.frame()
cor_data_Cluster_13<-cor(Cluster_13_X_sd, Cluster_13_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_13_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_13_DEGs<- colnames(Cluster_13_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_13_2 <- cor_data_Cluster_13[na.omit(match(Cluster_13_DEGs,rownames(cor_data_Cluster_13))),-na.omit(match(Cluster_13_DEGs,colnames(cor_data_Cluster_13)))] #na omit the ones that arent in my selected 100
  Cluster_13_sub_sample <- colnames(cor_data_Cluster_13_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_13_2 <- cor_data_Cluster_13_2[,na.omit(match(Cluster_13_sub_sample,colnames(cor_data_Cluster_13_2)))] #take all the rows and sample the columns
  bin_Cluster_13<-abs(cor_data_Cluster_13_2)
  bin_Cluster_13[which(bin_Cluster_13>sd(cor_data_Cluster_13_2))]<-1
  bin_Cluster_13[which(bin_Cluster_13!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_13<-bin_Cluster_13 %*% t(bin_Cluster_13) #adjacency matrix - 
  Cluster_13_entropy_result<- entropy(hyp_Cluster_13)
  Cluster_13_result_df <- rbind(Cluster_13_result_df, Cluster_13_entropy_result)
}
write.csv(Cluster_13_result_df, "path_to_file/Acute_cluster_13_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_14 ####
set.seed(3)
Cluster_14 <- spatial_obs[spatial_obs$cell_type == "14", ]
Cluster_14_X <- Acute_samples_X[na.omit(match(Cluster_14$X,rownames(Acute_samples_X))),]

sd.scores14<-apply(Cluster_14_X,2,sd) #2 = selects columns not rows
Cluster_14_X_sd<-Cluster_14_X[,which(sd.scores14>0)]


Cluster_14_result_df <- data.frame()
cor_data_Cluster_14<-cor(Cluster_14_X_sd, Cluster_14_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_14_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_14_DEGs<- colnames(Cluster_14_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_14_2 <- cor_data_Cluster_14[na.omit(match(Cluster_14_DEGs,rownames(cor_data_Cluster_14))),-na.omit(match(Cluster_14_DEGs,colnames(cor_data_Cluster_14)))] #na omit the ones that arent in my selected 100
  Cluster_14_sub_sample <- colnames(cor_data_Cluster_14_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_14_2 <- cor_data_Cluster_14_2[,na.omit(match(Cluster_14_sub_sample,colnames(cor_data_Cluster_14_2)))] #take all the rows and sample the columns
  bin_Cluster_14<-abs(cor_data_Cluster_14_2)
  bin_Cluster_14[which(bin_Cluster_14>sd(cor_data_Cluster_14_2))]<-1
  bin_Cluster_14[which(bin_Cluster_14!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_14<-bin_Cluster_14 %*% t(bin_Cluster_14) #adjacency matrix - 
  Cluster_14_entropy_result<- entropy(hyp_Cluster_14)
  Cluster_14_result_df <- rbind(Cluster_14_result_df, Cluster_14_entropy_result)
}
write.csv(Cluster_14_result_df, "path_to_file/Acute_cluster_14_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_15 ####
set.seed(3)
Cluster_15 <- spatial_obs[spatial_obs$cell_type == "15", ]
Cluster_15_X <- Acute_samples_X[na.omit(match(Cluster_15$X,rownames(Acute_samples_X))),]

sd.scores15<-apply(Cluster_15_X,2,sd) #2 = selects columns not rows
Cluster_15_X_sd<-Cluster_15_X[,which(sd.scores15>0)]


Cluster_15_result_df <- data.frame()
cor_data_Cluster_15<-cor(Cluster_15_X_sd, Cluster_15_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_15_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_15_DEGs<- colnames(Cluster_15_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_15_2 <- cor_data_Cluster_15[na.omit(match(Cluster_15_DEGs,rownames(cor_data_Cluster_15))),-na.omit(match(Cluster_15_DEGs,colnames(cor_data_Cluster_15)))] #na omit the ones that arent in my selected 100
  Cluster_15_sub_sample <- colnames(cor_data_Cluster_15_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_15_2 <- cor_data_Cluster_15_2[,na.omit(match(Cluster_15_sub_sample,colnames(cor_data_Cluster_15_2)))] #take all the rows and sample the columns
  bin_Cluster_15<-abs(cor_data_Cluster_15_2)
  bin_Cluster_15[which(bin_Cluster_15>sd(cor_data_Cluster_15_2))]<-1
  bin_Cluster_15[which(bin_Cluster_15!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_15<-bin_Cluster_15 %*% t(bin_Cluster_15) #adjacency matrix - 
  Cluster_15_entropy_result<- entropy(hyp_Cluster_15)
  Cluster_15_result_df <- rbind(Cluster_15_result_df, Cluster_15_entropy_result)
}
write.csv(Cluster_15_result_df, "path_to_file/Acute_cluster_15_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_16 ####
set.seed(3)
Cluster_16 <- spatial_obs[spatial_obs$cell_type == "16", ]
Cluster_16_X <- Acute_samples_X[na.omit(match(Cluster_16$X,rownames(Acute_samples_X))),]

sd.scores16<-apply(Cluster_16_X,2,sd) #2 = selects columns not rows
Cluster_16_X_sd<-Cluster_16_X[,which(sd.scores16>0)]


Cluster_16_result_df <- data.frame()
cor_data_Cluster_16<-cor(Cluster_16_X_sd, Cluster_16_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_16_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_16_DEGs<- colnames(Cluster_16_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_16_2 <- cor_data_Cluster_16[na.omit(match(Cluster_16_DEGs,rownames(cor_data_Cluster_16))),-na.omit(match(Cluster_16_DEGs,colnames(cor_data_Cluster_16)))] #na omit the ones that arent in my selected 100
  Cluster_16_sub_sample <- colnames(cor_data_Cluster_16_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_16_2 <- cor_data_Cluster_16_2[,na.omit(match(Cluster_16_sub_sample,colnames(cor_data_Cluster_16_2)))] #take all the rows and sample the columns
  bin_Cluster_16<-abs(cor_data_Cluster_16_2)
  bin_Cluster_16[which(bin_Cluster_16>sd(cor_data_Cluster_16_2))]<-1
  bin_Cluster_16[which(bin_Cluster_16!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_16<-bin_Cluster_16 %*% t(bin_Cluster_16) #adjacency matrix - 
  Cluster_16_entropy_result<- entropy(hyp_Cluster_16)
  Cluster_16_result_df <- rbind(Cluster_16_result_df, Cluster_16_entropy_result)
}
write.csv(Cluster_16_result_df, "path_to_file/Acute_cluster_16_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_17 ####
set.seed(3)

Cluster_17 <- spatial_obs[spatial_obs$cell_type == "17", ]
Cluster_17_X <- Acute_samples_X[na.omit(match(Cluster_17$X,rownames(Acute_samples_X))),]

sd.scores17<-apply(Cluster_17_X,2,sd) #2 = selects columns not rows
Cluster_17_X_sd<-Cluster_17_X[,which(sd.scores17>0)]

Cluster_17_result_df <- data.frame()
cor_data_Cluster_17<-cor(Cluster_17_X_sd, Cluster_17_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_17_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_17_DEGs<- colnames(Cluster_17_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_17_2 <- cor_data_Cluster_17[na.omit(match(Cluster_17_DEGs,rownames(cor_data_Cluster_17))),-na.omit(match(Cluster_17_DEGs,colnames(cor_data_Cluster_17)))] #na omit the ones that arent in my selected 100
  Cluster_17_sub_sample <- colnames(cor_data_Cluster_17_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_17_2 <- cor_data_Cluster_17_2[,na.omit(match(Cluster_17_sub_sample,colnames(cor_data_Cluster_17_2)))] #take all the rows and sample the columns
  bin_Cluster_17<-abs(cor_data_Cluster_17_2)
  bin_Cluster_17[which(bin_Cluster_17>sd(cor_data_Cluster_17_2))]<-1
  bin_Cluster_17[which(bin_Cluster_17!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_17<-bin_Cluster_17 %*% t(bin_Cluster_17) #adjacency matrix - 
  Cluster_17_entropy_result<- entropy(hyp_Cluster_17)
  Cluster_17_result_df <- rbind(Cluster_17_result_df, Cluster_17_entropy_result)
}
write.csv(Cluster_17_result_df, "path_to_file/Acute_cluster_17_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_18 #### 
set.seed(3)
Cluster_18 <- spatial_obs[spatial_obs$cell_type == "18", ]
Cluster_18_X <- Acute_samples_X[na.omit(match(Cluster_18$X,rownames(Acute_samples_X))),]


sd.scores18<-apply(Cluster_18_X,2,sd) #2 = selects columns not rows
Cluster_18_X_sd<-Cluster_18_X[,which(sd.scores18>0)]

Cluster_18_result_df <- data.frame()
cor_data_Cluster_18<-cor(Cluster_18_X_sd, Cluster_18_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_18_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_18_DEGs<- colnames(Cluster_18_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_18_2 <- cor_data_Cluster_18[na.omit(match(Cluster_18_DEGs,rownames(cor_data_Cluster_18))),-na.omit(match(Cluster_18_DEGs,colnames(cor_data_Cluster_18)))] #na omit the ones that arent in my selected 100
  Cluster_18_sub_sample <- colnames(cor_data_Cluster_18_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_18_2 <- cor_data_Cluster_18_2[,na.omit(match(Cluster_18_sub_sample,colnames(cor_data_Cluster_18_2)))] #take all the rows and sample the columns
  bin_Cluster_18<-abs(cor_data_Cluster_18_2)
  bin_Cluster_18[which(bin_Cluster_18>sd(cor_data_Cluster_18_2))]<-1
  bin_Cluster_18[which(bin_Cluster_18!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_18<-bin_Cluster_18 %*% t(bin_Cluster_18) #adjacency matrix - 
  Cluster_18_entropy_result<- entropy(hyp_Cluster_18)
  Cluster_18_result_df <- rbind(Cluster_18_result_df, Cluster_18_entropy_result)
}
write.csv(Cluster_18_result_df, "path_to_file/Acute_cluster_18_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_19 <90 SPOTS SO EXCLUDE####

#### Cluster_20 ####
set.seed(3)
Cluster_20 <- spatial_obs[spatial_obs$cell_type == "20", ]
Cluster_20_X <- Acute_samples_X[na.omit(match(Cluster_20$X,rownames(Acute_samples_X))),]


sd.scores20<-apply(Cluster_20_X,2,sd) #2 = selects columns not rows
Cluster_20_X_sd<-Cluster_20_X[,which(sd.scores20>0)]


Cluster_20_result_df <- data.frame()
cor_data_Cluster_20<-cor(Cluster_20_X_sd, Cluster_20_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_20_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_20_DEGs<- colnames(Cluster_20_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_20_2 <- cor_data_Cluster_20[na.omit(match(Cluster_20_DEGs,rownames(cor_data_Cluster_20))),-na.omit(match(Cluster_20_DEGs,colnames(cor_data_Cluster_20)))] #na omit the ones that arent in my selected 100
  Cluster_20_sub_sample <- colnames(cor_data_Cluster_20_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_20_2 <- cor_data_Cluster_20_2[,na.omit(match(Cluster_20_sub_sample,colnames(cor_data_Cluster_20_2)))] #take all the rows and sample the columns
  bin_Cluster_20<-abs(cor_data_Cluster_20_2)
  bin_Cluster_20[which(bin_Cluster_20>sd(cor_data_Cluster_20_2))]<-1
  bin_Cluster_20[which(bin_Cluster_20!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_20<-bin_Cluster_20 %*% t(bin_Cluster_20) #adjacency matrix - 
  Cluster_20_entropy_result<- entropy(hyp_Cluster_20)
  Cluster_20_result_df <- rbind(Cluster_20_result_df, Cluster_20_entropy_result)
}
write.csv(Cluster_20_result_df, "path_to_file/Acute_cluster_20_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_21 ####
set.seed(3)
Cluster_21 <- spatial_obs[spatial_obs$cell_type == "21", ]
Cluster_21_X <- Acute_samples_X[na.omit(match(Cluster_21$X,rownames(Acute_samples_X))),]

sd.scores21<-apply(Cluster_21_X,2,sd) #2 = selects columns not rows
Cluster_21_X_sd<-Cluster_21_X[,which(sd.scores21>0)]


Cluster_21_result_df <- data.frame()
cor_data_Cluster_21<-cor(Cluster_21_X_sd, Cluster_21_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_21_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_21_DEGs<- colnames(Cluster_21_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_21_2 <- cor_data_Cluster_21[na.omit(match(Cluster_21_DEGs,rownames(cor_data_Cluster_21))),-na.omit(match(Cluster_21_DEGs,colnames(cor_data_Cluster_21)))] #na omit the ones that arent in my selected 100
  Cluster_21_sub_sample <- colnames(cor_data_Cluster_21_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_21_2 <- cor_data_Cluster_21_2[,na.omit(match(Cluster_21_sub_sample,colnames(cor_data_Cluster_21_2)))] #take all the rows and sample the columns
  bin_Cluster_21<-abs(cor_data_Cluster_21_2)
  bin_Cluster_21[which(bin_Cluster_21>sd(cor_data_Cluster_21_2))]<-1
  bin_Cluster_21[which(bin_Cluster_21!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_21<-bin_Cluster_21 %*% t(bin_Cluster_21) #adjacency matrix - 
  Cluster_21_entropy_result<- entropy(hyp_Cluster_21)
  Cluster_21_result_df <- rbind(Cluster_21_result_df, Cluster_21_entropy_result)
}
write.csv(Cluster_21_result_df, "path_to_file/Acute_cluster_21_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_22 ####
set.seed(3)
Cluster_22 <- spatial_obs[spatial_obs$cell_type == "22", ]
Cluster_22_X <- Acute_samples_X[na.omit(match(Cluster_22$X,rownames(Acute_samples_X))),]

sd.scores22<-apply(Cluster_22_X,2,sd) #2 = selects columns not rows
Cluster_22_X_sd<-Cluster_22_X[,which(sd.scores22>0)]


Cluster_22_result_df <- data.frame()
cor_data_Cluster_22<-cor(Cluster_22_X_sd, Cluster_22_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_22_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_22_DEGs<- colnames(Cluster_22_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_22_2 <- cor_data_Cluster_22[na.omit(match(Cluster_22_DEGs,rownames(cor_data_Cluster_22))),-na.omit(match(Cluster_22_DEGs,colnames(cor_data_Cluster_22)))] #na omit the ones that arent in my selected 100
  Cluster_22_sub_sample <- colnames(cor_data_Cluster_22_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_22_2 <- cor_data_Cluster_22_2[,na.omit(match(Cluster_22_sub_sample,colnames(cor_data_Cluster_22_2)))] #take all the rows and sample the columns
  bin_Cluster_22<-abs(cor_data_Cluster_22_2)
  bin_Cluster_22[which(bin_Cluster_22>sd(cor_data_Cluster_22_2))]<-1
  bin_Cluster_22[which(bin_Cluster_22!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_22<-bin_Cluster_22 %*% t(bin_Cluster_22) #adjacency matrix - 
  Cluster_22_entropy_result<- entropy(hyp_Cluster_22)
  Cluster_22_result_df <- rbind(Cluster_22_result_df, Cluster_22_entropy_result)
}
write.csv(Cluster_22_result_df, "path_to_file/Acute_cluster_22_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_23 <90 SPOTS SO EXCLUDE####

#### Cluster_24 <90 SPOTS SO EXCLUDE####

#### Cluster_25 NOT PRESENT IN ACUTE ####

#### Cluster_26 NOT PRESENT IN ACUTE ####

#### Cluster_27 NOT PRESENT IN ACUTE ####

#### Cluster_28 <90 SPOTS SO EXCLUDE####

##### Combining ACUTE data - EXCLUDING CLUSTERS FOR <90 SPOTS #####
##create new df
All_acute_clusters_results_df <- data.frame()

Cluster_0_result_df <- read.csv("path_to_file/Acute_cluster_0_entropy_results_10000.csv",header=T, sep = ',')
Cluster_1_result_df <- read.csv("path_to_files/Acute_cluster_1_entropy_results_10000.csv",header=T, sep = ',')
Cluster_2_result_df <- read.csv("path_to_files/Acute_cluster_2_entropy_results_10000.csv",header=T, sep = ',')
Cluster_3_result_df <- read.csv("path_to_files/Acute_cluster_3_entropy_results_10000.csv",header=T, sep = ',')
Cluster_4_result_df <- read.csv("path_to_files/Acute_cluster_4_entropy_results_10000.csv",header=T, sep = ',')
Cluster_5_result_df <- read.csv("path_to_files/Acute_cluster_5_entropy_results_10000.csv",header=T, sep = ',')
Cluster_6_result_df <- read.csv("path_to_files/Acute_cluster_6_entropy_results_10000.csv",header=T, sep = ',')
Cluster_7_result_df <- read.csv("path_to_files/Acute_cluster_7_entropy_results_10000.csv",header=T, sep = ',')
Cluster_8_result_df <- read.csv("path_to_files/Acute_cluster_8_entropy_results_10000.csv",header=T, sep = ',')
Cluster_9_result_df <- read.csv("path_to_files/Acute_cluster_9_entropy_results_10000.csv",header=T, sep = ',')
Cluster_10_result_df <- read.csv("path_to_files/Acute_cluster_10_entropy_results_10000.csv",header=T, sep = ',')
Cluster_12_result_df <- read.csv("path_to_files/Acute_cluster_12_entropy_results_10000.csv",header=T, sep = ',')
Cluster_13_result_df <- read.csv("path_to_files/Acute_cluster_13_entropy_results_10000.csv",header=T, sep = ',')
Cluster_14_result_df <- read.csv("path_to_files/Acute_cluster_14_entropy_results_10000.csv",header=T, sep = ',')
Cluster_15_result_df <- read.csv("path_to_files/Acute_cluster_15_entropy_results_10000.csv",header=T, sep = ',')
Cluster_16_result_df <- read.csv("path_to_files/Acute_cluster_16_entropy_results_10000.csv",header=T, sep = ',')
Cluster_17_result_df <- read.csv("path_to_files/Acute_cluster_17_entropy_results_10000.csv",header=T, sep = ',')
Cluster_18_result_df <- read.csv("path_to_files/Acute_cluster_18_entropy_results_10000.csv",header=T, sep = ',')
Cluster_20_result_df <- read.csv("path_to_files/Acute_cluster_20_entropy_results_10000.csv",header=T, sep = ',')
Cluster_21_result_df <- read.csv("path_to_files/Acute_cluster_21_entropy_results_10000.csv",header=T, sep = ',')
Cluster_22_result_df <- read.csv("path_to_files/Acute_cluster_22_entropy_results_10000.csv",header=T, sep = ',')

# change first column name:
colnames(Cluster_0_result_df)[1] <- "Entropy"
colnames(Cluster_1_result_df)[1] <- "Entropy"
colnames(Cluster_2_result_df)[1] <- "Entropy"
colnames(Cluster_3_result_df)[1] <- "Entropy"
colnames(Cluster_4_result_df)[1] <- "Entropy"
colnames(Cluster_5_result_df)[1] <- "Entropy"
colnames(Cluster_6_result_df)[1] <- "Entropy"
colnames(Cluster_7_result_df)[1] <- "Entropy"
colnames(Cluster_8_result_df)[1] <- "Entropy"
colnames(Cluster_9_result_df)[1] <- "Entropy"
colnames(Cluster_10_result_df)[1] <- "Entropy"
colnames(Cluster_12_result_df)[1] <- "Entropy"
colnames(Cluster_13_result_df)[1] <- "Entropy"
colnames(Cluster_14_result_df)[1] <- "Entropy"
colnames(Cluster_15_result_df)[1] <- "Entropy"
colnames(Cluster_16_result_df)[1] <- "Entropy"
colnames(Cluster_17_result_df)[1] <- "Entropy"
colnames(Cluster_18_result_df)[1] <- "Entropy"
colnames(Cluster_20_result_df)[1] <- "Entropy"
colnames(Cluster_21_result_df)[1] <- "Entropy"
colnames(Cluster_22_result_df)[1] <- "Entropy"


##need to add columns and names to each df

All_acute_clusters_results_df <- rbind(Cluster_0_result_df,Cluster_1_result_df,Cluster_2_result_df,Cluster_3_result_df,Cluster_4_result_df,Cluster_5_result_df,Cluster_6_result_df,Cluster_7_result_df,Cluster_8_result_df,Cluster_9_result_df,Cluster_10_result_df,Cluster_12_result_df,
                                 Cluster_13_result_df,Cluster_14_result_df,Cluster_15_result_df,Cluster_16_result_df,Cluster_17_result_df,Cluster_18_result_df,Cluster_20_result_df,Cluster_21_result_df,Cluster_22_result_df)

All_acute_clusters_results_df$Iteration <- rep(1:10000, times = 21) # adds 1-1000, 29 times
#All_acute_clusters_results_df$Cluster <- rep(0:24,28, each = 1000) # add 0-27, 1000 times each

# Create a vector with 0:24 repeated 1000 times
zero_values <- rep("Adipocytes1", each = 10000)
one_values <- rep("FBs/SMCs1", each = 10000)
two_values <- rep("Keratinocytes/FBs1", each = 10000)
three_values<- rep("Keratinocytes1", each = 10000)
four_values <- rep("FBs1", each = 10000)
five_values <- rep("FBs/SMCs/Keratinocytes", each = 10000)
six_values<- rep("Keratinocytes2", each = 10000)
seven_values <- rep("FBs2", each = 10000)
eight_values <- rep("FBs3", each = 10000)
nine_values  <- rep("FBs4", each = 10000)
ten_values<- rep("FBs/SMCs2", each = 10000)
twelve_values <- rep("Keratinocytes3", each = 10000)
thirteen_values <- rep("SMCs", each = 10000)
fourteen_values <- rep("Keratinocytes4", each = 10000)
fifteen_values <- rep("Keratinocytes5", each = 10000)
sixteen_values<- rep("Adipocytes2", each = 10000)
seventeen_values <- rep("Keratinocytes6", each = 10000)
eighteen_values  <- rep("FBs5", each = 10000)
twenty_values <- rep("Keratinocytes7", each = 10000)
twentyone_values <- rep("Keratinocytes8", each = 10000)
twentytwo_values <- rep("Keratinocytes9", each = 10000)


# Combine both vectors
combined_values <- c(zero_values,one_values,two_values,three_values,four_values,five_values,six_values,seven_values,eight_values,
                     nine_values,ten_values,twelve_values,thirteen_values,fourteen_values,fifteen_values,
                     sixteen_values,seventeen_values,eighteen_values,twenty_values,twentyone_values,twentytwo_values)

# Assign the combined values to your data frame column
All_acute_clusters_results_df$Cluster <- combined_values




All_acute_clusters_results_df$Cluster<- as.factor(All_acute_clusters_results_df$Cluster)

write.csv(All_acute_clusters_results_df, "path_to_files/Final_acute_clusters_entropy_results_10000_Nov.csv", row.names=FALSE)

#All_acute_clusters_results_df <- read.csv("path_to_files/Final_acute_clusters_entropy_results_10000.csv")
All_acute_clusters_results_df$Cluster <- as.factor(All_acute_clusters_results_df$Cluster)

#assign cluster colours
cluster_colours<- c('Adipocytes1'='#F8766D',"FBs/SMCs1"= "#EF7F48","Keratinocytes/FBs1"= "#E48800","Keratinocytes1"= "#D69100","FBs1"= "#C69900","FBs/SMCs/Keratinocytes"= "#B4A000",
                    "Keratinocytes2"= "#9EA700","FBs2"="#83AD00","FBs3"= "#60B200","FBs4"= "#19B700","FBs/SMCs2"= "#00BB48","B-cells/FBs"= "#00BE6B","Keratinocytes3"= "#00C088","SMCs"= "#00C1A2",
                    "Keratinocytes4"="#00C0BA","Keratinocytes5"= "#00BECF","Adipocytes2"= "#00B9E1","Keratinocytes6"= "#00B3F1","FBs5"= "#00ABFD","Adipocytes3"= "#46A0FF","Keratinocytes7"= "#8794FF",
                    "Keratinocytes8"="#B086FF", "Keratinocytes9"="#CE79FF","Keratinocytes10"= "#E46DF6","Keratinocytes11"= "#F365E6","Keratinocytes/FBs2"= "#FC61D4","FBs6"= "#FF62BE",
                    "Keratinocytes/FBs3"="#FF67A6","Endothelial/Langerhans/Melanocytes/SMCs"="#FE6E8B")



tiff("path_to_file/Final_acute_clusters_entropy_boxplots_10000_iterations_named_clusters.tiff", height = 5000, width = 9000, res = 500)

#plot boxplot of acute cluster entropy values
p1 <- ggplot(data=All_acute_clusters_results_df, aes(y=factor(Cluster), x=Entropy)) +
  geom_boxplot(aes(color=Cluster))+
  coord_flip() +
  labs(title = "Entropy of Clusters in Acute Samples",
       x = "Entropy",
       y = "Cluster") +
  theme(
    axis.title = element_text(size = 14),  # Adjust axis title size
    axis.text = element_text(size = 12),   # Adjust axis text size
    legend.title = element_text(size = 14),# Adjust legend title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
p1 + scale_color_manual(values=c(cluster_colours)) + scale_x_continuous(breaks=seq(8.2, 9.2, 0.2), limits = c(8.2, 9.2))
dev.off()




### CHRONIC SAMPLES ####
Chronic_samples <- spatial_obs[spatial_obs$type == "Chronic", ]

Chronic_samples_X <- exp_data[na.omit(match(Chronic_samples$X,rownames(exp_data))),]

##### CHRONIC SAMPLES #####
#### Cluster_0 ####
set.seed(3)
Cluster_0 <- spatial_obs[spatial_obs$cell_type == "0", ]
Cluster_0_X <- Chronic_samples_X[na.omit(match(Cluster_0$X,rownames(Chronic_samples_X))),]


sd.scores0<-apply(Cluster_0_X,2,sd) #2 = selects columns not rows
Cluster_0_X_sd<-Cluster_0_X[,which(sd.scores0>0)]

Cluster_0_result_df <- data.frame()
cor_data_Cluster_0<-cor(Cluster_0_X_sd, Cluster_0_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_0_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_0_DEGs<- colnames(Cluster_0_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_0_2 <- cor_data_Cluster_0[na.omit(match(Cluster_0_DEGs,rownames(cor_data_Cluster_0))),-na.omit(match(Cluster_0_DEGs,colnames(cor_data_Cluster_0)))] #na omit the ones that arent in my selected 100
  Cluster_0_sub_sample <- colnames(cor_data_Cluster_0_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_0_2 <- cor_data_Cluster_0_2[,na.omit(match(Cluster_0_sub_sample,colnames(cor_data_Cluster_0_2)))] #take all the rows and sample the columns
  bin_Cluster_0<-abs(cor_data_Cluster_0_2)
  bin_Cluster_0[which(bin_Cluster_0>sd(cor_data_Cluster_0_2))]<-1
  bin_Cluster_0[which(bin_Cluster_0!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_0<-bin_Cluster_0 %*% t(bin_Cluster_0) #adjacency matrix - 
  Cluster_0_entropy_result<- entropy(hyp_Cluster_0)
  Cluster_0_result_df <- rbind(Cluster_0_result_df, Cluster_0_entropy_result)
}
write.csv(Cluster_0_result_df, "path_to_file/Chronic_cluster_0_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_1 ####
set.seed(3)
Cluster_1 <- spatial_obs[spatial_obs$cell_type == "1", ]
Cluster_1_X <- Chronic_samples_X[na.omit(match(Cluster_1$X,rownames(Chronic_samples_X))),]

sd.scores1<-apply(Cluster_1_X,2,sd) #2 = selects columns not rows
Cluster_1_X_sd<-Cluster_1_X[,which(sd.scores1>0)]


Cluster_1_result_df <- data.frame()
cor_data_Cluster_1<-cor(Cluster_1_X_sd, Cluster_1_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_1_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_1_DEGs<- colnames(Cluster_1_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_1_2 <- cor_data_Cluster_1[na.omit(match(Cluster_1_DEGs,rownames(cor_data_Cluster_1))),-na.omit(match(Cluster_1_DEGs,colnames(cor_data_Cluster_1)))] #na omit the ones that arent in my selected 100
  Cluster_1_sub_sample <- colnames(cor_data_Cluster_1_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_1_2 <- cor_data_Cluster_1_2[,na.omit(match(Cluster_1_sub_sample,colnames(cor_data_Cluster_1_2)))] #take all the rows and sample the columns
  bin_Cluster_1<-abs(cor_data_Cluster_1_2)
  bin_Cluster_1[which(bin_Cluster_1>sd(cor_data_Cluster_1_2))]<-1
  bin_Cluster_1[which(bin_Cluster_1!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_1<-bin_Cluster_1 %*% t(bin_Cluster_1) #adjacency matrix - 
  Cluster_1_entropy_result<- entropy(hyp_Cluster_1)
  Cluster_1_result_df <- rbind(Cluster_1_result_df, Cluster_1_entropy_result)
}
write.csv(Cluster_1_result_df, "path_to_file/Chronic_cluster_1_entropy_results_10000.csv", row.names=FALSE)


#### Cluster_2 ####
set.seed(3)
Cluster_2 <- spatial_obs[spatial_obs$cell_type == "2", ]
Cluster_2_X <- Chronic_samples_X[na.omit(match(Cluster_2$X,rownames(Chronic_samples_X))),]

sd.scores2<-apply(Cluster_2_X,2,sd) #2 = selects columns not rows
Cluster_2_X_sd<-Cluster_2_X[,which(sd.scores2>0)]

Cluster_2_result_df <- data.frame()
cor_data_Cluster_2<-cor(Cluster_2_X_sd, Cluster_2_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_2_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_2_DEGs<- colnames(Cluster_2_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_2_2 <- cor_data_Cluster_2[na.omit(match(Cluster_2_DEGs,rownames(cor_data_Cluster_2))),-na.omit(match(Cluster_2_DEGs,colnames(cor_data_Cluster_2)))] #na omit the ones that arent in my selected 100
  Cluster_2_sub_sample <- colnames(cor_data_Cluster_2_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_2_2 <- cor_data_Cluster_2_2[,na.omit(match(Cluster_2_sub_sample,colnames(cor_data_Cluster_2_2)))] #take all the rows and sample the columns
  bin_Cluster_2<-abs(cor_data_Cluster_2_2)
  bin_Cluster_2[which(bin_Cluster_2>sd(cor_data_Cluster_2_2))]<-1
  bin_Cluster_2[which(bin_Cluster_2!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_2<-bin_Cluster_2 %*% t(bin_Cluster_2) #adjacency matrix - 
  Cluster_2_entropy_result<- entropy(hyp_Cluster_2)
  Cluster_2_result_df <- rbind(Cluster_2_result_df, Cluster_2_entropy_result)
}
write.csv(Cluster_2_result_df, "path_to_file/Chronic_cluster_2_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_3 <90 SPOTS SO EXCLUDE####

#### Cluster_4 <90 SPOTS SO EXCLUDE####

#### Cluster_5 ####
set.seed(3)
Cluster_5 <- spatial_obs[spatial_obs$cell_type == "5", ]
Cluster_5_X <- Chronic_samples_X[na.omit(match(Cluster_5$X,rownames(Chronic_samples_X))),]

sd.scores5<-apply(Cluster_5_X,2,sd) #2 = selects columns not rows
Cluster_5_X_sd<-Cluster_5_X[,which(sd.scores5>0)]


Cluster_5_result_df <- data.frame()
cor_data_Cluster_5<-cor(Cluster_5_X_sd, Cluster_5_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_5_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_5_DEGs<- colnames(Cluster_5_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_5_2 <- cor_data_Cluster_5[na.omit(match(Cluster_5_DEGs,rownames(cor_data_Cluster_5))),-na.omit(match(Cluster_5_DEGs,colnames(cor_data_Cluster_5)))] #na omit the ones that arent in my selected 100
  Cluster_5_sub_sample <- colnames(cor_data_Cluster_5_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_5_2 <- cor_data_Cluster_5_2[,na.omit(match(Cluster_5_sub_sample,colnames(cor_data_Cluster_5_2)))] #take all the rows and sample the columns
  bin_Cluster_5<-abs(cor_data_Cluster_5_2)
  bin_Cluster_5[which(bin_Cluster_5>sd(cor_data_Cluster_5_2))]<-1
  bin_Cluster_5[which(bin_Cluster_5!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_5<-bin_Cluster_5 %*% t(bin_Cluster_5) #adjacency matrix - 
  Cluster_5_entropy_result<- entropy(hyp_Cluster_5)
  Cluster_5_result_df <- rbind(Cluster_5_result_df, Cluster_5_entropy_result)
}
write.csv(Cluster_5_result_df, "path_to_file/Chronic_cluster_5_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_6 ####
set.seed(3)
Cluster_6 <- spatial_obs[spatial_obs$cell_type == "6", ]
Cluster_6_X <- Chronic_samples_X[na.omit(match(Cluster_6$X,rownames(Chronic_samples_X))),]

sd.scores6<-apply(Cluster_6_X,2,sd) #2 = selects columns not rows
Cluster_6_X_sd<-Cluster_6_X[,which(sd.scores6>0)]


Cluster_6_result_df <- data.frame()
cor_data_Cluster_6<-cor(Cluster_6_X_sd, Cluster_6_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_6_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_6_DEGs<- colnames(Cluster_6_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_6_2 <- cor_data_Cluster_6[na.omit(match(Cluster_6_DEGs,rownames(cor_data_Cluster_6))),-na.omit(match(Cluster_6_DEGs,colnames(cor_data_Cluster_6)))] #na omit the ones that arent in my selected 100
  Cluster_6_sub_sample <- colnames(cor_data_Cluster_6_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_6_2 <- cor_data_Cluster_6_2[,na.omit(match(Cluster_6_sub_sample,colnames(cor_data_Cluster_6_2)))] #take all the rows and sample the columns
  bin_Cluster_6<-abs(cor_data_Cluster_6_2)
  bin_Cluster_6[which(bin_Cluster_6>sd(cor_data_Cluster_6_2))]<-1
  bin_Cluster_6[which(bin_Cluster_6!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_6<-bin_Cluster_6 %*% t(bin_Cluster_6) #adjacency matrix - 
  Cluster_6_entropy_result<- entropy(hyp_Cluster_6)
  Cluster_6_result_df <- rbind(Cluster_6_result_df, Cluster_6_entropy_result)
}
write.csv(Cluster_6_result_df, "path_to_file/Chronic_cluster_6_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_7 ####
set.seed(3)
Cluster_7 <- spatial_obs[spatial_obs$cell_type == "7", ]
Cluster_7_X <- Chronic_samples_X[na.omit(match(Cluster_7$X,rownames(Chronic_samples_X))),]

sd.scores7<-apply(Cluster_7_X,2,sd) #2 = selects columns not rows
Cluster_7_X_sd<-Cluster_7_X[,which(sd.scores7>0)]


Cluster_7_result_df <- data.frame()
cor_data_Cluster_7<-cor(Cluster_7_X_sd, Cluster_7_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_7_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_7_DEGs<- colnames(Cluster_7_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_7_2 <- cor_data_Cluster_7[na.omit(match(Cluster_7_DEGs,rownames(cor_data_Cluster_7))),-na.omit(match(Cluster_7_DEGs,colnames(cor_data_Cluster_7)))] #na omit the ones that arent in my selected 100
  Cluster_7_sub_sample <- colnames(cor_data_Cluster_7_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_7_2 <- cor_data_Cluster_7_2[,na.omit(match(Cluster_7_sub_sample,colnames(cor_data_Cluster_7_2)))] #take all the rows and sample the columns
  bin_Cluster_7<-abs(cor_data_Cluster_7_2)
  bin_Cluster_7[which(bin_Cluster_7>sd(cor_data_Cluster_7_2))]<-1
  bin_Cluster_7[which(bin_Cluster_7!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_7<-bin_Cluster_7 %*% t(bin_Cluster_7) #adjacency matrix - 
  Cluster_7_entropy_result<- entropy(hyp_Cluster_7)
  Cluster_7_result_df <- rbind(Cluster_7_result_df, Cluster_7_entropy_result)
}
write.csv(Cluster_7_result_df, "path_to_file/Chronic_cluster_7_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_8 <90 SPOTS SO EXCLUDE ####

#### Cluster_9 ####
set.seed(3)
Cluster_9 <- spatial_obs[spatial_obs$cell_type == "9", ]
Cluster_9_X <- Chronic_samples_X[na.omit(match(Cluster_9$X,rownames(Chronic_samples_X))),]

sd.scores9<-apply(Cluster_9_X,2,sd) #2 = selects columns not rows
Cluster_9_X_sd<-Cluster_9_X[,which(sd.scores9>0)]


Cluster_9_result_df <- data.frame()
cor_data_Cluster_9<-cor(Cluster_9_X_sd, Cluster_9_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_9_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_9_DEGs<- colnames(Cluster_9_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_9_2 <- cor_data_Cluster_9[na.omit(match(Cluster_9_DEGs,rownames(cor_data_Cluster_9))),-na.omit(match(Cluster_9_DEGs,colnames(cor_data_Cluster_9)))] #na omit the ones that arent in my selected 100
  Cluster_9_sub_sample <- colnames(cor_data_Cluster_9_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_9_2 <- cor_data_Cluster_9_2[,na.omit(match(Cluster_9_sub_sample,colnames(cor_data_Cluster_9_2)))] #take all the rows and sample the columns
  bin_Cluster_9<-abs(cor_data_Cluster_9_2)
  bin_Cluster_9[which(bin_Cluster_9>sd(cor_data_Cluster_9_2))]<-1
  bin_Cluster_9[which(bin_Cluster_9!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_9<-bin_Cluster_9 %*% t(bin_Cluster_9) #adjacency matrix - 
  Cluster_9_entropy_result<- entropy(hyp_Cluster_9)
  Cluster_9_result_df <- rbind(Cluster_9_result_df, Cluster_9_entropy_result)
}
write.csv(Cluster_9_result_df, "path_to_file/Chronic_cluster_9_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_10 ####
set.seed(3)

Cluster_10 <- spatial_obs[spatial_obs$cell_type == "10", ]
Cluster_10_X <- Chronic_samples_X[na.omit(match(Cluster_10$X,rownames(Chronic_samples_X))),]

sd.scores10<-apply(Cluster_10_X,2,sd) #2 = selects columns not rows
Cluster_10_X_sd<-Cluster_10_X[,which(sd.scores10>0)]


Cluster_10_result_df <- data.frame()
cor_data_Cluster_10<-cor(Cluster_10_X_sd, Cluster_10_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_10_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_10_DEGs<- colnames(Cluster_10_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_10_2 <- cor_data_Cluster_10[na.omit(match(Cluster_10_DEGs,rownames(cor_data_Cluster_10))),-na.omit(match(Cluster_10_DEGs,colnames(cor_data_Cluster_10)))] #na omit the ones that arent in my selected 100
  Cluster_10_sub_sample <- colnames(cor_data_Cluster_10_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_10_2 <- cor_data_Cluster_10_2[,na.omit(match(Cluster_10_sub_sample,colnames(cor_data_Cluster_10_2)))] #take all the rows and sample the columns
  bin_Cluster_10<-abs(cor_data_Cluster_10_2)
  bin_Cluster_10[which(bin_Cluster_10>sd(cor_data_Cluster_10_2))]<-1
  bin_Cluster_10[which(bin_Cluster_10!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_10<-bin_Cluster_10 %*% t(bin_Cluster_10) #adjacency matrix - 
  Cluster_10_entropy_result<- entropy(hyp_Cluster_10)
  Cluster_10_result_df <- rbind(Cluster_10_result_df, Cluster_10_entropy_result)
}
write.csv(Cluster_10_result_df, "path_to_file/Chronic_cluster_10_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_11 #### 
set.seed(3)
Cluster_11 <- spatial_obs[spatial_obs$cell_type == "11", ]
Cluster_11_X <- Chronic_samples_X[na.omit(match(Cluster_11$X,rownames(Chronic_samples_X))),]

sd.scores11<-apply(Cluster_11_X,2,sd) #2 = selects columns not rows
Cluster_11_X_sd<-Cluster_11_X[,which(sd.scores11>0)]


Cluster_11_result_df <- data.frame()
cor_data_Cluster_11<-cor(Cluster_11_X_sd, Cluster_11_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_11_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_11_DEGs<- colnames(Cluster_11_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_11_2 <- cor_data_Cluster_11[na.omit(match(Cluster_11_DEGs,rownames(cor_data_Cluster_11))),-na.omit(match(Cluster_11_DEGs,colnames(cor_data_Cluster_11)))] #na omit the ones that arent in my selected 100
  Cluster_11_sub_sample <- colnames(cor_data_Cluster_11_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_11_2 <- cor_data_Cluster_11_2[,na.omit(match(Cluster_11_sub_sample,colnames(cor_data_Cluster_11_2)))] #take all the rows and sample the columns
  bin_Cluster_11<-abs(cor_data_Cluster_11_2)
  bin_Cluster_11[which(bin_Cluster_11>sd(cor_data_Cluster_11_2))]<-1
  bin_Cluster_11[which(bin_Cluster_11!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_11<-bin_Cluster_11 %*% t(bin_Cluster_11) #adjacency matrix - 
  Cluster_11_entropy_result<- entropy(hyp_Cluster_11)
  Cluster_11_result_df <- rbind(Cluster_11_result_df, Cluster_11_entropy_result)
}
write.csv(Cluster_11_result_df, "path_to_file/Chronic_cluster_11_entropy_results_10000.csv", row.names=FALSE)


#### Cluster_12 ####
set.seed(3)
Cluster_12 <- spatial_obs[spatial_obs$cell_type == "12", ]
Cluster_12_X <- Chronic_samples_X[na.omit(match(Cluster_12$X,rownames(Chronic_samples_X))),]

sd.scores12<-apply(Cluster_12_X,2,sd) #2 = selects columns not rows
Cluster_12_X_sd<-Cluster_12_X[,which(sd.scores12>0)]


Cluster_12_result_df <- data.frame()
cor_data_Cluster_12<-cor(Cluster_12_X_sd, Cluster_12_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_12_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_12_DEGs<- colnames(Cluster_12_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_12_2 <- cor_data_Cluster_12[na.omit(match(Cluster_12_DEGs,rownames(cor_data_Cluster_12))),-na.omit(match(Cluster_12_DEGs,colnames(cor_data_Cluster_12)))] #na omit the ones that arent in my selected 100
  Cluster_12_sub_sample <- colnames(cor_data_Cluster_12_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_12_2 <- cor_data_Cluster_12_2[,na.omit(match(Cluster_12_sub_sample,colnames(cor_data_Cluster_12_2)))] #take all the rows and sample the columns
  bin_Cluster_12<-abs(cor_data_Cluster_12_2)
  bin_Cluster_12[which(bin_Cluster_12>sd(cor_data_Cluster_12_2))]<-1
  bin_Cluster_12[which(bin_Cluster_12!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_12<-bin_Cluster_12 %*% t(bin_Cluster_12) #adjacency matrix - 
  Cluster_12_entropy_result<- entropy(hyp_Cluster_12)
  Cluster_12_result_df <- rbind(Cluster_12_result_df, Cluster_12_entropy_result)
}
write.csv(Cluster_12_result_df, "path_to_file/Chronic_cluster_12_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_13 ####
set.seed(3)

Cluster_13 <- spatial_obs[spatial_obs$cell_type == "13", ]
Cluster_13_X <- Chronic_samples_X[na.omit(match(Cluster_13$X,rownames(Chronic_samples_X))),]

sd.scores13<-apply(Cluster_13_X,2,sd) #2 = selects columns not rows
Cluster_13_X_sd<-Cluster_13_X[,which(sd.scores13>0)]

Cluster_13_result_df <- data.frame()
cor_data_Cluster_13<-cor(Cluster_13_X_sd, Cluster_13_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_13_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_13_DEGs<- colnames(Cluster_13_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_13_2 <- cor_data_Cluster_13[na.omit(match(Cluster_13_DEGs,rownames(cor_data_Cluster_13))),-na.omit(match(Cluster_13_DEGs,colnames(cor_data_Cluster_13)))] #na omit the ones that arent in my selected 100
  Cluster_13_sub_sample <- colnames(cor_data_Cluster_13_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_13_2 <- cor_data_Cluster_13_2[,na.omit(match(Cluster_13_sub_sample,colnames(cor_data_Cluster_13_2)))] #take all the rows and sample the columns
  bin_Cluster_13<-abs(cor_data_Cluster_13_2)
  bin_Cluster_13[which(bin_Cluster_13>sd(cor_data_Cluster_13_2))]<-1
  bin_Cluster_13[which(bin_Cluster_13!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_13<-bin_Cluster_13 %*% t(bin_Cluster_13) #adjacency matrix - 
  Cluster_13_entropy_result<- entropy(hyp_Cluster_13)
  Cluster_13_result_df <- rbind(Cluster_13_result_df, Cluster_13_entropy_result)
}
write.csv(Cluster_13_result_df, "path_to_file/Chronic_cluster_13_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_14 ####
set.seed(3)
Cluster_14 <- spatial_obs[spatial_obs$cell_type == "14", ]
Cluster_14_X <- Chronic_samples_X[na.omit(match(Cluster_14$X,rownames(Chronic_samples_X))),]

sd.scores14<-apply(Cluster_14_X,2,sd) #2 = selects columns not rows
Cluster_14_X_sd<-Cluster_14_X[,which(sd.scores14>0)]


Cluster_14_result_df <- data.frame()
cor_data_Cluster_14<-cor(Cluster_14_X_sd, Cluster_14_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_14_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_14_DEGs<- colnames(Cluster_14_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_14_2 <- cor_data_Cluster_14[na.omit(match(Cluster_14_DEGs,rownames(cor_data_Cluster_14))),-na.omit(match(Cluster_14_DEGs,colnames(cor_data_Cluster_14)))] #na omit the ones that arent in my selected 100
  Cluster_14_sub_sample <- colnames(cor_data_Cluster_14_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_14_2 <- cor_data_Cluster_14_2[,na.omit(match(Cluster_14_sub_sample,colnames(cor_data_Cluster_14_2)))] #take all the rows and sample the columns
  bin_Cluster_14<-abs(cor_data_Cluster_14_2)
  bin_Cluster_14[which(bin_Cluster_14>sd(cor_data_Cluster_14_2))]<-1
  bin_Cluster_14[which(bin_Cluster_14!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_14<-bin_Cluster_14 %*% t(bin_Cluster_14) #adjacency matrix - 
  Cluster_14_entropy_result<- entropy(hyp_Cluster_14)
  Cluster_14_result_df <- rbind(Cluster_14_result_df, Cluster_14_entropy_result)
}
write.csv(Cluster_14_result_df, "path_to_file/Chronic_cluster_14_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_15 ####
set.seed(3)
Cluster_15 <- spatial_obs[spatial_obs$cell_type == "15", ]
Cluster_15_X <- Chronic_samples_X[na.omit(match(Cluster_15$X,rownames(Chronic_samples_X))),]

sd.scores15<-apply(Cluster_15_X,2,sd) #2 = selects columns not rows
Cluster_15_X_sd<-Cluster_15_X[,which(sd.scores15>0)]


Cluster_15_result_df <- data.frame()
cor_data_Cluster_15<-cor(Cluster_15_X_sd, Cluster_15_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_15_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_15_DEGs<- colnames(Cluster_15_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_15_2 <- cor_data_Cluster_15[na.omit(match(Cluster_15_DEGs,rownames(cor_data_Cluster_15))),-na.omit(match(Cluster_15_DEGs,colnames(cor_data_Cluster_15)))] #na omit the ones that arent in my selected 100
  Cluster_15_sub_sample <- colnames(cor_data_Cluster_15_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_15_2 <- cor_data_Cluster_15_2[,na.omit(match(Cluster_15_sub_sample,colnames(cor_data_Cluster_15_2)))] #take all the rows and sample the columns
  bin_Cluster_15<-abs(cor_data_Cluster_15_2)
  bin_Cluster_15[which(bin_Cluster_15>sd(cor_data_Cluster_15_2))]<-1
  bin_Cluster_15[which(bin_Cluster_15!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_15<-bin_Cluster_15 %*% t(bin_Cluster_15) #adjacency matrix - 
  Cluster_15_entropy_result<- entropy(hyp_Cluster_15)
  Cluster_15_result_df <- rbind(Cluster_15_result_df, Cluster_15_entropy_result)
}
write.csv(Cluster_15_result_df, "path_to_file/Chronic_cluster_15_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_16 <90 SPOTS SO EXCLUDE####

#### Cluster_17 <90 SPOTS SO EXCLUDE####

#### Cluster_18 <90 SPOTS SO EXCLUDE#### 

#### Cluster_19 ####
set.seed(3)
Cluster_19 <- spatial_obs[spatial_obs$cell_type == "19", ]
Cluster_19_X <- Chronic_samples_X[na.omit(match(Cluster_19$X,rownames(Chronic_samples_X))),]

sd.scores19<-apply(Cluster_19_X,2,sd) #2 = selects columns not rows
Cluster_19_X_sd<-Cluster_19_X[,which(sd.scores19>0)]

Cluster_19_result_df <- data.frame()
cor_data_Cluster_19<-cor(Cluster_19_X_sd, Cluster_19_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_19_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_19_DEGs<- colnames(Cluster_19_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_19_2 <- cor_data_Cluster_19[na.omit(match(Cluster_19_DEGs,rownames(cor_data_Cluster_19))),-na.omit(match(Cluster_19_DEGs,colnames(cor_data_Cluster_19)))] #na omit the ones that arent in my selected 100
  Cluster_19_sub_sample <- colnames(cor_data_Cluster_19_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_19_2 <- cor_data_Cluster_19_2[,na.omit(match(Cluster_19_sub_sample,colnames(cor_data_Cluster_19_2)))] #take all the rows and sample the columns
  bin_Cluster_19<-abs(cor_data_Cluster_19_2)
  bin_Cluster_19[which(bin_Cluster_19>sd(cor_data_Cluster_19_2))]<-1
  bin_Cluster_19[which(bin_Cluster_19!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_19<-bin_Cluster_19 %*% t(bin_Cluster_19) #adjacency matrix - 
  Cluster_19_entropy_result<- entropy(hyp_Cluster_19)
  Cluster_19_result_df <- rbind(Cluster_19_result_df, Cluster_19_entropy_result)
}
write.csv(Cluster_19_result_df, "path_to_file/Chronic_cluster_19_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_20 <90 SPOTS SO EXCLUDE####

#### Cluster_21 ####
set.seed(3)
Cluster_21 <- spatial_obs[spatial_obs$cell_type == "21", ]
Cluster_21_X <- Chronic_samples_X[na.omit(match(Cluster_21$X,rownames(Chronic_samples_X))),]

sd.scores21<-apply(Cluster_21_X,2,sd) #2 = selects columns not rows
Cluster_21_X_sd<-Cluster_21_X[,which(sd.scores21>0)]



Cluster_21_result_df <- data.frame()
cor_data_Cluster_21<-cor(Cluster_21_X_sd, Cluster_21_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_21_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_21_DEGs<- colnames(Cluster_21_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_21_2 <- cor_data_Cluster_21[na.omit(match(Cluster_21_DEGs,rownames(cor_data_Cluster_21))),-na.omit(match(Cluster_21_DEGs,colnames(cor_data_Cluster_21)))] #na omit the ones that arent in my selected 100
  Cluster_21_sub_sample <- colnames(cor_data_Cluster_21_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_21_2 <- cor_data_Cluster_21_2[,na.omit(match(Cluster_21_sub_sample,colnames(cor_data_Cluster_21_2)))] #take all the rows and sample the columns
  bin_Cluster_21<-abs(cor_data_Cluster_21_2)
  bin_Cluster_21[which(bin_Cluster_21>sd(cor_data_Cluster_21_2))]<-1
  bin_Cluster_21[which(bin_Cluster_21!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_21<-bin_Cluster_21 %*% t(bin_Cluster_21) #adjacency matrix - 
  Cluster_21_entropy_result<- entropy(hyp_Cluster_21)
  Cluster_21_result_df <- rbind(Cluster_21_result_df, Cluster_21_entropy_result)
}
write.csv(Cluster_21_result_df, "path_to_file/Chronic_cluster_21_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_22 <90 SPOTS SO EXCLUDE####

#### Cluster_23 <90 SPOTS SO EXCLUDE####

#### Cluster_24 ####
set.seed(3)

Cluster_24 <- spatial_obs[spatial_obs$cell_type == "24", ]
Cluster_24_X <- Chronic_samples_X[na.omit(match(Cluster_24$X,rownames(Chronic_samples_X))),]

sd.scores24<-apply(Cluster_24_X,2,sd) #2 = selects columns not rows
Cluster_24_X_sd<-Cluster_24_X[,which(sd.scores24>0)]


Cluster_24_result_df <- data.frame()
cor_data_Cluster_24<-cor(Cluster_24_X_sd, Cluster_24_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_24_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_24_DEGs<- colnames(Cluster_24_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_24_2 <- cor_data_Cluster_24[na.omit(match(Cluster_24_DEGs,rownames(cor_data_Cluster_24))),-na.omit(match(Cluster_24_DEGs,colnames(cor_data_Cluster_24)))] #na omit the ones that arent in my selected 100
  Cluster_24_sub_sample <- colnames(cor_data_Cluster_24_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_24_2 <- cor_data_Cluster_24_2[,na.omit(match(Cluster_24_sub_sample,colnames(cor_data_Cluster_24_2)))] #take all the rows and sample the columns
  bin_Cluster_24<-abs(cor_data_Cluster_24_2)
  bin_Cluster_24[which(bin_Cluster_24>sd(cor_data_Cluster_24_2))]<-1
  bin_Cluster_24[which(bin_Cluster_24!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_24<-bin_Cluster_24 %*% t(bin_Cluster_24) #adjacency matrix - 
  Cluster_24_entropy_result<- entropy(hyp_Cluster_24)
  Cluster_24_result_df <- rbind(Cluster_24_result_df, Cluster_24_entropy_result)
}
write.csv(Cluster_24_result_df, "path_to_file/Chronic_cluster_24_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_25 ####
set.seed(3)
Cluster_25 <- spatial_obs[spatial_obs$cell_type == "25", ]
Cluster_25_X <- Chronic_samples_X[na.omit(match(Cluster_25$X,rownames(Chronic_samples_X))),]

sd.scores25<-apply(Cluster_25_X,2,sd) #2 = selects columns not rows
Cluster_25_X_sd<-Cluster_25_X[,which(sd.scores25>0)]


Cluster_25_result_df <- data.frame()
cor_data_Cluster_25<-cor(Cluster_25_X_sd, Cluster_25_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_25_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_25_DEGs<- colnames(Cluster_25_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_25_2 <- cor_data_Cluster_25[na.omit(match(Cluster_25_DEGs,rownames(cor_data_Cluster_25))),-na.omit(match(Cluster_25_DEGs,colnames(cor_data_Cluster_25)))] #na omit the ones that arent in my selected 100
  Cluster_25_sub_sample <- colnames(cor_data_Cluster_25_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_25_2 <- cor_data_Cluster_25_2[,na.omit(match(Cluster_25_sub_sample,colnames(cor_data_Cluster_25_2)))] #take all the rows and sample the columns
  bin_Cluster_25<-abs(cor_data_Cluster_25_2)
  bin_Cluster_25[which(bin_Cluster_25>sd(cor_data_Cluster_25_2))]<-1
  bin_Cluster_25[which(bin_Cluster_25!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_25<-bin_Cluster_25 %*% t(bin_Cluster_25) #adjacency matrix - 
  Cluster_25_entropy_result<- entropy(hyp_Cluster_25)
  Cluster_25_result_df <- rbind(Cluster_25_result_df, Cluster_25_entropy_result)
}
write.csv(Cluster_25_result_df, "path_to_file/Chronic_cluster_25_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_26 ####
set.seed(3)
Cluster_26 <- spatial_obs[spatial_obs$cell_type == "26", ]
Cluster_26_X <- Chronic_samples_X[na.omit(match(Cluster_26$X,rownames(Chronic_samples_X))),]

sd.scores26<-apply(Cluster_26_X,2,sd) #2 = selects columns not rows
Cluster_26_X_sd<-Cluster_26_X[,which(sd.scores26>0)]

Cluster_26_result_df <- data.frame()
cor_data_Cluster_26<-cor(Cluster_26_X_sd, Cluster_26_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_26_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_26_DEGs<- colnames(Cluster_26_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_26_2 <- cor_data_Cluster_26[na.omit(match(Cluster_26_DEGs,rownames(cor_data_Cluster_26))),-na.omit(match(Cluster_26_DEGs,colnames(cor_data_Cluster_26)))] #na omit the ones that arent in my selected 100
  Cluster_26_sub_sample <- colnames(cor_data_Cluster_26_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_26_2 <- cor_data_Cluster_26_2[,na.omit(match(Cluster_26_sub_sample,colnames(cor_data_Cluster_26_2)))] #take all the rows and sample the columns
  bin_Cluster_26<-abs(cor_data_Cluster_26_2)
  bin_Cluster_26[which(bin_Cluster_26>sd(cor_data_Cluster_26_2))]<-1
  bin_Cluster_26[which(bin_Cluster_26!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_26<-bin_Cluster_26 %*% t(bin_Cluster_26) #adjacency matrix - 
  Cluster_26_entropy_result<- entropy(hyp_Cluster_26)
  Cluster_26_result_df <- rbind(Cluster_26_result_df, Cluster_26_entropy_result)
}
write.csv(Cluster_26_result_df, "path_to_file/Chronic_cluster_26_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_27  ####
set.seed(3)
Cluster_27 <- spatial_obs[spatial_obs$cell_type == "27", ]
Cluster_27_X <- Chronic_samples_X[na.omit(match(Cluster_27$X,rownames(Chronic_samples_X))),]

sd.scores27<-apply(Cluster_27_X,2,sd) #2 = selects columns not rows
Cluster_27_X_sd<-Cluster_27_X[,which(sd.scores27>0)]

Cluster_27_result_df <- data.frame()
cor_data_Cluster_27<-cor(Cluster_27_X_sd, Cluster_27_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_27_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_27_DEGs<- colnames(Cluster_27_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_27_2 <- cor_data_Cluster_27[na.omit(match(Cluster_27_DEGs,rownames(cor_data_Cluster_27))),-na.omit(match(Cluster_27_DEGs,colnames(cor_data_Cluster_27)))] #na omit the ones that arent in my selected 100
  Cluster_27_sub_sample <- colnames(cor_data_Cluster_27_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_27_2 <- cor_data_Cluster_27_2[,na.omit(match(Cluster_27_sub_sample,colnames(cor_data_Cluster_27_2)))] #take all the rows and sample the columns
  bin_Cluster_27<-abs(cor_data_Cluster_27_2)
  bin_Cluster_27[which(bin_Cluster_27>sd(cor_data_Cluster_27_2))]<-1
  bin_Cluster_27[which(bin_Cluster_27!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_27<-bin_Cluster_27 %*% t(bin_Cluster_27) #adjacency matrix - 
  Cluster_27_entropy_result<- entropy(hyp_Cluster_27)
  Cluster_27_result_df <- rbind(Cluster_27_result_df, Cluster_27_entropy_result)
}
write.csv(Cluster_27_result_df, "path_to_file/Chronic_cluster_27_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_28 <10 SPOTS SO EXCLUDE####

##### Combining Chronic data EXCLUDED CLUSTERS 8 AND 23 FOR <10 SPOTS #####
##create new df
All_Chronic_clusters_results_df <- data.frame()

Cluster_0_result_df <- read.csv("path_to_file/Chronic_cluster_0_entropy_results_10000.csv",header=T, sep = ',')
Cluster_1_result_df <- read.csv("path_to_file/Chronic_cluster_1_entropy_results_10000.csv",header=T, sep = ',')
Cluster_2_result_df <- read.csv("path_to_file/Chronic_cluster_2_entropy_results_10000.csv",header=T, sep = ',')
Cluster_5_result_df <- read.csv("path_to_file/Chronic_cluster_5_entropy_results_10000.csv",header=T, sep = ',')
Cluster_6_result_df <- read.csv("path_to_file/Chronic_cluster_6_entropy_results_10000.csv",header=T, sep = ',')
Cluster_7_result_df <- read.csv("path_to_file/Chronic_cluster_7_entropy_results_10000.csv",header=T, sep = ',')
Cluster_9_result_df <- read.csv("path_to_file/Chronic_cluster_9_entropy_results_10000.csv",header=T, sep = ',')
Cluster_10_result_df <- read.csv("path_to_file/Chronic_cluster_10_entropy_results_10000.csv",header=T, sep = ',')
Cluster_11_result_df <- read.csv("path_to_file/Chronic_cluster_11_entropy_results_10000.csv",header=T, sep = ',')
Cluster_12_result_df <- read.csv("path_to_file/Chronic_cluster_12_entropy_results_10000.csv",header=T, sep = ',')
Cluster_13_result_df <- read.csv("path_to_file/Chronic_cluster_13_entropy_results_10000.csv",header=T, sep = ',')
Cluster_14_result_df <- read.csv("path_to_file/Chronic_cluster_14_entropy_results_10000.csv",header=T, sep = ',')
Cluster_15_result_df <- read.csv("path_to_file/Chronic_cluster_15_entropy_results_10000.csv",header=T, sep = ',')
Cluster_19_result_df <- read.csv("path_to_file/Chronic_cluster_19_entropy_results_10000.csv",header=T, sep = ',')
Cluster_21_result_df <- read.csv("path_to_file/Chronic_cluster_21_entropy_results_10000.csv",header=T, sep = ',')
Cluster_24_result_df <- read.csv("path_to_file/Chronic_cluster_24_entropy_results_10000.csv",header=T, sep = ',')
Cluster_25_result_df <- read.csv("path_to_file/Chronic_cluster_25_entropy_results_10000.csv",header=T, sep = ',')
Cluster_26_result_df <- read.csv("path_to_file/Chronic_cluster_26_entropy_results_10000.csv",header=T, sep = ',')
Cluster_27_result_df <- read.csv("path_to_file/Chronic_cluster_27_entropy_results_10000.csv",header=T, sep = ',')

# change first column name:
colnames(Cluster_0_result_df)[1] <- "Entropy"
colnames(Cluster_1_result_df)[1] <- "Entropy"
colnames(Cluster_2_result_df)[1] <- "Entropy"
colnames(Cluster_5_result_df)[1] <- "Entropy"
colnames(Cluster_6_result_df)[1] <- "Entropy"
colnames(Cluster_7_result_df)[1] <- "Entropy"
colnames(Cluster_9_result_df)[1] <- "Entropy"
colnames(Cluster_10_result_df)[1] <- "Entropy"
colnames(Cluster_11_result_df)[1] <- "Entropy"
colnames(Cluster_12_result_df)[1] <- "Entropy"
colnames(Cluster_13_result_df)[1] <- "Entropy"
colnames(Cluster_14_result_df)[1] <- "Entropy"
colnames(Cluster_15_result_df)[1] <- "Entropy"
colnames(Cluster_19_result_df)[1] <- "Entropy"
colnames(Cluster_21_result_df)[1] <- "Entropy"
colnames(Cluster_24_result_df)[1] <- "Entropy"
colnames(Cluster_25_result_df)[1] <- "Entropy"
colnames(Cluster_26_result_df)[1] <- "Entropy"
colnames(Cluster_27_result_df)[1] <- "Entropy"

##need to add columns and names to each df

All_Chronic_clusters_results_df <- rbind(Cluster_0_result_df,Cluster_1_result_df,Cluster_2_result_df,Cluster_5_result_df,Cluster_6_result_df,Cluster_7_result_df,Cluster_9_result_df,Cluster_10_result_df,Cluster_11_result_df,Cluster_12_result_df,Cluster_13_result_df,Cluster_14_result_df,Cluster_15_result_df,
                                         Cluster_19_result_df,Cluster_21_result_df, Cluster_24_result_df,Cluster_25_result_df,Cluster_26_result_df,Cluster_27_result_df)

All_Chronic_clusters_results_df$Iteration <- rep(1:10000, times = 19) # adds 1-1000, 29 times

# Create a vector with 0:24 repeated 1000 times
first_values <- rep(0:2, each = 10000)

second_values <- rep(5:7, each = 10000)

third_values <- rep(9:15, each = 10000)

forth_values <- rep(19, each = 10000)

fifth_values <- rep(21, each = 10000)

sixth_values <- rep(24:27, each = 10000)

zero_values <- rep("Adipocytes1", each = 10000)
one_values <- rep("FBs/SMCs1", each = 10000)
two_values <- rep("Keratinocytes/FBs1", each = 10000)
five_values <- rep("FBs/SMCs/Keratinocytes", each = 10000)
six_values<- rep("Keratinocytes2", each = 10000)
seven_values <- rep("FBs2", each = 10000)
nine_values  <- rep("FBs4", each = 10000)
ten_values<- rep("FBs/SMCs2", each = 10000)
eleven_values <- rep("B-cells/FBs", each = 10000)
twelve_values <- rep("Keratinocytes3", each = 10000)
thirteen_values <- rep("SMCs", each = 10000)
fourteen_values <- rep("Keratinocytes4", each = 10000)
fifteen_values <- rep("Keratinocytes5", each = 10000)
nineteen_values<- rep("Adipocytes3", each = 10000)
twentyone_values <- rep("Keratinocytes8", each = 10000)
twentyfour_values <- rep("Keratinocytes11", each = 10000)
twentyfive_values <- rep("Keratinocytes/FBs2", each = 10000)
twentysix_values  <- rep("FBs6", each = 10000)
twentyseven_values<- rep("Keratinocytes/FBs3", each = 10000)

# Combine both vectors
combined_values <- c(zero_values, one_values,two_values,five_values,six_values,seven_values,
                     nine_values,ten_values,eleven_values,twelve_values,thirteen_values,fourteen_values,fifteen_values,
                     nineteen_values,twentyone_values,twentyfour_values,twentyfive_values,twentysix_values,twentyseven_values)

# Assign the combined values to your data frame column
All_Chronic_clusters_results_df$Cluster <- combined_values
# Combine both vectors
combined_values <- c(first_values, second_values, third_values, forth_values, fifth_values, sixth_values)


All_Chronic_clusters_results_df$Cluster<- as.character(All_Chronic_clusters_results_df$Cluster)

write.csv(All_Chronic_clusters_results_df, "path_to_file/Chronic_clusters_entropy_results/Final_chronic_clusters_entropy_results_10000_Nov.csv", row.names=FALSE)

All_Chronic_clusters_results_df <- read.csv("path_to_file/Final_chronic_clusters_entropy_results_10000.csv")
All_Chronic_clusters_results_df$Cluster <- as.factor(All_Chronic_clusters_results_df$Cluster)

# define colours
colours<- c('#F8766D', "#EF7F48", "#E48800", "#D69100", "#C69900", "#B4A000", "#9EA700","#83AD00", "#60B200", "#19B700", "#00BB48", "#00BE6B", "#00C088", "#00C1A2","#00C0BA", "#00BECF", "#00B9E1", "#00B3F1", "#00ABFD", "#46A0FF", "#8794FF","#B086FF", "#CE79FF", "#E46DF6", "#F365E6", "#FC61D4", "#FF62BE", "#FF67A6","#FE6E8B")

names(colours) <- c(0:28)

cluster_colours<- c('Adipocytes1'='#F8766D',"FBs/SMCs1"= "#EF7F48","Keratinocytes/FBs1"= "#E48800","Keratinocytes1"= "#D69100","FBs1"= "#C69900","FBs/SMCs/Keratinocytes"= "#B4A000",
                    "Keratinocytes2"= "#9EA700","FBs2"="#83AD00","FBs3"= "#60B200","FBs4"= "#19B700","FBs/SMCs2"= "#00BB48","B-cells/FBs"= "#00BE6B","Keratinocytes3"= "#00C088","SMCs"= "#00C1A2",
                    "Keratinocytes4"="#00C0BA","Keratinocytes5"= "#00BECF","Adipocytes2"= "#00B9E1","Keratinocytes6"= "#00B3F1","FBs5"= "#00ABFD","Adipocytes3"= "#46A0FF","Keratinocytes7"= "#8794FF",
                    "Keratinocytes8"="#B086FF", "Keratinocytes9"="#CE79FF","Keratinocytes10"= "#E46DF6","Keratinocytes11"= "#F365E6","Keratinocytes/FBs2"= "#FC61D4","FBs6"= "#FF62BE",
                    "Keratinocytes/FBs3"="#FF67A6","Endothelial/Langerhans/Melanocytes/SMCs"="#FE6E8B")



cluster_colours

tiff("path_to_file/Final_chronic_clusters_entropy_boxplots_10000_iterations_named_clusters.tiff", height = 5000, width = 9000, res = 500)
p2 <- ggplot(data=All_Chronic_clusters_results_df, aes(y=factor(Cluster), x=Entropy)) +
  geom_boxplot(aes(color=Cluster))+
  coord_flip() +
  labs(title = "Entropy of Clusters in Chronic Samples",
       x = "Entropy",
       y = "Cluster") +
  theme(
    axis.title = element_text(size = 14),  # Adjust axis title size
    axis.text = element_text(size = 12),   # Adjust axis text size
    legend.title = element_text(size = 14),# Adjust legend title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ))
p2 + scale_color_manual(values=c(cluster_colours)) + scale_x_continuous(breaks=seq(8.2, 9.2, 0.2), limits = c(8.2, 9.2))
dev.off()



### CONTROL SAMPLES ####

Control_samples <- spatial_obs[spatial_obs$type == "Control", ]

Control_samples_X <- exp_data[na.omit(match(Control_samples$X,rownames(exp_data))),]

##### CONTROL SAMPLES #####
#### Cluster_0 ####

set.seed(3)
Cluster_0 <- spatial_obs[spatial_obs$cell_type == "0", ]
Cluster_0_X <- Control_samples_X[na.omit(match(Cluster_0$X,rownames(Control_samples_X))),]


sd.scores0<-apply(Cluster_0_X,2,sd) #2 = selects columns not rows
Cluster_0_X_sd<-Cluster_0_X[,which(sd.scores0>0)]

Cluster_0_result_df <- data.frame()
cor_data_Cluster_0<-cor(Cluster_0_X_sd, Cluster_0_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_0_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_0_DEGs<- colnames(Cluster_0_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_0_2 <- cor_data_Cluster_0[na.omit(match(Cluster_0_DEGs,rownames(cor_data_Cluster_0))),-na.omit(match(Cluster_0_DEGs,colnames(cor_data_Cluster_0)))] #na omit the ones that arent in my selected 100
  Cluster_0_sub_sample <- colnames(cor_data_Cluster_0_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_0_2 <- cor_data_Cluster_0_2[,na.omit(match(Cluster_0_sub_sample,colnames(cor_data_Cluster_0_2)))] #take all the rows and sample the columns
  bin_Cluster_0<-abs(cor_data_Cluster_0_2)
  bin_Cluster_0[which(bin_Cluster_0>sd(cor_data_Cluster_0_2))]<-1
  bin_Cluster_0[which(bin_Cluster_0!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_0<-bin_Cluster_0 %*% t(bin_Cluster_0) #adjacency matrix - 
  Cluster_0_entropy_result<- entropy(hyp_Cluster_0)
  Cluster_0_result_df <- rbind(Cluster_0_result_df, Cluster_0_entropy_result)
}
write.csv(Cluster_0_result_df, "path_to_file/Control_cluster_0_entropy_results_10000.csv", row.names=FALSE)


#### Cluster_1 ####
set.seed(3)
Cluster_1 <- spatial_obs[spatial_obs$cell_type == "1", ]
Cluster_1_X <- Control_samples_X[na.omit(match(Cluster_1$X,rownames(Control_samples_X))),]

sd.scores1<-apply(Cluster_1_X,2,sd) #2 = selects columns not rows
Cluster_1_X_sd<-Cluster_1_X[,which(sd.scores1>0)]


Cluster_1_result_df <- data.frame()
cor_data_Cluster_1<-cor(Cluster_1_X_sd, Cluster_1_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_1_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_1_DEGs<- colnames(Cluster_1_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_1_2 <- cor_data_Cluster_1[na.omit(match(Cluster_1_DEGs,rownames(cor_data_Cluster_1))),-na.omit(match(Cluster_1_DEGs,colnames(cor_data_Cluster_1)))] #na omit the ones that arent in my selected 100
  Cluster_1_sub_sample <- colnames(cor_data_Cluster_1_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_1_2 <- cor_data_Cluster_1_2[,na.omit(match(Cluster_1_sub_sample,colnames(cor_data_Cluster_1_2)))] #take all the rows and sample the columns
  bin_Cluster_1<-abs(cor_data_Cluster_1_2)
  bin_Cluster_1[which(bin_Cluster_1>sd(cor_data_Cluster_1_2))]<-1
  bin_Cluster_1[which(bin_Cluster_1!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_1<-bin_Cluster_1 %*% t(bin_Cluster_1) #adjacency matrix - 
  Cluster_1_entropy_result<- entropy(hyp_Cluster_1)
  Cluster_1_result_df <- rbind(Cluster_1_result_df, Cluster_1_entropy_result)
}
write.csv(Cluster_1_result_df, "path_to_file/Control_cluster_1_entropy_results_10000.csv", row.names=FALSE)


#### Cluster_2 ####
set.seed(3)
Cluster_2 <- spatial_obs[spatial_obs$cell_type == "2", ]
Cluster_2_X <- Control_samples_X[na.omit(match(Cluster_2$X,rownames(Control_samples_X))),]

sd.scores2<-apply(Cluster_2_X,2,sd) #2 = selects columns not rows
Cluster_2_X_sd<-Cluster_2_X[,which(sd.scores2>0)]



Cluster_2_result_df <- data.frame()
cor_data_Cluster_2<-cor(Cluster_2_X_sd, Cluster_2_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_2_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_2_DEGs<- colnames(Cluster_2_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_2_2 <- cor_data_Cluster_2[na.omit(match(Cluster_2_DEGs,rownames(cor_data_Cluster_2))),-na.omit(match(Cluster_2_DEGs,colnames(cor_data_Cluster_2)))] #na omit the ones that arent in my selected 100
  Cluster_2_sub_sample <- colnames(cor_data_Cluster_2_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_2_2 <- cor_data_Cluster_2_2[,na.omit(match(Cluster_2_sub_sample,colnames(cor_data_Cluster_2_2)))] #take all the rows and sample the columns
  bin_Cluster_2<-abs(cor_data_Cluster_2_2)
  bin_Cluster_2[which(bin_Cluster_2>sd(cor_data_Cluster_2_2))]<-1
  bin_Cluster_2[which(bin_Cluster_2!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_2<-bin_Cluster_2 %*% t(bin_Cluster_2) #adjacency matrix - 
  Cluster_2_entropy_result<- entropy(hyp_Cluster_2)
  Cluster_2_result_df <- rbind(Cluster_2_result_df, Cluster_2_entropy_result)
}
write.csv(Cluster_2_result_df, "path_to_file/Control_cluster_2_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_3 ####
set.seed(3)
Cluster_3 <- spatial_obs[spatial_obs$cell_type == "3", ]
Cluster_3_X <- Control_samples_X[na.omit(match(Cluster_3$X,rownames(Control_samples_X))),]

sd.scores3<-apply(Cluster_3_X,2,sd) #2 = selects columns not rows
Cluster_3_X_sd<-Cluster_3_X[,which(sd.scores3>0)]


Cluster_3_result_df <- data.frame()
cor_data_Cluster_3<-cor(Cluster_3_X_sd, Cluster_3_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_3_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_3_DEGs<- colnames(Cluster_3_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_3_2 <- cor_data_Cluster_3[na.omit(match(Cluster_3_DEGs,rownames(cor_data_Cluster_3))),-na.omit(match(Cluster_3_DEGs,colnames(cor_data_Cluster_3)))] #na omit the ones that arent in my selected 100
  Cluster_3_sub_sample <- colnames(cor_data_Cluster_3_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_3_2 <- cor_data_Cluster_3_2[,na.omit(match(Cluster_3_sub_sample,colnames(cor_data_Cluster_3_2)))] #take all the rows and sample the columns
  bin_Cluster_3<-abs(cor_data_Cluster_3_2)
  bin_Cluster_3[which(bin_Cluster_3>sd(cor_data_Cluster_3_2))]<-1
  bin_Cluster_3[which(bin_Cluster_3!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_3<-bin_Cluster_3 %*% t(bin_Cluster_3) #adjacency matrix - 
  Cluster_3_entropy_result<- entropy(hyp_Cluster_3)
  Cluster_3_result_df <- rbind(Cluster_3_result_df, Cluster_3_entropy_result)
}
write.csv(Cluster_3_result_df, "path_to_file/Control_cluster_3_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_4 ####
set.seed(3)
Cluster_4 <- spatial_obs[spatial_obs$cell_type == "4", ]
Cluster_4_X <- Control_samples_X[na.omit(match(Cluster_4$X,rownames(Control_samples_X))),]

sd.scores4<-apply(Cluster_4_X,2,sd) #2 = selects columns not rows
Cluster_4_X_sd<-Cluster_4_X[,which(sd.scores4>0)]


Cluster_4_result_df <- data.frame()
cor_data_Cluster_4<-cor(Cluster_4_X_sd, Cluster_4_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_4_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_4_DEGs<- colnames(Cluster_4_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_4_2 <- cor_data_Cluster_4[na.omit(match(Cluster_4_DEGs,rownames(cor_data_Cluster_4))),-na.omit(match(Cluster_4_DEGs,colnames(cor_data_Cluster_4)))] #na omit the ones that arent in my selected 100
  Cluster_4_sub_sample <- colnames(cor_data_Cluster_4_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_4_2 <- cor_data_Cluster_4_2[,na.omit(match(Cluster_4_sub_sample,colnames(cor_data_Cluster_4_2)))] #take all the rows and sample the columns
  bin_Cluster_4<-abs(cor_data_Cluster_4_2)
  bin_Cluster_4[which(bin_Cluster_4>sd(cor_data_Cluster_4_2))]<-1
  bin_Cluster_4[which(bin_Cluster_4!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_4<-bin_Cluster_4 %*% t(bin_Cluster_4) #adjacency matrix - 
  Cluster_4_entropy_result<- entropy(hyp_Cluster_4)
  Cluster_4_result_df <- rbind(Cluster_4_result_df, Cluster_4_entropy_result)
}
write.csv(Cluster_4_result_df, "path_to_file/Control_cluster_4_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_5 <90 SPOTS SO EXCLUDE####

#### Cluster_6 <90 SPOTS SO EXCLUDE####

#### Cluster_7 ####
set.seed(3)
Cluster_7 <- spatial_obs[spatial_obs$cell_type == "7", ]
Cluster_7_X <- Control_samples_X[na.omit(match(Cluster_7$X,rownames(Control_samples_X))),]

sd.scores7<-apply(Cluster_7_X,2,sd) #2 = selects columns not rows
Cluster_7_X_sd<-Cluster_7_X[,which(sd.scores7>0)]


Cluster_7_result_df <- data.frame()
cor_data_Cluster_7<-cor(Cluster_7_X_sd, Cluster_7_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_7_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_7_DEGs<- colnames(Cluster_7_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_7_2 <- cor_data_Cluster_7[na.omit(match(Cluster_7_DEGs,rownames(cor_data_Cluster_7))),-na.omit(match(Cluster_7_DEGs,colnames(cor_data_Cluster_7)))] #na omit the ones that arent in my selected 100
  Cluster_7_sub_sample <- colnames(cor_data_Cluster_7_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_7_2 <- cor_data_Cluster_7_2[,na.omit(match(Cluster_7_sub_sample,colnames(cor_data_Cluster_7_2)))] #take all the rows and sample the columns
  bin_Cluster_7<-abs(cor_data_Cluster_7_2)
  bin_Cluster_7[which(bin_Cluster_7>sd(cor_data_Cluster_7_2))]<-1
  bin_Cluster_7[which(bin_Cluster_7!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_7<-bin_Cluster_7 %*% t(bin_Cluster_7) #adjacency matrix - 
  Cluster_7_entropy_result<- entropy(hyp_Cluster_7)
  Cluster_7_result_df <- rbind(Cluster_7_result_df, Cluster_7_entropy_result)
}
write.csv(Cluster_7_result_df, "path_to_file/Control_cluster_7_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_8 ####
set.seed(3)
Cluster_8 <- spatial_obs[spatial_obs$cell_type == "8", ]
Cluster_8_X <- Control_samples_X[na.omit(match(Cluster_8$X,rownames(Control_samples_X))),]

sd.scores8<-apply(Cluster_8_X,2,sd) #2 = selects columns not rows
Cluster_8_X_sd<-Cluster_8_X[,which(sd.scores8>0)]


Cluster_8_result_df <- data.frame()
cor_data_Cluster_8<-cor(Cluster_8_X_sd, Cluster_8_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_8_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_8_DEGs<- colnames(Cluster_8_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_8_2 <- cor_data_Cluster_8[na.omit(match(Cluster_8_DEGs,rownames(cor_data_Cluster_8))),-na.omit(match(Cluster_8_DEGs,colnames(cor_data_Cluster_8)))] #na omit the ones that arent in my selected 100
  Cluster_8_sub_sample <- colnames(cor_data_Cluster_8_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_8_2 <- cor_data_Cluster_8_2[,na.omit(match(Cluster_8_sub_sample,colnames(cor_data_Cluster_8_2)))] #take all the rows and sample the columns
  bin_Cluster_8<-abs(cor_data_Cluster_8_2)
  bin_Cluster_8[which(bin_Cluster_8>sd(cor_data_Cluster_8_2))]<-1
  bin_Cluster_8[which(bin_Cluster_8!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_8<-bin_Cluster_8 %*% t(bin_Cluster_8) #adjacency matrix - 
  Cluster_8_entropy_result<- entropy(hyp_Cluster_8)
  Cluster_8_result_df <- rbind(Cluster_8_result_df, Cluster_8_entropy_result)
}
write.csv(Cluster_8_result_df, "path_to_file/Control_cluster_8_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_9 <90 SPOTS SO EXCLUDE####

#### Cluster_10 <90 SPOTS SO EXCLUDE####

#### Cluster_11 <90 SPOTS SO EXCLUDE#### 

#### Cluster_12 <90 SPOTS SO EXCLUDE####

#### Cluster_13 ####
set.seed(3)

Cluster_13 <- spatial_obs[spatial_obs$cell_type == "13", ]
Cluster_13_X <- Control_samples_X[na.omit(match(Cluster_13$X,rownames(Control_samples_X))),]

sd.scores13<-apply(Cluster_13_X,2,sd) #2 = selects columns not rows
Cluster_13_X_sd<-Cluster_13_X[,which(sd.scores13>0)]

Cluster_13_result_df <- data.frame()
cor_data_Cluster_13<-cor(Cluster_13_X_sd, Cluster_13_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_13_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_13_DEGs<- colnames(Cluster_13_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_13_2 <- cor_data_Cluster_13[na.omit(match(Cluster_13_DEGs,rownames(cor_data_Cluster_13))),-na.omit(match(Cluster_13_DEGs,colnames(cor_data_Cluster_13)))] #na omit the ones that arent in my selected 100
  Cluster_13_sub_sample <- colnames(cor_data_Cluster_13_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_13_2 <- cor_data_Cluster_13_2[,na.omit(match(Cluster_13_sub_sample,colnames(cor_data_Cluster_13_2)))] #take all the rows and sample the columns
  bin_Cluster_13<-abs(cor_data_Cluster_13_2)
  bin_Cluster_13[which(bin_Cluster_13>sd(cor_data_Cluster_13_2))]<-1
  bin_Cluster_13[which(bin_Cluster_13!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_13<-bin_Cluster_13 %*% t(bin_Cluster_13) #adjacency matrix - 
  Cluster_13_entropy_result<- entropy(hyp_Cluster_13)
  Cluster_13_result_df <- rbind(Cluster_13_result_df, Cluster_13_entropy_result)
}
write.csv(Cluster_13_result_df, "path_to_file/Control_cluster_13_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_14 <90 SPOTS SO EXCLUDE####

#### Cluster_15 <90 SPOTS SO EXCLUDE####

#### Cluster_16 ####
set.seed(3)
Cluster_16 <- spatial_obs[spatial_obs$cell_type == "16", ]
Cluster_16_X <- Control_samples_X[na.omit(match(Cluster_16$X,rownames(Control_samples_X))),]

sd.scores16<-apply(Cluster_16_X,2,sd) #2 = selects columns not rows
Cluster_16_X_sd<-Cluster_16_X[,which(sd.scores16>0)]


Cluster_16_result_df <- data.frame()
cor_data_Cluster_16<-cor(Cluster_16_X_sd, Cluster_16_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_16_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_16_DEGs<- colnames(Cluster_16_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_16_2 <- cor_data_Cluster_16[na.omit(match(Cluster_16_DEGs,rownames(cor_data_Cluster_16))),-na.omit(match(Cluster_16_DEGs,colnames(cor_data_Cluster_16)))] #na omit the ones that arent in my selected 100
  Cluster_16_sub_sample <- colnames(cor_data_Cluster_16_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_16_2 <- cor_data_Cluster_16_2[,na.omit(match(Cluster_16_sub_sample,colnames(cor_data_Cluster_16_2)))] #take all the rows and sample the columns
  bin_Cluster_16<-abs(cor_data_Cluster_16_2)
  bin_Cluster_16[which(bin_Cluster_16>sd(cor_data_Cluster_16_2))]<-1
  bin_Cluster_16[which(bin_Cluster_16!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_16<-bin_Cluster_16 %*% t(bin_Cluster_16) #adjacency matrix - 
  Cluster_16_entropy_result<- entropy(hyp_Cluster_16)
  Cluster_16_result_df <- rbind(Cluster_16_result_df, Cluster_16_entropy_result)
}
write.csv(Cluster_16_result_df, "path_to_file/Control_cluster_16_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_17 ####
set.seed(3)

Cluster_17 <- spatial_obs[spatial_obs$cell_type == "17", ]
Cluster_17_X <- Control_samples_X[na.omit(match(Cluster_17$X,rownames(Control_samples_X))),]

sd.scores17<-apply(Cluster_17_X,2,sd) #2 = selects columns not rows
Cluster_17_X_sd<-Cluster_17_X[,which(sd.scores17>0)]

Cluster_17_result_df <- data.frame()
cor_data_Cluster_17<-cor(Cluster_17_X_sd, Cluster_17_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_17_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_17_DEGs<- colnames(Cluster_17_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_17_2 <- cor_data_Cluster_17[na.omit(match(Cluster_17_DEGs,rownames(cor_data_Cluster_17))),-na.omit(match(Cluster_17_DEGs,colnames(cor_data_Cluster_17)))] #na omit the ones that arent in my selected 100
  Cluster_17_sub_sample <- colnames(cor_data_Cluster_17_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_17_2 <- cor_data_Cluster_17_2[,na.omit(match(Cluster_17_sub_sample,colnames(cor_data_Cluster_17_2)))] #take all the rows and sample the columns
  bin_Cluster_17<-abs(cor_data_Cluster_17_2)
  bin_Cluster_17[which(bin_Cluster_17>sd(cor_data_Cluster_17_2))]<-1
  bin_Cluster_17[which(bin_Cluster_17!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_17<-bin_Cluster_17 %*% t(bin_Cluster_17) #adjacency matrix - 
  Cluster_17_entropy_result<- entropy(hyp_Cluster_17)
  Cluster_17_result_df <- rbind(Cluster_17_result_df, Cluster_17_entropy_result)
}
write.csv(Cluster_17_result_df, "path_to_file/Control_cluster_17_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_18 #### 
set.seed(3)
Cluster_18 <- spatial_obs[spatial_obs$cell_type == "18", ]
Cluster_18_X <- Control_samples_X[na.omit(match(Cluster_18$X,rownames(Control_samples_X))),]


sd.scores18<-apply(Cluster_18_X,2,sd) #2 = selects columns not rows
Cluster_18_X_sd<-Cluster_18_X[,which(sd.scores18>0)]

Cluster_18_result_df <- data.frame()
cor_data_Cluster_18<-cor(Cluster_18_X_sd, Cluster_18_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_18_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_18_DEGs<- colnames(Cluster_18_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_18_2 <- cor_data_Cluster_18[na.omit(match(Cluster_18_DEGs,rownames(cor_data_Cluster_18))),-na.omit(match(Cluster_18_DEGs,colnames(cor_data_Cluster_18)))] #na omit the ones that arent in my selected 100
  Cluster_18_sub_sample <- colnames(cor_data_Cluster_18_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_18_2 <- cor_data_Cluster_18_2[,na.omit(match(Cluster_18_sub_sample,colnames(cor_data_Cluster_18_2)))] #take all the rows and sample the columns
  bin_Cluster_18<-abs(cor_data_Cluster_18_2)
  bin_Cluster_18[which(bin_Cluster_18>sd(cor_data_Cluster_18_2))]<-1
  bin_Cluster_18[which(bin_Cluster_18!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_18<-bin_Cluster_18 %*% t(bin_Cluster_18) #adjacency matrix - 
  Cluster_18_entropy_result<- entropy(hyp_Cluster_18)
  Cluster_18_result_df <- rbind(Cluster_18_result_df, Cluster_18_entropy_result)
}
write.csv(Cluster_18_result_df, "path_to_file/Control_cluster_18_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_19 <90 SPOTS SO EXCLUDE####

#### Cluster_20 ####
set.seed(3)
Cluster_20 <- spatial_obs[spatial_obs$cell_type == "20", ]
Cluster_20_X <- Control_samples_X[na.omit(match(Cluster_20$X,rownames(Control_samples_X))),]


sd.scores20<-apply(Cluster_20_X,2,sd) #2 = selects columns not rows
Cluster_20_X_sd<-Cluster_20_X[,which(sd.scores20>0)]


Cluster_20_result_df <- data.frame()
cor_data_Cluster_20<-cor(Cluster_20_X_sd, Cluster_20_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_20_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_20_DEGs<- colnames(Cluster_20_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_20_2 <- cor_data_Cluster_20[na.omit(match(Cluster_20_DEGs,rownames(cor_data_Cluster_20))),-na.omit(match(Cluster_20_DEGs,colnames(cor_data_Cluster_20)))] #na omit the ones that arent in my selected 100
  Cluster_20_sub_sample <- colnames(cor_data_Cluster_20_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_20_2 <- cor_data_Cluster_20_2[,na.omit(match(Cluster_20_sub_sample,colnames(cor_data_Cluster_20_2)))] #take all the rows and sample the columns
  bin_Cluster_20<-abs(cor_data_Cluster_20_2)
  bin_Cluster_20[which(bin_Cluster_20>sd(cor_data_Cluster_20_2))]<-1
  bin_Cluster_20[which(bin_Cluster_20!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_20<-bin_Cluster_20 %*% t(bin_Cluster_20) #adjacency matrix - 
  Cluster_20_entropy_result<- entropy(hyp_Cluster_20)
  Cluster_20_result_df <- rbind(Cluster_20_result_df, Cluster_20_entropy_result)
}
write.csv(Cluster_20_result_df, "path_to_file/Control_cluster_20_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_21 <90 SPOTS SO EXCLUDE####

#### Cluster_22 <90 SPOTS SO EXCLUDE####

#### Cluster_23 ####
set.seed(3)

Cluster_23 <- spatial_obs[spatial_obs$cell_type == "23", ]
Cluster_23_X <- Control_samples_X[na.omit(match(Cluster_23$X,rownames(Control_samples_X))),]

sd.scores23<-apply(Cluster_23_X,2,sd) #2 = selects columns not rows
Cluster_23_X_sd<-Cluster_23_X[,which(sd.scores23>0)]


Cluster_23_result_df <- data.frame()
cor_data_Cluster_23<-cor(Cluster_23_X_sd, Cluster_23_X_sd) #all genes x all genes
num_columns <- ncol(Cluster_23_X_sd)
for (i in 1:10000) {  # Outer loop (runs 1000 times)
  print(i)
  Cluster_23_DEGs<- colnames(Cluster_23_X_sd)[sample(1:num_columns,size=100)]
  cor_data_Cluster_23_2 <- cor_data_Cluster_23[na.omit(match(Cluster_23_DEGs,rownames(cor_data_Cluster_23))),-na.omit(match(Cluster_23_DEGs,colnames(cor_data_Cluster_23)))] #na omit the ones that arent in my selected 100
  Cluster_23_sub_sample <- colnames(cor_data_Cluster_23_2)[sample(1:(num_columns-100),size=1000)] 
  cor_data_Cluster_23_2 <- cor_data_Cluster_23_2[,na.omit(match(Cluster_23_sub_sample,colnames(cor_data_Cluster_23_2)))] #take all the rows and sample the columns
  bin_Cluster_23<-abs(cor_data_Cluster_23_2)
  bin_Cluster_23[which(bin_Cluster_23>sd(cor_data_Cluster_23_2))]<-1
  bin_Cluster_23[which(bin_Cluster_23!=1)]<-0 # this is the hypergraph incidence matrix
  hyp_Cluster_23<-bin_Cluster_23 %*% t(bin_Cluster_23) #adjacency matrix - 
  Cluster_23_entropy_result<- entropy(hyp_Cluster_23)
  Cluster_23_result_df <- rbind(Cluster_23_result_df, Cluster_23_entropy_result)
}
write.csv(Cluster_23_result_df, "path_to_file/Control_cluster_23_entropy_results_10000.csv", row.names=FALSE)

#### Cluster_24 <90 SPOTS SO EXCLUDE####

#### Cluster_25 <90 SPOTS SO EXCLUDE####

#### Cluster_26 <90 SPOTS SO EXCLUDE####

#### Cluster_27 <90 SPOTS SO EXCLUDE####

#### Cluster_28 <90 SPOTS SO EXCLUDE####

##### Combining Control data EXCLUDING CLUSTERS 6, 11, 15, 24, 25, 26, 27 FOR <90 SPOTS#####
##create new df
All_Control_clusters_results_df <- data.frame()

Cluster_0_result_df <- read.csv("path_to_file/Control_cluster_0_entropy_results_10000.csv",header=T, sep = ',')
Cluster_1_result_df <- read.csv("path_to_file/Control_cluster_1_entropy_results_10000.csv",header=T, sep = ',')
Cluster_2_result_df <- read.csv("path_to_file/Control_cluster_2_entropy_results_10000.csv",header=T, sep = ',')
Cluster_3_result_df <- read.csv("path_to_file/Control_cluster_3_entropy_results_10000.csv",header=T, sep = ',')
Cluster_4_result_df <- read.csv("path_to_file/Control_cluster_4_entropy_results_10000.csv",header=T, sep = ',')
Cluster_7_result_df <- read.csv("path_to_file/Control_cluster_7_entropy_results_10000.csv",header=T, sep = ',')
Cluster_8_result_df <- read.csv("path_to_file/Control_cluster_8_entropy_results_10000.csv",header=T, sep = ',')
Cluster_13_result_df <- read.csv("path_to_file/Control_cluster_13_entropy_results_10000.csv",header=T, sep = ',')
Cluster_16_result_df <- read.csv("path_to_file/Control_cluster_16_entropy_results_10000.csv",header=T, sep = ',')
Cluster_17_result_df <- read.csv("path_to_file/Control_cluster_17_entropy_results_10000.csv",header=T, sep = ',')
Cluster_18_result_df <- read.csv("path_to_file/Control_cluster_18_entropy_results_10000.csv",header=T, sep = ',')
Cluster_20_result_df <- read.csv("path_to_file/Control_cluster_20_entropy_results_10000.csv",header=T, sep = ',')
Cluster_23_result_df <- read.csv("path_to_file/Control_cluster_23_entropy_results_10000.csv",header=T, sep = ',')

# change first column name:
colnames(Cluster_0_result_df)[1] <- "Entropy"
colnames(Cluster_1_result_df)[1] <- "Entropy"
colnames(Cluster_2_result_df)[1] <- "Entropy"
colnames(Cluster_3_result_df)[1] <- "Entropy"
colnames(Cluster_4_result_df)[1] <- "Entropy"
colnames(Cluster_7_result_df)[1] <- "Entropy"
colnames(Cluster_8_result_df)[1] <- "Entropy"
colnames(Cluster_13_result_df)[1] <- "Entropy"
colnames(Cluster_16_result_df)[1] <- "Entropy"
colnames(Cluster_17_result_df)[1] <- "Entropy"
colnames(Cluster_18_result_df)[1] <- "Entropy"
colnames(Cluster_20_result_df)[1] <- "Entropy"
colnames(Cluster_23_result_df)[1] <- "Entropy"


##need to add columns and names to each df

All_Control_clusters_results_df <- rbind(Cluster_0_result_df,Cluster_1_result_df,Cluster_2_result_df,Cluster_3_result_df,Cluster_4_result_df, Cluster_7_result_df,Cluster_8_result_df,
                                         Cluster_13_result_df,Cluster_16_result_df,Cluster_17_result_df,Cluster_18_result_df,Cluster_20_result_df,Cluster_23_result_df)

All_Control_clusters_results_df$Iteration <- rep(1:10000, times = 13) # adds 1-1000, 29 times
#All_Control_clusters_results_df$Cluster <- rep(0:5, each = 10000) # add 0-27, 1000 times each

# Create a vector with 0:24 repeated 1000 times
first_values <- rep(0:4, each = 10000)

second_values <- rep(7:8, each = 10000)

third_values <- rep(13, each = 10000)

fourth_values <- rep(16:18, each = 10000)

fifth_values <- rep(20, each = 10000)

sixth_values <- rep(23, each = 10000)

# Combine both vectors
combined_values <- c(first_values, second_values, third_values, fourth_values, fifth_values, sixth_values)


# Create a vector with 0:24 repeated 1000 times
zero_values <- rep("Adipocytes1", each = 10000)
one_values <- rep("FBs/SMCs1", each = 10000)
two_values <- rep("Keratinocytes/FBs1", each = 10000)
three_values<- rep("Keratinocytes1", each = 10000)
four_values <- rep("FBs1", each = 10000)
seven_values <- rep("FBs2", each = 10000)
eight_values <- rep("FBs3", each = 10000)
thirteen_values <- rep("SMCs", each = 10000)
sixteen_values<- rep("Adipocytes2", each = 10000)
seventeen_values <- rep("Keratinocytes6", each = 10000)
eighteen_values  <- rep("FBs5", each = 10000)
twenty_values <- rep("Keratinocytes7", each = 10000)
twentythree_values <- rep("Keratinocytes10", each = 10000)



# Combine both vectors
combined_values <- c(zero_values, one_values,two_values,three_values,four_values,seven_values,eight_values,
                     thirteen_values,sixteen_values,seventeen_values,eighteen_values,twenty_values,twentythree_values)

# Assign the combined values to your data frame column
All_Control_clusters_results_df$Cluster <- combined_values



#All_Control_clusters_results_df$Cluster<- as.character(All_Control_clusters_results_df$Cluster)

write.csv(All_Control_clusters_results_df, "path_to_file/Final_control_clusters_entropy_results_10000_Nov.csv", row.names=FALSE)

All_Control_clusters_results_df <- read.csv("path_to_file/Final_control_clusters_entropy_results_10000.csv")
All_Control_clusters_results_df$Cluster <- as.factor(All_Control_clusters_results_df$Cluster)


# define colours
colours<- c('#F8766D', "#EF7F48", "#E48800", "#D69100", "#C69900", "#B4A000", "#9EA700","#83AD00", "#60B200", "#19B700", "#00BB48", "#00BE6B", "#00C088", "#00C1A2","#00C0BA", "#00BECF", "#00B9E1", "#00B3F1", "#00ABFD", "#46A0FF", "#8794FF","#B086FF", "#CE79FF", "#E46DF6", "#F365E6", "#FC61D4", "#FF62BE", "#FF67A6","#FE6E8B")

names(colours) <- c(0:28)
cluster_colours<- c('Adipocytes1'='#F8766D',"FBs/SMCs1"= "#EF7F48","Keratinocytes/FBs1"= "#E48800","Keratinocytes1"= "#D69100","FBs1"= "#C69900","FBs/SMCs/Keratinocytes"= "#B4A000",
                    "Keratinocytes2"= "#9EA700","FBs2"="#83AD00","FBs3"= "#60B200","FBs4"= "#19B700","FBs/SMCs2"= "#00BB48","B-cells/FBs"= "#00BE6B","Keratinocytes3"= "#00C088","SMCs"= "#00C1A2",
                    "Keratinocytes4"="#00C0BA","Keratinocytes5"= "#00BECF","Adipocytes2"= "#00B9E1","Keratinocytes6"= "#00B3F1","FBs5"= "#00ABFD","Adipocytes3"= "#46A0FF","Keratinocytes7"= "#8794FF",
                    "Keratinocytes8"="#B086FF", "Keratinocytes9"="#CE79FF","Keratinocytes10"= "#E46DF6","Keratinocytes11"= "#F365E6","Keratinocytes/FBs2"= "#FC61D4","FBs6"= "#FF62BE",
                    "Keratinocytes/FBs3"="#FF67A6","Endothelial/Langerhans/Melanocytes/SMCs"="#FE6E8B")



tiff("path_to_file/Final_control_clusters_entropy_boxplots_10000_iterations_named_clusters.tiff", height = 5000, width = 9000, res = 500)
p1 <- ggplot(data=All_Control_clusters_results_df, aes(y=factor(Cluster), x=Entropy)) +
  geom_boxplot(aes(color=Cluster))+
  coord_flip() +
  labs(title = "Entropy of Clusters in Control Samples",
       x = "Entropy",
       y = "Cluster") +
  theme(
    axis.title = element_text(size = 14),  # Adjust axis title size
    axis.text = element_text(size = 12),   # Adjust axis text size
    legend.title = element_text(size = 14),# Adjust legend title size
    legend.text = element_text(size = 12),  # Adjust legend text size
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 ))
p1 + scale_color_manual(values=c(cluster_colours)) + scale_x_continuous(breaks=seq(8.2, 9.2, 0.2), limits = c(8.2, 9.2))
dev.off()



##### MERGE ALL THE DATASETS TOGETHER #####
acute_values <- rep("Acute", each = 210000)
chronic_values <- rep("Chronic", each = 190000)
control_values <- rep("Control", each = 130000)

All_acute_clusters_results_df$Pathology <- acute_values
All_Chronic_clusters_results_df$Pathology <- chronic_values
All_Control_clusters_results_df$Pathology <- control_values

All_data <- rbind(All_acute_clusters_results_df, All_Chronic_clusters_results_df, All_Control_clusters_results_df)
All_data$Cluster<- as.character(All_data$Cluster)



tiff("path_to_file/Final_control_clusters_entropy_boxplots_10000_iterations.tiff", height = 5000, width = 9000, res = 500)
p1 <- ggplot(data=All_data, aes(y=factor(Cluster, levels=0:28), x=Entropy)) +
  geom_boxplot(aes(color=Cluster))+
  coord_flip() +
  labs(title = "Entropy of Clusters in Control Samples",
       x = "Entropy",
       y = "Cluster") +
  theme(
    axis.title = element_text(size = 14),  # Adjust axis title size
    axis.text = element_text(size = 12),   # Adjust axis text size
    legend.title = element_text(size = 14),# Adjust legend title size
    legend.text = element_text(size = 12)  # Adjust legend text size
  )
p1 + scale_color_manual(values=c(colours)) + scale_x_continuous(breaks=seq(8.2, 9.2, 0.2), limits = c(8.2, 9.2))
dev.off()






