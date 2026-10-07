library(data.table)
library(dplyr)
library(chemometrics)
library(vegan)
library(pcaPP)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # Windows路径

# 十类药物 (使用人数>20, 包括Statins)
med_cat10 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f",
               "med_hbp1_f","med_hbp2_f","med_hbp3_6_f","med_hbp4_f","med_hbp5_f",
               "med_lip1_f")
cov_traits <- c("sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg")



beta_diversity_list_mp4_g <- readRDS(paste0("results/indices/beta_diversity_mp4_genus.rds"))
# # Bray-Curtis dissimilarity结果默认是dist对象，可以转换为矩阵查看
# matrix_bray <- as.matrix(beta_diversity_list_mp4_g[[1]])
# # 查看前5x5矩阵
# print(matrix_bray[1:5, 1:5])
#### PERMANOVA (MP4, at the genus level) (用药与否的二分类表型组间的genus水平的Bray-Curtis dissimilarity差异) ----
dat <- beta_diversity_list_mp4_g[[1]]
# id <- names(dat) # 矩阵列名 (ID)
# med <- factor(med, levels = c("0","1"))
cov <- beta_diversity_list_mp4_g[[2]]
cov <- cov[,-1]
# 分类校正变量必须转换为因子
for (i in c(med_cat10,cov_traits)) {
  cov[[i]] <- factor(cov[[i]]) 
}

permanova_results_all <- data.frame()
for (i in 1:length(med_cat10)){
  
  print(paste0(i," out of ",length(med_cat10)," || ",med_cat10[i]," (多因素，校正其他药物) || ",Sys.time()))
  ## 多因素PERMANOVA (校正其他药物) 注意：主要分析元素要放在最后
  cov_order <- c(med_cat10[-i],med_cat10[i]) # 移动med_cat10[i]元素到末尾
  f2 <- as.formula(paste0("dat ~ ",paste(cov_order,collapse = " + "))) # 构建多因素分析公式
  set.seed(123)  # 固定随机种子
  adonis_result2 <- adonis2(f2, data = cov, permutations = 999, parallel = 8)
  # print(adonis_result2)
  adonis_result2$n <- nrow(cov)
  # 计算调整R方
  adonis_result2$adjusted_R2 <- RsquareAdj(adonis_result2$R2, adonis_result2$n, adonis_result2$Df)
  # 计算调整R方
  adonis_result2$cov <- row.names(adonis_result2)
  adonis_result2$trait <- med_cat10[i]
  adonis_result2$adj <- "multivariable"
  adonis_result2$formula <- paste0("dat ~ ",paste(cov_order,collapse = " + "))
  permanova_results_all <- rbind(permanova_results_all,adonis_result2)
  ## 多因素PERMANOVA (校正其他药物) 注意：主要分析元素要放在最后
}
openxlsx::write.xlsx(permanova_results_all, paste0("results/permanova/medication10_permanova_beta_diversity_mp4_genus.xlsx"))
#### PERMANOVA (MP4, at the genus level) (用药与否的二分类表型组间的genus水平的Bray-Curtis dissimilarity差异) ####
