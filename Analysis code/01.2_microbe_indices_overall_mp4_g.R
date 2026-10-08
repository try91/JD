library(data.table)
library(dplyr)
library(chemometrics)
library(vegan)
library(pcaPP)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/analyte_measurements_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/gut_microbial_composition_function_pathway_profiles_dat_20261006.txt", header = TRUE)
# 菌群MP4填补0值丰度 #
{
  # 菌群2014菌群 (分类和连续)
  # 丰度>0.0001, 出现率>10%的微生物 (物种和属)
  mp4_s_names <- colnames(micro_dat)[3:361]
  mp4_g_names <- colnames(micro_dat)[721:912]
  
  mp4_names <- c(mp4_s_names,mp4_g_names)
  
  for (col_name in mp4_names) {
    col_name_bin <- paste0(col_name,"_bin")  # 构建新变量名
    micro_dat[[col_name_bin]] <- ifelse(micro_dat[[col_name]] == 0, 0, 1)  # 创建新变量并赋值
  }
  
  for (col_name in mp4_names) {
    col_data <- micro_dat[[col_name]]  # 获取当前列的数据
    non_zero <- col_data[col_data!= 0]  # 找出该列的非零值
    
    if (length(non_zero) > 0) {
      min_non_zero <- min(non_zero,na.rm = TRUE)  # 获取非零最小值
      col_data[col_data == 0] <- min_non_zero / 2  # 将0值替换为非零最小值除以2
    }
    
    col_name_zero <- paste0(col_name,"_zero")  # 构建新变量名
    micro_dat[[col_name_zero]] <- col_data # 将处理后的列数据赋值回数据框新建的"填补0值后的列" (用于log转换；填充的数据，不能用来计算diversity和richness，也不能用来permanova)
  }
}
# 菌群MP4填补0值丰度 #


phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  right_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp0_3"


#### 变量整理 ----
# 菌群2014菌群 (分类和连续)
# 丰度>0.0001, 出现率>10%的微生物 (物种和属)
mp4_s_names <- colnames(micro_dat)[3:361]
mp4_g_names <- colnames(micro_dat)[721:912]
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp4_s_names_short <- mp4_s_names[!grepl("_GGB",mp4_s_names)] # 排除未分类的菌属（GGB）, 未分类菌种（SGB）先保留
mp4_g_names_short <- mp4_g_names[!grepl("_GGB",mp4_g_names)] # 排除未分类的菌属（GGB）
# 排除未分类的菌属（GGB）和菌种（SGB） #

# 转换后的菌的名称
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)
mp4_g_zero <- paste0(mp4_g_names,"_zero") # 菌群MP4填补0值丰度 (属层面)


# 2010污染物 (连续)
edc_traits <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                "BPA","BPS","BPF",
                "TCC","TCS")
edc_traits2 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS")
edc_traits3 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP")
edc_traits3_q2 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP") # 检出率>50%的PAE6
edc_traits3_q4 <- c("MEHP","MECPP","MEHHP","MEP") # 检出率>75%的PAE4
edc_traits4 <- c("BPA","BPS","BPF")
edc_traits4_q2 <- c("BPA") # 检出率>50%的BP1
edc_traits5 <- c("TCC","TCS")
edc_traits6 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                 "BPA","TCC","TCS") # 检出率>50%
edc_traits7 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                 "MEHP","MECPP","MEHHP","MEP") # 检出率>75%
edc_traits8 <- c("MnBP","MCPP","MBzP",
                 "BPS","BPF") # 检出率<50%

edc_traits_log10 <- paste0(edc_traits,"_log10")
edc_traits2_log10 <- paste0(edc_traits2,"_log10")
edc_traits3_log10 <- paste0(edc_traits3,"_log10")
edc_traits3_q2_log10 <- paste0(edc_traits3_q2,"_log10")
edc_traits3_q4_log10 <- paste0(edc_traits3_q4,"_log10")
edc_traits4_log10 <- paste0(edc_traits4,"_log10")
edc_traits4_q2_log10 <- paste0(edc_traits4_q2,"_log10")
edc_traits5_log10 <- paste0(edc_traits5,"_log10")
edc_traits6_log10 <- paste0(edc_traits6,"_log10")
edc_traits7_log10 <- paste0(edc_traits7,"_log10")
edc_traits8_log10 <- paste0(edc_traits8,"_log10")
# EDC INDEX 变量名
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")
edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")


# 2021、2014死亡和新发表型 (分类)
phy_incident_cat <- c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014")
phy_incident_time <- c("timecvd_1021","timecvd_1014","timeckd_1014","timedm_1014")
phy_censor_cat <- c("censorall_1021","censorall_1014")
phy_censor_time <- c("timeall_1021","timeall_1014")
# 2014、2010表型 (分类)
phy_out_cat <- c("cvd_f","ckd_f","dm_f")
phy_traits_cat <- c("cvd_b","ckd_b","dm_b","as_imt_f","as_imt_b","hpt_f","hpt_b","nafld_f","nafld_b",
                    "ob_f","ob_b","abob_f","abob_b","dyslip_f","dyslip_b","hua_f","hua_b","ir_f","ir_b","mets_f","mets_b")
# 2014、2010表型 (连续)
phy_traits_cont <- c("bmi_f","bmi_b","wc_f","wc_b","hc_f","hc_b","whr_f","whr_b","height_f","height_b","weight_f","weight_b",
                     "hdl_f","hdl_b","ldl_f","ldl_b","apoa_f","apoa_b","apob_f","apob_b","chol_f","chol_b","tg_f","tg_b","nonhdl_f","nonhdl_b",
                     "alt_f","alt_b","ast_f","ast_b","ggt_f","ggt_b","scr_f","scr_b","egfr_f","egfr_b","ua_f","ua_b","bia_f","bia_b",
                     "glu0_f","glu0_b","glu120_f","glu120_b","vhba1c_f","vhba1c_b","ins0_f","ins0_b","ins120_f","ins120_b","homair_f","homair_b","homab_f","homab_b",
                     
                     "sbp_f","sbp_b","dbp_f","dbp_b","pr_f","pr_b",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "wbc_f","wbc_b","crp_f",
                     "plt_f","plt_b","hgb_f","hgb_b","eos_f","lym_f","mon_f","neu_f",
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f")
# 2014药物 (分类)
# 十类药物 (使用人数>20, 包括Statins)
med_cat10 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f",
               "med_hbp1_f","med_hbp2_f","med_hbp3_6_f","med_hbp4_f","med_hbp5_f",
               "med_lip1_f")
# 六类与菌群显著相关药物 (Sulfonylureas, Biguanides, Thiazolidinediones, AGIs, ARBs, Calcium antagonists) + Statins (MP4数据)
med_cat7 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f", 
              "med_hbp1_f","med_hbp4_f", 
              "med_lip1_f")
#### 变量整理 ####


#### 整体人群的菌群(MP4, at the genus level)指标(richness, shannon index, and uniqueness) ----
## 构建计算Richness & Shannon的函数 (MP4, genus) ##
RS <- function(DAT){
  pr <- DAT[,mp4_g_names] # 行为个体样本，列为微生物的数据格式 (使用未填补0值的数据)
  pr <- as.matrix(pr)
  row.names(pr) <- DAT$id14_15
  
  ### Richness (number of observed taxa) (R ‘vegan’ package (version 2.6-4))
  richness <- specnumber(pr, MARGIN = 1) # MARGIN=1表示按行计算
  richness <- as.data.frame(richness)
  ### Shannon diversity (within-sample, alpha-diversity) (R ‘vegan’ package (version 2.6-4))
  shannon <- diversity(pr, index = "shannon", MARGIN = 1) # MARGIN=1表示按行计算
  shannon <- as.data.frame(shannon)
  
  # 整合Richness & Shannon
  dat.rs <- cbind(richness,shannon)
  
  return(dat.rs)
}
## 构建计算Distance的函数 (MP4, genus) ##
Distance <- function(DAT){
  pr <- DAT[,mp4_g_names] # 行为个体样本，列为微生物的数据格式 (使用未填补0值的数据)
  pr <- as.matrix(pr)
  row.names(pr) <- DAT$id14_15
  
  pr.zero <- DAT[,mp4_g_zero] # 行为个体样本，列为微生物的数据格式 (使用填补0值后的数据)
  pr.zero <- as.matrix(pr.zero)
  row.names(pr.zero) <- DAT$id14_15
  pr.clr <- clr(pr.zero) # centered log-ratio-transformed species abundances
  
  ### 4 uniqueness indices using different distance metrics to capture different aspects of within-sample microbial variations: 
  ## Bray–Curtis dissimilarity (‘vegdist’ function in the R package ‘vegan’ (version 2.6-4, method=‘bray’))
  d.bray <- vegdist(pr, method = "bray")
  ## Jaccard distance (‘vegdist’ function in the R package ‘vegan’ (version 2.6-4, method=‘jaccard’))
  d.jaccard <- vegdist(pr, method = "jaccard")
  ## Aitchison distance (‘vegdist’ function in the R package ‘vegan’ (version 2.6-4, method=‘euclidean’ with centered log-ratio-transformed species abundances))
  d.aitchison <- vegdist(pr.clr, method = "euclidean")
  ## Kendall’s tau coefficient (‘cor.fk’ function in the R package ‘pcaPP’ (version 2.0-5))
  d.kendalls <- 1 - cor.fk(t(pr.clr)) # 将相关系数转换为相异性 (0=完全正相关，2=完全负相关)
  # Spearman's correlation ('cor' function in the R package 'stats' (version 4.3.1))
  d.spearman <- 1 - cor(t(pr.clr),method = "spearman") # 将相关系数转换为相异性 (0=完全正相关，2=完全负相关)
  
  # 整合5种距离指标
  dat.dist <- list(Bray = d.bray,
                   Jaccard = d.jaccard,
                   Aitchison = d.aitchison,
                   Kendall = d.kendalls,
                   Spearman = d.spearman
  )
  
  return(dat.dist)
}
## 构建计算Uniqueness的函数 (MP4, genus) ##
Uniqueness <- function(DAT){
  d <- as.matrix(DAT)
  u <- matrix(NA,nrow(d),1)
  u <- as.data.frame(u)
  rownames(u) <- rownames(d)
  colnames(u) <- "Uniqueness"
  for(x in 1:nrow(d)){
    a <- as.numeric(d[x,-x])
    u[x,1] <- min(a)
  }
  return(u)
}


indices_results <- data.frame()
for (i in c("phy_edc_temp0_3")){
  print(paste0(i," || MP4 Genus || ",Sys.time()))
  
  dat_temp <- phy_edc_dat
  
  # 计算 Richness & Shannon
  rs_temp <- RS(dat_temp)
  
  # 计算 Distance
  dis_temp <- Distance(dat_temp)
  
  # 计算 Uniqueness
  uniq_temp <- lapply(dis_temp,Uniqueness)
  uniq_temp <- do.call(cbind,uniq_temp)
  colnames(uniq_temp) <- names(dis_temp)
  uniq_temp <- as.data.frame(uniq_temp)
  
  # 整合Richness & Shannon and Uniqueness
  dat_temp_results <- cbind(rs_temp,uniq_temp)
  colnames(dat_temp_results) <- paste0(colnames(dat_temp_results))
  
  dat_temp_results$id14_15 <- row.names(dat_temp_results)
  
  if (ncol(indices_results) == 0) {  # 首次合并时为空
    indices_results <- dat_temp_results
  } else {                           # 后续合并
    indices_results <- left_join(indices_results, dat_temp_results, by="id14_15")
  }
 
}

### 结果数据 ###
openxlsx::write.xlsx(indices_results, paste0("results/indices/overall_indices7_mp4_genus.xlsx"))
#### 整体人群的菌群(MP4, at the genus level)指标(richness, shannon index, and uniqueness) ####
