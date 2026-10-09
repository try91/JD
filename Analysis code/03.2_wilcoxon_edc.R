library(data.table)
library(dplyr)
library(ppcor)
library(effsize)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
phenotype_dat$age_g_b <- ifelse(phenotype_dat$age_b < 60, 0, 1)
phenotype_dat$high_fruveg_f <- ifelse(phenotype_dat$high_fruveg == 99, NA, phenotype_dat$high_fruveg)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")

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

mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)


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


# 分类变量名
cat_traits <- c("age_g_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg_f")
# EDC变量名
edc_traits <- c(edc_traits_log10)


#### 构建Wilcoxon分析函数 ----
# 创建一个函数，用于计算分组变量之间的显著性差异
wilcoxon_test <- function(CONT, CAT, DAT) {
  
  DAT[[CAT]] <- ifelse(DAT[[CAT]] == 99, NA, DAT[[CAT]])
  
  DAT <- DAT[!is.na(DAT[[CAT]]),]
  # 首先，确保CAT列是因子类型
  DAT[[CAT]] <- factor(DAT[[CAT]])

  # Mann-Whitney U (Wilcoxon秩和检验)
  wilcoxon_test <- wilcox.test(DAT[[CONT]] ~ DAT[[CAT]], data = DAT, alternative = "two.sided")
  wilcoxon_test2 <- wilcox.test(DAT[[CONT]] ~ DAT[[CAT]], data = DAT, alternative = "greater") # 检测 group0 的中位数是否大于 group1 的中位数
  wilcoxon_test3 <- wilcox.test(DAT[[CONT]] ~ DAT[[CAT]], data = DAT, alternative = "less") # 检测 group0 的中位数是否小于 group1 的中位数
  cliff_delta <- cliff.delta(DAT[[CONT]][DAT[[CAT]] == 0], DAT[[CONT]][DAT[[CAT]] == 1])  # 效应量（Cliff's Delta），输出δ值及置信区间 （当 Cliff's Delta > 0 时，表示第一组有更高的倾向获得更大的值）
  
  temp_wilcoxon <- data.frame(
    exp_name = CONT,
    out_name = CAT,
    method = "wilcoxon",
    median_group0 = median(DAT[[CONT]][DAT[[CAT]] == 0], na.rm = TRUE),
    median_group1 = median(DAT[[CONT]][DAT[[CAT]] == 1], na.rm = TRUE),
    p = wilcoxon_test$p.value,
    p_group0_gt_group1 = wilcoxon_test2$p.value,
    p_group1_gt_group0 = wilcoxon_test3$p.value,
    cliff.delta = cliff_delta$estimate,
    cliff.delta_lower = cliff_delta[["conf.int"]][["lower"]],
    cliff.delta_upper = cliff_delta[["conf.int"]][["upper"]],
    sample = sample_name,
    n = length(na.omit(DAT[[CONT]]))
  )
  
  return(temp_wilcoxon)
  
}
#### 构建Wilcoxon分析函数 ####


#### EDC ----
wilcoxon_results_list <- list()
for (i in c("phy_edc_temp")) {
  # i <- "phy_edc_temp"
  
  sample_name <- i
  
  # 由于cat变量包含大于二分类的变量，筛选仅含两种唯一值的变量名（忽略NA）
  temp <- phy_edc_dat[,cat_traits]
  two_value_vars <- names(temp)[sapply(temp, function(x) length(unique(na.omit(x))) == 2)]
  
  #### Wilcoxon rank-sum test (二分类结局表型组间的indices(overall)差异) ----
  wilcoxon_results_edc <- data.frame()
  for (k in edc_traits) {
    print(paste0(sample_name," |Wilcoxon EDC| ",which(edc_traits == k)," out of ",length(edc_traits),": ",k))
    
    for (m in two_value_vars) {
      
      ## 分析样本选取 (排除缺失的项) ##
      cols <- c(k, m)
      phy_edc_dat_temp <- phy_edc_dat[,cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      
      result <- wilcoxon_test(k, m, phy_edc_dat_temp)
      wilcoxon_results_edc <- rbind(wilcoxon_results_edc, result)
    }
  }
  wilcoxon_results_list[[paste0("wilcoxon_results_edc_(",sample_name,")")]] <- wilcoxon_results_edc
}
### 保存单个数据 ###
saveRDS(wilcoxon_results_list, paste0("results/correlations/wilcoxon/wilcoxon_results_edc.rds"))
#### EDC ####
