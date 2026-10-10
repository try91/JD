library(data.table)
library(dplyr)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders


#### 变量整理 ----
# # 菌群2014菌群 (分类和连续)
# # 丰度>0.0001, 出现率>10%的微生物 (物种和属)
# mp4_s_names <- colnames(micro_dat)[3:361]
# mp4_g_names <- colnames(micro_dat)[721:912]
# # 排除未分类的菌属（GGB）和菌种（SGB） #
# mp4_s_names_short <- mp4_s_names[!grepl("_GGB",mp4_s_names)] # 排除未分类的菌属（GGB）, 未分类菌种（SGB）先保留
# mp4_g_names_short <- mp4_g_names[!grepl("_GGB",mp4_g_names)] # 排除未分类的菌属（GGB）
# # 排除未分类的菌属（GGB）和菌种（SGB） #
# 
# # 转换后的菌的名称
# mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
# mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
# 
# mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
# mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)


# 2010污染物 (连续)
edc_traits <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                "TCC","TCS",
                "BPA","BPS","BPF")
edc_traits2 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS")
edc_traits3 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP")
edc_traits3_q2 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP") # 检出率>50%的PAE6
edc_traits3_q4 <- c("MEHP","MECPP","MEHHP","MEP") # 检出率>75%的PAE4
edc_traits4 <- c("BPA","BPS","BPF")
edc_traits4_q2 <- c("BPA") # 检出率>50%的BP1
edc_traits5 <- c("TCC","TCS")
edc_traits6 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                 "TCC","TCS","BPA") # 检出率>50%
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
phy_out_cat <- c("dm_f","ckd_f","cvd_f")
phy_traits_cat <- c("dm_f","ckd_f","cvd_f","dm_b","ckd_b","cvd_b",
                    "ob_f","ob_b","abob_f","abob_b","ir_f","ir_b","dyslip_f","dyslip_b",
                    "mets_f","mets_b","nafld_f","nafld_b","hua_f","hua_b","hpt_f","hpt_b","as_imt_f","as_imt_b")
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


#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ----
### 读取分类组间EDC差异数据
wilcoxon_results <- readRDS("results/correlations/wilcoxon/wilcoxon_results_edc.rds")
wilcoxon_results <- wilcoxon_results[[1]]
wilcoxon_results <- wilcoxon_results[wilcoxon_results$exp_name %in% edc_traits_log10 & wilcoxon_results$out_name %in% c("age_g_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg_f"),]
wilcoxon_results$exp_name <- factor(wilcoxon_results$exp_name,
                                    levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                      "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                      "TCC","TCS",
                                                      "BPA","BPS","BPF"),"_log10"),
                                    labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                               "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                               "TCC","TCS",
                                               "BPA","BPS","BPF"))
wilcoxon_results$out_name <- factor(wilcoxon_results$out_name, 
                                    levels = c("age_g_b","sex_b_rev","high_edu_b","smk1_b","drk1_b","paactive3_g_b","high_fruveg_f"), 
                                    labels = c("Older age (≥60 y)","Men","Higher educational attainment","Current smoking","Current drinking","Sufficient physical activity","Adequate fruits and vegetables intake"))
wilcoxon_results <- wilcoxon_results %>%
  arrange(out_name,exp_name)
# 以每个outcome表型为单位，校正exposure
wilcoxon_results <- wilcoxon_results %>%
  group_by(out_name) %>%  # 按outcome分组，校正exposure
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


df <- wilcoxon_results
# 当 Cliff's Delta > 0 时，表示第一组有更高的倾向获得更大的值。分析时第一组是0，第二组是1，此处统一方向加一个负号，表示Cliff's Delta > 0 时，表示第二组有更高的倾向获得更大的值。
df$cliff.delta <- -df$cliff.delta

df <- df[,c("out_name","exp_name","p","p_adj_bh","cliff.delta")]


colnames(df) <- c("Factor","Analyte","P value","BH-adjusted P","Effect size (Cliff's delta, two-sided Wilcoxon rank-sum test)")
openxlsx::write.xlsx(df,"tables/(stable06)_cov_edc_wilcoxon.xlsx")
#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ####

#### EDC-EDC 相关性 (Partial Spearman) 填补后 ----
### 读取EDC-EDC关联性数据
p_spearman_results <- readRDS("results/correlations/spearman/partial_spearman_results_edc-edc_mp4_out.rds")
p_spearman_results <- p_spearman_results[[1]]
p_spearman_results <- p_spearman_results[p_spearman_results$exp_name %in% edc_traits_log10 & p_spearman_results$out_name %in% edc_traits_log10,]

p_spearman_results$exp_name <- factor(p_spearman_results$exp_name, 
                                      levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                        "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                        "TCC","TCS",
                                                        "BPA","BPS","BPF"),"_log10"),
                                      labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                 "TCC","TCS",
                                                 "BPA","BPS","BPF"))
p_spearman_results$out_name <- factor(p_spearman_results$out_name,
                                      levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                        "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                        "TCC","TCS",
                                                        "BPA","BPS","BPF"),"_log10"),
                                      labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                 "TCC","TCS",
                                                 "BPA","BPS","BPF"))


p_spearman_results <- p_spearman_results %>%
  arrange(exp_name,out_name)
# 以每个outcome表型为单位，校正exposure
p_spearman_results <- p_spearman_results %>%
  group_by(out_name) %>%  # 按outcome分组，校正exposure
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


df <- p_spearman_results
df <- df[df$exp_name != df$out_name,]
df <- df[,c("exp_name","out_name","estimate","p","p_adj_bh")]
# 核心：无向组合去重 #
df$key <- apply(df[, c("exp_name", "out_name")], 1, function(x) {
  paste(sort(x), collapse = "_")  # 排序后拼接成唯一键
})

# 保留每个唯一组合的第一行，删除重复的 (a,b) (b,a) #
df <- df[!duplicated(df$key), ]

# 删掉辅助列（可选） #
df$key <- NULL


colnames(df) <- c("Analyte 1","Analyte 2","Spearman's rho (two-sided)","P value","BH-adjusted P")
openxlsx::write.xlsx(df,"tables/(stable07)_edc_edc_correaltion.xlsx")
#### EDC-EDC 相关性 (Partial Spearman) 填补后 ####
