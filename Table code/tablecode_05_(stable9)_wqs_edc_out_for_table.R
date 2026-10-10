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


#### 数据处理 (for stable) ----
### 设置纳入图片的暴露和结局
## 暴露：EDC
exposure <- c("EDC_14","PFAS","PAE6","TC","BP1")
# 标准化名称
exposure_label <- c("EDCs (14)","PFAS (5)","PAEs (6)","Antimicrobials (2)","Bisphenols (1)")
exp_dat <- as.data.frame(exposure)

## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
outcome <- c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021",
             "dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b",
             "bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
             "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
             "alt_b","ast_b","ggt_b","bia_b",
             "egfr_b","scr_b","ua_b",
             "glu0_b","glu120_b","vhba1c_b",
             "ins0_b","ins120_b","homair_b","homab_b",
             "sbp_b","dbp_b","pr_b",
             "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
             "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
             "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")

outcome_final <- outcome
out_dat <- as.data.frame(outcome_final)
# 标准化名称
outcome_label <- c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                   "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT",
                   "BMI","Height","Weight","WHR","WC","HC",
                   "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                   "ALT","AST","GGT","Bile acid",
                   "eGFR","Serum creatinine","UA",
                   "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                   "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                   "SBP","DBP","PR",
                   "FT3","FT4","TSH","TPOAb","TgAb",
                   "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                   "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")

## 构建有完整exposure和outcome的数据框 ##
out_dat1 <- out_dat
out_dat1$exposure <- exposure[1]
out_dat2 <- out_dat
out_dat2$exposure <- exposure[2]
out_dat3 <- out_dat
out_dat3$exposure <- exposure[3]
out_dat4 <- out_dat
out_dat4$exposure <- exposure[4]
out_dat5 <- out_dat
out_dat5$exposure <- exposure[5]
dat_exp_out_all <- rbind(out_dat1,out_dat2,out_dat3,out_dat4,out_dat5)
colnames(dat_exp_out_all) <- c("out", "exp")
## 构建有完整exposure和outcome的数据框 ##


# 设置人群名称
sample_name1 <- "phy_edc_temp" # 用总人群

# 读取数据 (总人群)
results_wqs_bin <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_bin_(",sample_name1,").xlsx"))
# 读取数据 (总人群)
results_wqs_cont <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_cont_(",sample_name1,").xlsx"))

results_wqs <- rbind(results_wqs_bin, results_wqs_cont) %>%
  filter(exp %in% exposure & out %in% outcome_final & cov == "wqs" & type == "Q2")

results_wqs_pos <- results_wqs[results_wqs$direction == "pos",] # 正向权重wqs分析结果
unique(results_wqs_pos$exp) # 缺BP1
unique(results_wqs_pos$out)
results_wqs_pos <- left_join(dat_exp_out_all,results_wqs_pos,by=c("exp", "out"))
results_wqs_pos$p <- ifelse(is.na(results_wqs_pos$p), 1, results_wqs_pos$p)

results_wqs_neg <- results_wqs[results_wqs$direction == "neg",] # 负向权重wqs分析结果
unique(results_wqs_neg$exp) # 缺BP1
unique(results_wqs_neg$out) # 缺whr_b
results_wqs_neg <- left_join(dat_exp_out_all,results_wqs_neg,by=c("exp", "out"))
results_wqs_neg$p <- ifelse(is.na(results_wqs_neg$p), 1, results_wqs_neg$p)

# FDR 校正
# 以每个exposure表型为单位，校正outcome
dat_wqs_pos <- results_wqs_pos %>%
  group_by(exp) %>%  # 按exposure分组，校正outcome
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()

dat_wqs_neg <- results_wqs_neg %>%
  group_by(exp) %>%  # 按exposure分组，校正outcome
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


# EDC 权重数据整理(正)
dat_wqs_pos$numb_pos_weight <- rowSums(dat_wqs_pos[, 13:26] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
dat_wqs_pos$numb_neg_weight <- rowSums(dat_wqs_pos[, 13:26] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量

dat_wqs_pos$weight_threshold_pos <- 1 / dat_wqs_pos$numb_pos_weight
dat_wqs_pos$weight_threshold_neg <- 1 / dat_wqs_pos$numb_neg_weight

dat_wqs_pos_edc14 <- dat_wqs_pos[dat_wqs_pos$exp == "EDC_14",]
dat_wqs_pos_edc14_long <- tidyr::gather(dat_wqs_pos_edc14, edc, weight, 13:26, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！


# EDC 权重数据整理(负)
dat_wqs_neg$numb_pos_weight <- rowSums(dat_wqs_neg[, 13:26] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
dat_wqs_neg$numb_neg_weight <- rowSums(dat_wqs_neg[, 13:26] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量

dat_wqs_neg$weight_threshold_pos <- 1 / dat_wqs_neg$numb_pos_weight
dat_wqs_neg$weight_threshold_neg <- 1 / dat_wqs_neg$numb_neg_weight

dat_wqs_neg_edc14 <- dat_wqs_neg[dat_wqs_neg$exp == "EDC_14",]
dat_wqs_neg_edc14_long <- tidyr::gather(dat_wqs_neg_edc14, edc, weight, 13:26, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！


# WQS positive model #
dat_wqs_pos$exp <- factor(dat_wqs_pos$exp, levels = exposure, labels = exposure_label)
dat_wqs_pos$out <- factor(dat_wqs_pos$out, levels = outcome_final, labels = outcome_label)
dat_wqs_pos <- dat_wqs_pos %>% arrange(out,exp)

dat_wqs_pos$lci <- dat_wqs_pos$estimate - 1.96*dat_wqs_pos$se
dat_wqs_pos$uci <- dat_wqs_pos$estimate + 1.96*dat_wqs_pos$se

dat_wqs_pos$hror <- exp(dat_wqs_pos$estimate)
dat_wqs_pos$lci.hror <- exp(dat_wqs_pos$lci)
dat_wqs_pos$uci.hror <- exp(dat_wqs_pos$uci)

dat_wqs_pos$text_beta <- paste0(sprintf("%.3f", dat_wqs_pos$estimate)," (",sprintf("%.3f", dat_wqs_pos$lci),", ",sprintf("%.3f", dat_wqs_pos$uci),")")
dat_wqs_pos$text_hror <- paste0(sprintf("%.3f", dat_wqs_pos$hror)," (",sprintf("%.3f", dat_wqs_pos$lci.hror),", ",sprintf("%.3f", dat_wqs_pos$uci.hror),")")

dat_wqs_pos$text_beta <- ifelse(!dat_wqs_pos$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                      "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_wqs_pos$text_beta, "/")
dat_wqs_pos$text_hror <- ifelse(dat_wqs_pos$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                     "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_wqs_pos$text_hror, "/")

dat_wqs_pos <- dat_wqs_pos[,c("exp","out","text_hror","text_beta","z","p","p_adj_bh",
                              paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                       "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                                       "TCC","TCS",
                                       "BPA"),"_log10"))]

# 第8到21列：NA → -10000
dat_wqs_pos[, 8:21][is.na(dat_wqs_pos[, 8:21])] <- -10000
dat_wqs_pos[, 8:21] <- round(dat_wqs_pos[, 8:21], 3)
dat_wqs_pos[, 8:21] <- lapply(dat_wqs_pos[, 8:21], as.character)
dat_wqs_pos[dat_wqs_pos == -10000] <- "/"


# WQS negative model #
dat_wqs_neg$exp <- factor(dat_wqs_neg$exp, levels = exposure, labels = exposure_label)
dat_wqs_neg$out <- factor(dat_wqs_neg$out, levels = outcome_final, labels = outcome_label)
dat_wqs_neg <- dat_wqs_neg %>% arrange(out,exp)

dat_wqs_neg$lci <- dat_wqs_neg$estimate - 1.96*dat_wqs_neg$se
dat_wqs_neg$uci <- dat_wqs_neg$estimate + 1.96*dat_wqs_neg$se

dat_wqs_neg$hror <- exp(dat_wqs_neg$estimate)
dat_wqs_neg$lci.hror <- exp(dat_wqs_neg$lci)
dat_wqs_neg$uci.hror <- exp(dat_wqs_neg$uci)

dat_wqs_neg$text_beta <- paste0(sprintf("%.3f", dat_wqs_neg$estimate)," (",sprintf("%.3f", dat_wqs_neg$lci),", ",sprintf("%.3f", dat_wqs_neg$uci),")")
dat_wqs_neg$text_hror <- paste0(sprintf("%.3f", dat_wqs_neg$hror)," (",sprintf("%.3f", dat_wqs_neg$lci.hror),", ",sprintf("%.3f", dat_wqs_neg$uci.hror),")")

dat_wqs_neg$text_beta <- ifelse(!dat_wqs_neg$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                        "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_wqs_neg$text_beta, "/")
dat_wqs_neg$text_hror <- ifelse(dat_wqs_neg$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                       "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_wqs_neg$text_hror, "/")

dat_wqs_neg <- dat_wqs_neg[,c("exp","out","text_hror","text_beta","z","p","p_adj_bh",
                              paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                       "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                                       "TCC","TCS",
                                       "BPA"),"_log10"))]

# 第8到21列：NA → -10000
dat_wqs_neg[, 8:21][is.na(dat_wqs_neg[, 8:21])] <- -10000
dat_wqs_neg[, 8:21] <- round(dat_wqs_neg[, 8:21], 3)
dat_wqs_neg[, 8:21] <- lapply(dat_wqs_neg[, 8:21], as.character)
dat_wqs_neg[dat_wqs_neg == -10000] <- "/"
#### 数据处理 (for stable) ####


colnames(dat_wqs_pos) <- c("Analyte exposure in mixtures","Clinical biomarkers and outcomes","OR (95% CI)","Beta (95% CI)","Z value","P value","BH-adjusted P",
                           "PFOS","PFOA","PFNA","PFDA","PFHxS",
                           "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                           "TCC","TCS",
                           "BPA")
dat_wqs_pos <- dat_wqs_pos[!dat_wqs_pos$`Clinical biomarkers and outcomes` %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"),]

openxlsx::write.xlsx(dat_wqs_pos,paste0("tables/(stable09_1)_wqs_edc_out_(positive_model).xlsx"))


colnames(dat_wqs_neg) <- c("Analyte exposure in mixtures","Clinical biomarkers and outcomes","OR (95% CI)","Beta (95% CI)","Z value","P value","BH-adjusted P",
                           "PFOS","PFOA","PFNA","PFDA","PFHxS",
                           "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                           "TCC","TCS",
                           "BPA")
dat_wqs_neg <- dat_wqs_neg[!dat_wqs_neg$`Clinical biomarkers and outcomes` %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"),]

openxlsx::write.xlsx(dat_wqs_neg,paste0("tables/(stable09_2)_wqs_edc_out_(negative_model).xlsx"))
