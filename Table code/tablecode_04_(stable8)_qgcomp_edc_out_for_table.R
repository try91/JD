library(data.table)
library(dplyr)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

#### 变量整理 ----
# 菌群2014菌群 (分类和连续)
# 丰度>0.0001, 出现率>10%的微生物 (物种和属)
mp4_s_names <- read.table("jiading/sourceDataTaxon/mpa4/species_names_mp4_10%.txt")
mp4_s_names <- mp4_s_names[,1]
mp4_g_names <- read.table("jiading/sourceDataTaxon/mpa4/genus_names_mp4_10%.txt")
mp4_g_names <- mp4_g_names[,1]

# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
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

# 2021、2014死亡和新发表型 (分类)
phy_incident_cat <- c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014")    # 新发 cvd, ckd, dm 去除基线 case (只做EDC对outcome，不做cvd_incident_1421)
phy_incident_time <- c("timecvd_1021","timecvd_1014","timeckd_1014","timedm_1014")
phy_censor_cat <- c("censorall_1021","censorall_1014")
phy_censor_time <- c("timeall_1021","timeall_1014")
# 2014、2010表型 (分类)
phy_out_cat <- c("cvd_f","ckd_f","dm_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014或2021新发case，横断面数据)
phy_traits_cat <- c("cvd_b","ckd_b","dm_b","as_imt_f","as_imt_b","hpt_f","hpt_b","nafld_f","nafld_b",
                    "ob_f","ob_b","abob_f","abob_b","dyslip_f","dyslip_b","hua_f","hua_b","ir_f","ir_b","mets_f","mets_b",
                    
                    # "smk1_f","drk1_f","paactive3_g_f",
                    
                    "sitduration_f","sitduration_b","sleeptg_f",
                    "dm_treat_f","dm_treat_b","hpt_treat_f","hpt_treat_b","hpl_treat_f","hpl_treat_b",
                    "diet_score_g_f","high_fruveg_f","low_ssb_f","low_meat_f","high_fish_f")
# 2014、2010表型 (连续)
phy_traits_cont <- c("bmi_f","bmi_b","wc_f","wc_b","hc_f","hc_b","whr_f","whr_b","height_f","height_b","weight_f","weight_b",
                     "hdl_f","hdl_b","ldl_f","ldl_b","apoa_f","apoa_b","apob_f","apob_b","chol_f","chol_b","tg_f","tg_b","nonhdl_f","nonhdl_b",
                     "alt_f","alt_b","ast_f","ast_b","ggt_f","ggt_b","scr_f","scr_b","egfr_f","egfr_b","acr_f","acr_b","ua_f","ua_b","bia_f","bia_b",
                     "glu0_f","glu0_b","glu120_f","glu120_b","vhba1c_f","vhba1c_b","ins0_f","ins0_b","ins120_f","ins120_b","homair_f","homair_b","homab_f","homab_b",
                     # "dmduration_f", "dmduration_b",
                     "sbp_f","sbp_b","dbp_f","dbp_b","pr_f","pr_b",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "wbc_f","wbc_b","crp_f",
                     "plt_f","plt_b","hgb_f","hgb_b","eos_f","lym_f","mon_f","neu_f",
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f",
                     
                     "sleept_f","sittimet_f","sittimet_b","sum_met_f","sum_met_b",
                     "alco_f","alco_b","diet_score_f")
# 2014药物 (分类)
# 二十类(所有)药物
med_cat20 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f","med_dm5_f","med_dm6_f","med_dm7_f",
               "med_hbp1_f","med_hbp2_f","med_hbp3_6_f","med_hbp4_f","med_hbp5_f",
               "med_lip1_f","med_lip2_f","med_lip3_f",
               "med_ua1_f","med_ua2_f",
               "med_thy1_f","med_thy2_f",
               "med_oth_f")
# 十类药物 (使用人数>20, 包括Statins)
med_cat10 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f",
               "med_hbp1_f","med_hbp2_f","med_hbp3_6_f","med_hbp4_f","med_hbp5_f",
               "med_lip1_f")
# 六类与菌群显著相关药物 (Sulfonylureas, Biguanides, Thiazolidinediones, AGIs, ARBs, Calcium antagonists) + Statins (MP4数据)
med_cat7 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f", 
              "med_hbp1_f","med_hbp4_f", 
              "med_lip1_f") 
# 五类与菌群显著相关药物 (Biguanides, Thiazolidinediones, AGIs, ARBs, Calcium antagonists) + Statins (MP3数据)
med_cat6 <- c("med_dm2_f","med_dm3_f","med_dm4_f",
              "med_hbp1_f","med_hbp4_f",
              "med_lip1_f")
# 汇总的所有10类、7类和6类药物
med_all <- c("med_all10","med_all7","med_all6")
#### 变量整理 ####


# # 总人群人群
# sample_name1 <- "phy_edc_temp" # 用总人群

# # 男性 (总人群, n=3930)
# sample_name1 <- "phy_edc_temp0_1"

# # 女性 (总人群, n=6394)
# sample_name1 <- "phy_edc_temp0_2"

for (sample_name1 in c("phy_edc_temp","phy_edc_temp0_1","phy_edc_temp0_2")) {
  
  #### 数据处理 (for stable) ----
  ### 设置纳入图片的暴露和结局
  ## 暴露：EDC
  exposure <- c("EDC_14","EDC_pos","EDC_neg","PFAS","PAE_6","TC","BP_1")
  # 标准化名称
  exposure_label <- c("EDCs (14)","EDCs (positive weight)","EDCs (negative weight)","PFAS (5)","PAEs (6)","Antimicrobials (2)","Bisphenols (1)")
  
  ## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
  outcome <- c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021",
               "dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b",
               "bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
               "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
               "alt_b","ast_b","ggt_b","bia_b",
               "egfr_b","scr_b","ua_b", #"acr_b",
               "glu0_b","glu120_b","vhba1c_b",
               "ins0_b","ins120_b","homair_b","homab_b",
               "sbp_b","dbp_b","pr_b", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
               "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
               "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
               "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")
  outcome_final <- outcome
  # 标准化名称
  outcome_label <- c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                     "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT",
                     "BMI","Height","Weight","WHR","WC","HC",
                     "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                     "ALT","AST","GGT","Bile acid",
                     "eGFR","Serum creatinine","UA", #"UACR",
                     "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                     "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                     "SBP","DBP","PR", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
                     "FT3","FT4","TSH","TPOAb","TgAb",
                     "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                     "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")
  
  # 读取数据
  results_qg_all1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",sample_name1,")_20260702.xlsx"))
  results_qg_edc1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(",sample_name1,")_20260702.xlsx"))
  
  results_qg <- rbind(results_qg_all1, results_qg_edc1) %>%
    filter(exp %in% exposure & out %in% outcome_final)
  results_qg <- results_qg[!(results_qg$out %in% c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014") & results_qg$method == "qgcomp bin"),] # 删除"cvd_incident_1021","ckd_incident_1014","dm_incident_1014"的logistic结果
  colnames(results_qg)[9:22] <- gsub("_log10","",colnames(results_qg)[9:22])
  
  results_qg$z <- results_qg$estimate/results_qg$se
  
  # FDR 校正
  # 以每个exposure表型为单位，校正outcome
  dat <- results_qg %>%
    group_by(exp) %>%  # 按exposure分组，校正outcome
    mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
    ungroup()
  
  # EDC 权重数据整理
  dat$numb_pos_weight <- rowSums(dat[, 9:22] > 0, na.rm = TRUE)# 计算每一行第9列到第22列中正数的数量
  dat$numb_neg_weight <- rowSums(dat[, 9:22] < 0, na.rm = TRUE)# 计算每一行第9列到第22列中负数的数量
  
  dat$weight_threshold_pos <- 1 / dat$numb_pos_weight
  dat$weight_threshold_neg <- 1 / dat$numb_neg_weight
  
  dat_edc <- dat[dat$exp == "EDC_14",]
  dat_edc_19_long_qgcomp <- tidyr::gather(dat_edc, edc, weight, 9:22, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
  
  dat_qgcomp <- dat
  dat_qgcomp$exp <- factor(dat_qgcomp$exp, levels = exposure, labels = exposure_label)
  dat_qgcomp$out <- factor(dat_qgcomp$out, levels = outcome_final, labels = outcome_label)
  dat_qgcomp <- dat_qgcomp %>% arrange(out,exp)
  
  
  dat_qgcomp$lci <- dat_qgcomp$estimate - 1.96*dat_qgcomp$se
  dat_qgcomp$uci <- dat_qgcomp$estimate + 1.96*dat_qgcomp$se
  
  dat_qgcomp$hror <- exp(dat_qgcomp$estimate)
  dat_qgcomp$lci.hror <- exp(dat_qgcomp$lci)
  dat_qgcomp$uci.hror <- exp(dat_qgcomp$uci)
  
  dat_qgcomp$text_beta <- paste0(sprintf("%.3f", dat_qgcomp$estimate)," (",sprintf("%.3f", dat_qgcomp$lci),", ",sprintf("%.3f", dat_qgcomp$uci),")")
  dat_qgcomp$text_hror <- paste0(sprintf("%.3f", dat_qgcomp$hror)," (",sprintf("%.3f", dat_qgcomp$lci.hror),", ",sprintf("%.3f", dat_qgcomp$uci.hror),")")
  
  dat_qgcomp$text_beta <- ifelse(!dat_qgcomp$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                        "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_qgcomp$text_beta, "/")
  dat_qgcomp$text_hror <- ifelse(dat_qgcomp$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",
                                                       "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT"), dat_qgcomp$text_hror, "/")
  
  dat_qgcomp <- dat_qgcomp[,c("exp","out","text_hror","text_beta","z","p","p_adj_bh",
                              "PFOS","PFOA","PFNA","PFDA","PFHxS",
                              "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                              "TCC","TCS",
                              "BPA")]
  
  # 第8到21列：NA → -10000
  dat_qgcomp[, 8:21][is.na(dat_qgcomp[, 8:21])] <- -10000
  dat_qgcomp[, 8:21] <- round(dat_qgcomp[, 8:21], 3)
  dat_qgcomp[, 8:21] <- lapply(dat_qgcomp[, 8:21], as.character)
  dat_qgcomp[dat_qgcomp == -10000] <- "/"
  #### 数据处理 (for stable) ####
  
  colnames(dat_qgcomp)[1:7] <- c("Analyte exposure in mixtures","Clinical biomarkers and outcomes","OR (95% CI)","Beta (95% CI)","Z value","P value","BH-adjusted P")
  # dat_qgcomp <- dat_qgcomp[dat_qgcomp$`Analyte exposure in mixtures` %in% c("EDCs (14)", "PFAS (5)", "PAEs (6)", "Antimicrobials (2)", "Bisphenols (1)"),]
  dat_qgcomp <- dat_qgcomp[!dat_qgcomp$`Clinical biomarkers and outcomes` %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"),]
  
  openxlsx::write.xlsx(dat_qgcomp,paste0("tables/qgcomp_edc_out_(",sample_name1,")_20260714.xlsx"))
}
