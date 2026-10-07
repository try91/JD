library(data.table)
library(dplyr)
library(ggplot2)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

#### 变量整理 ----
# 菌群2014菌群 (分类和连续)
# 丰度>0.0001, 出现率>10%的微生物 (物种和属)
mp4_s_names <- read.table("jiading/sourceDataTaxon/mpa4/species_names_mp4_10%.txt")
mp4_s_names <- mp4_s_names[,1]
mp4_g_names <- read.table("jiading/sourceDataTaxon/mpa4/genus_names_mp4_10%.txt")
mp4_g_names <- mp4_g_names[,1]
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp4_s_names_short <- mp4_s_names[!grepl("_GGB",mp4_s_names)] # 排除未分类的菌属（GGB）, 未分类菌种（SGB）先保留
mp4_g_names_short <- mp4_g_names[!grepl("_GGB",mp4_g_names)] # 排除未分类的菌属（GGB）
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp3_s_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_10%.txt")
mp3_s_names <- mp3_s_names[,1]
mp3_g_names <- read.table("jiading/sourceDataTaxon/mpa3/genus_names_mp3_10%.txt")
mp3_g_names <- mp3_g_names[,1]
# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)
mp4_g_zero <- paste0(mp4_g_names,"_zero") # 菌群MP4填补0值丰度 (属层面)

mp3_s_bin <- paste0(mp3_s_names,"_bin") # 菌群MP3出现与否的分类变量 (物种层面)
mp3_s_log10 <- paste0(mp3_s_names,"_log10") # 菌群MP3丰度的log10转换 (物种层面)
mp3_s_zero <- paste0(mp3_s_names,"_zero") # 菌群MP3填补0值丰度 (物种层面)

mp3_g_bin <- paste0(mp3_g_names,"_bin") # 菌群MP3出现与否的分类变量 (属层面)
mp3_g_log10 <- paste0(mp3_g_names,"_log10") # 菌群MP3丰度的log10转换 (属层面)
mp3_g_zero <- paste0(mp3_g_names,"_zero") # 菌群MP3填补0值丰度 (属层面)
# # mp3中用于构建ma的菌的名称
# mp3_ma_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_ma.txt")
# mp3_ma_names <- mp3_ma_names$V1
# mp3_ma_names_log10 <- paste0(mp3_ma_names,"_log10")
# # microbial age (MA)
# # mean(phy_edc$MA,na.rm = TRUE)

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
phy_out_cat <- c("cvd_f","ckd_f","dm_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014新发case，横断面数据)
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

#### EDC INDEX 变量名 ----
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")
edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### EDC INDEX 变量名 ####

#### 结局变量名汇总 ----
## 结局 横断面研究结果
select_out_cat <- c("dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b")
# 标准化名称
select_out_cat_labels <- c("Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT")
#### 结局变量名汇总 ####


# 总人群结果
for (select_sample in c("phy_edc_temp")) {
  # select_sample <- "phy_edc_temp"
  
  #### 数据处理 ----
  ## 读取COX & Logistic分析结果
  cox_results_edc_incident_all <- readxl::read_xlsx("results/cox/cox_results_edc_incident_20260702.xlsx")
  cox_results_edc_index_incident_all <- readxl::read_xlsx("results/cox/cox_results_edc_index_incident_20260702.xlsx")
  
  logistic_results_edc_incident_all <- readxl::read_xlsx("results/glm/logistic_results_edc_incident_20260702.xlsx")
  logistic_results_edc_index_incident_all <- readxl::read_xlsx("results/glm/logistic_results_edc_index_incident_20260702.xlsx")
  
  ## 读取Qgcomp分析结果
  qg_results_q2 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(phy_edc_temp)_20260702.xlsx"))
  qg_results_q2_pn <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(phy_edc_temp)_20260702.xlsx"))
  
  ## 整合分析结果
  cox_results <- rbind(cox_results_edc_incident_all,cox_results_edc_index_incident_all)
  cox_results <- cox_results[cox_results$sample == select_sample,]
  cox_results <- cox_results[cox_results$adjust == "adj",]
  cox_results <- cox_results[cox_results$outcome %in% phy_incident_cat,]
  cox_results <- cox_results[cox_results$rowname %in% edc_index_b_keep,]
  cox_results <- cox_results[,c(1,3:13)]
  colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
  out_cox <- unique(cox_results$outcome) # 提取结局变量
  
  logistic_results <- rbind(logistic_results_edc_incident_all,logistic_results_edc_index_incident_all)
  logistic_results <- logistic_results[logistic_results$sample == select_sample,]
  logistic_results <- logistic_results[logistic_results$adjust == "adj",]
  logistic_results <- logistic_results[logistic_results$outcome %in% select_out_cat,]
  logistic_results <- logistic_results[logistic_results$rowname %in% edc_index_b_keep,]
  colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
  out_logistic <- unique(logistic_results$outcome) # 提取结局变量
  
  
  qg_results_q2 <- rbind(qg_results_q2,qg_results_q2_pn)
  qg_results_q2 <- qg_results_q2[qg_results_q2$sample == select_sample,]
  qg_results_q2 <- qg_results_q2[,c(1:8)]
  qg_results_q2$z <- qg_results_q2$estimate / qg_results_q2$se 
  colnames(qg_results_q2)[1:2] <- c("rowname","outcome")
  qg_results_q2$outcome <- ifelse((qg_results_q2$outcome %in% c(phy_incident_cat,phy_censor_cat)) & (qg_results_q2$method == "qgcomp bin"), paste0(qg_results_q2$outcome,"_logistic"), qg_results_q2$outcome)
  qg_results_q2 <- qg_results_q2[qg_results_q2$outcome %in% c(phy_incident_cat,select_out_cat),]
  qg_results_q2 <- qg_results_q2[qg_results_q2$rowname %in% c("EDC_14", "PFAS", "PAE_6", "BP_1", "TC"),]
  
  
  dat_result1 <- rbind(cox_results,logistic_results)
  dat_result1_1 <- dat_result1[dat_result1$rowname %in% c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_1 <- dat_result1_1 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_2 <- dat_result1[dat_result1$rowname %in% c("edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_2 <- dat_result1_2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  # 以每个OUT表型为单位进行校正
  qg_results_q2 <- qg_results_q2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_1$method_type <- "edc_count_index"
  dat_result1_2$method_type <- "edc_score"
  qg_results_q2$method_type <- "qgcomp"
  
  
  table(dat_result1_1$rowname,dat_result1_1$outcome)
  table(dat_result1_2$rowname,dat_result1_2$outcome)
  table(qg_results_q2$rowname,qg_results_q2$outcome)
  dat_result_all <- bind_rows(dat_result1_1,dat_result1_2,qg_results_q2)
  
  
  dat_result_all$exposure_type <- ifelse(dat_result_all$rowname %in% c("EDC_14","edc_count2_edc14_b","edc_score_edc14_b"), "EDCs (14)",
                                         ifelse(dat_result_all$rowname %in% c("PFAS","edc_count2_pfas_b","edc_score_pfas_b"), "PFAS (5)",
                                                ifelse(dat_result_all$rowname %in% c("PAE_6","edc_count2_pae6_b","edc_score_pae6_b"), "PAEs (6)",
                                                       ifelse(dat_result_all$rowname %in% c("TC","edc_count2_tc_b","edc_score_tc_b"), "Antimicrobials (2)", "Bisphenols (1)"))))
  dat_result_all$exposure_type <- factor(dat_result_all$exposure_type, levels = c("EDCs (14)", "PFAS (5)", "PAEs (6)", "Antimicrobials (2)", "Bisphenols (1)"))
  
  dat_result_all$rowname <- factor(dat_result_all$rowname, levels = c("EDC_14", "PFAS", "PAE_6", "TC", "BP_1",
                                                                      "edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b",
                                                                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),
                                   labels = c("qg-comp (14 EDCs)", "qg-comp (PFAS)", "qg-comp (PAEs)", "qg-comp (antimicrobials)", "qg-comp (bisphenols)",
                                              "EDC Scoremedian (14 EDCs)", "EDC Scoremedian (PFAS)", "EDC Scoremedian (PAEs)", "EDC Scoremedian (antimicrobials)", "EDC Scoremedian (bisphenols)",
                                              "EDC Scorequartile (14 EDCs)", "EDC Scorequartile (PFAS)", "EDC Scorequartile (PAEs)", "EDC Scorequartile (antimicrobials)", "EDC Scorequartile (bisphenols)"))
  
  dat_result_all$outcome <- factor(dat_result_all$outcome,
                                   levels = c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021",select_out_cat),
                                   labels = c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",select_out_cat_labels))
  
  dat_result_all$method_type <- factor(dat_result_all$method_type, 
                                       levels = c("qgcomp", "edc_count_index", "edc_score"),
                                       labels = c("Mixture effectqg-comp","EDC Scoremedian","EDC Scorequartile"))
  
  dat_result_all <- dat_result_all %>%
    arrange(outcome,exposure_type,method_type)
  
  # 计算95%CI
  dat_result_all$z.lci <- dat_result_all$z - 1.96
  dat_result_all$z.uci <- dat_result_all$z + 1.96
  
  dat_result_all$lci <- dat_result_all$estimate - 1.96*dat_result_all$se
  dat_result_all$uci <- dat_result_all$estimate + 1.96*dat_result_all$se
  
  dat_result_all$hror <- exp(dat_result_all$estimate)
  dat_result_all$lci.hror <- exp(dat_result_all$lci)
  dat_result_all$uci.hror <- exp(dat_result_all$uci)
  
  ### 如果p<0.05需要确保uci<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_result_all[dat_result_all$p < 0.05 & (dat_result_all$uci.hror >= 0.9995 & dat_result_all$uci.hror < 1),]
    dat_result_all$uci.hror <- ifelse(dat_result_all$p < 0.05 & (dat_result_all$uci.hror >= 0.9995 & dat_result_all$uci.hror < 1), 0.999, dat_result_all$uci.hror)
  }
  
  dat_result_all$text_beta <- paste0(sprintf("%.3f", dat_result_all$estimate)," (",sprintf("%.3f", dat_result_all$lci),", ",sprintf("%.3f", dat_result_all$uci),")")
  dat_result_all$text_hror <- paste0(sprintf("%.3f", dat_result_all$hror)," (",sprintf("%.3f", dat_result_all$lci.hror),", ",sprintf("%.3f", dat_result_all$uci.hror),")")
  dat_result_all$text_hr <- ifelse(dat_result_all$outcome %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"), dat_result_all$text_hror, "/")
  dat_result_all$text_or <- ifelse(!dat_result_all$outcome %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"), dat_result_all$text_hror, "/")
  
  dat_result_all <- dat_result_all[,c("exposure_type","method_type","outcome","text_hr","text_or","z","p","p_adj_bh")]
  #### 数据处理 ####
  
  colnames(dat_result_all) <- c("Classes of analytes","Exposure types","Outcomes","HR (95% CI)","OR (95% CI)","Z value","P value","BH-adjusted P")
  openxlsx::write.xlsx(dat_result_all,paste0("tables/compare_qgcomp_edc_score_(",select_sample,")_20260702.xlsx"))
}

# 性别特异性人群结果
for (select_sample in c("phy_edc_temp0_1","phy_edc_temp0_2")) {
  # select_sample <- "phy_edc_temp0_1"
  
  #### 数据处理 ----
  ## 读取COX & Logistic分析结果
  cox_results_edc_incident_all <- readxl::read_xlsx("results/cox/sex_specific_cox_results_edc_incident_20260702.xlsx")
  cox_results_edc_index_incident_all <- readxl::read_xlsx("results/cox/sex_specific_cox_results_edc_index_incident_20260702.xlsx")
  
  logistic_results_edc_incident_all <- readxl::read_xlsx("results/glm/sex_specific_logistic_results_edc_incident_20260702.xlsx")
  logistic_results_edc_index_incident_all <- readxl::read_xlsx("results/glm/sex_specific_logistic_results_edc_index_incident_20260702.xlsx")
  
  ## 读取Qgcomp分析结果
  qg_results_q2 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",select_sample,")_20260702.xlsx"))
  qg_results_q2_pn <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(",select_sample,")_20260702.xlsx"))
  
  
  
  
  ## 整合分析结果
  cox_results <- rbind(cox_results_edc_incident_all,cox_results_edc_index_incident_all)
  cox_results <- cox_results[cox_results$sample == select_sample,]
  cox_results <- cox_results[cox_results$adjust == "adj",]
  cox_results <- cox_results[cox_results$outcome %in% phy_incident_cat,]
  cox_results <- cox_results[cox_results$rowname %in% edc_index_b_keep,]
  cox_results <- cox_results[,c(1,3:13)]
  colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
  out_cox <- unique(cox_results$outcome) # 提取结局变量
  
  logistic_results <- rbind(logistic_results_edc_incident_all,logistic_results_edc_index_incident_all)
  logistic_results <- logistic_results[logistic_results$sample == select_sample,]
  logistic_results <- logistic_results[logistic_results$adjust == "adj",]
  logistic_results <- logistic_results[logistic_results$outcome %in% select_out_cat,]
  logistic_results <- logistic_results[logistic_results$rowname %in% edc_index_b_keep,]
  colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
  out_logistic <- unique(logistic_results$outcome) # 提取结局变量
  
  
  qg_results_q2 <- rbind(qg_results_q2,qg_results_q2_pn)
  qg_results_q2 <- qg_results_q2[qg_results_q2$sample == select_sample,]
  qg_results_q2 <- qg_results_q2[,c(1:8)]
  qg_results_q2$z <- qg_results_q2$estimate / qg_results_q2$se 
  colnames(qg_results_q2)[1:2] <- c("rowname","outcome")
  qg_results_q2$outcome <- ifelse((qg_results_q2$outcome %in% c(phy_incident_cat,phy_censor_cat)) & (qg_results_q2$method == "qgcomp bin"), paste0(qg_results_q2$outcome,"_logistic"), qg_results_q2$outcome)
  qg_results_q2 <- qg_results_q2[qg_results_q2$outcome %in% c(phy_incident_cat,select_out_cat),]
  qg_results_q2 <- qg_results_q2[qg_results_q2$rowname %in% c("EDC_14", "PFAS", "PAE_6", "BP_1", "TC"),]
  
  
  dat_result1 <- rbind(cox_results,logistic_results)
  dat_result1_1 <- dat_result1[dat_result1$rowname %in% c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_1 <- dat_result1_1 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_2 <- dat_result1[dat_result1$rowname %in% c("edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_2 <- dat_result1_2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  # 以每个OUT表型为单位进行校正
  qg_results_q2 <- qg_results_q2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_1$method_type <- "edc_count_index"
  dat_result1_2$method_type <- "edc_score"
  qg_results_q2$method_type <- "qgcomp"
  
  
  table(dat_result1_1$rowname,dat_result1_1$outcome)
  table(dat_result1_2$rowname,dat_result1_2$outcome)
  table(qg_results_q2$rowname,qg_results_q2$outcome)
  dat_result_all <- bind_rows(dat_result1_1,dat_result1_2,qg_results_q2)
  
  
  dat_result_all$exposure_type <- ifelse(dat_result_all$rowname %in% c("EDC_14","edc_count2_edc14_b","edc_score_edc14_b"), "EDCs (14)",
                                         ifelse(dat_result_all$rowname %in% c("PFAS","edc_count2_pfas_b","edc_score_pfas_b"), "PFAS (5)",
                                                ifelse(dat_result_all$rowname %in% c("PAE_6","edc_count2_pae6_b","edc_score_pae6_b"), "PAEs (6)",
                                                       ifelse(dat_result_all$rowname %in% c("TC","edc_count2_tc_b","edc_score_tc_b"), "Antimicrobials (2)", "Bisphenols (1)"))))
  dat_result_all$exposure_type <- factor(dat_result_all$exposure_type, levels = c("EDCs (14)", "PFAS (5)", "PAEs (6)", "Antimicrobials (2)", "Bisphenols (1)"))
  
  dat_result_all$rowname <- factor(dat_result_all$rowname, levels = c("EDC_14", "PFAS", "PAE_6", "TC", "BP_1",
                                                                      "edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b",
                                                                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),
                                   labels = c("qg-comp (14 EDCs)", "qg-comp (PFAS)", "qg-comp (PAEs)", "qg-comp (antimicrobials)", "qg-comp (bisphenols)",
                                              "EDC Scoremedian (14 EDCs)", "EDC Scoremedian (PFAS)", "EDC Scoremedian (PAEs)", "EDC Scoremedian (antimicrobials)", "EDC Scoremedian (bisphenols)",
                                              "EDC Scorequartile (14 EDCs)", "EDC Scorequartile (PFAS)", "EDC Scorequartile (PAEs)", "EDC Scorequartile (antimicrobials)", "EDC Scorequartile (bisphenols)"))
  
  dat_result_all$outcome <- factor(dat_result_all$outcome,
                                   levels = c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021",select_out_cat),
                                   labels = c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",select_out_cat_labels))
  
  dat_result_all$method_type <- factor(dat_result_all$method_type, 
                                       levels = c("qgcomp", "edc_count_index", "edc_score"),
                                       labels = c("Mixture effectqg-comp","EDC Scoremedian","EDC Scorequartile"))
  
  dat_result_all <- dat_result_all %>%
    arrange(outcome,exposure_type,method_type)
  
  # 计算95%CI
  dat_result_all$z.lci <- dat_result_all$z - 1.96
  dat_result_all$z.uci <- dat_result_all$z + 1.96
  
  dat_result_all$lci <- dat_result_all$estimate - 1.96*dat_result_all$se
  dat_result_all$uci <- dat_result_all$estimate + 1.96*dat_result_all$se
  
  dat_result_all$hror <- exp(dat_result_all$estimate)
  dat_result_all$lci.hror <- exp(dat_result_all$lci)
  dat_result_all$uci.hror <- exp(dat_result_all$uci)
  
  ### 如果p>0.05需要确保lci<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_result_all[dat_result_all$p >= 0.05 & (dat_result_all$lci.hror >= 0.9995 & dat_result_all$lci.hror < 1),]
    dat_result_all$lci.hror <- ifelse(dat_result_all$p >= 0.05 & (dat_result_all$lci.hror >= 0.9995 & dat_result_all$lci.hror < 1), 0.999, dat_result_all$lci.hror)
  }
  
  dat_result_all$text_beta <- paste0(sprintf("%.3f", dat_result_all$estimate)," (",sprintf("%.3f", dat_result_all$lci),", ",sprintf("%.3f", dat_result_all$uci),")")
  dat_result_all$text_hror <- paste0(sprintf("%.3f", dat_result_all$hror)," (",sprintf("%.3f", dat_result_all$lci.hror),", ",sprintf("%.3f", dat_result_all$uci.hror),")")
  dat_result_all$text_hr <- ifelse(dat_result_all$outcome %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"), dat_result_all$text_hror, "/")
  dat_result_all$text_or <- ifelse(!dat_result_all$outcome %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"), dat_result_all$text_hror, "/")
  
  dat_result_all <- dat_result_all[,c("exposure_type","method_type","outcome","text_hr","text_or","z","p","p_adj_bh")]
  #### 数据处理 ####
  
  colnames(dat_result_all) <- c("Classes of analytes","Exposure types","Outcomes","HR (95% CI)","OR (95% CI)","Z value","P value","BH-adjusted P")
  openxlsx::write.xlsx(dat_result_all,paste0("tables/compare_qgcomp_edc_score_(",select_sample,")_20260702.xlsx"))
}
