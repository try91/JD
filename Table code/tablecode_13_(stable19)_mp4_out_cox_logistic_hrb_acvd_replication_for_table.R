library(data.table)
library(dplyr)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)


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


#### 构建函数把species_name转换为可以作图的标准化名称 ----
trans_mp4_s_names <- function(STRING){
  
  string_temp <- gsub("s__","",STRING)
  string_temp <- gsub("_log10","",string_temp)
  # # 删除所有"_" (一般情况下第一个"_"分隔属和种)
  # string_temp <- gsub("_"," ",string_temp)
  # 删除"_" (一般情况下第一个"_"分隔属和种)
  string_temp <- sub("_"," ",string_temp) # 使用 sub()替换第一个匹配项。
  # 删除特殊位置的"_"
  string_temp <- gsub("_bacterium"," bacterium",string_temp)
  string_temp <- gsub("bacterium_","bacterium ",string_temp)
  string_temp <- gsub("_unclassified"," unclassified",string_temp)
  string_temp <- gsub("unclassified_","unclassified ",string_temp)
  string_temp <- gsub("copri_","copri ",string_temp)
  string_temp <- gsub("Cibionibacter_","Cibionibacter ",string_temp)
  string_temp <- gsub("Marseille_","Marseille ",string_temp)
  string_temp <- gsub("oral_taxon_","oral taxon ",string_temp)
  # # 添加特殊位置的"_"
  # string_temp <- gsub(" Family","_Family",string_temp)
  # string_temp <- gsub(" phylum","_phylum",string_temp)
  # string_temp <- gsub("Candidatus ","Candidatus_",string_temp)
  # 处理sp
  string_temp <- gsub("sp_","sp. ",string_temp)
  string_temp <- gsub("_sp."," sp.",string_temp)
  return(string_temp)
}
#### 构建函数把species_name转换为可以作图的标准化名称 ####

#### 读取 (JD HRB ACVD) Replication结果 ----
## JD结果 ##
{
  ### 整理 mp4-out 分析结果 ###
  {
    #### 读取MP4-Outcomes COX & logistic分析结果 ----
    cox_results <- readxl::read_xlsx("results/cox/cox_results_mp4_incident.xlsx")
    cox_results <- cox_results[cox_results$adjust == "adj" & cox_results$outcome == "cvd_incident_1421",]
    cox_results <- cox_results[,c(1,3:13)]
    colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
    out_cox <- unique(cox_results$outcome) # 提取结局变量
    exp_cox <- unique(cox_results$exposure) # 提取暴露变量
    
    logistic_results <- readxl::read_xlsx("results/glm/logistic_results_mp4_incident.xlsx")
    logistic_results <- logistic_results[logistic_results$adjust == "adj" & logistic_results$outcome %in% c("ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"),]
    logistic_results <- logistic_results[logistic_results$adjust == "adj",]
    colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
    out_logistic <- unique(logistic_results$outcome) # 提取结局变量
    exp_logistic <- unique(logistic_results$exposure) # 提取暴露变量
    
    cox_logistic_results <- rbind(cox_results,logistic_results)
    cox_logistic_results <- cox_logistic_results[cox_logistic_results$rowname %in% c(mp4_s_log10),]
    # 以每个OUT表型为单位进行校正（以outcome为组，校正每个菌）
    cox_logistic_results <- cox_logistic_results %>%
      group_by(outcome) %>%  # 按outcome分组
      mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
      ungroup()
    cox_logistic_short <- cox_logistic_results[,c("rowname","outcome","estimate","se","p","p_adj_bh")]
    colnames(cox_logistic_short) <- c("exp_name","out_name","estimate_cox_logistic","se_cox_logistic","p","p_adj_bh")
    
    cox_logistic_short$out_name <- factor(cox_logistic_short$out_name, levels = c("cvd_incident_1421","ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"))
    cox_logistic_short <- cox_logistic_short %>%
      arrange(exp_name,out_name)
    
    # effect size 转换为 HR/OR
    cox_logistic_short$exp_estimate <- exp(cox_logistic_short$estimate_cox_logistic)
    
    # cox 和 logistic 显著的菌-Outcome pairs
    pair_mp4_out_sig <- cox_logistic_short[cox_logistic_short$p < 0.05,] # (p显著)
    unique(pair_mp4_out_sig$exp_name) # 127个菌
    pair_mp4_cvd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "cvd_incident_1421",] # 26个CVD显著相关的菌种
    pair_mp4_cvd_sig <- pair_mp4_cvd_sig[,"exp_name"]
    pair_mp4_cvd_sig$cvd_related <- 1
    pair_mp4_ckd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "ckd_incident_1014_no_self_report",] # 84个CKD显著相关的菌种
    pair_mp4_ckd_sig <- pair_mp4_ckd_sig[,"exp_name"]
    pair_mp4_ckd_sig$ckd_related <- 1
    pair_mp4_dm_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "dm_incident_1014_no_self_report",] # 35个DM显著相关的菌种
    pair_mp4_dm_sig <- pair_mp4_dm_sig[,"exp_name"]
    pair_mp4_dm_sig$dm_related <- 1
    #### 读取MP4-Outcomes COX & logistic分析结果 ####
    
    
    mp4_s <- data.frame(mp4_s_log10)
    mp4_s_matched <- left_join(mp4_s, pair_mp4_dm_sig, by=c("mp4_s_log10" = "exp_name")) %>%
      left_join(pair_mp4_ckd_sig, by=c("mp4_s_log10" = "exp_name")) %>%
      left_join(pair_mp4_cvd_sig, by=c("mp4_s_log10" = "exp_name"))
    # 保留第2到11列中，至少有一个不是NA的行
    mp4_s_cleaned <- mp4_s_matched[rowSums(!is.na(mp4_s_matched[, 2:4])) > 0, ]
    
    
    cox_logistic_short_dm_selected <- cox_logistic_short[cox_logistic_short$exp_name %in% mp4_s_cleaned$mp4_s_log10 & cox_logistic_short$out_name == "dm_incident_1014_no_self_report",]
    colnames(cox_logistic_short_dm_selected) <- c("exp_name","out_name","estimate_dm","se_dm","p_dm","p_adj_bh_dm","hr/or_dm")
    cox_logistic_short_dm_selected <- cox_logistic_short_dm_selected[,c("exp_name","estimate_dm","se_dm","p_dm")]
    
    cox_logistic_short_ckd_selected <- cox_logistic_short[cox_logistic_short$exp_name %in% mp4_s_cleaned$mp4_s_log10 & cox_logistic_short$out_name == "ckd_incident_1014_no_self_report",]
    colnames(cox_logistic_short_ckd_selected) <- c("exp_name","out_name","estimate_ckd","se_ckd","p_ckd","p_adj_bh_ckd","hr/or_ckd")
    cox_logistic_short_ckd_selected <- cox_logistic_short_ckd_selected[,c("exp_name","estimate_ckd","se_ckd","p_ckd")]
    
    cox_logistic_short_cvd_selected <- cox_logistic_short[cox_logistic_short$exp_name %in% mp4_s_cleaned$mp4_s_log10 & cox_logistic_short$out_name == "cvd_incident_1421",]
    colnames(cox_logistic_short_cvd_selected) <- c("exp_name","out_name","estimate_cvd","se_cvd","p_cvd","p_adj_bh_cvd","hr/or_cvd")
    cox_logistic_short_cvd_selected <- cox_logistic_short_cvd_selected[,c("exp_name","estimate_cvd","se_cvd","p_cvd")]
    
    
    mp4_s_matched2 <- left_join(mp4_s, cox_logistic_short_dm_selected, by=c("mp4_s_log10" = "exp_name")) %>%
      left_join(cox_logistic_short_ckd_selected, by=c("mp4_s_log10" = "exp_name")) %>%
      left_join(cox_logistic_short_cvd_selected, by=c("mp4_s_log10" = "exp_name"))
    # 保留第2到20列中，至少有一个不是NA的行
    mp4_s_cleaned2 <- mp4_s_matched2[rowSums(!is.na(mp4_s_matched2[, 2:10])) > 0, ]
    
    mp4_s_cleaned2$mp4_s_log10 <- gsub("_log10","",mp4_s_cleaned2$mp4_s_log10)
    colnames(mp4_s_cleaned2)[1] <- "species"
    
    mp4_s_for_rep <- mp4_s_cleaned2
    row.names(mp4_s_for_rep) <- NULL
    mp4_s_for_rep[is.na(mp4_s_for_rep)] <- 0
  }
  results_jd <- mp4_s_for_rep
  ### 整理 mp4-out 分析结果 ###
  
  dm_jd <- results_jd[,c("species","estimate_dm","se_dm","p_dm")]
  colnames(dm_jd) <- c("species","estimate_dm_jd","se_dm_jd","p_dm_jd")
  dm_jd$estimate_dm_jd <- as.numeric(dm_jd$estimate_dm_jd)
  dm_jd$p_dm_jd <- as.numeric(dm_jd$p_dm_jd)
  nrow(dm_jd[dm_jd$p_dm_jd < 0.05,]) # 35 dm-related sig species in JD (no self-report cases)
  
  ckd_jd <- results_jd[,c("species","estimate_ckd","se_ckd","p_ckd")]
  colnames(ckd_jd) <- c("species","estimate_ckd_jd","se_ckd_jd","p_ckd_jd")
  ckd_jd$estimate_ckd_jd <- as.numeric(ckd_jd$estimate_ckd_jd)
  ckd_jd$p_ckd_jd <- as.numeric(ckd_jd$p_ckd_jd)
  nrow(ckd_jd[ckd_jd$p_ckd_jd < 0.05,]) # 84 ckd-related sig species in JD
  
  cvd_jd <- results_jd[,c("species","estimate_cvd","se_cvd","p_cvd")]
  colnames(cvd_jd) <- c("species","estimate_cvd_jd","se_cvd_jd","p_cvd_jd")
  cvd_jd$estimate_cvd_jd <- as.numeric(cvd_jd$estimate_cvd_jd)
  cvd_jd$p_cvd_jd <- as.numeric(cvd_jd$p_cvd_jd)
  nrow(cvd_jd[cvd_jd$p_cvd_jd < 0.05,]) # 26 cvd-related sig species in JD
}

## HEB结果 ##
{
  dm_hrb <- readxl::read_xlsx("raw_data/validation_DM_n127.xlsx")
  ckd_hrb <- readxl::read_xlsx("raw_data/validation_CKD_n127.xlsx")
}

## JieZ ACVD结果 ##
{
  cvd_gd <- readxl::read_xlsx("raw_data/validation_ASCVD_n127.xlsx")
}
#### 读取 (JD HRB ACVD) Replication结果 ####

#### (JD HRB ACVD) Replication结果处理 ----
# DM结果处理 #
{
  results_dm <- left_join(dm_jd, dm_hrb, by = "species") %>%
    filter(p_dm_jd < 0.05) %>%
    arrange(desc(estimate_dm_jd))
  
  results_dm$outcome <- "Diabetes"
  results_dm <- results_dm[,c("species","outcome","estimate_dm_hrb","se_dm_hrb","p_dm_hrb","p_fdr_dm_hrb")]
  
  colnames(results_dm) <- c("exp_name","out_name","estimate","se","p","p_adj_bh")
}

# CKD结果处理 #
{
  results_ckd <- left_join(ckd_jd, ckd_hrb, by = "species") %>%
    filter(p_ckd_jd < 0.05) %>%
    arrange(desc(estimate_ckd_jd))
  
  results_ckd$outcome <- "CKD"
  results_ckd <- results_ckd[,c("species","outcome","estimate_ckd_hrb","se_ckd_hrb","p_ckd_hrb","p_fdr_ckd_hrb")]
  
  colnames(results_ckd) <- c("exp_name","out_name","estimate","se","p","p_adj_bh")
}

# ASCVD结果处理 #
{
  results_acvd <- left_join(cvd_jd, cvd_gd, by = "species") %>%
    filter(p_cvd_jd < 0.05) %>%
    arrange(desc(estimate_cvd_jd))
  
  results_acvd$outcome <- "CVD"
  results_acvd <- results_acvd[,c("species","outcome","estimate_cvd_gd","se_cvd_gd","p_cvd_gd","p_fdr_cvd_gd")]
  
  colnames(results_acvd) <- c("exp_name","out_name","estimate","se","p","p_adj_bh")
}

# 合并显著的DM, CKD, CVD验证结果
results_all <- rbind(results_dm, results_ckd, results_acvd)
#### (JD HRB ACVD) Replication结果处理 ----

results_all$lci <- results_all$estimate - 1.96*results_all$se
results_all$uci <- results_all$estimate + 1.96*results_all$se

results_all$or <- exp(results_all$estimate)
results_all$lci.or <- exp(results_all$lci)
results_all$uci.or <- exp(results_all$uci)

# 如果p<0.05需要确保UCI<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999
test <- results_all[results_all$p < 0.05 & (results_all$uci.or >= 0.9995 & results_all$uci.or < 1),]
results_all$uci.or <- ifelse(results_all$p < 0.05 & (results_all$uci.or >= 0.9995 & results_all$uci.or < 1), 0.999, results_all$uci.or)

results_all$text_or <- paste0(sprintf("%.3f", results_all$or)," (",sprintf("%.3f", results_all$lci.or),", ",sprintf("%.3f", results_all$uci.or),")")

results_all$exp_name <- factor(results_all$exp_name, 
                               levels = mp4_s_names,
                               labels = trans_mp4_s_names(mp4_s_names))
results_all$out_name <- factor(results_all$out_name, 
                               levels = c("Diabetes","CKD","CVD"))
results_all <- results_all%>%
  arrange(out_name,exp_name)

dat_for_stable <- results_all[,c("exp_name","out_name","text_or","p","p_adj_bh")]

colnames(dat_for_stable) <- c("Gut microbial species","Outcomes","OR (95% CI)","P value","BH-adjusted P")
openxlsx::write.xlsx(dat_for_stable,"tables/(stable19)_mp4_out_cox_logistic_hrb_acvd_replication.xlsx")
