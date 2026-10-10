library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

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
    
  # 筛选JD和HRB方向一致的结果
  results_dm_same_direction <- results_dm[results_dm$estimate_dm_jd * results_dm$estimate_dm_hrb > 0,]
  results_dm_same_direction_hrb_sig_fdr <- results_dm_same_direction[results_dm_same_direction$p_fdr_dm_hrb < 0.2,]
  
  results_dm_estimate <- results_dm[,c("species","estimate_dm_jd","estimate_dm_hrb")]
  colnames(results_dm_estimate) <- c("species","JD","HRB")
  
  results_dm_se <- results_dm[,c("species","se_dm_jd","se_dm_hrb")]
  colnames(results_dm_se) <- c("species","JD","HRB")
  
  results_dm_p <- results_dm[,c("species","p_dm_jd","p_dm_hrb")]
  colnames(results_dm_p) <- c("species","JD","HRB")
  
  # 转换成long data
  results_dm_estimate_long <- results_dm_estimate %>%
    pivot_longer(
      cols = -species,                   
      names_to = "study",              
      values_to = "estimate"           
    )
  results_dm_se_long <- results_dm_se %>%
    pivot_longer(
      cols = -species,                   
      names_to = "study",              
      values_to = "se"          
    )
  results_dm_p_long <- results_dm_p %>%
    pivot_longer(
      cols = -species,                  
      names_to = "study",              
      values_to = "p"          
    )
  results_dm_long <- left_join(results_dm_estimate_long, results_dm_se_long, by=c("species","study")) %>%
    left_join(results_dm_p_long, by=c("species","study"))
}

# CKD结果处理 #
{
  results_ckd <- left_join(ckd_jd, ckd_hrb, by = "species") %>%
    filter(p_ckd_jd < 0.05) %>%
    arrange(desc(estimate_ckd_jd))
    
  # 筛选JD和HRB方向一致的结果
  results_ckd_same_direction <- results_ckd[results_ckd$estimate_ckd_jd * results_ckd$estimate_ckd_hrb > 0,]
  results_ckd_same_direction_hrb_sig_fdr <- results_ckd_same_direction[results_ckd_same_direction$p_fdr_ckd_hrb < 0.2,]
  
  results_ckd_estimate <- results_ckd[,c("species","estimate_ckd_jd","estimate_ckd_hrb")]
  colnames(results_ckd_estimate) <- c("species","JD","HRB")
  
  results_ckd_se <- results_ckd[,c("species","se_ckd_jd","se_ckd_hrb")]
  colnames(results_ckd_se) <- c("species","JD","HRB")
  
  results_ckd_p <- results_ckd[,c("species","p_ckd_jd","p_ckd_hrb")]
  colnames(results_ckd_p) <- c("species","JD","HRB")
  # 转换成long data
  results_ckd_estimate_long <- results_ckd_estimate %>%
    pivot_longer(
      cols = -species,                   
      names_to = "study",              
      values_to = "estimate"           
    )
  results_ckd_se_long <- results_ckd_se %>%
    pivot_longer(
      cols = -species,                   
      names_to = "study",              
      values_to = "se"           
    )
  results_ckd_p_long <- results_ckd_p %>%
    pivot_longer(
      cols = -species,                   
      names_to = "study",             
      values_to = "p"           
    )
  results_ckd_long <- left_join(results_ckd_estimate_long, results_ckd_se_long, by=c("species","study")) %>%
    left_join(results_ckd_p_long, by=c("species","study"))
}

# ASCVD结果处理 #
{
  results_acvd <- left_join(cvd_jd, cvd_gd, by = "species") %>%
    filter(p_cvd_jd < 0.05) %>%
    arrange(desc(estimate_cvd_jd))
    
  # 筛选JD和GD方向一致的结果
  results_acvd_same_direction <- results_acvd[results_acvd$estimate_cvd_jd * results_acvd$estimate_cvd_gd > 0,]
  results_acvd_same_direction_acvd_sig_fdr <- results_acvd_same_direction[results_acvd_same_direction$p_fdr_cvd_gd < 0.2,]
  
  results_acvd_estimate <- results_acvd[,c("species","estimate_cvd_jd","estimate_cvd_gd")]
  colnames(results_acvd_estimate) <- c("species","JD","GD")
  
  results_acvd_se <- results_acvd[,c("species","se_cvd_jd","se_cvd_gd")]
  colnames(results_acvd_se) <- c("species","JD","GD")
  
  results_acvd_p <- results_acvd[,c("species","p_cvd_jd","p_cvd_gd")]
  colnames(results_acvd_p) <- c("species","JD","GD")
  # 转换成long data
  results_acvd_estimate_long <- results_acvd_estimate %>%
    pivot_longer(
      cols = -species,                  
      names_to = "study",             
      values_to = "estimate"          
    )
  results_acvd_se_long <- results_acvd_se %>%
    pivot_longer(
      cols = -species,                  
      names_to = "study",             
      values_to = "se"           
    )
  results_acvd_p_long <- results_acvd_p %>%
    pivot_longer(
      cols = -species,                  
      names_to = "study",             
      values_to = "p"          
    )
  results_acvd_long <- left_join(results_acvd_estimate_long, results_acvd_se_long, by=c("species","study")) %>%
    left_join(results_acvd_p_long, by=c("species","study"))
}
#### (JD HRB ACVD) Replication结果处理 ----

#### 保留验证队列方向一致的结果 ----
# for fig 6a #
results_dm_same_direction_for_fig6a <- results_dm_same_direction
results_ckd_same_direction_for_fig6a <- results_ckd_same_direction
results_acvd_same_direction_for_fig6a <- results_acvd_same_direction

results_dm_same_direction_for_fig6a$species <- paste0(results_dm_same_direction_for_fig6a$species,"_log10")
results_ckd_same_direction_for_fig6a$species <- paste0(results_ckd_same_direction_for_fig6a$species,"_log10")
results_acvd_same_direction_for_fig6a$species <- paste0(results_acvd_same_direction_for_fig6a$species,"_log10")

openxlsx::write.xlsx(results_dm_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_dm.xlsx")
openxlsx::write.xlsx(results_ckd_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_ckd.xlsx")
openxlsx::write.xlsx(results_acvd_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_acvd.xlsx")

# 验证队列中方向一致且FDR-P<0.2的结果
results_dm_same_direction_hrb_sig <- results_dm_same_direction_hrb_sig_fdr
results_ckd_same_direction_hrb_sig <- results_ckd_same_direction_hrb_sig_fdr
results_acvd_same_direction_acvd_sig <- results_acvd_same_direction_acvd_sig_fdr

results_dm_same_direction_hrb_sig$species <- paste0(results_dm_same_direction_hrb_sig$species,"_log10")
results_ckd_same_direction_hrb_sig$species <- paste0(results_ckd_same_direction_hrb_sig$species,"_log10")
results_acvd_same_direction_acvd_sig$species <- paste0(results_acvd_same_direction_acvd_sig$species,"_log10")

openxlsx::write.xlsx(results_dm_same_direction_hrb_sig, "results/replication/replication_same_direction/replication_same_direction_sig_dm.xlsx")
openxlsx::write.xlsx(results_ckd_same_direction_hrb_sig, "results/replication/replication_same_direction/replication_same_direction_sig_ckd.xlsx")
openxlsx::write.xlsx(results_acvd_same_direction_acvd_sig, "results/replication/replication_same_direction/replication_same_direction_sig_acvd.xlsx")
#### 保留验证队列方向一致的结果 ####

#### 作图Forest (JD HRB) ----
# Logistic分析结果(MP4-DM CKD)作图 #
# 构建绘制森林图的函数
forest_function <- function(DAT, TITLE){
  
  f <- ggplot(data=DAT, aes(x=estimate, y=species)) +
    
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, 
                      group = study,              # 关键1：按study分组
                      color = errorbar_color),    # 关键2：误差线颜色映射到errorbar_color变量
                  width=0, cex=0.9, position = position_dodge(0.75)) +
    
    geom_point(aes(group = study,            # 关键1：按study分组
                   color = point_color),     # 关键2：点颜色映射到point_color变量
               size = 3.8,
                   position = position_dodge(0.75)) +
    
    geom_text(aes(label = text, group = study),
              hjust = -0.5,
              size = 6.6,
              position = position_dodge(width = 0.9)) +
    
    # 将颜色、图例标题、整合到 scale_colour_manual
    scale_colour_manual(
      name = "Cohort",  # legend title
      
      # 设置颜色
      values = c("jd" = "#8e44ad",
                 "rep" = "#2980b9"),
      # 直接在scale_color_manual中设置标签
      labels = c("jd" = "JD_sub", # "Discovery: JD_sub"
                 "rep" = "Validation"), # "Validation: (a) HRB_sub1 (diabetes), (b) HRB_sub2 (CKD), (c) ACVD_2017 (CVD)"
      # 通过breaks参数明确指定图例顺序
      breaks = c("jd", "rep")
    ) +
    
    scale_x_continuous(
      # labels = function(x) round(exp(x), 1),  # 显示为 exp(x)，保留1位小数 (不保留小数点后最后一位的0)
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.1)
    ) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(title = TITLE,
         x="OR/HR (95%CI)") +
    
    theme_classic() +
    
    theme(
      # plot.margin = margin(1, 1, 1, 1, "cm"),  # 四周各 1cm 边距
      plot.title = element_text(size = 25),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 25, margin = margin(t = 5, r = 0, b = 0, l = 0)),
      
      panel.spacing.x = unit(8, "mm"),  # 横向间距
      panel.spacing.y = unit(8, "mm"), # 纵向间距
      
      axis.ticks = element_line(size = 0.8, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.15, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      
      axis.text.x = element_text(size = 25, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 23, colour = "black"),
      
      # legend.position = "bottom",
      # legend.title = element_text(size = 15, colour = "black"),
      legend.key.height = unit(10, "mm"),
      legend.key.width = unit(10, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 25, colour = "black")   # 调整legend文本大小
    )
  
  return(f)
}

# DM数据整理for plot #
{
  dat_dm_for_forest <- results_dm_long[results_dm_long$species %in% results_dm_same_direction_hrb_sig_fdr$species,]
  dat_dm_for_forest$hror <- exp(dat_dm_for_forest$estimate)
  dat_dm_for_forest$lci <- exp(dat_dm_for_forest$estimate - 1.96*dat_dm_for_forest$se)
  dat_dm_for_forest$uci <- exp(dat_dm_for_forest$estimate + 1.96*dat_dm_for_forest$se)
  
  ### 如果p<0.05需要确保UCI<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_dm_for_forest[dat_dm_for_forest$p < 0.05 & (dat_dm_for_forest$uci >= 0.9995 & dat_dm_for_forest$uci < 1),]
    dat_dm_for_forest$uci <- ifelse(dat_dm_for_forest$p < 0.05 & (dat_dm_for_forest$uci >= 0.9995 & dat_dm_for_forest$uci < 1), 0.999, dat_dm_for_forest$uci)
  }
  
  dat_dm_for_forest$text <- paste0(sprintf("%.3f", dat_dm_for_forest$hror)," (",sprintf("%.3f", dat_dm_for_forest$lci),", ",sprintf("%.3f", dat_dm_for_forest$uci),")") 
  
  dat_dm_for_forest$point_color <- ifelse(dat_dm_for_forest$study == "JD", "jd", "rep")
  dat_dm_for_forest$errorbar_color <- ifelse(dat_dm_for_forest$study == "JD", "jd", "rep")
}
dat_dm_for_forest$species <- factor(dat_dm_for_forest$species, 
                                    levels = rev(unique(dat_dm_for_forest$species)), 
                                    labels = rev(trans_mp4_s_names(unique(dat_dm_for_forest$species))))
f1 <- forest_function(dat_dm_for_forest, "")

# CKD数据整理for plot #
{
  dat_ckd_for_forest <- results_ckd_long[results_ckd_long$species %in% results_ckd_same_direction_hrb_sig_fdr$species,]
  dat_ckd_for_forest$hror <- exp(dat_ckd_for_forest$estimate)
  dat_ckd_for_forest$lci <- exp(dat_ckd_for_forest$estimate - 1.96*dat_ckd_for_forest$se)
  dat_ckd_for_forest$uci <- exp(dat_ckd_for_forest$estimate + 1.96*dat_ckd_for_forest$se)
  
  ### 如果p<0.05需要确保UCI<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_ckd_for_forest[dat_ckd_for_forest$p < 0.05 & (dat_ckd_for_forest$uci >= 0.9995 & dat_ckd_for_forest$uci < 1),]
    dat_ckd_for_forest$uci <- ifelse(dat_ckd_for_forest$p < 0.05 & (dat_ckd_for_forest$uci >= 0.9995 & dat_ckd_for_forest$uci < 1), 0.999, dat_ckd_for_forest$uci)
  }
  
  dat_ckd_for_forest$text <- paste0(sprintf("%.3f", dat_ckd_for_forest$hror)," (",sprintf("%.3f", dat_ckd_for_forest$lci),", ",sprintf("%.3f", dat_ckd_for_forest$uci),")") 
  
  dat_ckd_for_forest$point_color <- ifelse(dat_ckd_for_forest$study == "JD", "jd", "rep")
  dat_ckd_for_forest$errorbar_color <- ifelse(dat_ckd_for_forest$study == "JD", "jd", "rep")
}
dat_ckd_for_forest$species <- factor(dat_ckd_for_forest$species, 
                                     levels = rev(unique(dat_ckd_for_forest$species)),
                                     labels = rev(trans_mp4_s_names(unique(dat_ckd_for_forest$species))))
f2 <- forest_function(dat_ckd_for_forest, "")

# CVD数据整理for plot #
{
  dat_acvd_for_forest <- results_acvd_long[results_acvd_long$species %in% results_acvd_same_direction_acvd_sig_fdr$species,]
  dat_acvd_for_forest$hror <- exp(dat_acvd_for_forest$estimate)
  dat_acvd_for_forest$lci <- exp(dat_acvd_for_forest$estimate - 1.96*dat_acvd_for_forest$se)
  dat_acvd_for_forest$uci <- exp(dat_acvd_for_forest$estimate + 1.96*dat_acvd_for_forest$se)
  dat_acvd_for_forest$text <- paste0(sprintf("%.3f", dat_acvd_for_forest$hror)," (",sprintf("%.3f", dat_acvd_for_forest$lci),", ",sprintf("%.3f", dat_acvd_for_forest$uci),")") 
  
  dat_acvd_for_forest$point_color <- ifelse(dat_acvd_for_forest$study == "JD", "jd", "rep")
  dat_acvd_for_forest$errorbar_color <- ifelse(dat_acvd_for_forest$study == "JD", "jd", "rep")
}
dat_acvd_for_forest$species <- factor(dat_acvd_for_forest$species, 
                                    levels = rev(unique(dat_acvd_for_forest$species)), 
                                    labels = rev(trans_mp4_s_names(unique(dat_acvd_for_forest$species))))
f3 <- forest_function(dat_acvd_for_forest, "")


f1_2 <- cowplot::plot_grid(f1, NULL,
                           nrow = 2,
                           rel_heights = c(1, 0.25),
                           align = "v")
f3_2 <- cowplot::plot_grid(f3, NULL,
                           nrow = 2,
                           rel_heights = c(1, 2.6),
                           align = "v")

f <- cowplot::plot_grid(f1_2, f2, f3_2,
                        ncol = 3,
                        rel_widths = c(1, 1, 0.8))

ggsave(paste0("figures/main_figures/(fig5bcd)_forest_plot_replication_dm-ckd-cvd.pdf"),
       f, width = 34, height = 20, limitsize = FALSE)
#### 作图Forest (JD HRB) ####


#### 添加菌的门信息 (作图数据处理) ----
{
  dat_dm_for_forest <- dat_dm_for_forest[dat_dm_for_forest$study == "JD",]
  dat_ckd_for_forest <- dat_ckd_for_forest[dat_ckd_for_forest$study == "JD",]
  dat_acvd_for_forest <- dat_acvd_for_forest[dat_acvd_for_forest$study == "JD",]
  
  #### 读取clade name数据 ----
  ## 读取clade name ##
  clade_name_mp4_s <- read.table("raw_data/JD.mp4.n4491_clade_name.txt", header = TRUE)
  clade_name_mp4_s$species <- trans_mp4_s_names(clade_name_mp4_s$species)
  clade_name_mp4_s <- clade_name_mp4_s[,c("phylum","species")]
  table(clade_name_mp4_s$species)
  table(clade_name_mp4_s$phylum)
  #### 读取clade name数据 ####
  
  dat_dm_for_phy_heatmap <- data.frame(species = unique(dat_dm_for_forest$species)) 
  dat_dm_for_phy_heatmap <- left_join(dat_dm_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_dm_for_phy_heatmap$y_axis <- 1
  
  dat_ckd_for_phy_heatmap <- data.frame(species = unique(dat_ckd_for_forest$species)) 
  dat_ckd_for_phy_heatmap <- left_join(dat_ckd_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_ckd_for_phy_heatmap$y_axis <- 1
  
  dat_acvd_for_phy_heatmap <- data.frame(species = unique(dat_acvd_for_forest$species)) 
  dat_acvd_for_phy_heatmap <- left_join(dat_acvd_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_acvd_for_phy_heatmap$y_axis <- 1
  
  
  dat_dm_for_phy_heatmap$species <- factor(dat_dm_for_phy_heatmap$species,
                                           levels = c(dat_dm_for_phy_heatmap$species))
  dat_ckd_for_phy_heatmap$species <- factor(dat_ckd_for_phy_heatmap$species,
                                            levels = c(dat_ckd_for_phy_heatmap$species))
  dat_acvd_for_phy_heatmap$species <- factor(dat_acvd_for_phy_heatmap$species,
                                             levels = c(dat_acvd_for_phy_heatmap$species))
  
  phylum_colors <- c(
    "p__Actinobacteria"  = "#cf6a87",
    "p__Bacteroidetes"   = "#f19066",
    "p__Candidatus_Saccharibacteria" = "#f5cd79",
    "p__Firmicutes"      = "#2ecc71",
    "p__Proteobacteria"  = "#34ace0"
  )
}
# 添加菌的门信息 (作图)
f_heatmap_phylum_dm <- ggplot(dat_dm_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_dm, filename=paste0("figures/main_figures/(fig5b)_forest_mp4_out_phylum_dm.pdf"), width = 14, height = 3.8, limitsize = FALSE)


f_heatmap_phylum_ckd1 <- ggplot(dat_ckd_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_ckd1, filename=paste0("figures/main_figures/(fig5c)_forest_mp4_out_phylum_ckd.pdf"), width = 14, height = 3.8, limitsize = FALSE)


f_heatmap_phylum_cvd <- ggplot(dat_acvd_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_cvd, filename=paste0("figures/main_figures/(fig5d)_forest_mp4_out_phylum_cvd.pdf"), width = 14, height = 3.8, limitsize = FALSE)
#### 添加菌的门信息 (作图数据处理) ####
