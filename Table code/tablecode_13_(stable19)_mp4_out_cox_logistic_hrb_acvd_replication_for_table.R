library(data.table)
library(dplyr)
library(tidyr)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

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

# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)
mp4_g_zero <- paste0(mp4_g_names,"_zero") # 菌群MP4填补0值丰度 (属层面)

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
                     
                     "sleept_f","sittimet_f","sittimet_b","sum_met_f","sum_met_b", "sum_met_work_f","sum_met_work_b", "sum_met_all_f","sum_met_all_b",
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

#### 读取 (JD HRB) Replication结果 ----
## JD结果 ##
{
  results_jd <- 
    results_jd <- 
    
    results_jd <- read.csv("results/replication/mp4_s_for_replication/disease_interation_related_mp4_species_20260728.csv")
  
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
  dm_hrb <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "T2D")
  dm_hrb <- dm_hrb[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(dm_hrb) <- c("species","estimate_dm_hrb","se_dm_hrb","p_dm_hrb","p_fdr_dm_hrb")
  dm_hrb <- dm_hrb %>%
    mutate(p_adj_bh_dm_hrb = p.adjust(p_dm_hrb, method = "BH"))
  
  ckd_hrb <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "CKD")
  ckd_hrb <- ckd_hrb[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(ckd_hrb) <- c("species","estimate_ckd_hrb","se_ckd_hrb","p_ckd_hrb","p_fdr_ckd_hrb")
  ckd_hrb <- ckd_hrb %>%
    mutate(p_adj_bh_ckd_hrb = p.adjust(p_ckd_hrb, method = "BH"))
}

## JieZ ACVD结果 ##
{
  cvd_gd <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "CVD") # 之前由分组变量设置了“Control”和“ASCVD”默认ASCVD作为ref，导致所有结果反向，已更正
  cvd_gd <- cvd_gd[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(cvd_gd) <- c("species","estimate_cvd_gd","se_cvd_gd","p_cvd_gd","p_fdr_cvd_gd")
  cvd_gd <- cvd_gd %>%
    mutate(p_adj_bh_cvd_gd = p.adjust(p_cvd_gd, method = "BH"))
}
#### 读取 (JD HRB) Replication结果 ####

#### Replication结果处理 ----
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
#### Replication结果处理 ----

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
openxlsx::write.xlsx(dat_for_stable,"tables/mp4_out_cox_logistic_hrb_acvd_replication_20260802.xlsx")
