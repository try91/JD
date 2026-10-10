library(data.table)
library(dplyr)
library(qgcomp)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")

# 亚组人群 #
{
  # 人群0_1（男性）
  phy_edc_temp0_1 <- phy_edc_dat %>%
    filter(sex_b_rev == 1)
  # 人群0_2（女性）
  phy_edc_temp0_2 <- phy_edc_dat %>%
    filter(sex_b_rev == 0)
  # # 有2014肠道菌群数据人群
  # phy_edc_temp0_3 <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  #   right_join(micro_dat, by = "ID")
  phy_edc_temp_list <- list()
  phy_edc_temp_list[["phy_edc_temp"]] <- phy_edc_dat
  phy_edc_temp_list[["phy_edc_temp0_1"]] <- phy_edc_temp0_1
  phy_edc_temp_list[["phy_edc_temp0_2"]] <- phy_edc_temp0_2
  # phy_edc_temp_list[["phy_edc_temp0_3"]] <- phy_edc_temp0_3
  }

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


#### 性别特异性qgcomp分析 (EDC14 Q2) (q=2, bayes=FALSE) ----
# 构建分析新发病的qgcomp模型 (q=2, 没有bayes这个参数) (校正："age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
qg_analysis_cox_q2 <- function(DAT,OUT,TIME,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE"){
    edc_group <- edc_traits3_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "PAE_4"){   # 检出率>75%的PAE4
    edc_group <- edc_traits3_q4_log10
  }else if(EDC == "BP"){
    edc_group <- edc_traits4_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC"){
    edc_group <- edc_traits5_log10
  }else if(EDC == "EDC_19"){
    edc_group <- edc_traits_log10
  }else if(EDC == "EDC_14"){  # 检出率>50%的EDC14
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_9"){   # 检出率>75%的EDC9
    edc_group <- edc_traits7_log10
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_pos2"){
    results_qg_temp_pos <- results_qg[results_qg$exp == "EDC_pos" & results_qg$out == OUT,c(9:22)]
    results_qg_temp_pos <- results_qg_temp_pos[,which(!is.na(results_qg_temp_pos[1,]))]
    results_qg_temp_neg <- results_qg[results_qg$exp == "EDC_neg" & results_qg$out == OUT,c(9:22)]
    results_qg_temp_neg <- results_qg_temp_neg[,which(!is.na(results_qg_temp_neg[1,]))]
    results_qg_temp_neg_short <- results_qg_temp_neg[,which(results_qg_temp_neg[1,]>0)]
    results_qg_temp_neg_short <- results_qg_temp_neg_short[,which(results_qg_temp_neg_short[1,] == max(abs(results_qg_temp_neg_short[1,])))]
    edc_group <- c(colnames(results_qg_temp_pos),colnames(results_qg_temp_neg_short))
  }else if(EDC == "EDC_neg2"){
    results_qg_temp_neg <- results_qg[results_qg$exp == "EDC_neg" & results_qg$out == OUT,c(9:22)]
    results_qg_temp_neg <- results_qg_temp_neg[,which(!is.na(results_qg_temp_neg[1,]))]
    results_qg_temp_neg_short <- results_qg_temp_neg[,which(results_qg_temp_neg[1,]>0)]
    results_qg_temp_neg_short <- results_qg_temp_neg_short[,which(results_qg_temp_neg_short[1,] == max(abs(results_qg_temp_neg_short[1,])))]
    if(length(results_qg_temp_neg_short) > 0){
      edc_group <- colnames(results_qg_temp_neg)[which(colnames(results_qg_temp_neg) != colnames(results_qg_temp_neg_short))]
    }else{
      edc_group <- colnames(results_qg_temp_neg)
    }
  }
  
  # 汇总变量名
  cols <- c(OUT,TIME,edc_group,cov_traits)
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  colnames(QGCOMP)[2] <- 't'
  QGCOMP$y <- as.numeric(as.character(QGCOMP$y)) # COX 分析结局是连续变量
  QGCOMP$high_fruveg <- factor(QGCOMP$high_fruveg)
  
  toxic_chems <- edc_group
  
  results = qgcomp.cox.noboot(Surv(t, y)~., expnms=toxic_chems, 
                              data = QGCOMP, 
                              q = 2)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|cox")]] <- results
  return(temp_results)
  
}
# 构建分析分类表型的qgcomp模型 (q=2, bayes=FALSE) (校正："age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
qg_analysis_bin_q2 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE"){
    edc_group <- edc_traits3_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "PAE_4"){   # 检出率>75%的PAE4
    edc_group <- edc_traits3_q4_log10
  }else if(EDC == "BP"){
    edc_group <- edc_traits4_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC"){
    edc_group <- edc_traits5_log10
  }else if(EDC == "EDC_19"){
    edc_group <- edc_traits_log10
  }else if(EDC == "EDC_14"){  # 检出率>50%的EDC14
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_9"){   # 检出率>75%的EDC9
    edc_group <- edc_traits7_log10
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$y <- factor(QGCOMP$y) # 分类变量结局转化为因子
  QGCOMP$high_fruveg <- factor(QGCOMP$high_fruveg)
  
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=binomial(), 
                          q = 2)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|bin")]] <- results
  return(temp_results)
  
}
# 构建分析连续表型的qgcomp模型 (q=2, bayes=FALSE) (校正："age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
qg_analysis_cont_q2 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE"){
    edc_group <- edc_traits3_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "PAE_4"){   # 检出率>75%的PAE4
    edc_group <- edc_traits3_q4_log10
  }else if(EDC == "BP"){
    edc_group <- edc_traits4_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC"){
    edc_group <- edc_traits5_log10
  }else if(EDC == "EDC_19"){
    edc_group <- edc_traits_log10
  }else if(EDC == "EDC_14"){  # 检出率>50%的EDC14
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_9"){   # 检出率>75%的EDC9
    edc_group <- edc_traits7_log10
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:22)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$high_fruveg <- factor(QGCOMP$high_fruveg)
  
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=gaussian(),
                          q = 2)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|cont")]] <- results
  return(temp_results)
  
}


# 新发结局 (生存)
out_traits_cox <- c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014")
time_traits_cox <- c("timecvd_1021","timecvd_1014","timeckd_1014","timedm_1014")

# 结局变量 (分类)
out_traits_bin <- c("cvd_b","ckd_b","dm_b","as_imt_b","hpt_b","nafld_b","ob_b","abob_b","dyslip_b","hua_b","ir_b","mets_b")

# 指标变量 (连续)
out_traits_cont <- c("bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
                     "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
                     "alt_b","ast_b","ggt_b","bia_b",
                     "egfr_b","scr_b","ua_b",
                     "glu0_b","glu120_b","vhba1c_b",
                     "ins0_b","ins120_b","homair_b","homab_b",
                     "sbp_b","dbp_b","pr_b",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")

# 校正变量
cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")

### 原始分析 (校正："age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
qgcomp_results_list <- list()
for (i in c("phy_edc_temp0_1","phy_edc_temp0_2")) {
  # i <- "phy_edc_temp0_1"
  phy_edc_dat <- phy_edc_temp_list[[i]] # 提取subgroup
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
  phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
  # 分类协变量转换为因子
  
  #### survival outcome ----
  results_qg_cox <- list()
  for (j in 1:length(out_traits_cox)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_cox),", survival outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
    
    # result_1 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"EDC_19")
    result_2 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"EDC_14") # 检出率>50%的EDC14
    result_3 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"EDC_9") # 检出率>75%的EDC9
    result_4 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"PFAS")
    result_5 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"PAE_6") # 检出率>50%的PAE6
    result_6 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"BP_1") # 检出率>50%的BP1
    result_7 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"TC")
    
    results_qg_cox <- do.call(c, list(results_qg_cox, result_2, result_3, result_4, result_5, result_6, result_7))
  }
  #### survival outcome ####
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
    
    # result_1 <- qg_analysis_bin2_q2(phy_edc_dat,out_traits_bin[j],"EDC_19")
    result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_14") # 检出率>50%的EDC14
    result_3 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%的EDC9
    result_4 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS")
    result_5 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE_6") # 检出率>50%的PAE6
    result_6 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"BP_1") # 检出率>50%的BP1
    result_7 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"TC")
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_2, result_3, result_4, result_5, result_6, result_7))
  }
  #### binary outcome ####
  
  #### continuous outcome ----
  results_qg_cont <- list()
  for (j in 1:length(out_traits_cont)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_cont),", continuous outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
    
    # result_1 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_19")
    result_2 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_14") # 检出率>50%的EDC14
    result_3 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_9") # 检出率>75%的EDC9
    result_4 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"PFAS")
    result_5 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"PAE_6") # 检出率>50%的PAE6
    result_6 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"BP_1") # 检出率>50%的BP1
    result_7 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"TC")
    
    results_qg_cont <- do.call(c, list(results_qg_cont, result_2, result_3, result_4, result_5, result_6, result_7))
  }
  #### continuous outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
  for (k in 1:length(results_qg_cox)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_cox)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_cox)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_cox)[k])),
                    estimate = results_qg_cox[[k]][["psi"]][["psi1"]],
                    se = results_qg_cox[[k]][["psi"]][["psi1"]]/results_qg_cox[[k]][["zstat"]],
                    p = results_qg_cox[[k]][["pval"]],
                    sample = sample_name,
                    n = length(results_qg_cox[[k]][["fit"]][["y"]]),
                    t(results_qg_cox[[k]][["pos.weights"]]),
                    t(results_qg_cox[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  for (k in 1:length(results_qg_bin)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_bin)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_bin)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_bin)[k])),
                    estimate = results_qg_bin[[k]][["psi"]][["psi1"]],
                    se = results_qg_bin[[k]][["psi"]][["psi1"]]/results_qg_bin[[k]][["zstat"]][2],
                    p = results_qg_bin[[k]][["pval"]][2],
                    sample = sample_name,
                    n = length(results_qg_bin[[k]][["fit"]][["y"]]),
                    t(results_qg_bin[[k]][["pos.weights"]]),
                    t(results_qg_bin[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  for (k in 1:length(results_qg_cont)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_cont)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_cont)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_cont)[k])),
                    estimate = results_qg_cont[[k]][["psi"]][["psi1"]],
                    se = results_qg_cont[[k]][["psi"]][["psi1"]]/results_qg_cont[[k]][["tstat"]][2],
                    p = results_qg_cont[[k]][["pval"]][2],
                    sample = sample_name,
                    n = length(results_qg_cont[[k]][["fit"]][["y"]]),
                    t(results_qg_cont[[k]][["pos.weights"]]),
                    t(results_qg_cont[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",sample_name,").xlsx"))
}
#### 性别特异性qgcomp分析 (EDC14 Q2) (q=2, bayes=FALSE) ####

#### 性别特异性qgcomp分析EDC结果，分别筛选正负EDC (EDC14 Q2) (q=2, bayes=FALSE) ----
### 基于原始分析，提取同方向EDC的进一步分析 (校正："age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
for (i in c("phy_edc_temp0_1","phy_edc_temp0_2")) {
  # i <- "phy_edc_temp0_1"
  phy_edc_dat <- phy_edc_temp_list[[i]] # 提取subgroup
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
  phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
  # 分类协变量转换为因子
  
  # 读取上一步qgcomp分析结果
  results_qg <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",sample_name,").xlsx"))
  results_qg <- results_qg[results_qg$exp == "EDC_14",]
  
  #### survival outcome ----
  results_qg_cox <- list()
  for (j in 1:length(out_traits_cox)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_cox),", survival outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") # "high_fruveg"包括unknown组, "high_fruveg_f"不包括unknown组
    
    result_1 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"EDC_pos")
    result_2 <- qg_analysis_cox_q2(phy_edc_dat,out_traits_cox[j],time_traits_cox[j],"EDC_neg")
    
    results_qg_cox <- do.call(c, list(results_qg_cox, result_1, result_2))
  }
  #### survival outcome ####
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") # "high_fruveg"包括unknown组, "high_fruveg_f"不包括unknown组
    
    result_1 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_pos")
    result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_neg")
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_1, result_2))
  }
  #### binary outcome ####
  
  #### continuous outcome ----
  results_qg_cont <- list()
  for (j in 1:length(out_traits_cont)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_cont),", continuous outcome"))
    
    cov_traits <- c("age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") # "high_fruveg"包括unknown组, "high_fruveg_f"不包括unknown组
    
    result_1 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_pos")
    result_2 <- qg_analysis_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_neg")
    
    results_qg_cont <- do.call(c, list(results_qg_cont, result_1, result_2))
  }
  #### continuous outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
  for (k in 1:length(results_qg_cox)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_cox)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_cox)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_cox)[k])),
                    estimate = results_qg_cox[[k]][["psi"]][["psi1"]],
                    se = results_qg_cox[[k]][["psi"]][["psi1"]]/results_qg_cox[[k]][["zstat"]],
                    p = results_qg_cox[[k]][["pval"]],
                    sample = sample_name,
                    n = length(results_qg_cox[[k]][["fit"]][["y"]]),
                    t(results_qg_cox[[k]][["pos.weights"]]),
                    t(results_qg_cox[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  for (k in 1:length(results_qg_bin)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_bin)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_bin)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_bin)[k])),
                    estimate = results_qg_bin[[k]][["psi"]][["psi1"]],
                    se = results_qg_bin[[k]][["psi"]][["psi1"]]/results_qg_bin[[k]][["zstat"]][2],
                    p = results_qg_bin[[k]][["pval"]][2],
                    sample = sample_name,
                    n = length(results_qg_bin[[k]][["fit"]][["y"]]),
                    t(results_qg_bin[[k]][["pos.weights"]]),
                    t(results_qg_bin[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  for (k in 1:length(results_qg_cont)) {
    r <- data.frame(exp = sub("\\|.*$", "", names(results_qg_cont)[k]),
                    out = sub("^.*\\|(.*?)\\|.*$", "\\1", names(results_qg_cont)[k]),
                    method = paste0("qgcomp ",sub(".*\\|(.*)$", "\\1", names(results_qg_cont)[k])),
                    estimate = results_qg_cont[[k]][["psi"]][["psi1"]],
                    se = results_qg_cont[[k]][["psi"]][["psi1"]]/results_qg_cont[[k]][["tstat"]][2],
                    p = results_qg_cont[[k]][["pval"]][2],
                    sample = sample_name,
                    n = length(results_qg_cont[[k]][["fit"]][["y"]]),
                    t(results_qg_cont[[k]][["pos.weights"]]),
                    t(results_qg_cont[[k]][["neg.weights"]]*-1))
    sum_results_qg_all <- bind_rows(sum_results_qg_all, r)
    
  }
  
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(",sample_name,").xlsx"))
}
#### 性别特异性qgcomp分析EDC结果，分别筛选正负EDC (EDC14 Q2) (q=2, bayes=FALSE) ####
