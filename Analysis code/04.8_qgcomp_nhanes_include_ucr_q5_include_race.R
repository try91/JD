library(data.table)
library(dplyr)
library(qgcomp) 

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

nhanes_dat <- read.csv("raw_data/nhanes_dat_2003-2018_for_analysis.csv")

#### Subgroup for analysis ----
# 13 EDCs (PFAS_5, PAEs_6, TC_1, BP_1); 2011-2016 (3c) #
{
  cols_edc <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                "TCS",
                "BPA")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat1 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]),]
  table(nhanes_dat1$release)
}
# PAEs_6; 2003-2018 (8c) #
{
  cols_edc <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat5 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]),]
  table(nhanes_dat5$release)
}
#### Subgroup for analysis ####

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
                "BPA","BPS","BPF",
                "TCC","TCS")
edc_traits2 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS")
edc_traits3 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP")
edc_traits3_q2 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP") # 检出率>50%的PAE6
edc_traits3_q4 <- c("MEHP","MECPP","MEHHP","MEP") # 检出率>75%的PAE4
edc_traits4 <- c("BPA","BPS","BPF")
edc_traits4_q2 <- c("BPA") # 检出率>50%的BP1
edc_traits5 <- c("TCC","TCS")
edc_traits5_q2 <- c("TCS")
edc_traits6 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                 "BPA",
                 "TCS") # 检出率>50%
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
edc_traits5_q2_log10 <- paste0(edc_traits5_q2,"_log10")
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


#### qgcomp分析 (EDC13 Q2) (q=2, bayes=FALSE) ----
# 构建分析分类表型的qgcomp模型 (q=2, bayes=FALSE) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
qg_analysis_bin_q2 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){    # 检出率>50%的TC1
    edc_group <- edc_traits5_q2_log10
  }else if(EDC == "EDC_13"){  # 检出率>50%的EDC13
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:21)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:21)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$y <- factor(QGCOMP$y) # 分类变量结局转化为因子
  # QGCOMP$SDMVPSU <- factor(QGCOMP$SDMVPSU) # 校正聚类变量转化为因子
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=binomial(), q = 2)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|bin")]] <- results
  return(temp_results)
  
}
# 构建分析分类表型的qgcomp模型 (q=4, bayes=FALSE) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
qg_analysis_bin_q4 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){    # 检出率>50%的TC1
    edc_group <- edc_traits5_q2_log10
  }else if(EDC == "EDC_13"){  # 检出率>50%的EDC13
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:21)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:21)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$y <- factor(QGCOMP$y) # 分类变量结局转化为因子
  # QGCOMP$SDMVPSU <- factor(QGCOMP$SDMVPSU) # 校正聚类变量转化为因子
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=binomial(), q = 4)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|bin")]] <- results
  return(temp_results)
  
}


# 结局变量 (分类)
out_traits_bin <- c("dm","ckd","cvd")
# 校正变量
cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")


### 原始分析 (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle") ###
qgcomp_results_list <- list()
for (i in c("nhanes_dat1_1114")) {
  # i <- "nhanes_dat1_1114"
  
  phy_edc_dat <- nhanes_dat1[nhanes_dat1$release %in% c("2011-2012","2013-2014"),]
  phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/2 # 两轮数据除以2
  cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    if(all(is.na(phy_edc_dat[[out_traits_bin[j]]]))){
      cat("第", out_traits_bin[j], "个元素全为空，跳过\n")
      next  # 进入下一次循环
    }
    
    # result_1 <- qg_analysis_bin2_q2(phy_edc_dat,out_traits_bin[j],"EDC_19")
    result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_13") # 检出率>50%的EDC13
    # result_3 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%的EDC9
    result_4 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>50%的PFAS5
    result_5 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE_6") # 检出率>50%的PAE6
    result_6 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"BP_1") # 检出率>50%的BP1
    result_7 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"TC_1") # 检出率>50%的TC1
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_2, result_4, result_5, result_6, result_7))
  }
  #### binary outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
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
  
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name,")_NOsamplingweights.xlsx"))
  
}
#### qgcomp分析 (EDC13 Q2) (q=2, bayes=FALSE) ####

#### 读取qgcomp分析EDC结果，分别筛选正负EDC (EDC13 Q2) (q=2, bayes=FALSE) ----
### 基于原始分析，提取同方向EDC的进一步分析 (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle") ###
for (i in c("nhanes_dat1_1114")) {
  
  # i <- "nhanes_dat1_1114"
  
    phy_edc_dat <- nhanes_dat1[nhanes_dat1$release %in% c("2011-2012","2013-2014"),]
    phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/2 # 两轮数据除以2
    cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  # 读取上一步qgcomp分析结果
  results_qg <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name,")_NOsamplingweights.xlsx"))
  results_qg <- results_qg[results_qg$exp == "EDC_13",]
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    if(all(is.na(phy_edc_dat[[out_traits_bin[j]]]))){
      cat("第", out_traits_bin[j], "个元素全为空，跳过\n")
      next  # 进入下一次循环
    }
    
    result_1 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_pos")
    result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_neg")
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_1, result_2))
  }
  #### binary outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
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
  
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(edc+-)_(",sample_name,")_NOsamplingweights.xlsx"))
}
#### 读取qgcomp分析EDC结果，分别筛选正负EDC (EDC13 Q2) (q=2, bayes=FALSE) ####


#### qgcomp分析 (PAE_6 Q2) (q=2, bayes=FALSE) ----
# 构建分析分类表型的qgcomp模型 (q=2, bayes=FALSE) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
qg_analysis_bin_q2 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){    # 检出率>50%的TC1
    edc_group <- edc_traits5_q2_log10
  }else if(EDC == "EDC_13"){  # 检出率>50%的EDC13
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_11"){  # 检出率>50%的EDC11
    edc_group <- paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                          "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP"),"_log10")
  }else if(EDC == "EDC_8"){  # 检出率>50%的EDC11
    edc_group <- paste0(c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                          "TCS","BPA"),"_log10")
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:14)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:14)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB")
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$y <- factor(QGCOMP$y) # 分类变量结局转化为因子
  # QGCOMP$SDMVPSU <- factor(QGCOMP$SDMVPSU) # 校正聚类变量转化为因子
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=binomial(), q = 2)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|bin")]] <- results
  return(temp_results)
  
}
# 构建分析分类表型的qgcomp模型 (q=4, bayes=FALSE) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
qg_analysis_bin_q4 <- function(DAT,OUT,EDC){
  if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){   # 检出率>50%的PAE6
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){    # 检出率>50%的BP1
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){    # 检出率>50%的TC1
    edc_group <- edc_traits5_q2_log10
  }else if(EDC == "EDC_13"){  # 检出率>50%的EDC13
    edc_group <- edc_traits6_log10
  }else if(EDC == "EDC_11"){  # 检出率>50%的EDC11
    edc_group <- paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                          "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP"),"_log10")
  }else if(EDC == "EDC_8"){  # 检出率>50%的EDC11
    edc_group <- paste0(c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                          "TCS","BPA"),"_log10")
  }else if(EDC == "EDC_pos"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:14)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]>0)]
    edc_group <- colnames(results_qg_temp)
  }else if(EDC == "EDC_neg"){
    results_qg_temp <- results_qg[results_qg$out == OUT,c(9:14)]
    results_qg_temp <- results_qg_temp[,which(results_qg_temp[1,]<0)]
    edc_group <- colnames(results_qg_temp)
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB")
  # 提取所需变量
  QGCOMP <- DAT[,cols]
  QGCOMP <- na.omit(QGCOMP)
  
  colnames(QGCOMP)[1] <- 'y'
  QGCOMP$y <- factor(QGCOMP$y) # 分类变量结局转化为因子
  # QGCOMP$SDMVPSU <- factor(QGCOMP$SDMVPSU) # 校正聚类变量转化为因子
  toxic_chems <- edc_group
  
  results = qgcomp.noboot(y~., expnms=toxic_chems,
                          data = QGCOMP, family=binomial(), q = 4)
  
  temp_results <- list()
  temp_results[[paste0(EDC,"|",OUT,"|bin")]] <- results
  return(temp_results)
  
}


# 结局变量 (分类)
out_traits_bin <- c("dm","ckd","cvd")
# 校正变量
cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")


### 原始分析 (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle") ###
qgcomp_results_list <- list()
for (i in c("nhanes_dat5_0318")) {
  # i <- "nhanes"
  
  phy_edc_dat <- nhanes_dat5
  phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/8 # 八轮数据除以8
  cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    # result_1 <- qg_analysis_bin2_q2(phy_edc_dat,out_traits_bin[j],"EDC_19")
    # result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_8") # 检出率>50%的EDC13
    # result_3 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%的EDC9
    # result_4 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>50%的PFAS5
    result_5 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE_6") # 检出率>50%的PAE6
    # result_6 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"BP_1") # 检出率>50%的BP1
    # result_7 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"TC_1") # 检出率>50%的TC1
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_5))
  }
  #### binary outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
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
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name,")_NOsamplingweights.xlsx"))
  
}
#### qgcomp分析 (PAE_6 Q2) (q=2, bayes=FALSE) ####

#### 读取qgcomp分析EDC结果，分别筛选正负EDC (PAE_6 Q2) (q=2, bayes=FALSE) ----
### 基于原始分析，提取同方向EDC的进一步分析 (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle") ###
for (i in c("nhanes_dat5_0318")) {
  # i <- "nhanes"
  
  phy_edc_dat <- nhanes_dat5
  phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/8 # 八轮数据除以8
  cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  # 读取上一步qgcomp分析结果
  results_qg <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name,")_NOsamplingweights.xlsx"))
  results_qg <- results_qg[results_qg$exp == "PAE_6",]
  
  #### binary outcome ----
  results_qg_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # 打印当前系统时间
    print(Sys.time())
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcome"))
    
    result_1 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_pos")
    result_2 <- qg_analysis_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_neg")
    
    results_qg_bin <- do.call(c, list(results_qg_bin, result_1, result_2))
  }
  #### binary outcome ####
  
  ### 汇总结果格 ###
  ## 总体人群 ##
  sum_results_qg_all <- data.frame()
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
  openxlsx::write.xlsx(sum_results_qg_all, paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(edc+-)_(",sample_name,")_NOsamplingweights.xlsx"))
}
#### 读取qgcomp分析EDC结果，分别筛选正负EDC (PAE_6 Q2) (q=2, bayes=FALSE) ####
