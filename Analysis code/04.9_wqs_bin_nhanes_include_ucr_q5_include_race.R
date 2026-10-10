library(data.table)
library(dplyr)
library(gWQS)
set.seed(20241108)  # 设置随机数种子，确保结果可重复

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
  nhanes_dat1 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]), ]
  table(nhanes_dat1$release)
}
# PAEs_6; 2003-2018 (8c) #
{
  cols_edc <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat5 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]), ]
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


# nhanes_dat1 (13EDC) #
#### WQS ----
# 构建分析分类表型的WQS模型(EDC 2分类) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
wqs_bin_q2 <- function(DAT,OUT,EDC){
  # out <- out_traits_bin[1]
  # edc <- "PFAS"
  if(EDC == "EDC_13"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){
    edc_group <- edc_traits5_q2_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age + race_eth + sex_rev + smk1 + drk1 + high_edu + mvpa_g + kcal_q4 + ucr_q5 + release_cycle"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    # weights = WQS$WTMEC2YR_COMB,
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       # weights = WQS$WTMEC2YR_COMB,
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}
# 构建分析分类表型的WQS模型(EDC 4分类) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
wqs_bin_q4 <- function(DAT,OUT,EDC){
  # out <- out_traits_bin[1]
  # edc <- "PFAS"
  if(EDC == "EDC_13"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){
    edc_group <- edc_traits5_q2_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age + race_eth + sex_rev + smk1 + drk1 + high_edu + mvpa_g + kcal_q4 + ucr_q5 + release_cycle"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    # weights = WQS$WTMEC2YR_COMB,
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       # weights = WQS$WTMEC2YR_COMB,
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })

  return(temp_results)
  
}


# 结局变量 (分类)
out_traits_bin <- c("dm","ckd","cvd")
# 校正变量
cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")


### 分析 ###
options(future.globals.maxSize = 1024 * 1024 * 1024 * 4)  # 提高 future.globals.maxSize 的阈值（设为4GB）
for (i in c("nhanes_dat1_1114")) {
  # i <- "nhanes_dat1_1114"
  
  phy_edc_dat <- nhanes_dat1[nhanes_dat1$release %in% c("2011-2012","2013-2014"),]
  phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/2 # 两轮数据除以2
  cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  #### 分类结局 ----
  results_wqs_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # j<-1
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcomes"))
    
    print(paste0(Sys.time(),"; EDC_13"))
    result_1 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_13") # 检出率>50%
    print(paste0(Sys.time(),"; PFAS"))
    result_2 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>50%
    print(paste0(Sys.time(),"; PAE_6"))
    result_3 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE_6") # 检出率>50%
    print(paste0(Sys.time(),"; BP_1"))
    result_4 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"BP_1") # 检出率>50%
    print(paste0(Sys.time(),"; TC_1"))
    result_5 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"TC_1") # 检出率>50%
    
    # print(paste0(Sys.time(),"; EDC_9"))
    # result_6 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%
    # print(paste0(Sys.time(),"; PFAS"))
    # result_7 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>75%
    # print(paste0(Sys.time(),"; PAE_4"))
    # result_8 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PAE_4") # 检出率>75%
    
    results_wqs_bin <- do.call(c, list(results_wqs_bin, result_1, result_2, result_3, result_4, result_5))
  }
  
  ## 把结果汇总成表格 ##
  wqs_results_bin_table <- data.frame()
  for (k in 1:length(results_wqs_bin)) {
    
    temp_sum <- summary(results_wqs_bin[[k]])
    temp_sum2 <- temp_sum[["coefficients"]]
    
    strings <- strsplit(names(results_wqs_bin)[k], "\\|")[[1]]
    
    r <- data.frame(f1 = as.character(results_wqs_bin[[k]][["formula"]][3]),
                    cov = row.names(temp_sum2),
                    exp = strings[1],
                    type = strings[2],
                    out = strings[3],
                    method = strings[5],
                    direction = strings[4],
                    sample = sample_name,
                    temp_sum2,
                    t(results_wqs_bin[[k]][["final_weights"]][["mean_weight"]]))
    
    weight_names <- as.character(results_wqs_bin[[k]][["final_weights"]][["mix_name"]])
    colnames(r) <- c("formula","cov","exp","type","out","method","direction","sample","estimate","se","z","p", gsub("_log10", "", weight_names))
    
    wqs_results_bin_table <- bind_rows(wqs_results_bin_table, r)
  }
  #### 分类结局 ####
  
  openxlsx::write.xlsx(wqs_results_bin_table, paste0("results/correlations/wqs/wqs_results_bin_(",sample_name,")_NOsamplingweights.xlsx"))
}

#### WQS ####


# nhanes_dat5 (PAE_6) #
#### WQS ----
# 构建分析分类表型的WQS模型(EDC 2分类) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
wqs_bin_q2 <- function(DAT,OUT,EDC){
  # out <- out_traits_bin[1]
  # edc <- "PFAS"
  if(EDC == "EDC_13"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){
    edc_group <- edc_traits5_q2_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age + race_eth + sex_rev + smk1 + drk1 + high_edu + mvpa_g + kcal_q4 + ucr_q5 + release_cycle"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    # weights = WQS$WTMEC2YR_COMB,
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       # weights = WQS$WTMEC2YR_COMB,
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}
# 构建分析分类表型的WQS模型(EDC 4分类) (校正："age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
wqs_bin_q4 <- function(DAT,OUT,EDC){
  # out <- out_traits_bin[1]
  # edc <- "PFAS"
  if(EDC == "EDC_13"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE_6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP_1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC_1"){
    edc_group <- edc_traits5_q2_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits,"WTMEC2YR_COMB","SDMVPSU","SDMVSTRA")
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age + race_eth + sex_rev + smk1 + drk1 + high_edu + mvpa_g + kcal_q4 + ucr_q5 + release_cycle"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    # weights = WQS$WTMEC2YR_COMB,
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       # weights = WQS$WTMEC2YR_COMB,
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}


# 结局变量 (分类)
out_traits_bin <- c("dm","ckd","cvd")
# 校正变量
cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")


### 分析 ###
options(future.globals.maxSize = 1024 * 1024 * 1024 * 4)  # 提高 future.globals.maxSize 的阈值（设为4GB）
for (i in c("nhanes_dat5_0318")) {
  # i <- "nhanes_dat5_0318"
  
  phy_edc_dat <- nhanes_dat5
  phy_edc_dat$WTMEC2YR_COMB <- phy_edc_dat$WTMEC2YR/8 # 八轮数据除以8
  cov_traits <- c("age","race_eth","sex_rev","smk1","drk1","high_edu","mvpa_g","kcal_q4","ucr_q5","release_cycle")
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_rev <- factor(phy_edc_dat$sex_rev) # 0/1（女/男）
  phy_edc_dat$race_eth <- factor(phy_edc_dat$race_eth)
  phy_edc_dat$smk1 <- factor(phy_edc_dat$smk1)
  phy_edc_dat$drk1 <- factor(phy_edc_dat$drk1)
  phy_edc_dat$high_edu <- factor(phy_edc_dat$high_edu)
  phy_edc_dat$mvpa_g <- factor(phy_edc_dat$mvpa_g)
  phy_edc_dat$kcal_q4 <- factor(phy_edc_dat$kcal_q4)
  phy_edc_dat$ucr_q5 <- factor(phy_edc_dat$ucr_q5)
  phy_edc_dat$release_cycle <- factor(phy_edc_dat$release_cycle)
  # 分类协变量转换为因子
  
  #### 分类结局 ----
  results_wqs_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # j<-1
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcomes"))
    
    # print(paste0(Sys.time(),"; EDC_13"))
    # result_1 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_13") # 检出率>50%
    # print(paste0(Sys.time(),"; PFAS"))
    # result_2 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>50%
    print(paste0(Sys.time(),"; PAE_6"))
    result_3 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE_6") # 检出率>50%
    # print(paste0(Sys.time(),"; BP_1"))
    # result_4 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"BP_1") # 检出率>50%
    # print(paste0(Sys.time(),"; TC_1"))
    # result_5 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"TC_1") # 检出率>50%
    
    # print(paste0(Sys.time(),"; EDC_9"))
    # result_6 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%
    # print(paste0(Sys.time(),"; PFAS"))
    # result_7 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>75%
    # print(paste0(Sys.time(),"; PAE_4"))
    # result_8 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PAE_4") # 检出率>75%
    
    results_wqs_bin <- do.call(c, list(results_wqs_bin, result_3))
  }
  
  ## 把结果汇总成表格 ##
  wqs_results_bin_table <- data.frame()
  for (k in 1:length(results_wqs_bin)) {
    
    temp_sum <- summary(results_wqs_bin[[k]])
    temp_sum2 <- temp_sum[["coefficients"]]
    
    strings <- strsplit(names(results_wqs_bin)[k], "\\|")[[1]]
    
    r <- data.frame(f1 = as.character(results_wqs_bin[[k]][["formula"]][3]),
                    cov = row.names(temp_sum2),
                    exp = strings[1],
                    type = strings[2],
                    out = strings[3],
                    method = strings[5],
                    direction = strings[4],
                    sample = sample_name,
                    temp_sum2,
                    t(results_wqs_bin[[k]][["final_weights"]][["mean_weight"]]))
    
    weight_names <- as.character(results_wqs_bin[[k]][["final_weights"]][["mix_name"]])
    colnames(r) <- c("formula","cov","exp","type","out","method","direction","sample","estimate","se","z","p", gsub("_log10", "", weight_names))
    
    wqs_results_bin_table <- bind_rows(wqs_results_bin_table, r)
  }
  #### 分类结局 ####
  
  openxlsx::write.xlsx(wqs_results_bin_table, paste0("results/correlations/wqs/wqs_results_bin_(",sample_name,")_NOsamplingweights.xlsx"))
}
#### WQS ####
