library(data.table)
library(dplyr)
library(gWQS)
set.seed(20241108)  # 设置随机数种子，确保结果可重复

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp"  # 提取subgroup的名称

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


#### WQS ----
# 构建分析分类表型的WQS模型(EDC 2分类) (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
wqs_bin_q2 <- function(DAT,OUT,EDC){
  # OUT <- out_traits_bin[1]
  # EDC <- "PFAS"
  if(EDC == "EDC_14"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC"){
    edc_group <- edc_traits5_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  WQS$high_fruveg <- factor(WQS$high_fruveg) # 必须在循环中去除missing变量后才能转换为因子，不然因子中元素数量为0会报错
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(edc,"|",out,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}
# 构建分析连续表型的WQS模型(EDC 2分类) (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
wqs_cont_q2 <- function(DAT,OUT,EDC){
  # OUT <- out_traits[12]
  # EDC <- "PFAS"
  if(EDC == "EDC_14"){
    edc_group <- edc_traits6_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE6"){
    edc_group <- edc_traits3_q2_log10
  }else if(EDC == "BP1"){
    edc_group <- edc_traits4_q2_log10
  }else if(EDC == "TC"){
    edc_group <- edc_traits5_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS$high_fruveg <- factor(WQS$high_fruveg) # 必须在循环中去除missing变量后才能转换为因子，不然因子中元素数量为0会报错
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "gaussian",
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|pos|cont")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 2, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "gaussian",
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q2|",OUT,"|neg|cont")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}
# 构建分析分类表型的WQS模型(EDC 4分类) (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
wqs_bin_q4 <- function(DAT,OUT,EDC){
  # OUT <- out_traits_bin[1]
  # EDC <- "PFAS"
  if(EDC == "EDC_9"){
    edc_group <- edc_traits7_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE4"){
    edc_group <- edc_traits3_q4_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS[[1]] <- factor(WQS[[1]])
  WQS$high_fruveg <- factor(WQS$high_fruveg) # 必须在循环中去除missing变量后才能转换为因子，不然因子中元素数量为0会报错
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "binomial",
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|pos|bin")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "binomial",
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|neg|bin")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })

  return(temp_results)
  
}
# 构建分析连续表型的WQS模型(EDC 4分类) (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
wqs_cont_q4 <- function(DAT,OUT,EDC){
  # OUT <- out_traits[12]
  # EDC <- "PFAS"
  if(EDC == "EDC_9"){
    edc_group <- edc_traits7_log10
  }else if(EDC == "PFAS"){
    edc_group <- edc_traits2_log10
  }else if(EDC == "PAE4"){
    edc_group <- edc_traits3_q4_log10
  }
  
  # 汇总变量名
  cols <- c(OUT,edc_group,cov_traits)
  # 提取所需变量
  WQS <- DAT[,cols]
  WQS <- na.omit(WQS)
  WQS$high_fruveg <- factor(WQS$high_fruveg) # 必须在循环中去除missing变量后才能转换为因子，不然因子中元素数量为0会报错
  
  # # 1000人小样本测试
  # WQS <- WQS[sample(nrow(WQS), 1000), ]
  # # 1000人小样本测试
  toxic_chems <- edc_group
  # 重置临时列表
  temp_results <- list()
  
  # 构建分析模型
  f1 <- formula(paste0(OUT," ~ wqs + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg"))
  
  # 正向权重
  tryCatch({
    results <- gwqs(f1, mix_name = toxic_chems,
                    data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                    b1_pos = TRUE, b_constr = TRUE, family = "gaussian",
                    seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|pos|cont")]] <- results
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|pos|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  # 负向权重
  tryCatch({
    results_ne <- gwqs(f1, mix_name = toxic_chems,
                       data = WQS, q = 4, validation = 0.6, b = 10, # Number of bootstrap samples used in parameter estimation
                       b1_pos = FALSE, b_constr = TRUE, family = "gaussian",
                       seed = 20241108, plots = TRUE, tables = TRUE)
    
    temp_results[[paste0(EDC,"|Q4|",OUT,"|neg|cont")]] <- results_ne
  }, error = function(e) {
    # temp_results[[paste0(EDC,"|",OUT,"|neg|error")]] <- paste("Encounter an error")
    # 打印实际的错误内容
    print(paste("Encounter an error:", e$message))
  })
  
  return(temp_results)
  
}



# 结局变量 (分类)
out_traits_bin <- c("dm_incident_1014","ckd_incident_1014","cvd_incident_1021","cvd_incident_1014",
                    "dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b")
# 指标变量 (连续)
out_traits_cont <- c("bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
                     "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
                     "alt_b","ast_b","ggt_b","bia_b",
                     "egfr_b","scr_b","ua_b",
                     "glu0_b","glu120_b","vhba1c_b",
                     "ins0_b","ins120_b","homair_b","homab_b",
                     "sbp_b","dbp_b","p_b",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")

# 校正变量
cov_traits <- c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")


### 分析 ###
options(future.globals.maxSize = 1024 * 1024 * 1024 * 4)  # 提高 future.globals.maxSize 的阈值（设为4GB）
for (i in c("phy_edc_temp")) {
  # i <- "phy_edc_temp"
  
  sample_name <- i # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
  phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
  phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
  # 分类协变量转换为因子
  
  
  #### 分类结局 ----
  results_wqs_bin <- list()
  for (j in 1:length(out_traits_bin)) {
    # j<-1
    # 打印目前进度
    print(paste0(sample_name," || ",j," out of ",length(out_traits_bin),", binary outcomes"))
    
    print(paste0(Sys.time(),"; EDC_14"))
    result_1 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"EDC_14") # 检出率>50%
    print(paste0(Sys.time(),"; PFAS"))
    result_2 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>50%
    print(paste0(Sys.time(),"; PAE6"))
    result_3 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"PAE6") # 检出率>50%
    print(paste0(Sys.time(),"; BP1"))
    result_4 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"BP1") # 检出率>50%
    print(paste0(Sys.time(),"; TC"))
    result_5 <- wqs_bin_q2(phy_edc_dat,out_traits_bin[j],"TC") # 检出率>50%
    
    # print(paste0(Sys.time(),"; EDC_9"))
    # result_6 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"EDC_9") # 检出率>75%
    # print(paste0(Sys.time(),"; PFAS"))
    # result_7 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PFAS") # 检出率>75%
    # print(paste0(Sys.time(),"; PAE4"))
    # result_8 <- wqs_bin_q4(phy_edc_dat,out_traits_bin[j],"PAE4") # 检出率>75%
    
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
    colnames(r) <- c("formula","cov","exp","type","out","method","direction","sample","estimate","se","z","p",weight_names)
    
    wqs_results_bin_table <- bind_rows(wqs_results_bin_table, r)
  }
  #### 分类结局 ####
  
  
  # #### 连续指标 ----
  # results_wqs_cont <- list()
  # for (j in 1:length(out_traits_cont)) {
  #   # j<-1
  #   # 打印目前进度
  #   print(paste0(sample_name," || ",j," out of ",length(out_traits_cont),", continue outcomes"))
  #   
  #   print(paste0(Sys.time(),"; EDC_14"))
  #   result_1 <- wqs_cont_q2(phy_edc_dat,out_traits_cont[j],"EDC_14") # 检出率>50%
  #   print(paste0(Sys.time(),"; PFAS"))
  #   result_2 <- wqs_cont_q2(phy_edc_dat,out_traits_cont[j],"PFAS") # 检出率>50%
  #   print(paste0(Sys.time(),"; PAE6"))
  #   result_3 <- wqs_cont_q2(phy_edc_dat,out_traits_cont[j],"PAE6") # 检出率>50%
  #   print(paste0(Sys.time(),"; BP1"))
  #   result_4 <- wqs_cont_q2(phy_edc_dat,out_traits_cont[j],"BP1") # 检出率>50%
  #   print(paste0(Sys.time(),"; TC"))
  #   result_5 <- wqs_cont_q2(phy_edc_dat,out_traits_cont[j],"TC") # 检出率>50%
  #   
  #   # print(paste0(Sys.time(),"; EDC_9"))
  #   # result_6 <- wqs_cont_q4(phy_edc_dat,out_traits_cont[j],"EDC_9") # 检出率>75%
  #   # print(paste0(Sys.time(),"; PFAS"))
  #   # result_7 <- wqs_cont_q4(phy_edc_dat,out_traits_cont[j],"PFAS") # 检出率>75%
  #   # print(paste0(Sys.time(),"; PAE4"))
  #   # result_8 <- wqs_cont_q4(phy_edc_dat,out_traits_cont[j],"PAE4") # 检出率>75%
  #   
  #   results_wqs_cont <- do.call(c, list(results_wqs_cont, result_1, result_2, result_3, result_4, result_5, result_6, result_7, result_8))
  # }
  # 
  # ## 把结果汇总成表格 ##
  # wqs_results_cont_table <- data.frame()
  # for (k in 1:length(results_wqs_cont)) {
  #   
  #   temp_sum <- summary(results_wqs_cont[[k]])
  #   temp_sum2 <- temp_sum[["coefficients"]]
  #   
  #   strings <- strsplit(names(results_wqs_cont)[k], "\\|")[[1]]
  #   
  #   r <- data.frame(f1 = as.character(results_wqs_cont[[k]][["formula"]][3]),
  #                   cov = row.names(temp_sum2),
  #                   exp = strings[1],
  #                   type = strings[2],
  #                   out = strings[3],
  #                   method = strings[5],
  #                   direction = strings[4],
  #                   sample = sample_name,
  #                   temp_sum2,
  #                   t(results_wqs_cont[[k]][["final_weights"]][["mean_weight"]]))
  #   
  #   weight_names <- as.character(results_wqs_cont[[k]][["final_weights"]][["mix_name"]])
  #   colnames(r) <- c("formula","cov","exp","type","out","method","direction","sample","estimate","se","z","p",weight_names)
  #   
  #   wqs_results_cont_table <- bind_rows(wqs_results_cont_table, r)
  # }
  # #### 连续指标 ####
  
  openxlsx::write.xlsx(wqs_results_bin_table, paste0("results/correlations/wqs/wqs_results_bin_(",sample_name,").xlsx"))
  # openxlsx::write.xlsx(wqs_results_cont_table, paste0("results/correlations/wqs/wqs_results_cont_(",sample_name,").xlsx"))
}
#### WQS ####
