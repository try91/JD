library(data.table)
library(dplyr)
library(survival)

setwd("your_file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/analyte_measurements_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/gut_microbial_composition_function_pathway_profiles_dat_20261006.txt", header = TRUE)

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


#### COX EDC与incidence CVD (2010~2021), CKD (2010~2014), DM (2010~2014) 关系 ----
### 构建分析EDC与新发病的cox模型 (不校正) ###
cox_unadj <- function(DAT, SAMPLE, EDC, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, paste0(EDC,"_log10"))
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 使用构建的公式拟合Cox模型
  formula_log10 <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", EDC, "_log10"))
  
  ### 分析模型
  cox_fit_log10 <- coxph(formula_log10, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",EDC,"_log10","|",CENSOR,"|unadj")]] <- summary(cox_fit_log10)
  
  return(cox_results)
}
### 构建分析EDC与新发病的cox模型 (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
cox_adj <- function(DAT, SAMPLE, EDC, COV, MODEL, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, paste0(EDC,"_log10"), COV)
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_cox$high_fruveg <- factor(phy_edc_cox$high_fruveg)
  
  ### 使用构建的公式拟合Cox模型 (使用 as.formula 和 paste0 来动态引用列名)
  formula_log10 <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", EDC, "_log10", MODEL))
  
  ### 分析模型
  cox_fit_log10 <- coxph(formula_log10, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",EDC,"_log10","|",CENSOR,"|adj")]] <- summary(cox_fit_log10)
  
  return(cox_results)
}
### 构建分析EDC INDEX与新发病的cox模型 (不校正) ###
cox_edc_index_unadj <- function(DAT, SAMPLE, EDC_INDEX, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, EDC_INDEX)
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 使用构建的公式拟合Cox模型
  formula <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", EDC_INDEX))
  
  ### 分析模型
  cox_fit <- coxph(formula, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",EDC_INDEX,"|",CENSOR,"|unadj")]] <- summary(cox_fit)
  
  return(cox_results)
}
### 构建分析EDC INDEX与新发病的cox模型 (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
cox_edc_index_adj <- function(DAT, SAMPLE, EDC_INDEX, COV, MODEL, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, EDC_INDEX, COV)
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_cox$high_fruveg <- factor(phy_edc_cox$high_fruveg)
  
  ### 使用构建的公式拟合Cox模型 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", EDC_INDEX, MODEL))
  
  ### 分析模型
  cox_fit <- coxph(formula, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",EDC_INDEX,"|",CENSOR,"|adj")]] <- summary(cox_fit)
  
  return(cox_results)
}


### 使用as.formula和paste0构建公式 ###
covariate <- list()
covariate[[1]] <- paste0(" + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg")


### COX分析 ###
cox_results_edc_incident_list <- list()
cox_results_edc_index_incident_list <- list()
for (i in c("phy_edc_temp")) { # 不同亚组
  # i <- "phy_edc_temp"
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
  phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
  phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
  
  # 单个EDC浓度分析
  for (j in edc_traits) { # EDC
    # j <- edc_traits[1]
    print(paste0("COX: ",sample_name," || ",j,": ",which(edc_traits == j)," out of ",length(edc_traits)))
    
    cox_results_cvd21_unadj <- cox_unadj(phy_edc_dat, sample_name, j, "timecvd_1021", "cvd_incident_1021")
    cox_results_cvd21_adj <- cox_adj(phy_edc_dat, sample_name, j, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timecvd_1021", "cvd_incident_1021")
    
    cox_results_cvd14_unadj <- cox_unadj(phy_edc_dat, sample_name, j, "timecvd_1014", "cvd_incident_1014")
    cox_results_cvd14_adj <- cox_adj(phy_edc_dat, sample_name, j, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timecvd_1014", "cvd_incident_1014")
    
    cox_results_ckd15_unadj <- cox_unadj(phy_edc_dat, sample_name, j, "timeckd_1014", "ckd_incident_1014")
    cox_results_ckd15_adj <- cox_adj(phy_edc_dat, sample_name, j, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timeckd_1014", "ckd_incident_1014")
    
    cox_results_dm15_unadj<- cox_unadj(phy_edc_dat, sample_name, j, "timedm_1014", "dm_incident_1014")
    cox_results_dm15_adj <- cox_adj(phy_edc_dat, sample_name, j, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timedm_1014", "dm_incident_1014")
    
    cox_results_edc_incident_list <- do.call(c, list(cox_results_edc_incident_list, 
                                                     cox_results_cvd21_unadj, cox_results_cvd21_adj, 
                                                     cox_results_cvd14_unadj, cox_results_cvd14_adj, 
                                                     cox_results_ckd15_unadj, cox_results_ckd15_adj, 
                                                     cox_results_dm15_unadj, cox_results_dm15_adj))
    
  }
  
  # EDC INDEX分析
  for (k in c(edc_index_b_keep)) { # EDC_INDEX
    
    print(paste0("COX: ",sample_name," || ",k,": ",which(c(edc_index_b_keep) == k)," out of ",length(c(edc_index_b_keep))))
    
    cox_results_edc_index_cvd21_unadj <- cox_edc_index_unadj(phy_edc_dat, sample_name, k, "timecvd_1021", "cvd_incident_1021")
    cox_results_edc_index_cvd21_adj <- cox_edc_index_adj(phy_edc_dat, sample_name, k, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timecvd_1021", "cvd_incident_1021")
    
    cox_results_edc_index_cvd14_unadj <- cox_edc_index_unadj(phy_edc_dat, sample_name, k, "timecvd_1014", "cvd_incident_1014")
    cox_results_edc_index_cvd14_adj <- cox_edc_index_adj(phy_edc_dat, sample_name, k, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timecvd_1014", "cvd_incident_1014")
    
    cox_results_edc_index_ckd15_unadj <- cox_edc_index_unadj(phy_edc_dat, sample_name, k, "timeckd_1014", "ckd_incident_1014")
    cox_results_edc_index_ckd15_adj <- cox_edc_index_adj(phy_edc_dat, sample_name, k, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timeckd_1014", "ckd_incident_1014")
    
    cox_results_edc_index_dm15_unadj<- cox_edc_index_unadj(phy_edc_dat, sample_name, k, "timedm_1014", "dm_incident_1014")
    cox_results_edc_index_dm15_adj <- cox_edc_index_adj(phy_edc_dat, sample_name, k, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], "timedm_1014", "dm_incident_1014")
    
    cox_results_edc_index_incident_list <- do.call(c, list(cox_results_edc_index_incident_list, 
                                                           cox_results_edc_index_cvd21_unadj, cox_results_edc_index_cvd21_adj, 
                                                           cox_results_edc_index_cvd14_unadj, cox_results_edc_index_cvd14_adj, 
                                                           cox_results_edc_index_ckd15_unadj, cox_results_edc_index_ckd15_adj, 
                                                           cox_results_edc_index_dm15_unadj, cox_results_edc_index_dm15_adj))
    
  }
  
}
# 保存原始COX分析数据 #
saveRDS(cox_results_edc_incident_list, paste0("results/cox/cox_results_edc_incident.rds"))
saveRDS(cox_results_edc_index_incident_list, paste0("results/cox/cox_results_edc_index_incident.rds"))


### 整理COX分析数据 ###
# 读取原始COX分析数据 #
cox_results_edc_incident_list <- readRDS("results/cox/cox_results_edc_incident.rds")
cox_results_edc_incident_all <- data.frame()
for (i in 1:length(cox_results_edc_incident_list)){
  
  list_name <- names(cox_results_edc_incident_list[i])
  
  sample <- sub("\\|.*", "", list_name)
  exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
  edc <- gsub("_log10", "", exp)
  if(grepl("_log10",exp)){
    edc_type <- "log10"
  }
  out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
  adj <- sub(".*\\|", "", list_name)
  
  cox_results_edc_incident_list_temp <- cox_results_edc_incident_list[[i]]
  # cox_results_temp <- data.frame()
  
  cox_results_edc_incident_temp <- data.frame(cox_results_edc_incident_list_temp[["coefficients"]])
  cox_results_edc_incident_temp$rowname <- row.names(cox_results_edc_incident_temp)
  
  cox_results_edc_incident_temp$exposure <- edc
  cox_results_edc_incident_temp$edc_type <- edc_type
  cox_results_edc_incident_temp$outcome <- out
  cox_results_edc_incident_temp$method <- "cox"
  cox_results_edc_incident_temp$adjust <- adj
  cox_results_edc_incident_temp$sample <- sample
  cox_results_edc_incident_temp$n <- cox_results_edc_incident_list_temp[["n"]]
  row.names(cox_results_edc_incident_temp) <-NULL
  
  cox_results_edc_incident_all <- rbind(cox_results_edc_incident_all, cox_results_edc_incident_temp)
}

cox_results_edc_index_incident_list <- readRDS("results/cox/cox_results_edc_index_incident.rds")
cox_results_edc_index_incident_all <- data.frame()
for (i in 1:length(cox_results_edc_index_incident_list)){
  
  list_name <- names(cox_results_edc_index_incident_list[i])
  
  sample <- sub("\\|.*", "", list_name)
  exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
  edc <- gsub("_log10", "", exp)
  edc_type <- "edc_index"
  
  out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
  adj <- sub(".*\\|", "", list_name)
  
  cox_results_edc_index_incident_list_temp <- cox_results_edc_index_incident_list[[i]]
  # cox_results_temp <- data.frame()
  
  cox_results_edc_index_incident_temp <- data.frame(cox_results_edc_index_incident_list_temp[["coefficients"]])
  cox_results_edc_index_incident_temp$rowname <- row.names(cox_results_edc_index_incident_temp)
  
  cox_results_edc_index_incident_temp$exposure <- edc
  cox_results_edc_index_incident_temp$edc_type <- edc_type
  cox_results_edc_index_incident_temp$outcome <- out
  cox_results_edc_index_incident_temp$method <- "cox"
  cox_results_edc_index_incident_temp$adjust <- adj
  cox_results_edc_index_incident_temp$sample <- sample
  cox_results_edc_index_incident_temp$n <- cox_results_edc_index_incident_list_temp[["n"]]
  row.names(cox_results_edc_index_incident_temp) <-NULL
  
  cox_results_edc_index_incident_all <- rbind(cox_results_edc_index_incident_all, cox_results_edc_index_incident_temp)
}

# 保存汇总cox结果数据 #
openxlsx::write.xlsx(cox_results_edc_incident_all,"results/cox/cox_results_edc_incident.xlsx")
openxlsx::write.xlsx(cox_results_edc_index_incident_all,"results/cox/cox_results_edc_index_incident.xlsx")
#### COX EDC与incidence CVD (2010~2021), CKD (2010~2014), DM (2010~2014) 关系 ####

#### Logistic regression EDC与二分类结局表型 ----
### 构建分析EDC与二分类表型的logistic模型 (不校正) ###
logistic_unadj <- function(DAT, SAMPLE, EDC, OUT){
  
  cols <- c(OUT, paste0(EDC,"_log10"), paste0(EDC,"_detected"))
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 使用构建的公式拟合Logistic模型
  formula_log10 <- as.formula(paste0(OUT, " ~ ", EDC, "_log10"))
  
  ### 分析模型
  glm_log10 <- glm(formula_log10, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",EDC,"_log10","|",OUT,"|unadj")]] <- summary(glm_log10)
  
  return(logistic_results)
}
### 构建分析EDC与二分类表型的logistic模型 (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
logistic_adj <- function(DAT, SAMPLE, EDC, COV, MODEL, OUT){
  
  cols <- c(OUT, paste0(EDC,"_log10"), paste0(EDC,"_detected"), COV)
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_logistic$high_fruveg <- factor(phy_edc_logistic$high_fruveg)
  
  ### 使用构建的公式拟合Cox模型 (使用 as.formula 和 paste0 来动态引用列名)
  formula_log10 <- as.formula(paste0(OUT, " ~ ", EDC, "_log10", MODEL))
  
  ### 分析模型
  glm_log10 <- glm(formula_log10, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",EDC,"_log10","|",OUT,"|adj")]] <- summary(glm_log10)
  
  return(logistic_results)
}
### 构建分析EDC INDEX与二分类表型的logistic模型 (不校正) ###
logistic_edc_index_unadj <- function(DAT, SAMPLE, EDC_INDEX, OUT){
  
  cols <- c(OUT, EDC_INDEX)
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 使用构建的公式拟合Logistic模型
  formula <- as.formula(paste0(OUT, " ~ ", EDC_INDEX))
  
  ### 分析模型
  glm <- glm(formula, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",EDC_INDEX,"|",OUT,"|unadj")]] <- summary(glm)
  
  return(logistic_results)
}
### 构建分析EDC INDEX与二分类表型的logistic模型 (校正："age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") ###
logistic_edc_index_adj <- function(DAT, SAMPLE, EDC_INDEX, COV, MODEL, OUT){
  
  cols <- c(OUT, EDC_INDEX, COV)
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_logistic$high_fruveg <- factor(phy_edc_logistic$high_fruveg)
  
  ### 使用构建的公式拟合Cox模型 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0(OUT, " ~ ", EDC_INDEX, MODEL))
  
  ### 分析模型
  glm <- glm(formula, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",EDC_INDEX,"|",OUT,"|adj")]] <- summary(glm)
  
  return(logistic_results)
}


### 使用as.formula和paste0构建公式 ###
covariate <- list()
covariate[[1]] <- paste0(" + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg")


### Logistic分析 ###
logistic_results_edc_incident_all <- data.frame()
logistic_results_edc_index_incident_all <- data.frame()
for (i in c("phy_edc_temp")) { # 不同亚组
  # i <- "phy_edc_temp"
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
  phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
  phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
  
  # 单个EDC浓度分析
  for (j in edc_traits) { # EDC
    print(paste0("Logistic: ",sample_name," || ",j,": ",which(edc_traits == j)," out of ",length(edc_traits)))
    
    for (k in c("cvd_b","ckd_b","dm_b","as_imt_b","hpt_b","nafld_b","ob_b","abob_b","dyslip_b","hua_b","ir_b","mets_b")) {
      
      logistic_results_edc_incident_unadj <- logistic_unadj(phy_edc_dat, sample_name, j, k)
      logistic_results_edc_incident_adj <- logistic_adj(phy_edc_dat, sample_name, j, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], k)
      
      ### 整理Logistic分析数据 ###
      # unadj #
      {
        for (a in 1:length(logistic_results_edc_incident_unadj)){
          list_name <- names(logistic_results_edc_incident_unadj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10|_quantile|_detected", "", exp)
          if(grepl("_log10",exp)){
            edc_type <- "log10"
          }else if(grepl("_quantile",exp)){
            edc_type <- "quantile"
          }else{
            edc_type <- "detected"
          }
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          logistic_results_edc_incident_list_temp <- logistic_results_edc_incident_unadj[[a]]
          
          logistic_results_edc_incident_temp <- data.frame(logistic_results_edc_incident_list_temp[["coefficients"]])
          logistic_results_edc_incident_temp$rowname <- row.names(logistic_results_edc_incident_temp)
          
          logistic_results_edc_incident_temp$exposure <- edc
          logistic_results_edc_incident_temp$edc_type <- edc_type
          logistic_results_edc_incident_temp$outcome <- out
          logistic_results_edc_incident_temp$method <- "logistic"
          logistic_results_edc_incident_temp$adjust <- adj
          logistic_results_edc_incident_temp$sample <- sample
          logistic_results_edc_incident_temp$n <- logistic_results_edc_incident_list_temp[["df.null"]] + 1
          colnames(logistic_results_edc_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_edc_incident_temp) <-NULL
          
          logistic_results_edc_incident_all <- rbind(logistic_results_edc_incident_all, logistic_results_edc_incident_temp)
        }
      }
      # adj #
      {
        for (a in 1:length(logistic_results_edc_incident_adj)){
          list_name <- names(logistic_results_edc_incident_adj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10|_quantile|_detected", "", exp)
          if(grepl("_log10",exp)){
            edc_type <- "log10"
          }else if(grepl("_quantile",exp)){
            edc_type <- "quantile"
          }else{
            edc_type <- "detected"
          }
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          logistic_results_edc_incident_list_temp <- logistic_results_edc_incident_adj[[a]]
          
          logistic_results_edc_incident_temp <- data.frame(logistic_results_edc_incident_list_temp[["coefficients"]])
          logistic_results_edc_incident_temp$rowname <- row.names(logistic_results_edc_incident_temp)
          
          logistic_results_edc_incident_temp$exposure <- edc
          logistic_results_edc_incident_temp$edc_type <- edc_type
          logistic_results_edc_incident_temp$outcome <- out
          logistic_results_edc_incident_temp$method <- "logistic"
          logistic_results_edc_incident_temp$adjust <- adj
          logistic_results_edc_incident_temp$sample <- sample
          logistic_results_edc_incident_temp$n <- logistic_results_edc_incident_list_temp[["df.null"]] + 1
          colnames(logistic_results_edc_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_edc_incident_temp) <-NULL
          
          logistic_results_edc_incident_all <- rbind(logistic_results_edc_incident_all, logistic_results_edc_incident_temp)
        }
      }
    }
  }
  
  # EDC INDEX分析
  for (l in c(edc_index_b_keep)) { # EDC_INDEX
    print(paste0("Logistic: ",sample_name," || ",l,": ",which(c(edc_index_b_keep) == l)," out of ",length(c(edc_index_b_keep))))
    
    for (m in c("cvd_b","ckd_b","dm_b","as_imt_b","hpt_b","nafld_b","ob_b","abob_b","dyslip_b","hua_b","ir_b","mets_b")) {
      
      logistic_results_edc_index_incident_unadj <- logistic_edc_index_unadj(phy_edc_dat, sample_name, l, m)
      logistic_results_edc_index_incident_adj <- logistic_edc_index_adj(phy_edc_dat, sample_name, l, c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg"), covariate[[1]], m)
      
      ### 整理Logistic分析数据 ###
      # unadj #
      {
        for (a in 1:length(logistic_results_edc_index_incident_unadj)){
          list_name <- names(logistic_results_edc_index_incident_unadj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10|_quantile|_detected", "", exp)
          edc_type <- "edc_index"
          
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          logistic_results_edc_index_incident_list_temp <- logistic_results_edc_index_incident_unadj[[a]]
          
          logistic_results_edc_index_incident_temp <- data.frame(logistic_results_edc_index_incident_list_temp[["coefficients"]])
          logistic_results_edc_index_incident_temp$rowname <- row.names(logistic_results_edc_index_incident_temp)
          
          logistic_results_edc_index_incident_temp$exposure <- edc
          logistic_results_edc_index_incident_temp$edc_type <- edc_type
          logistic_results_edc_index_incident_temp$outcome <- out
          logistic_results_edc_index_incident_temp$method <- "logistic"
          logistic_results_edc_index_incident_temp$adjust <- adj
          logistic_results_edc_index_incident_temp$sample <- sample
          logistic_results_edc_index_incident_temp$n <- logistic_results_edc_index_incident_list_temp[["df.null"]]
          colnames(logistic_results_edc_index_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_edc_index_incident_temp) <-NULL
          
          logistic_results_edc_index_incident_all <- rbind(logistic_results_edc_index_incident_all, logistic_results_edc_index_incident_temp)
        }
      }
      # adj #
      {
        for (a in 1:length(logistic_results_edc_index_incident_adj)){
          list_name <- names(logistic_results_edc_index_incident_adj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10|_quantile|_detected", "", exp)
          edc_type <- "edc_index"
          
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          logistic_results_edc_index_incident_list_temp <- logistic_results_edc_index_incident_adj[[a]]
          
          logistic_results_edc_index_incident_temp <- data.frame(logistic_results_edc_index_incident_list_temp[["coefficients"]])
          logistic_results_edc_index_incident_temp$rowname <- row.names(logistic_results_edc_index_incident_temp)
          
          logistic_results_edc_index_incident_temp$exposure <- edc
          logistic_results_edc_index_incident_temp$edc_type <- edc_type
          logistic_results_edc_index_incident_temp$outcome <- out
          logistic_results_edc_index_incident_temp$method <- "logistic"
          logistic_results_edc_index_incident_temp$adjust <- adj
          logistic_results_edc_index_incident_temp$sample <- sample
          logistic_results_edc_index_incident_temp$n <- logistic_results_edc_index_incident_list_temp[["df.null"]] + 1
          colnames(logistic_results_edc_index_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_edc_index_incident_temp) <-NULL
          
          logistic_results_edc_index_incident_all <- rbind(logistic_results_edc_index_incident_all, logistic_results_edc_index_incident_temp)
        }
      }
    }
  }
}
# 保存汇总Logistic结果数据 #
openxlsx::write.xlsx(logistic_results_edc_incident_all,"results/glm/logistic_results_edc_incident.xlsx")
openxlsx::write.xlsx(logistic_results_edc_index_incident_all,"results/glm/logistic_results_edc_index_incident.xlsx")
#### Logistic regression EDC与二分类结局表型 ####
