library(data.table)
library(dplyr)
library(survival)

setwd("your_file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/analyte_measurements_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/gut_microbial_composition_function_pathway_profiles_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  right_join(micro_dat, by = "ID")

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

#### COX MP4与incidence CVD (2014~2021), CKD (2010~2014), DM (2010~2014) 关系 ----
### 构建分析MP4 log10与新发病的cox模型 (不校正) ###
cox_mp4_unadj <- function(DAT, SAMPLE, MP4, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, MP4)
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 构建公式 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", MP4))
  
  ### 分析模型
  cox_fit <- coxph(formula, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",MP4,"|",CENSOR,"|unadj")]] <- summary(cox_fit)
  
  return(cox_results)
}
### 构建分析MP4 log10与新发病的cox模型 (校正："age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7") ###
cox_mp4_adj <- function(DAT, SAMPLE, MP4, COV, MODEL, TIME, CENSOR){
  
  cols <- c(TIME, CENSOR, MP4, COV)
  phy_edc_cox <- DAT[,cols]
  phy_edc_cox <- na.omit(phy_edc_cox)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_cox$high_fruveg <- factor(phy_edc_cox$high_fruveg)
  
  ### 构建公式 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ ", MP4, MODEL))
  
  ### 分析模型
  cox_fit <- coxph(formula, data = phy_edc_cox)
  
  # 保存模型摘要
  cox_results <- list()
  cox_results[[paste0(SAMPLE,"|",MP4,"|",CENSOR,"|adj")]] <- summary(cox_fit)
  
  return(cox_results)
}

### 使用as.formula和paste0构建公式 ###
covariate <- list()
covariate[[2]] <- paste0(" + age_f + sex_b_rev + smk1_f + drk1_f + high_edu_b + paactive3_g_f + high_fruveg + med_all7")

### COX分析 ###
cox_results_mp4_incident_list <- list()
for (i in c("phy_edc_temp0_3")) { # 不同亚组
  # i <- "phy_edc_temp0_3"
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
  phy_edc_dat$smk1_f <- factor(phy_edc_dat$smk1_f)
  phy_edc_dat$drk1_f <- factor(phy_edc_dat$drk1_f)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_f <- factor(phy_edc_dat$paactive3_g_f)
  phy_edc_dat$med_all7 <- factor(phy_edc_dat$med_all7)
  
  # MP4分析
  for (k in c(mp4_s_log10)) { # MP4
    
    print(paste0("COX: ",sample_name," || ",k,": ",which(c(mp4_s_log10) == k)," out of ",length(c(mp4_s_log10))))
    
    cox_results_mp4_cvd_unadj <- cox_mp4_unadj(phy_edc_dat, sample_name, k, "timecvd_1421", "cvd_incident_1421")
    cox_results_mp4_cvd_adj <- cox_mp4_adj(phy_edc_dat, sample_name, k, c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7"), covariate[[2]], "timecvd_1421", "cvd_incident_1421")
    
    cox_results_mp4_incident_list <- do.call(c, list(cox_results_mp4_incident_list, 
                                                     cox_results_mp4_cvd_unadj, cox_results_mp4_cvd_adj))
  }
  
}
# 保存原始COX分析数据 #
saveRDS(cox_results_mp4_incident_list, paste0("results/cox/cox_results_mp4_incident.rds"))

### 整理COX分析数据 ###
# 读取原始COX分析数据 #
cox_results_mp4_incident_list <- readRDS("results/cox/cox_results_mp4_incident.rds")
cox_results_mp4_incident_all <- data.frame()
for (i in 1:length(cox_results_mp4_incident_list)){
  
  list_name <- names(cox_results_mp4_incident_list[i])
  
  sample <- sub("\\|.*", "", list_name)
  exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
  edc <- gsub("_log10", "", exp)
  edc_type <- "mp4"
  
  out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
  adj <- sub(".*\\|", "", list_name)
  
  cox_results_mp4_incident_list_temp <- cox_results_mp4_incident_list[[i]]
  # cox_results_temp <- data.frame()
  
  cox_results_mp4_incident_temp <- data.frame(cox_results_mp4_incident_list_temp[["coefficients"]])
  cox_results_mp4_incident_temp$rowname <- row.names(cox_results_mp4_incident_temp)
  
  cox_results_mp4_incident_temp$exposure <- edc
  cox_results_mp4_incident_temp$edc_type <- edc_type
  cox_results_mp4_incident_temp$outcome <- out
  cox_results_mp4_incident_temp$method <- "cox"
  cox_results_mp4_incident_temp$adjust <- adj
  cox_results_mp4_incident_temp$sample <- sample
  cox_results_mp4_incident_temp$n <- cox_results_mp4_incident_list_temp[["n"]]
  row.names(cox_results_mp4_incident_temp) <-NULL
  
  cox_results_mp4_incident_all <- rbind(cox_results_mp4_incident_all, cox_results_mp4_incident_temp)
}

# 保存汇总cox结果数据 #
openxlsx::write.xlsx(cox_results_mp4_incident_all,"results/cox/cox_results_mp4_incident.xlsx")
#### COX MP4与incidence CVD (2014~2021), CKD (2010~2014), DM (2010~2014) 关系 ####

#### Logistic regression MP4与二分类结局表型 ----
### 构建分析MP4 log10与二分类表型的logistic模型 (不校正) ###
logistic_mp4_unadj <- function(DAT, SAMPLE, MP4, OUT){
  
  cols <- c(OUT, MP4)
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 构建公式 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0(OUT, " ~ ", MP4))
  
  ### 分析模型
  glm <- glm(formula, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",MP4,"|",OUT,"|unadj")]] <- summary(glm)
  
  return(logistic_results)
}
### 构建分析MP4 log10与二分类表型的logistic模型 (校正："age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7") ###
logistic_mp4_adj <- function(DAT, SAMPLE, MP4, COV, MODEL, OUT){
  
  cols <- c(OUT, MP4, COV)
  phy_edc_logistic <- DAT[,cols]
  phy_edc_logistic <- na.omit(phy_edc_logistic)
  
  ### 把水果蔬菜变量转换为因子
  phy_edc_logistic$high_fruveg <- factor(phy_edc_logistic$high_fruveg)
  
  ### 构建公式 (使用 as.formula 和 paste0 来动态引用列名)
  formula <- as.formula(paste0(OUT, " ~ ", MP4, MODEL))
  
  ### 分析模型
  glm <- glm(formula, family = binomial, data = phy_edc_logistic)
  
  # 保存模型摘要
  logistic_results <- list()
  logistic_results[[paste0(SAMPLE,"|",MP4,"|",OUT,"|adj")]] <- summary(glm)
  
  return(logistic_results)
}

### 使用as.formula和paste0构建公式 ###
covariate <- list()
covariate[[2]] <- paste0(" + age_f + sex_b_rev + smk1_f + drk1_f + high_edu_b + paactive3_g_f + high_fruveg + med_all7")

### Logistic分析 ###
logistic_results_mp4_incident_all <- data.frame()
for (i in c("phy_edc_temp0_3")) { # 不同亚组
  # i <- "phy_edc_temp0_3"
  
  sample_name <- i  # 提取subgroup的名称
  
  # 分类协变量转换为因子
  phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
  phy_edc_dat$smk1_f <- factor(phy_edc_dat$smk1_f)
  phy_edc_dat$drk1_f <- factor(phy_edc_dat$drk1_f)
  phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
  phy_edc_dat$paactive3_g_f <- factor(phy_edc_dat$paactive3_g_f)
  phy_edc_dat$med_all7 <- factor(phy_edc_dat$med_all7)
  
  # MP4分析
  for (l in c(mp4_s_log10)) { # MP4
    print(paste0("Logistic: ",sample_name," || ",l,": ",which(c(mp4_s_log10) == l)," out of ",length(c(mp4_s_log10))))
    
    for (m in c("ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report", phy_traits_cat)) {
      
      logistic_results_mp4_incident_unadj <- logistic_mp4_unadj(phy_edc_dat, sample_name, l, m)
      logistic_results_mp4_incident_adj <- logistic_mp4_adj(phy_edc_dat, sample_name, l, c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7"), covariate[[2]], m)
      
      ### 整理Logistic分析数据 ###
      # unadj #
      {
        for (a in 1:length(logistic_results_mp4_incident_unadj)){
          list_name <- names(logistic_results_mp4_incident_unadj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10", "", exp)
          edc_type <- "mp4"
          
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          
          logistic_results_mp4_incident_list_temp <- logistic_results_mp4_incident_unadj[[a]]
          
          logistic_results_mp4_incident_temp <- data.frame(logistic_results_mp4_incident_list_temp[["coefficients"]])
          logistic_results_mp4_incident_temp$rowname <- row.names(logistic_results_mp4_incident_temp)
          
          logistic_results_mp4_incident_temp$exposure <- edc
          logistic_results_mp4_incident_temp$edc_type <- edc_type
          logistic_results_mp4_incident_temp$outcome <- out
          logistic_results_mp4_incident_temp$method <- "logistic"
          logistic_results_mp4_incident_temp$adjust <- adj
          logistic_results_mp4_incident_temp$sample <- sample
          logistic_results_mp4_incident_temp$n <- logistic_results_mp4_incident_list_temp[["df.null"]] + 1
          colnames(logistic_results_mp4_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_mp4_incident_temp) <-NULL
          
          logistic_results_mp4_incident_all <- rbind(logistic_results_mp4_incident_all, logistic_results_mp4_incident_temp)
        }
      }
      # adj #
      {
        for (a in 1:length(logistic_results_mp4_incident_adj)){
          list_name <- names(logistic_results_mp4_incident_adj[a])
          
          sample <- sub("\\|.*", "", list_name)
          exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
          edc <- gsub("_log10", "", exp)
          edc_type <- "mp4"
          
          out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
          adj <- sub(".*\\|", "", list_name)
          
          
          logistic_results_mp4_incident_list_temp <- logistic_results_mp4_incident_adj[[a]]
          
          logistic_results_mp4_incident_temp <- data.frame(logistic_results_mp4_incident_list_temp[["coefficients"]])
          logistic_results_mp4_incident_temp$rowname <- row.names(logistic_results_mp4_incident_temp)
          
          logistic_results_mp4_incident_temp$exposure <- edc
          logistic_results_mp4_incident_temp$edc_type <- edc_type
          logistic_results_mp4_incident_temp$outcome <- out
          logistic_results_mp4_incident_temp$method <- "logistic"
          logistic_results_mp4_incident_temp$adjust <- adj
          logistic_results_mp4_incident_temp$sample <- sample
          logistic_results_mp4_incident_temp$n <- logistic_results_mp4_incident_list_temp[["df.null"]] + 1
          colnames(logistic_results_mp4_incident_temp)[c(1,2,3,4)] <- c("coef","se.coef.","t or z","p")
          row.names(logistic_results_mp4_incident_temp) <-NULL
          
          logistic_results_mp4_incident_all <- rbind(logistic_results_mp4_incident_all, logistic_results_mp4_incident_temp)
        }
      }
      
    }
  }
}
# 保存汇总Logistic结果数据 #
openxlsx::write.xlsx(logistic_results_mp4_incident_all,"results/glm/logistic_results_mp4_incident.xlsx")
#### Logistic regression MP4与二分类结局表型 ####
