library(data.table)
library(dplyr)
# library(ggplot2)
# library(reshape2)
# library(ppcor)
library(survival)
library(rms)

setwd("your_file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/analyte_measurements_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/gut_microbial_composition_function_pathway_profiles_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")

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


#### RCS (EDC INDEX & 19 EDC与CVD, CKD, DM结局表型) ----
rcs_results_list <- list()
for (i in c("phy_edc_temp")) { # 不同亚组
  # i <- "phy_edc_temp"
  
  sample_name <- i  # 提取subgroup的名称
  
  COV <- c("age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg") # 设置COV校正变量
  
  
  for (j in c(edc_traits_log10)) { # 19 EDC
    EDC <- j # 设置EDC
    print(paste0("RCS: ",sample_name," || ",j,": ",which(c(edc_traits_log10) == j)," out of ",length(c(edc_traits_log10))))
    
    for (k in c(3,4)) { # 不同节点选择
      KNOT <- k # 设置KNOT
      
      for (l in c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021")) { # COX模型结局
        
        if(l == "cvd_incident_1021"){
          CENSOR <- "cvd_incident_1021"
          TIME <- "timecvd_1021"
        }else if(l == "cvd_incident_1014"){
          CENSOR <- "cvd_incident_1014"
          TIME <- "timecvd_1014"
        }else if(l == "ckd_incident_1014"){
          CENSOR <- "ckd_incident_1014"
          TIME <- "timeckd_1014"
        }else{
          CENSOR <- "dm_incident_1014"
          TIME <- "timedm_1014"
        }
        
        # 如果变量是常数，从数据集中删除它们
        cols <- c(TIME, CENSOR, EDC, COV)
        phy_edc_cox <- phy_edc_dat[,cols]
        phy_edc_cox <- na.omit(phy_edc_cox)
        
        # 分类协变量转换为因子
        phy_edc_cox$sex_b_rev <- factor(phy_edc_cox$sex_b_rev) # 0/1（女/男）
        phy_edc_cox$smk1_b <- factor(phy_edc_cox$smk1_b)
        phy_edc_cox$drk1_b <- factor(phy_edc_cox$drk1_b)
        phy_edc_cox$high_edu_b <- factor(phy_edc_cox$high_edu_b)
        phy_edc_cox$paactive3_g_b <- factor(phy_edc_cox$paactive3_g_b)
        phy_edc_cox$high_fruveg <- factor(phy_edc_cox$high_fruveg)
        
        colnames(phy_edc_cox)[3] <- "EDC" # 把所有EDC统一名称，方便后续作图
        # 对数据进行打包，整理
        dd <- datadist(phy_edc_cox) #为后续程序设定数据环境
        options(datadist='dd') #为后续程序设定数据环境
        ### 使用构建的公式拟合Cox模型 (使用 as.formula 和 paste0 来动态引用列名)
        formula <- as.formula(paste0("Surv(", TIME, ", ", CENSOR, ") ~ rcs(EDC,",KNOT,") + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg"))
        # 拟合Cox比例风险模型，使用rcs函数添加限制性立方样条
        rcs_fit <- cph(formula, data=phy_edc_cox)
        # 保存模型
        rcs_results_list[[paste0(sample_name,"|",EDC,"|",CENSOR,"|adj|cox|",KNOT)]] <- rcs_fit
      }
    }
  }
}
# 保存原始RCS分析数据 #
saveRDS(rcs_results_list, paste0("results/cox/rcs_results.rds"))

#### RCS (EDC INDEX & 19 EDC与CVD, CKD, DM结局表型) ####
