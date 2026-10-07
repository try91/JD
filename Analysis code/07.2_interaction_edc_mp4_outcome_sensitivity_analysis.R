library(data.table)
library(dplyr)
library(chemometrics)
library(ggplot2)
library(survival)
library(interactionR)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))

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

# 2021、2014死亡和新发表型 (分类)
phy_incident_cat <- c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014")    # 新发 cvd, ckd, dm 去除基线 case (只做EDC对outcome，不做cvd_incident_1421)
phy_incident_time <- c("timecvd_1021","timecvd_1014","timeckd_1014","timedm_1014")
phy_censor_cat <- c("censorall_1021","censorall_1014")
phy_censor_time <- c("timeall_1021","timeall_1014")
# 2014、2010表型 (分类)
phy_out_cat <- c("cvd_f","ckd_f","dm_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014新发case，横断面数据)
phy_traits_cat <- c("cvd_b","ckd_b","dm_b","as_imt_f","as_imt_b","hpt_f","hpt_b","nafld_f","nafld_b",
                    "ob_f","ob_b","abob_f","abob_b","dyslip_f","dyslip_b","hua_f","hua_b","ir_f","ir_b","mets_f","mets_b")
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
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f")
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

#### 读取EDC INDEX数据 ----
edc_count_b <- read.table("jiading/sourceDataEDCs/index/edc_count_2010_20260313.txt", header = TRUE)
edc_score_b <- read.table("jiading/sourceDataEDCs/index/edc_score_2010_20260313.txt", header = TRUE)

edc_count_f <- read.table("jiading/sourceDataEDCs/index/edc_count_2014_20260313.txt", header = TRUE)
edc_score_f <- read.table("jiading/sourceDataEDCs/index/edc_score_2014_20260313.txt", header = TRUE)

# EDC INDEX 变量名
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")
edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### 读取EDC INDEX数据 ####

#### Link EDC INDEX ----
phy_edc_dat <- phy_edc_temp_list[["phy_edc_temp0_3"]] # 提取subgroup
sample_name <- "phy_edc_temp0_3"  # 提取subgroup的名称

phy_edc_dat <- left_join(phy_edc_dat,edc_count_f,by="ID") %>%
  left_join(edc_score_f,by="ID") %>%
  left_join(edc_count_b,by="ID") %>%
  left_join(edc_score_b,by="ID")

# 分类协变量转换为因子 #
phy_edc_dat$sex_b_rev <- factor(phy_edc_dat$sex_b_rev) # 0/1（女/男）
phy_edc_dat$smk1_b <- factor(phy_edc_dat$smk1_b)
phy_edc_dat$smk1_f <- factor(phy_edc_dat$smk1_f)
phy_edc_dat$drk1_b <- factor(phy_edc_dat$drk1_b)
phy_edc_dat$drk1_f <- factor(phy_edc_dat$drk1_f)
phy_edc_dat$high_edu_b <- factor(phy_edc_dat$high_edu_b)
phy_edc_dat$paactive3_g_b <- factor(phy_edc_dat$paactive3_g_b)
phy_edc_dat$paactive3_g_f <- factor(phy_edc_dat$paactive3_g_f)
phy_edc_dat$high_fruveg <- factor(phy_edc_dat$high_fruveg) ### 把水果蔬菜变量转换为因子
phy_edc_dat$med_all7 <- factor(phy_edc_dat$med_all7)
#### Link EDC INDEX ####

#### 对筛选得到的interaction分析结果进行敏感性分析 (有相乘交互的结果) ----
results_multiplicative_interaction_sig_keep <- readxl::read_xlsx("results/cox/interaction/multi_interaction_sig_20260728.xlsx")
unique(results_multiplicative_interaction_sig_keep$exposure)
unique(results_multiplicative_interaction_sig_keep$mediator)
unique(results_multiplicative_interaction_sig_keep$outcome)
unique(results_multiplicative_interaction_sig_keep$keep_exp_med_out) %>% sort()


dat_mp4_detected_rate <- data.frame()
for (i in unique(results_multiplicative_interaction_sig_keep$mediator)) {
  
  i_bin <- paste0(gsub("_log10","",i),"_bin")
  temp <- data.frame(species = i,
                     detect_rate = nrow(phy_edc_dat[phy_edc_dat[[i_bin]] == 1,])/nrow(phy_edc_dat))
  
  dat_mp4_detected_rate <- rbind(dat_mp4_detected_rate, temp)
}

mp4_tertile <- dat_mp4_detected_rate[dat_mp4_detected_rate$detect_rate >= 0.6667,] # 可以三分位的菌
mp4_tertile_names <- mp4_tertile$species
mp4_tertile_not <- dat_mp4_detected_rate[dat_mp4_detected_rate$detect_rate < 0.6667,] # 不可以三分位的菌
mp4_tertile_not_names <- mp4_tertile_not$species
#### 对筛选得到的interaction分析结果进行敏感性分析 (有相乘交互的结果) ####


#### interaction分析 (EDC-outcome/菌分组) ----
pair_mp4_out_for_interaction <- results_multiplicative_interaction_sig_keep[,c("exposure","mediator","outcome")]
colnames(pair_mp4_out_for_interaction) <- c("exp_name","mediator","out_name")
pair_mp4_out_for_interaction$exp_med_out <- paste0(pair_mp4_out_for_interaction$exp_name,"|",pair_mp4_out_for_interaction$mediator,"|",pair_mp4_out_for_interaction$out_name)

# interaction 分析
exp_var <- unique(pair_mp4_out_for_interaction$exp_name)
med_var <- unique(pair_mp4_out_for_interaction$mediator)
out_var <- unique(pair_mp4_out_for_interaction$out_name)
exp_med_out <- unique(pair_mp4_out_for_interaction$exp_med_out) # 86对潜在的EDC-菌-Outcome pairs (有相乘交互的结果)

results_all <- data.frame()
for (x in out_var) {
  for (i in exp_var) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," || EDC (",which(exp_var == i)," out of ",length(exp_var),"): ",i,", Outcome (",which(out_var == x)," out of ",length(out_var),"): ",x))
    for (j in med_var) {
      
      # x <- out_var[1]
      # i <- exp_var[1]
      # j <- med_var[1]
      if(!paste0(i,"|",j,"|",x) %in% exp_med_out){
        next
      }
      
      ## 分析样本选取 (由于各步骤间样本量要统一，因此在计算前首先选择各自的样本，排除结局缺失的项) ##
      if(x == "cvd_incident_1021"){
        outcome <- "cvd_incident_1021"
        cols <- c("cvd_incident_1021", "timecvd_1021", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else if(x == "ckd_incident_1014"){
        outcome <- "ckd_incident_1014"
        cols <- c("ckd_incident_1014", "timeckd_1014", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else{
        outcome <- "dm_incident_1014"
        cols <- c("dm_incident_1014", "timedm_1014", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      colnames(phy_edc_dat_temp)[c(1,2,3,4)] <- c("CENSOR","TIME","EDC","MP4")
      
      ## 根据菌群丰度进行三分类 ##
      {
        ## 根据菌群丰度进行三分类 ##
        phy_edc_dat_temp$mp4_group <- NULL
        
        if(j %in% mp4_tertile_names){
          # 求三分位数（2个分割点：33%、67%）
          mp4_s_q3 <- quantile(phy_edc_dat_temp$MP4, probs = c(1/3, 2/3))
          phy_edc_dat_temp$mp4_group <- ifelse(phy_edc_dat_temp$MP4 <= mp4_s_q3[1], 0, 
                                               ifelse(phy_edc_dat_temp$MP4 <= mp4_s_q3[2], 1, 2)) # 0为低菌群丰度组，1为中菌群丰度组，2为高菌群丰度组
          
        }else if(j %in% mp4_tertile_not_names){
          mp4_s_min <- min(phy_edc_dat_temp$MP4)
          mp4_s_median <- median(phy_edc_dat_temp[phy_edc_dat_temp$MP4 != mp4_s_min,]$MP4)
          
          phy_edc_dat_temp$mp4_group <- ifelse(phy_edc_dat_temp$MP4 == mp4_s_min, 0, 
                                               ifelse(phy_edc_dat_temp$MP4 <= mp4_s_median, 1, 2)) # 0为低菌群丰度组，1为中菌群丰度组，2为高菌群丰度组
          
        }
        # table(phy_edc_dat_temp$mp4_group)
        # by(phy_edc_dat_temp$MP4, phy_edc_dat_temp$mp4_group, summary)
        
        
        ## 根据EDC浓度进行二分类 ##
        if(i %in% edc_traits8_log10){ # 检出率<50%, detected or not 分组
          edc_name <- gsub("_log10","",i)
          phy_edc_dat_temp$edc_group <- phy_edc_dat_temp[[paste0(edc_name,"_detected")]] # 0为EDC未检出组，1为EDC检出组
        }else{ # 检出率>50% or EDC index, 中位数分组
          edc_median <- median(phy_edc_dat_temp$EDC)
          phy_edc_dat_temp$edc_group <- ifelse(phy_edc_dat_temp$EDC <= edc_median, 0, 1) # 0为低EDC浓度组，1为高EDC浓度组
        }
        
        # 分组EDC-低菌群丰度亚组分析 (direction: "1_2")
        {
          phy_edc_dat_temp_low <- phy_edc_dat_temp[phy_edc_dat_temp$mp4_group == 0,]
          # EDC分组转换为因子
          phy_edc_dat_temp_low$edc_group <- factor(phy_edc_dat_temp_low$edc_group)
          
          cox_fit_low <- coxph(Surv(TIME, CENSOR) ~ edc_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_low)
          cox_fit_low_summary <- summary(cox_fit_low)
          
          cox_results_low_temp <- data.frame(cox_fit_low_summary[["coefficients"]])
          cox_results_low_temp$rowname <- row.names(cox_results_low_temp)
          cox_results_low_temp$exposure <- i
          cox_results_low_temp$mediator <- j
          cox_results_low_temp$outcome <- outcome
          cox_results_low_temp$method <- "cox"
          cox_results_low_temp$adjust <- "adj"
          cox_results_low_temp$sample <- sample_name
          cox_results_low_temp$n <- nrow(phy_edc_dat_temp_low)
          cox_results_low_temp$group <- "mp4_low"
          cox_results_low_temp$direction <- "1_2"
          row.names(cox_results_low_temp) <-NULL
        }
        # 分组EDC-中菌群丰度亚组分析 (direction: "1_2")
        {
          phy_edc_dat_temp_middle <- phy_edc_dat_temp[phy_edc_dat_temp$mp4_group == 1,]
          # EDC分组转换为因子
          phy_edc_dat_temp_middle$edc_group <- factor(phy_edc_dat_temp_middle$edc_group)
          
          cox_fit_middle <- coxph(Surv(TIME, CENSOR) ~ edc_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_middle)
          cox_fit_middle_summary <- summary(cox_fit_middle)
          
          cox_results_middle_temp <- data.frame(cox_fit_middle_summary[["coefficients"]])
          cox_results_middle_temp$rowname <- row.names(cox_results_middle_temp)
          cox_results_middle_temp$exposure <- i
          cox_results_middle_temp$mediator <- j
          cox_results_middle_temp$outcome <- outcome
          cox_results_middle_temp$method <- "cox"
          cox_results_middle_temp$adjust <- "adj"
          cox_results_middle_temp$sample <- sample_name
          cox_results_middle_temp$n <- nrow(phy_edc_dat_temp_middle)
          cox_results_middle_temp$group <- "mp4_middle"
          cox_results_middle_temp$direction <- "1_2"
          row.names(cox_results_middle_temp) <-NULL
        }
        # 分组EDC-高菌群丰度亚组分析 (direction: "1_2")
        {
          phy_edc_dat_temp_high <- phy_edc_dat_temp[phy_edc_dat_temp$mp4_group == 2,]
          # EDC分组转换为因子
          phy_edc_dat_temp_high$edc_group <- factor(phy_edc_dat_temp_high$edc_group)
          
          cox_fit_high <- coxph(Surv(TIME, CENSOR) ~ edc_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_high)
          cox_fit_high_summary <- summary(cox_fit_high)
          
          cox_results_high_temp <- data.frame(cox_fit_high_summary[["coefficients"]])
          cox_results_high_temp$rowname <- row.names(cox_results_high_temp)
          cox_results_high_temp$exposure <- i
          cox_results_high_temp$mediator <- j
          cox_results_high_temp$outcome <- outcome
          cox_results_high_temp$method <- "cox"
          cox_results_high_temp$adjust <- "adj"
          cox_results_high_temp$sample <- sample_name
          cox_results_high_temp$n <- nrow(phy_edc_dat_temp_high)
          cox_results_high_temp$group <- "mp4_high"
          cox_results_high_temp$direction <- "1_2"
          row.names(cox_results_high_temp) <-NULL
        }
        # 分组EDC-设置交互项分析 (direction: "1_2")
        {
          # 菌群分组转换为因子
          phy_edc_dat_temp$mp4_group <- factor(phy_edc_dat_temp$mp4_group)
          # EDC分组转换为因子
          phy_edc_dat_temp$edc_group <- factor(phy_edc_dat_temp$edc_group)
          
          cox_fit_interaction <- coxph(Surv(TIME, CENSOR) ~ edc_group + mp4_group + edc_group*mp4_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp)
          cox_fit_interaction_summary <- summary(cox_fit_interaction)
          
          cox_results_interaction_temp <- data.frame(cox_fit_interaction_summary[["coefficients"]])
          cox_results_interaction_temp$rowname <- row.names(cox_results_interaction_temp)
          cox_results_interaction_temp$exposure <- i
          cox_results_interaction_temp$mediator <- j
          cox_results_interaction_temp$outcome <- outcome
          cox_results_interaction_temp$method <- "cox"
          cox_results_interaction_temp$adjust <- "adj"
          cox_results_interaction_temp$sample <- sample_name
          cox_results_interaction_temp$n <- nrow(phy_edc_dat_temp)
          cox_results_interaction_temp$group <- "interaction"
          cox_results_interaction_temp$direction <- "1_2"
          row.names(cox_results_interaction_temp) <-NULL
        }
        results_all <- bind_rows(results_all, cox_results_low_temp, cox_results_middle_temp, cox_results_high_temp, cox_results_interaction_temp)
        
      }
    }
  }
}
# results_all$exp.lci <- ifelse(results_all$direction != "1_4", exp(results_all$coef - 1.96*results_all$se.coef.), results_all$exp.lci) # additive interaction结果包含OR和β，不转换
# results_all$exp.uci <- ifelse(results_all$direction != "1_4", exp(results_all$coef + 1.96*results_all$se.coef.), results_all$exp.uci) # additive interaction结果包含OR和β，不转换


openxlsx::write.xlsx(results_all,"results/cox/interaction/cox_results_edc_mp4_interaction_sensitivity_20260728.xlsx")
#### interaction分析 (EDC-outcome/菌分组) ####


#### 整理interaction分析结果 ----
results_all <- readxl::read_xlsx("results/cox/interaction/cox_results_edc_mp4_interaction_sensitivity_20260728.xlsx")
unique(results_all$mediator)
# 挑选暴露结局变量
results_short <- results_all[results_all$rowname %in% c("edc_group1",
                                                        "mp4_group1","mp4_group2",
                                                        "edc_group1:mp4_group1","edc_group1:mp4_group2") &
                               results_all$mediator %in% c(mp4_s_log10),]
results_short <- results_short[(results_short$rowname %in% c("edc_group1") & results_short$group != "interaction") | 
                                 results_short$rowname %in% c("edc_group1:mp4_group1","edc_group1:mp4_group2"),]
results_short$OUTCOME <- ifelse(results_short$outcome %in% c("cvd_incident_1021","cvd_incident_1421"), "CVD",
                                ifelse(results_short$outcome %in% c("ckd_incident_1014"), "CKD", "DM"))

results_short$exp_med_out <- paste0(results_short$exposure,"|",results_short$mediator,"|",results_short$outcome)



## 多层次验证 ##
# 第一层：EDC_group1和MP4_group1以及MP4_group2都有显著交互作用 (检验都为同向)
results_short_level1 <- results_short %>%
  group_by(exposure, mediator, outcome) %>%
  filter(
    any(rowname == "edc_group1:mp4_group1" & `Pr...z..` < 0.05) &
      any(rowname == "edc_group1:mp4_group2" & `Pr...z..` < 0.05)
  ) %>%
  ungroup()
results_short_level1 <- results_short_level1[results_short_level1$rowname == "edc_group1:mp4_group2",]
results_short_level1$level <- 1

# 第二层：EDC_group1和MP4_group1以及MP4_group2其中一个有显著交互作用
results_short_level2 <- results_short %>%
  group_by(exposure, mediator, outcome) %>%
  filter(
    (any(rowname == "edc_group1:mp4_group1" & `Pr...z..` < 0.05) &
      any(rowname == "edc_group1:mp4_group2" & `Pr...z..` >= 0.05)) | 
      (any(rowname == "edc_group1:mp4_group1" & `Pr...z..` >= 0.05) &
         any(rowname == "edc_group1:mp4_group2" & `Pr...z..` < 0.05))
  ) %>%
  ungroup()
results_short_level2 <- results_short_level2[results_short_level2$rowname %in% c("edc_group1:mp4_group1","edc_group1:mp4_group2"),]
results_short_level2 <- results_short_level2[results_short_level2$Pr...z.. < 0.05,]
results_short_level2$level <- 2

# 第三层：EDC_group1和MP4_group1以及MP4_group2都没有显著交互作用 (同正或同负)
results_short_level3 <- results_short[!results_short$exp_med_out %in% c(unique(results_short_level1$exp_med_out),unique(results_short_level2$exp_med_out)),] %>%
  group_by(exposure, mediator, outcome) %>%
  # 判断这两行 coef 是否同号
  filter(
    {
      coef1 <- coef[rowname == "edc_group1:mp4_group1"]
      coef2 <- coef[rowname == "edc_group1:mp4_group2"]
      coef1 * coef2 > 0  # 同号则乘积 > 0
    }
  ) %>%
  ungroup()
results_short_level3 <- results_short_level3[results_short_level3$rowname == "edc_group1:mp4_group2",]
results_short_level3$level <- 3

# 第四层：EDC_group1和MP4_group1以及MP4_group2都没有显著交互作用 (符号不同)
results_short_level4 <- results_short[!results_short$exp_med_out %in% c(unique(results_short_level1$exp_med_out),unique(results_short_level2$exp_med_out)),] %>%
  group_by(exposure, mediator, outcome) %>%
  # 判断这两行 coef 是否同号
  filter(
    {
      coef1 <- coef[rowname == "edc_group1:mp4_group1"]
      coef2 <- coef[rowname == "edc_group1:mp4_group2"]
      coef1 * coef2 < 0  # 不同号则乘积 < 0
    }
  ) %>%
  ungroup()
results_short_level4 <- results_short_level4[results_short_level4$rowname %in% c("edc_group1:mp4_group1","edc_group1:mp4_group2"),]
results_short_level4$level <- 4

results_short_level_all <- rbind(results_short_level1, results_short_level2, results_short_level3)
results_short_level_all <- results_short_level_all[,c("coef","Pr...z..","rowname","exposure","mediator","outcome","exp_med_out","level")]

openxlsx::write.xlsx(results_short_level_all,"results/cox/interaction/multi_interaction_sig_sensitivity_for_validate_20260728.xlsx")
#### 整理interaction分析结果 ####


#### 比较"EDC/菌二分类结果"和"EDC二分类/菌三分类结果" ----
results_both_interaction_sig_keep <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig_20260728.xlsx")
# 因为二分类分析中为了保证菌群丰度对于结局是正相关，部分菌群调整了方向以高丰度作为reference，因此这部分菌的相乘交互相乘交互效应值需要调整方向
results_both_interaction_sig_keep1 <- results_both_interaction_sig_keep[results_both_interaction_sig_keep$rowname %in% c("EDC high:GM low"),]
results_both_interaction_sig_keep1$multiplicative_scale_coef <- log(results_both_interaction_sig_keep1$exp.coef.)
results_both_interaction_sig_keep1$multiplicative_scale_coef_trans <- -results_both_interaction_sig_keep1$multiplicative_scale_coef
# 二分类中未调整reference的菌，相乘交互效应值不需要调整方向
results_both_interaction_sig_keep2 <- results_both_interaction_sig_keep[results_both_interaction_sig_keep$rowname %in% c("EDC high:GM high"),]
results_both_interaction_sig_keep2$multiplicative_scale_coef <- log(results_both_interaction_sig_keep2$exp.coef.)
results_both_interaction_sig_keep2$multiplicative_scale_coef_trans <- results_both_interaction_sig_keep2$multiplicative_scale_coef
# 合并两个方向的相乘交互结果
results_both_interaction_sig_keep <- rbind(results_both_interaction_sig_keep1, results_both_interaction_sig_keep2)
results_both_interaction_sig_keep$outcome <- gsub("_no_self_report","",results_both_interaction_sig_keep$outcome)
results_both_interaction_sig_keep$keep_exp_med_out <- gsub("_no_self_report","",results_both_interaction_sig_keep$keep_exp_med_out)


results_short_level_all <- readxl::read_xlsx("results/cox/interaction/multi_interaction_sig_sensitivity_for_validate_20260728.xlsx")
results_short_level_all <- results_short_level_all[,c("coef","exp_med_out","level")]


dat_validate <- inner_join(results_both_interaction_sig_keep, results_short_level_all, by=c("keep_exp_med_out"="exp_med_out")) %>%
  left_join(dat_mp4_detected_rate, by=c("mediator"="species"))
dat_validate <- dat_validate[,c("exposure","mediator","outcome","keep_exp_med_out","multiplicative_scale_coef_trans","coef","level","detect_rate")]
unique(dat_validate$mediator)

# dat_validate <- dat_validate[dat_validate$mediator %in% mp4_s_log10_short,] # 只保留了有菌属名称的结果
# unique(dat_validate$mediator)

dat_validate_same_direction <- dat_validate[dat_validate$multiplicative_scale_coef_trans * dat_validate$coef > 0,]
dat_validate_same_direction_sig <- dat_validate_same_direction[dat_validate_same_direction$level %in% c(1,2),]

openxlsx::write.xlsx(dat_validate_same_direction_sig,"results/cox/interaction/both_interaction_sig_validated_resutls_20260728.xlsx")
#### 比较"EDC/菌二分类结果"和"EDC二分类/菌三分类结果" ----
