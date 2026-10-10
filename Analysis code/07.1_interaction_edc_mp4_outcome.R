library(data.table)
library(dplyr)
library(survival)
library(interactionR)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  right_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp0_3"
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


#### 筛选菌-Outcome pairs ----
## 菌群-outcome数据 ##
cox_results <- readxl::read_xlsx("results/cox/cox_results_mp4_incident.xlsx")
cox_results <- cox_results[cox_results$adjust == "adj",]
cox_results <- cox_results[,c(1,3:13)]
colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
out_cox <- unique(cox_results$outcome) # 提取结局变量
exp_cox <- unique(cox_results$exposure) # 提取暴露变量

logistic_results <- readxl::read_xlsx("results/glm/logistic_results_mp4_incident.xlsx")
logistic_results <- logistic_results[logistic_results$adjust == "adj",]
colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
out_logistic <- unique(logistic_results$outcome) # 提取结局变量
exp_logistic <- unique(logistic_results$exposure) # 提取暴露变量

cox_logistic_results <- rbind(cox_results,logistic_results)
cox_logistic_results <- cox_logistic_results[cox_logistic_results$outcome %in% c("cvd_incident_1421","ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"),]
cox_logistic_results <- cox_logistic_results[cox_logistic_results$rowname %in% c(mp4_s_log10),]
# 以每个OUT表型为单位进行校正（以outcome为组，校正每个菌）
cox_logistic_results <- cox_logistic_results %>%
  group_by(outcome) %>%  # 按outcome分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
cox_logistic_short <- cox_logistic_results[,c("rowname","outcome","estimate","p","p_adj_bh")]
colnames(cox_logistic_short)[c(1:5)] <- c("exp_name","out_name","estimate_cox_logistic","p","p_adj_bh")

cox_logistic_short$out_name <- factor(cox_logistic_short$out_name, levels = c("cvd_incident_1421","ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"))
cox_logistic_short <- cox_logistic_short %>%
  arrange(exp_name,out_name)

# effect size 转换为 HR/OR
cox_logistic_short$exp_estimate <- exp(cox_logistic_short$estimate_cox_logistic)

# cox 和 logistic 显著的菌-Outcome pairs
pair_mp4_out_sig <- cox_logistic_short[cox_logistic_short$p < 0.05,] # (p显著)
## 分别提取三个结局结果 ##
pair_mp4_cvd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "cvd_incident_1421" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 26个CVD显著相关的菌种
pair_mp4_ckd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "ckd_incident_1014_no_self_report" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 84个CKD显著相关的菌种
pair_mp4_dm_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "dm_incident_1014_no_self_report" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 35个DM显著相关的菌种
#### 筛选菌-Outcome pairs ####


#### interaction分析 (EDC-outcome/菌分组) ----
pair_mp4_out_for_interaction <- rbind(pair_mp4_cvd_sig, pair_mp4_ckd_sig, pair_mp4_dm_sig)
pair_mp4_out_for_interaction$med_out <- paste0(pair_mp4_out_for_interaction$exp_name,"|",pair_mp4_out_for_interaction$out_name)

# interaction 分析
exp_var <- c(edc_index_f_keep,edc_traits6_log10) # 2014年样本中计算得到的EDC Scores和单个EDC
med_var <- unique(pair_mp4_out_for_interaction$exp_name)
out_var <- unique(pair_mp4_out_for_interaction$out_name)
med_out <- unique(pair_mp4_out_for_interaction$med_out) # 145对潜在的菌-Outcome pairs(菌-outcome显著,p<0.05)

results_all <- data.frame()
for (x in out_var) {
  for (i in exp_var) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," || EDC (",which(exp_var == i)," out of ",length(exp_var),"): ",i,", Outcome (",which(out_var == x)," out of ",length(out_var),"): ",x))
    for (j in med_var) {
      
      # x <- out_var[1]
      # i <- exp_var[26]
      # j <- med_var[1]
      if(!paste0(j,"|",x) %in% med_out){
        next
      }
      
      ## 分析样本选取 (由于各步骤间样本量要统一，因此在计算前首先选择各自的样本，排除结局缺失的项) ##
      if(x == "cvd_incident_1421"){
        outcome <- "cvd_incident_1021"
        cols <- c("cvd_incident_1021", "timecvd_1021", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else if(x == "ckd_incident_1014_no_self_report"){
        outcome <- "ckd_incident_1014"
        cols <- c("ckd_incident_1014", "timeckd_1014", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else{
        outcome <- "dm_incident_1014_no_self_report"
        cols <- c("dm_incident_1014", "timedm_1014", i, j, "age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }
      phy_edc_dat_temp <- phy_edc_dat[,cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      colnames(phy_edc_dat_temp)[c(1,2,3,4)] <- c("CENSOR","TIME","EDC","MP4")
      
      ## 根据菌群丰度进行二分类 ##
      {
        ## 根据菌群丰度进行二分类 ##
        mp4_s_median <- median(phy_edc_dat_temp$MP4)
        phy_edc_dat_temp$mp4_group <- ifelse(phy_edc_dat_temp$MP4 <= mp4_s_median, 0, 1) # 0为低菌群丰度组，1为高菌群丰度组
        
        ## 根据EDC浓度进行二分类 ##
        if(i %in% edc_traits8_log10){ # 检出率<50%, detected or not 分组
          edc_name <- gsub("_log10","",i)
          phy_edc_dat_temp$edc_group <- phy_edc_dat_temp[[paste0(edc_name,"_detected")]] # 0为EDC未检出组，1为EDC检出组
        }else{ # 检出率>50% or EDC index, 中位数分组
          edc_median <- median(phy_edc_dat_temp$EDC)
          phy_edc_dat_temp$edc_group <- ifelse(phy_edc_dat_temp$EDC <= edc_median, 0, 1) # 0为低EDC浓度组，1为高EDC浓度组
        }
        
        # 连续EDC-低菌群丰度亚组分析 (direction: "1_1")
        {
          phy_edc_dat_temp_low <- phy_edc_dat_temp[phy_edc_dat_temp$MP4 <= mp4_s_median,]
          
          cox_fit_low <- coxph(Surv(TIME, CENSOR) ~ EDC + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_low)
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
          cox_results_low_temp$direction <- "1_1"
          row.names(cox_results_low_temp) <-NULL
        }
        # 连续EDC-高菌群丰度亚组分析 (direction: "1_1")
        {
          phy_edc_dat_temp_high <- phy_edc_dat_temp[phy_edc_dat_temp$MP4 > mp4_s_median,]
          
          cox_fit_high <- coxph(Surv(TIME, CENSOR) ~ EDC + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_high)
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
          cox_results_high_temp$direction <- "1_1"
          row.names(cox_results_high_temp) <-NULL
        }
        # 连续EDC-设置交互项分析
        {
          # 菌群分组转换为因子
          phy_edc_dat_temp$mp4_group <- factor(phy_edc_dat_temp$mp4_group)
          
          cox_fit_interaction <- coxph(Surv(TIME, CENSOR) ~ EDC + mp4_group + EDC*mp4_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp)
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
          cox_results_interaction_temp$direction <- "1_1"
          row.names(cox_results_interaction_temp) <-NULL
        }
        results_all <- bind_rows(results_all, cox_results_low_temp, cox_results_high_temp, cox_results_interaction_temp)
        
        # 分组EDC-低菌群丰度亚组分析 (direction: "1_2")
        {
          phy_edc_dat_temp_low <- phy_edc_dat_temp[phy_edc_dat_temp$MP4 <= mp4_s_median,]
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
        # 分组EDC-高菌群丰度亚组分析 (direction: "1_2")
        {
          phy_edc_dat_temp_high <- phy_edc_dat_temp[phy_edc_dat_temp$MP4 > mp4_s_median,]
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
        results_all <- bind_rows(results_all, cox_results_low_temp, cox_results_high_temp, cox_results_interaction_temp)
        
        # recode (edc_group or mp4_group)
        {
          # if one of the exposures is protective - such that the stratum with the lowest risk becomes the new reference category
          phy_edc_dat_temp$mp4_group <- as.numeric(as.character(phy_edc_dat_temp$mp4_group))
          phy_edc_dat_temp$mp4_group_recode <- 1 - phy_edc_dat_temp$mp4_group
          
          phy_edc_dat_temp$edc_group <- as.numeric(as.character(phy_edc_dat_temp$edc_group))
          phy_edc_dat_temp$edc_group_recode <- 1 - phy_edc_dat_temp$edc_group
          # if one of the exposures is protective - such that the stratum with the lowest risk becomes the new reference category
          
          for (a in c(1:4)) { # 尝试四种组合
            if(a == 1){
              EDC_test <- "edc_group"
              MP4_test <- "mp4_group"
              edc_group_recode_check <- "not"
              mp4_group_recode_check <- "not"
            }else if(a == 2){
              EDC_test <- "edc_group_recode"
              MP4_test <- "mp4_group"
              edc_group_recode_check <- "recode"
              mp4_group_recode_check <- "not"
            }else if(a == 3){
              EDC_test <- "edc_group"
              MP4_test <- "mp4_group_recode"
              edc_group_recode_check <- "not"
              mp4_group_recode_check <- "recode"
            }else if(a == 4){
              EDC_test <- "edc_group_recode"
              MP4_test <- "mp4_group_recode"
              edc_group_recode_check <- "recode"
              mp4_group_recode_check <- "recode"
            }
            
            cols <- c("CENSOR","TIME",EDC_test,MP4_test,"age_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7")
            phy_edc_dat_temp_recode <- phy_edc_dat_temp[,cols]
            colnames(phy_edc_dat_temp_recode)[c(3,4)] <- c("edc_group_recode","mp4_group_recode")
            
            # 菌群分组转换为因子
            phy_edc_dat_temp_recode$mp4_group_recode <- factor(phy_edc_dat_temp_recode$mp4_group_recode)
            # EDC分组转换为因子
            phy_edc_dat_temp_recode$edc_group_recode <- factor(phy_edc_dat_temp_recode$edc_group_recode)
            
            
            # recode分组EDC-低菌群丰度亚组分析 (direction: "1_3")
            {
              phy_edc_dat_temp_recode_low <- phy_edc_dat_temp_recode[phy_edc_dat_temp_recode$mp4_group_recode == 0,]
              
              recode_fit_low <- coxph(Surv(TIME, CENSOR) ~ edc_group_recode + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_recode_low)
              recode_fit_low_summary <- summary(recode_fit_low)
              
              recode_fit_low_temp <- data.frame(recode_fit_low_summary[["coefficients"]])
              recode_fit_low_temp$rowname <- row.names(recode_fit_low_temp)
              recode_fit_low_temp$exposure <- i
              recode_fit_low_temp$mediator <- j
              recode_fit_low_temp$outcome <- outcome
              recode_fit_low_temp$method <- "cox"
              recode_fit_low_temp$adjust <- "adj"
              recode_fit_low_temp$sample <- sample_name
              recode_fit_low_temp$n <- nrow(phy_edc_dat_temp_recode_low)
              recode_fit_low_temp$group <- "mp4_recode_low"
              recode_fit_low_temp$direction <- "1_3"
              recode_fit_low_temp$edc_group_recode_check <- edc_group_recode_check
              recode_fit_low_temp$mp4_group_recode_check <- mp4_group_recode_check
              row.names(recode_fit_low_temp) <-NULL
            }
            # recode分组EDC-高菌群丰度亚组分析 (direction: "1_3")
            {
              phy_edc_dat_temp_recode_high <- phy_edc_dat_temp_recode[phy_edc_dat_temp_recode$mp4_group_recode == 1,]
              
              recode_fit_high <- coxph(Surv(TIME, CENSOR) ~ edc_group_recode + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_recode_high)
              recode_fit_high_summary <- summary(recode_fit_high)
              
              recode_fit_high_temp <- data.frame(recode_fit_high_summary[["coefficients"]])
              recode_fit_high_temp$rowname <- row.names(recode_fit_high_temp)
              recode_fit_high_temp$exposure <- i
              recode_fit_high_temp$mediator <- j
              recode_fit_high_temp$outcome <- outcome
              recode_fit_high_temp$method <- "cox"
              recode_fit_high_temp$adjust <- "adj"
              recode_fit_high_temp$sample <- sample_name
              recode_fit_high_temp$n <- nrow(phy_edc_dat_temp_recode_high)
              recode_fit_high_temp$group <- "mp4_recode_high"
              recode_fit_high_temp$direction <- "1_3"
              recode_fit_high_temp$edc_group_recode_check <- edc_group_recode_check
              recode_fit_high_temp$mp4_group_recode_check <- mp4_group_recode_check
              row.names(recode_fit_high_temp) <-NULL
            }
            # recode分组EDC-设置交互项分析 (direction: "1_3")
            {
              recode_fit <- coxph(Surv(TIME, CENSOR) ~ edc_group_recode + mp4_group_recode + edc_group_recode*mp4_group_recode + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg + med_all7, data = phy_edc_dat_temp_recode)
              recode_fit_summary <- summary(recode_fit)
              
              recode_fit_temp <- data.frame(recode_fit_summary[["coefficients"]])
              recode_fit_temp$rowname <- row.names(recode_fit_temp)
              recode_fit_temp$exposure <- i
              recode_fit_temp$mediator <- j
              recode_fit_temp$outcome <- outcome
              recode_fit_temp$method <- "cox"
              recode_fit_temp$adjust <- "adj"
              recode_fit_temp$sample <- sample_name
              recode_fit_temp$n <- nrow(phy_edc_dat_temp_recode)
              recode_fit_temp$group <- "interaction_recode"
              recode_fit_temp$direction <- "1_3"
              recode_fit_temp$edc_group_recode_check <- edc_group_recode_check
              recode_fit_temp$mp4_group_recode_check <- mp4_group_recode_check
              row.names(recode_fit_temp) <-NULL
            }
            
            if(recode_fit_temp[recode_fit_temp$rowname == "edc_group_recode1","coef"] > 0 & recode_fit_temp[recode_fit_temp$rowname == "mp4_group_recode1","coef"] > 0){
              break # 如果筛选到edc和mp4都为危害因素的变量组合，则跳出循环
            }
            
          }
        }
        results_all <- bind_rows(results_all, recode_fit_low_temp, recode_fit_high_temp, recode_fit_temp)
        
        # additive interaction分析
        tryCatch(expr = {
          # additive interaction analysis (direction: "1_4")
          {
            aia_edc_mp4_group <- interactionR(recode_fit,
                                              exposure_names = c("edc_group_recode", "mp4_group_recode"),
                                              ci.type = "mover", ci.level = 0.95,
                                              em = FALSE, recode = FALSE)
            
            additive_interaction_temp <- aia_edc_mp4_group$dframe
            
            colnames(additive_interaction_temp) <- c("rowname","exp.coef.","exp.lci","exp.uci","Pr...z..")
            # additive_interaction_temp$rowname <- "edc_group_recode:mp4_group_recode"
            additive_interaction_temp$exposure <- i
            additive_interaction_temp$mediator <- j
            additive_interaction_temp$outcome <- outcome
            additive_interaction_temp$method <- "cox"
            additive_interaction_temp$adjust <- "adj"
            additive_interaction_temp$sample <- sample_name
            additive_interaction_temp$n <- nrow(phy_edc_dat_temp)
            additive_interaction_temp$group <- "additive_interaction"
            additive_interaction_temp$direction <- "1_4"
            additive_interaction_temp$edc_group_recode_check <- edc_group_recode_check
            additive_interaction_temp$mp4_group_recode_check <- mp4_group_recode_check
          }
          results_all <- bind_rows(results_all, additive_interaction_temp)
          # additive interaction analysis
        },
        error = function(e) {
          message("错误信息: ", e$message)
          return(NA)  # 返回默认值
        })
        
      }
    }
  }
}
results_all$exp.lci <- ifelse(results_all$direction != "1_4", exp(results_all$coef - 1.96*results_all$se.coef.), results_all$exp.lci) # additive interaction结果包含OR和β，不转换
results_all$exp.uci <- ifelse(results_all$direction != "1_4", exp(results_all$coef + 1.96*results_all$se.coef.), results_all$exp.uci) # additive interaction结果包含OR和β，不转换


openxlsx::write.xlsx(results_all,"results/cox/interaction/cox_results_edc_mp4_interaction_all.xlsx")
#### interaction分析 (EDC-outcome/菌分组) ####


#### 整理interaction分析结果 ----
results_all <- readxl::read_xlsx("results/cox/interaction/cox_results_edc_mp4_interaction_all.xlsx")

# 挑选暴露结局变量
results_short <- results_all[results_all$rowname %in% c("EDC","mp4_group1","EDC:mp4_group1",
                                                        "edc_group1","edc_group1:mp4_group1",
                                                        "edc_group_recode1","mp4_group_recode1","edc_group_recode1:mp4_group_recode1",
                                                        "Multiplicative scale","RERI","AP","SI",
                                                        "MP4","edc_group1","MP4:edc_group1",
                                                        "mp4_group1","mp4_group1:edc_group1") &
                               results_all$mediator %in% c(mp4_s_log10),]
results_short$OUTCOME <- ifelse(results_short$outcome %in% c("cvd_incident_1021"), "CVD",
                                ifelse(results_short$outcome %in% c("ckd_incident_1014"), "CKD", "DM"))

results_short_test1 <- results_short[results_short$OUTCOME == "DM" & results_short$mediator == "s__Flavonifractor_plautii_log10" & results_short$exposure == "edc_score_edc14_f",]
results_short_test2 <- results_short[results_short$OUTCOME == "CVD" & results_short$mediator == "s__Dysosmobacter_sp_BX15_log10" & results_short$exposure == "edc_score_pae6_f",]

# interactionR包P值结果不准确通过95%CI确定显著性
results_short$sig_flag <- "not_sig"
results_short$sig_flag <- ifelse((results_short$direction != "1_4") & (results_short$Pr...z.. < 0.05), "sig", results_short$sig_flag)
results_short$sig_flag <- ifelse((results_short$rowname == "Multiplicative scale") & ((results_short$exp.lci > 1 & results_short$exp.uci > 1) | (results_short$exp.lci < 1 & results_short$exp.uci < 1)), "sig", results_short$sig_flag)
results_short$sig_flag <- ifelse((results_short$rowname == "RERI") & ((results_short$exp.lci > 0 & results_short$exp.uci > 0) | (results_short$exp.lci < 0 & results_short$exp.uci < 0)), "sig", results_short$sig_flag)
results_short$sig_flag <- ifelse((results_short$rowname == "AP") & ((results_short$exp.lci > 0 & results_short$exp.uci > 0) | (results_short$exp.lci < 0 & results_short$exp.uci < 0)), "sig", results_short$sig_flag)
results_short$sig_flag <- ifelse((results_short$rowname == "SI") & ((results_short$exp.lci > 1 & results_short$exp.uci > 1) | (results_short$exp.lci < 1 & results_short$exp.uci < 1)), "sig", results_short$sig_flag)
results_short$sig_flag <- ifelse((results_short$rowname == "SI") & (is.na(results_short$exp.coef.)), "not", results_short$sig_flag)
# 以每个OUT表型为单位进行校正
results_short_for_fdr <- results_short[results_short$rowname %in% c("EDC:mp4_group1","edc_group1:mp4_group1","edc_group_recode1:mp4_group_recode1",
                                                                    "MP4:edc_group1","mp4_group1:edc_group1"),]
results_short_for_fdr <- results_short_for_fdr %>%
  group_by(OUTCOME, direction) %>%  # 按OUTCOME中同一个方向分组
  mutate(p_adj_bh = p.adjust(`Pr...z..`, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
results_short_for_fdr <- results_short_for_fdr[,c("rowname","exposure","mediator","outcome","direction","p_adj_bh")]
results_short <- left_join(results_short, results_short_for_fdr, by=c("rowname","exposure","mediator","outcome","direction"))

# edc_group项标准化命名 (direction 1_1 1_2)
{
  # EDC group: edc_group1代表高EDC group
  results_short$rowname <- ifelse((results_short$direction != "1_3"), gsub("edc_group1", "EDC high", results_short$rowname), results_short$rowname)
  
  # MP4 group: mp4_group1代表高MP4 group
  results_short$rowname <- ifelse((results_short$direction != "1_3"), gsub("mp4_group1", "GM high", results_short$rowname), results_short$rowname)
}
# recode项标准化命名 (direction 1_3)
{
  # recode后的EDC group: edc_group_recode1代表转换前的低EDC group
  results_short$rowname <- ifelse((results_short$edc_group_recode_check == "recode" & results_short$direction == "1_3"), gsub("edc_group_recode1", "EDC low", results_short$rowname), results_short$rowname)
  results_short$rowname <- ifelse((results_short$edc_group_recode_check == "not" & results_short$direction == "1_3"), gsub("edc_group_recode1", "EDC high", results_short$rowname), results_short$rowname)
  
  # recode后的MP4 group: mp4_group_recode1代表转换前的低MP4 group
  results_short$rowname <- ifelse((results_short$mp4_group_recode_check == "recode" & results_short$direction == "1_3"), gsub("mp4_group_recode1", "GM low", results_short$rowname), results_short$rowname)
  results_short$rowname <- ifelse((results_short$mp4_group_recode_check == "not" & results_short$direction == "1_3"), gsub("mp4_group_recode1", "GM high", results_short$rowname), results_short$rowname)
  
  
  # recode后的MP4 group: mp4_recode_low分组代表转换前的高MP4 group, mp4_recode_high分组代表转换前的低MP4 group
  results_short$group <- ifelse((results_short$mp4_group_recode_check == "recode" & results_short$direction == "1_3" & results_short$group == "mp4_recode_low"), "mp4_high", results_short$group)
  results_short$group <- ifelse((results_short$mp4_group_recode_check == "recode" & results_short$direction == "1_3" & results_short$group == "mp4_recode_high"), "mp4_low", results_short$group)
  results_short$group <- ifelse((results_short$mp4_group_recode_check == "not" & results_short$direction == "1_3" & results_short$group == "mp4_recode_low"), "mp4_low", results_short$group)
  results_short$group <- ifelse((results_short$mp4_group_recode_check == "not" & results_short$direction == "1_3" & results_short$group == "mp4_recode_high"), "mp4_high", results_short$group)
}


# 保留recode group结果，保留interaction分析结果
results_short <- results_short[results_short$direction %in% c("1_3","1_4"),]
# results_short <- results_short[!(results_short$direction %in% c("1_3") & results_short$group == "interaction_recode"),]
results_short <- results_short[,c("rowname","exposure","mediator","outcome","group","n","exp.coef.","exp.lci","exp.uci","Pr...z..","sig_flag","p_adj_bh","direction","OUTCOME","edc_group_recode_check","mp4_group_recode_check")]

results_short$exposure <- factor(results_short$exposure, levels = c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                                                                    "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f",
                                                                    edc_traits_log10))
results_short$mediator <- factor(results_short$mediator, levels = c(mp4_s_log10))
results_short$OUTCOME <- factor(results_short$OUTCOME, levels = c("CVD","CKD","DM"))


### 交互分析结果筛选流程 ###
## 第一步: 分别筛选相乘和相加交互显著的结果 ##
results_short$keep_exp_med_out <- paste0(results_short$exposure,"|",results_short$mediator,"|",results_short$outcome)
# 保留"相乘交互作用"显著的组(Multiplicative scale P<0.05)
results_multiplicative_interaction_sig <- results_short %>%
  group_by(exposure, mediator, OUTCOME) %>%
  # 直接判断：分组内是否存在 乘法尺度+显著 的行
  filter(any(rowname == "Multiplicative scale" & sig_flag == "sig")) %>%
  ungroup()
results_multiplicative_interaction_sig <- results_multiplicative_interaction_sig %>%
  arrange(OUTCOME,exposure,mediator,direction)

# 保留"相加交互作用"显著的组(以RERI显著为标准)
results_additive_interaction_sig <- results_short %>%
  group_by(exposure, mediator, OUTCOME) %>%
  # 直接判断：分组内是否存在 加法尺度+显著 的行
  filter(any(rowname == "RERI" & sig_flag == "sig")) %>%
  ungroup()
results_additive_interaction_sig <- results_additive_interaction_sig %>%
  arrange(OUTCOME,exposure,mediator,direction)

## 第二步: 根据外来物质有害理论，我们只保留EDC高对于结局有危害效应的结果 ##
# 筛选菌群高or低组中，EDC对于结局至少有一个显著危害效应的组
results_multiplicative_interaction_sig_keep <- results_multiplicative_interaction_sig %>%
  group_by(exposure, mediator, OUTCOME) %>%
  filter(any(rowname == "EDC high")) %>%
  # 给每个大组打标记：组内是否在 direction=1_3 下满足条件
  mutate(
    keep_group = any(
      (direction == "1_3" & group != "interaction_recode") & sig_flag == "sig" & exp.coef. > 1,
      na.rm = TRUE
    )
  ) %>%
  # 只保留符合条件的大组（所有direction都保留）
  filter(keep_group == TRUE) %>%
  ungroup() %>%
  # 清理辅助列，取消分组
  dplyr::select(!keep_group)
# 筛选有相乘交互的组
uniq_dat1 <- unique(results_multiplicative_interaction_sig_keep$keep_exp_med_out)


# results_multiplicative_interaction_sig_notkeep <- results_multiplicative_interaction_sig %>%
#   group_by(exposure, mediator, OUTCOME) %>%
#   filter(any(rowname == "EDC low")) %>%
#   ungroup()
  
 
results_additive_interaction_sig_keep <- results_additive_interaction_sig %>%
  group_by(exposure, mediator, OUTCOME) %>%
  filter(any(rowname == "EDC high")) %>%
  # 给每个大组打标记：组内是否在 direction=1_3 下满足条件
  mutate(
    keep_group = any(
      (direction == "1_3" & group != "interaction_recode") & sig_flag == "sig" & exp.coef. > 1,
      na.rm = TRUE
    )
  ) %>%
  # 只保留符合条件的大组（所有direction都保留）
  filter(keep_group == TRUE) %>%
  ungroup() %>%
  # 清理辅助列，取消分组
  dplyr::select(!keep_group)
# 筛选有相加交互的组
uniq_dat2 <- unique(results_additive_interaction_sig_keep$keep_exp_med_out)


uniq_dat3 <- uniq_dat1[uniq_dat1 %in% uniq_dat2]
results_both_interaction_sig_keep <- results_multiplicative_interaction_sig[results_multiplicative_interaction_sig$keep_exp_med_out %in% uniq_dat3,]


openxlsx::write.xlsx(results_both_interaction_sig_keep,"results/cox/interaction/both_interaction_sig.xlsx")
openxlsx::write.xlsx(results_multiplicative_interaction_sig_keep,"results/cox/interaction/multi_interaction_sig.xlsx")
openxlsx::write.xlsx(results_additive_interaction_sig_keep,"results/cox/interaction/add_interaction_sig.xlsx")
#### 整理interaction分析结果 ####
