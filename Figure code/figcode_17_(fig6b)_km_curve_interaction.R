library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(survival)
library(survminer)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  right_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp0_3"

### 分类协变量转换为因子 ###
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
### 分类协变量转换为因子 ###


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


#### 读取+处理interaction分析结果 ----
# "相乘交互作用"显著(Multiplicative scale P<0.05), 同时"相加交互作用"显著(以RERI显著为标准)的组
results_interaction_sig1 <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig.xlsx")
# "相乘交互作用"显著的组(Multiplicative scale P<0.05)
results_interaction_sig2 <- readxl::read_xlsx("results/cox/interaction/multi_interaction_sig.xlsx")
results_interaction_sig2 <- results_interaction_sig2[!results_interaction_sig2$keep_exp_med_out %in% results_interaction_sig1$keep_exp_med_out,] # 排除相乘和相加交互同时显著的结果
# "相加交互作用"显著的组(以RERI显著为标准)
results_interaction_sig3 <- readxl::read_xlsx("results/cox/interaction/add_interaction_sig.xlsx")
results_interaction_sig3 <- results_interaction_sig3[!results_interaction_sig3$keep_exp_med_out %in% results_interaction_sig1$keep_exp_med_out,] # 排除相乘和相加交互同时显著的结果

results_interaction_sig1 <- results_interaction_sig1[results_interaction_sig1$mediator %in% mp4_s_log10_short,]


results_interaction_sig1_dm <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "DM",]
results_interaction_sig1_dm$hr_95ci <- sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_dm$exp.coef., results_interaction_sig1_dm$exp.lci, results_interaction_sig1_dm$exp.uci)
results_interaction_sig1_dm$text <- paste0(results_interaction_sig1_dm$rowname,": ",results_interaction_sig1_dm$hr_95ci)
results_interaction_sig1_dm$text <- ifelse(results_interaction_sig1_dm$text == "SI: NA (NA, NA)", "", results_interaction_sig1_dm$text)
results_interaction_sig1_dm$p_text <- ifelse(results_interaction_sig1_dm$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_dm$Pr...z..,3)))
results_interaction_sig1_dm <- results_interaction_sig1_dm[results_interaction_sig1_dm$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_dm <- results_interaction_sig1_dm[!(results_interaction_sig1_dm$rowname %in% c("EDC high") & results_interaction_sig1_dm$group %in% c("interaction_recode")),]

results_interaction_sig1_ckd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CKD",]
results_interaction_sig1_ckd$hr_95ci <- sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_ckd$exp.coef., results_interaction_sig1_ckd$exp.lci, results_interaction_sig1_ckd$exp.uci)
results_interaction_sig1_ckd$text <- paste0(results_interaction_sig1_ckd$rowname,": ",results_interaction_sig1_ckd$hr_95ci)
results_interaction_sig1_ckd$text <- ifelse(results_interaction_sig1_ckd$text == "SI: NA (NA, NA)", "", results_interaction_sig1_ckd$text)
results_interaction_sig1_ckd$p_text <- ifelse(results_interaction_sig1_ckd$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_ckd$Pr...z..,3)))
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[!(results_interaction_sig1_ckd$rowname %in% c("EDC high") & results_interaction_sig1_ckd$group %in% c("interaction_recode")),]

results_interaction_sig1_cvd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CVD",]
results_interaction_sig1_cvd$hr_95ci <-sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_cvd$exp.coef., results_interaction_sig1_cvd$exp.lci, results_interaction_sig1_cvd$exp.uci)
results_interaction_sig1_cvd$text <- paste0(results_interaction_sig1_cvd$rowname,": ",results_interaction_sig1_cvd$hr_95ci)
results_interaction_sig1_cvd$text <- ifelse(results_interaction_sig1_cvd$text == "SI: NA (NA, NA)", NA, results_interaction_sig1_cvd$text)
results_interaction_sig1_cvd$p_text <- ifelse(results_interaction_sig1_cvd$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_cvd$Pr...z..,3)))
results_interaction_sig1_cvd <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_cvd <- results_interaction_sig1_cvd[!(results_interaction_sig1_cvd$rowname %in% c("EDC high") & results_interaction_sig1_cvd$group %in% c("interaction_recode")),]



results_interaction_sig1_dm$point_color <- ifelse(results_interaction_sig1_dm$rowname == "EDC high" & results_interaction_sig1_dm$group == "mp4_low", "cat1", 
                                                  ifelse(results_interaction_sig1_dm$rowname == "EDC high" & results_interaction_sig1_dm$group == "mp4_high", "cat2", 
                                                         ifelse(results_interaction_sig1_dm$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_dm$point_color <- factor(results_interaction_sig1_dm$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_dm <- results_interaction_sig1_dm %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_dm1 <- results_interaction_sig1_dm[results_interaction_sig1_dm$direction == "1_3",]
results_interaction_sig1_dm2 <- results_interaction_sig1_dm[results_interaction_sig1_dm$direction == "1_4", "text"]
colnames(results_interaction_sig1_dm2) <- "additive_text"
results_interaction_sig1_dm3 <- cbind(results_interaction_sig1_dm1,results_interaction_sig1_dm2)
results_interaction_sig1_dm3$hr_95ci <- paste0("HR: ",results_interaction_sig1_dm3$hr_95ci)


results_interaction_sig1_ckd$point_color <- ifelse(results_interaction_sig1_ckd$rowname == "EDC high" & results_interaction_sig1_ckd$group == "mp4_low", "cat1", 
                                                   ifelse(results_interaction_sig1_ckd$rowname == "EDC high" & results_interaction_sig1_ckd$group == "mp4_high", "cat2", 
                                                          ifelse(results_interaction_sig1_ckd$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_ckd$point_color <- factor(results_interaction_sig1_ckd$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_ckd <- results_interaction_sig1_ckd %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_ckd1 <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$direction == "1_3",] 
results_interaction_sig1_ckd2 <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$direction == "1_4", "text"]
colnames(results_interaction_sig1_ckd2) <- "additive_text"
results_interaction_sig1_ckd3 <- cbind(results_interaction_sig1_ckd1,results_interaction_sig1_ckd2)
results_interaction_sig1_ckd3$hr_95ci <- paste0("HR: ",results_interaction_sig1_ckd3$hr_95ci)


results_interaction_sig1_cvd$point_color <- ifelse(results_interaction_sig1_cvd$rowname == "EDC high" & results_interaction_sig1_cvd$group == "mp4_low", "cat1", 
                                                   ifelse(results_interaction_sig1_cvd$rowname == "EDC high" & results_interaction_sig1_cvd$group == "mp4_high", "cat2", 
                                                          ifelse(results_interaction_sig1_cvd$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_cvd$point_color <- factor(results_interaction_sig1_cvd$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_cvd <- results_interaction_sig1_cvd %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_cvd1 <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$direction == "1_3",] 
results_interaction_sig1_cvd2 <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$direction == "1_4", "text"]
colnames(results_interaction_sig1_cvd2) <- "additive_text"
results_interaction_sig1_cvd3 <- cbind(results_interaction_sig1_cvd1,results_interaction_sig1_cvd2)
results_interaction_sig1_cvd3$hr_95ci <- paste0("HR: ",results_interaction_sig1_cvd3$hr_95ci)
#### 读取+处理interaction分析结果 ####


############################################# Kaplan-Meier 生存曲线  #############################################
#### KM survival curve 数据处理 (CVD) ----
# 进一步筛选MP4 group中EDC high vs EDC low有显著的pairs
results_interaction_sig1 <- results_interaction_sig1 %>%
  group_by(exposure,mediator,OUTCOME) %>%
  slice(1:2) %>%
  ungroup
results_interaction_sig1_keep <- results_interaction_sig1[results_interaction_sig1$sig_flag =="sig",] %>%
  distinct(exposure,mediator,OUTCOME,.keep_all = TRUE)

select_exp_med <- results_interaction_sig1_keep
select_exp_med <- select_exp_med[select_exp_med$exposure %in% c(edc_index_f_keep, edc_traits6_log10),]
#### KM survival curve 数据处理 (CVD) ####

#### Kaplan-Meier survival curve ----
select_exp_med <- select_exp_med[select_exp_med$OUTCOME == "CVD",]

plot_list_cvd <- list()
plot_number_at_risk_cvd <- list()
n_case_subgroup <- data.frame()
for (i in 4) {
  # i <- 4
  ## 分析样本选取 (由于各步骤间样本量要统一，因此在计算前首先选择各自的样本，排除结局缺失的项) ##
  if(select_exp_med$OUTCOME[i] == "CVD"){
    outcome <- "Incident CVD (2010-2021)"
    cols <- c("cvd_incident_1021", "timecvd_1021", select_exp_med$exposure[i],select_exp_med$mediator[i])
  }else if(select_exp_med$OUTCOME[i] == "CKD"){
    outcome <- "Incident CKD (2010-2014)"
    cols <- c("ckd_incident_1014", "timeckd_1014", select_exp_med$exposure[i],select_exp_med$mediator[i])
  }else{
    outcome <- "Incident DM (2010-2014)"
    cols <- c("dm_incident_1014", "timedm_1014", select_exp_med$exposure[i],select_exp_med$mediator[i])
  }
  
  phy_edc_dat_temp <- phy_edc_dat[,cols]
  phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
  colnames(phy_edc_dat_temp)[c(1,2,3,4)] <- c("CENSOR","TIME","EDC","MP4")
  
  ## 根据EDC浓度进行二分类 ##
  edc_median <- median(phy_edc_dat_temp$EDC)
  phy_edc_dat_temp$edc_group <- ifelse(phy_edc_dat_temp$EDC <= edc_median, 0, 1) # 0为低EDC浓度组，1为高EDC浓度组
  phy_edc_dat_temp$edc_group <- factor(phy_edc_dat_temp$edc_group)
  
  ## 根据菌群丰度进行二分类 ##
  mp4_s_median <- median(phy_edc_dat_temp$MP4)
  phy_edc_dat_temp$mp4_group <- ifelse(phy_edc_dat_temp$MP4 <= mp4_s_median, 0, 1) # 0为低菌群丰度组，1为高菌群丰度组
  phy_edc_dat_temp$mp4_group <- factor(phy_edc_dat_temp$mp4_group)
  
  # table(phy_edc_dat_temp$edc_group, phy_edc_dat_temp$mp4_group)
  # 统计EDC group和MP4 group里有多少case
  dat_n_case_subgroup <- phy_edc_dat_temp %>%
    group_by(edc_group,mp4_group) %>%
    mutate(n_total = n(),
           n_case = sum(CENSOR == 1, na.rm = TRUE),  # CENSOR=1的数量
           n_case_0y = sum(TIME >= 0, na.rm = TRUE),  # 第0年，未发病人群
           n_case_1y = sum(TIME >= 1, na.rm = TRUE),  # 第1年，未发病人群
           n_case_2y = sum(TIME >= 2, na.rm = TRUE),  # 第2年，未发病人群
           n_case_3y = sum(TIME >= 3, na.rm = TRUE),  # 第3年，未发病人群
           n_case_4y = sum(TIME >= 4, na.rm = TRUE),  # 第4年，未发病人群
           n_case_5y = sum(TIME >= 5, na.rm = TRUE),  # 第5年，未发病人群
           n_case_6y = sum(TIME >= 6, na.rm = TRUE),  # 第6年，未发病人群
           n_case_7y = sum(TIME >= 7, na.rm = TRUE),  # 第7年，未发病人群
           n_case_8y = sum(TIME >= 8, na.rm = TRUE),  # 第8年，未发病人群
           n_case_9y = sum(TIME >= 9, na.rm = TRUE),  # 第9年，未发病人群
           n_case_10y = sum(TIME >= 10, na.rm = TRUE),  # 第10年，未发病人群
           n_case_11y = sum(TIME >= 11, na.rm = TRUE),  # 第11年，未发病人群
    ) %>% 
    ungroup() %>%
    distinct(edc_group,mp4_group,n_total,n_case,
             n_case_0y,n_case_1y,n_case_2y,n_case_3y,n_case_4y,n_case_5y,n_case_6y,n_case_7y,n_case_8y,n_case_9y,n_case_10y,n_case_11y) %>%
    arrange(mp4_group,desc(edc_group))
  dat_n_case_subgroup$exp <- select_exp_med$exposure[i]
  dat_n_case_subgroup$med <- select_exp_med$mediator[i]
  dat_n_case_subgroup$out <- select_exp_med$OUTCOME[i]
  n_case_subgroup <- rbind(n_case_subgroup,dat_n_case_subgroup)
  # 统计EDC group和MP4 group里有多少case
  
  phy_edc_dat_temp$group <- ifelse(phy_edc_dat_temp$mp4_group == 0 & phy_edc_dat_temp$edc_group == 0, "Low species abundance-Low EDC",
                                   ifelse(phy_edc_dat_temp$mp4_group == 0 & phy_edc_dat_temp$edc_group == 1, "Low species abundance-High EDC",
                                          ifelse(phy_edc_dat_temp$mp4_group == 1 & phy_edc_dat_temp$edc_group == 0, "High species abundance-Low EDC", "High species abundance-High EDC")))
  phy_edc_dat_temp$group <- factor(phy_edc_dat_temp$group, levels = c("Low species abundance-High EDC", "Low species abundance-Low EDC", "High species abundance-High EDC", "High species abundance-Low EDC"))
  
  # KM 曲线
  {
    fit <- survfit(Surv(TIME,CENSOR) ~ group, data = phy_edc_dat_temp)
    f1 <- ggsurvplot(fit,                     
                     data = phy_edc_dat_temp,
                     fun = "event",
                     conf.int = FALSE,         # If TRUE, plots confidence interval
                     palette = c("#b5be26", "#78e08f", "#5e2a2a", "#e55039"), 
                     xlim = c(0, max(phy_edc_dat_temp$TIME)),         
                     xlab = paste0("Time in years"),   
                     ggtheme = theme_light(),
                     ncensor.plot.height = 0.5,
                     size = 0.4,             # 设置生存曲线的粗细
                     censor.size = 2,        # 调整删失点的大小
                     censor.shape = "|"      # 调整删失点的形状，竖线通常更醒目
    )
    
    # # 单独提取Log-rank检验结果 #
    # phy_edc_dat1 <- phy_edc_dat_temp[phy_edc_dat_temp$group %in% c("Low species abundance-High EDC","Low species abundance-Low EDC"),]
    # surv_diff1 <- survdiff(Surv(TIME,CENSOR) ~ group, data = phy_edc_dat1) # 黄绿色 vs. 绿色
    # 
    # phy_edc_dat2 <- phy_edc_dat_temp[phy_edc_dat_temp$group %in% c("High species abundance-High EDC","High species abundance-Low EDC"),]
    # surv_diff2 <- survdiff(Surv(TIME,CENSOR) ~ group, data = phy_edc_dat2) # 深红 vs. 红色
    # 
    # phy_edc_dat3 <- phy_edc_dat_temp[phy_edc_dat_temp$group %in% c("Low species abundance-High EDC","High species abundance-High EDC"),]
    # surv_diff3 <- survdiff(Surv(TIME,CENSOR) ~ group, data = phy_edc_dat3) # 黄绿色 vs.深红
    # 
    # phy_edc_dat4 <- phy_edc_dat_temp[phy_edc_dat_temp$group %in% c("Low species abundance-High EDC","High species abundance-Low EDC"),]
    # surv_diff4 <- survdiff(Surv(TIME,CENSOR) ~ group, data = phy_edc_dat4) # 黄绿色 vs.红色
    
        f1$plot <- f1$plot + 
      
      labs(title = outcome,
           subtitle = paste0("GM: ", gsub("_log10","",select_exp_med$mediator[i]), "\n",
                             "EDC: ", gsub("_log10","",select_exp_med$exposure[i])),
           y="Cumulative incidence of CVD (%)") +
      
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, size = 0.6), # 面板区域边框
        axis.line = element_line(size = 0.4, color = "black"), # 统一调整坐标轴粗细
        
        plot.title = element_text(size = 14, colour = "black"),
        plot.subtitle = element_text(size = 13, colour = "black"),
        
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)),
        
        
        axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
        axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
        
        axis.text.x = element_text(size = 12, colour = "black"), # 调整x轴文字
        axis.text.y = element_text(size = 12, colour = "black"),
        
        legend.title = element_blank(),
        legend.text = element_text(size = 10, colour = "black", margin = margin(t = 1, r = 0, b = 1, l = 0, unit = "mm")), # 调整legend文本大小
        # 图例键大小
        legend.key.height = unit(2, "mm"),
        legend.key.width = unit(2, "mm")  
      ) + 
      
      guides(
        # 假设图例基于颜色 (color)，如果是线型 (linetype) 也需同样设置
        color = guide_legend(nrow = 4, byrow = TRUE)
      )+
      # 关键：expand 控制轴线与数据两端空白
      scale_x_continuous(limits = c(0,12),            # x轴范围：0到12
                         breaks = c(0:11),            # 刻度位置：0,3,6,9,12
                         expand = expansion(mult = c(0.02, 0.02))) +  # x 轴左右留白 2%
      scale_y_continuous(labels = function(y) y*100,
                         expand = expansion(mult = c(0.02, 0.25)))     # y 轴上下留白 2%
    
  }
  
  # Number at risk
  {
    dat_n_case_subgroup_long <- tidyr::gather(dat_n_case_subgroup, year, number_at_risk, n_case_0y:n_case_11y, na.rm=TRUE, factor_key=TRUE)
    dat_n_case_subgroup_long$group <- ifelse(dat_n_case_subgroup_long$mp4_group == 0 & dat_n_case_subgroup_long$edc_group == 0, "Low species abundance-Low EDC",
                                             ifelse(dat_n_case_subgroup_long$mp4_group == 0 & dat_n_case_subgroup_long$edc_group == 1, "Low species abundance-High EDC",
                                                    ifelse(dat_n_case_subgroup_long$mp4_group == 1 & dat_n_case_subgroup_long$edc_group == 0, "High species abundance-Low EDC", "High species abundance-High EDC")))
    dat_n_case_subgroup_long$group <- factor(dat_n_case_subgroup_long$group, levels = rev(c("Low species abundance-High EDC", "Low species abundance-Low EDC", "High species abundance-High EDC", "High species abundance-Low EDC")))
    
    
    f2 <- ggplot(dat_n_case_subgroup_long, aes(x = year, y = group)) +
      # 添加text标签
      geom_text(aes(label = number_at_risk),
                size = 2.8,
                show.legend = FALSE ) +  # 可选：避免文本标签在图例中产生重复项
      
      labs(x = "Time in years", y = "") +
      
      ggtitle(paste0("Number at risk | ",select_exp_med$exposure[i],"|",select_exp_med$mediator[i])) +
      
      # 设置主题
      theme_classic() +
      
      theme(
        plot.margin = margin(5, 5, 5, 5), # t,r,b,l
        
        axis.line = element_line(size = 0.4, color = "black"), # 统一调整坐标轴粗细
        
        plot.title = element_text(size = 14, colour = "black"),
        plot.subtitle = element_text(size = 13, colour = "black"),
        
        axis.title.y = element_text(size = 12, colour = "black", margin = margin(t = 0, r = 2, b = 0, l = 0)),
        axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 2, r = 0, b = 0, l = 0)),
        
        
        axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
        axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
        
        axis.text.x = element_text(size = 10, colour = "black"), # 调整x轴文字
        axis.text.y = element_text(size = 12, colour = "black")
      ) +
      # 关键：expand 控制轴线与数据两端空白
      scale_x_discrete(label = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
                       expand = expansion(mult = c(0.05, 0.05))) +  # x 轴左右留白 1%
      scale_y_discrete(expand = expansion(mult = c(0.15, 0.15)))    # y 轴上下留白 1%
    
  }
  
  if(select_exp_med$OUTCOME[i] == "CVD"){
    plot_list_cvd[[paste0(select_exp_med$exposure[i], "|", select_exp_med$mediator[i], "|", outcome)]] <- f1
    plot_number_at_risk_cvd[[paste0(select_exp_med$exposure[i], "|", select_exp_med$mediator[i], "|", outcome)]] <- f2
  }
}


combined_plot1 <- arrange_ggsurvplots(plot_list_cvd,
                                      print = FALSE,
                                      ncol = 1, nrow = 1) # 指定排列的列数和行数
ggsave(paste0("figures/main_figures/(fig6b)_km_curve_cvd.pdf"),
       combined_plot1, width = 8, height = 4, limitsize = FALSE)


ggsave(paste0("figures/main_figures/(fig6b)_km_curve_cvd_number_at_risk.pdf"),
       plot_number_at_risk_cvd[[1]], width = 7, height = 1.5, limitsize = FALSE)
#### Kaplan-Meier survival curve ####
