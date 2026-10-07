library(data.table)
library(dplyr)
library(chemometrics)
library(ggplot2)
library(survival)
library(interactionR)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))

#### 构建函数把species_name转换为可以作图的标准化名称 ----
trans_mp4_s_names <- function(STRING){
  
  string_temp <- gsub("s__","",STRING)
  string_temp <- gsub("_log10","",string_temp)
  # # 删除所有"_" (一般情况下第一个"_"分隔属和种)
  # string_temp <- gsub("_"," ",string_temp)
  # 删除"_" (一般情况下第一个"_"分隔属和种)
  string_temp <- sub("_"," ",string_temp) # 使用 sub()替换第一个匹配项。
  # 删除特殊位置的"_"
  string_temp <- gsub("_bacterium"," bacterium",string_temp)
  string_temp <- gsub("bacterium_","bacterium ",string_temp)
  string_temp <- gsub("_unclassified"," unclassified",string_temp)
  string_temp <- gsub("unclassified_","unclassified ",string_temp)
  string_temp <- gsub("copri_","copri ",string_temp)
  string_temp <- gsub("Cibionibacter_","Cibionibacter ",string_temp)
  string_temp <- gsub("Marseille_","Marseille ",string_temp)
  string_temp <- gsub("oral_taxon_","oral taxon ",string_temp)
  # # 添加特殊位置的"_"
  # string_temp <- gsub(" Family","_Family",string_temp)
  # string_temp <- gsub(" phylum","_phylum",string_temp)
  # string_temp <- gsub("Candidatus ","Candidatus_",string_temp)
  # 处理sp
  string_temp <- gsub("sp_","sp. ",string_temp)
  string_temp <- gsub("_sp."," sp.",string_temp)
  return(string_temp)
}
#### 构建函数把species_name转换为可以作图的标准化名称 ####

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

#### 筛选菌-Outcome pairs ----
## 菌群-outcome数据 ##
cox_results <- readxl::read_xlsx("results/cox/cox_results_mp4_incident_20260728.xlsx")
cox_results <- cox_results[cox_results$adjust == "adj",]
cox_results <- cox_results[,c(1,3:13)]
colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
out_cox <- unique(cox_results$outcome) # 提取结局变量
exp_cox <- unique(cox_results$exposure) # 提取暴露变量

logistic_results <- readxl::read_xlsx("results/glm/logistic_results_mp4_incident_20260728.xlsx")
logistic_results <- logistic_results[logistic_results$adjust == "adj",]
colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
out_logistic <- unique(logistic_results$outcome) # 提取结局变量
exp_logistic <- unique(logistic_results$exposure) # 提取暴露变量

cox_logistic_results <- rbind(cox_results,logistic_results)
cox_logistic_results <- cox_logistic_results[cox_logistic_results$outcome %in% c("cvd_incident_1421","ckd_incident_1014","dm_incident_1014_no_self_report"),]
cox_logistic_results <- cox_logistic_results[cox_logistic_results$rowname %in% c(mp4_s_log10),]
# 以每个OUT表型为单位进行校正（以outcome为组，校正每个菌）
cox_logistic_results <- cox_logistic_results %>%
  group_by(outcome) %>%  # 按outcome分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
cox_logistic_short <- cox_logistic_results[,c("rowname","outcome","estimate","p","p_adj_bh")]
colnames(cox_logistic_short)[c(1:5)] <- c("exp_name","out_name","estimate_cox_logistic","p","p_adj_bh")

cox_logistic_short$out_name <- factor(cox_logistic_short$out_name, levels = c("cvd_incident_1421","ckd_incident_1014","dm_incident_1014_no_self_report"))
cox_logistic_short <- cox_logistic_short %>%
  arrange(exp_name,out_name)

# effect size 转换为 HR/OR
cox_logistic_short$exp_estimate <- exp(cox_logistic_short$estimate_cox_logistic)

# cox 和 logistic 显著的菌-Outcome pairs
pair_mp4_out_sig <- cox_logistic_short[cox_logistic_short$p < 0.05,] # (p显著)
## 分别提取三个结局结果 ##
pair_mp4_cvd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "cvd_incident_1421" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 26个CVD显著相关的菌种
pair_mp4_ckd_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "ckd_incident_1014" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 84个CKD显著相关的菌种
pair_mp4_dm_sig <- pair_mp4_out_sig[pair_mp4_out_sig$out_name == "dm_incident_1014_no_self_report" & pair_mp4_out_sig$exp_name %in% c(mp4_s_log10),] # 35个DM显著相关的菌种
#### 筛选菌-Outcome pairs ####


#### 二分类EDC和二分类GM分别对结局的关联 ----
pair_mp4_out_for_interaction <- rbind(pair_mp4_cvd_sig, pair_mp4_ckd_sig, pair_mp4_dm_sig)
pair_mp4_out_for_interaction$med_out <- paste0(pair_mp4_out_for_interaction$exp_name,"|",pair_mp4_out_for_interaction$out_name)

# interaction 分析
exp_var <- c(edc_index_f_keep,edc_traits6_log10) # 2014年样本中计算得到的EDC Scores和单个EDC
med_var <- unique(pair_mp4_out_for_interaction$exp_name)
out_var <- unique(pair_mp4_out_for_interaction$out_name)
med_out <- unique(pair_mp4_out_for_interaction$med_out) # 145对潜在的菌-Outcome pairs(菌-outcome显著,p<0.05)

results_binary_exp <- data.frame()
for (x in out_var) {
  for (i in exp_var) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," || EDC (",which(exp_var == i)," out of ",length(exp_var),"): ",i,", Outcome (",which(out_var == x)," out of ",length(out_var),"): ",x))
    for (j in med_var) {
      
      # x <- out_var[2]
      # i <- exp_var[1]
      # j <- med_var[1]
      if(!paste0(j,"|",x) %in% med_out){
        next
      }
      
      ## 分析样本选取 (由于各步骤间样本量要统一，因此在计算前首先选择各自的样本，排除结局缺失的项) ##
      if(x == "cvd_incident_1421"){
        outcome <- "cvd_incident_1421"
        cols <- c("cvd_incident_1421", "timecvd_1421", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else if(x == "ckd_incident_1014"){
        outcome <- "ckd_incident_1014"
        cols <- c("ckd_incident_1014", "timeckd_1014", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else{
        outcome <- "dm_incident_1014_no_self_report"
        cols <- c("dm_incident_1014_no_self_report", "timedm_1014", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      colnames(phy_edc_dat_temp)[c(1,2,3,4)] <- c("CENSOR","TIME","EDC","MP4")
      
      
      ## 根据菌群丰度进行二分类 ##
      mp4_s_median <- median(phy_edc_dat_temp$MP4)
      phy_edc_dat_temp$mp4_group <- ifelse(phy_edc_dat_temp$MP4 <= mp4_s_median, 0, 1) # 0为低菌群丰度组，1为高菌群丰度组
      phy_edc_dat_temp$mp4_group <- factor(phy_edc_dat_temp$mp4_group) # 菌群分组转换为因子
      # logistic 分析 #
      if(outcome == "cvd_incident_1421"){
        cox_fit_interaction <- coxph(Surv(TIME, CENSOR) ~ mp4_group + age_f + sex_b_rev + smk1_f + drk1_f + high_edu_b + paactive3_g_f + high_fruveg + med_all7, data = phy_edc_dat_temp)
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
        cox_results_interaction_temp$direction <- "2_1"
        row.names(cox_results_interaction_temp) <-NULL
        
        results_binary_exp <- bind_rows(results_binary_exp, cox_results_interaction_temp)
      }else{
        logistic_fit_interaction <- glm(CENSOR ~ mp4_group + age_f + sex_b_rev + smk1_f + drk1_f + high_edu_b + paactive3_g_f + high_fruveg + med_all7, family = binomial, data = phy_edc_dat_temp)
        logistic_fit_interaction_summary <- summary(logistic_fit_interaction)
        
        logistic_fit_interaction_temp <- data.frame(logistic_fit_interaction_summary[["coefficients"]])
        logistic_fit_interaction_temp$rowname <- row.names(logistic_fit_interaction_temp)
        logistic_fit_interaction_temp$exposure <- i
        logistic_fit_interaction_temp$mediator <- j
        logistic_fit_interaction_temp$outcome <- outcome
        logistic_fit_interaction_temp$method <- "logistic"
        logistic_fit_interaction_temp$adjust <- "adj"
        logistic_fit_interaction_temp$sample <- sample_name
        logistic_fit_interaction_temp$n <- nrow(phy_edc_dat_temp)
        logistic_fit_interaction_temp$group <- "interaction"
        logistic_fit_interaction_temp$direction <- "2_1"
        colnames(logistic_fit_interaction_temp)[c(1,2,3)] <- c("coef","se.coef.","z")
        row.names(logistic_fit_interaction_temp) <-NULL
        
        results_binary_exp <- bind_rows(results_binary_exp, logistic_fit_interaction_temp)
      }
      
      
      
      ## 分析样本选取 (由于各步骤间样本量要统一，因此在计算前首先选择各自的样本，排除结局缺失的项) ##
      if(x == "cvd_incident_1421"){
        outcome <- "cvd_incident_1021"
        cols <- c("cvd_incident_1021", "timecvd_1021", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else if(x == "ckd_incident_1014"){
        outcome <- "ckd_incident_1014"
        cols <- c("ckd_incident_1014", "timeckd_1014", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }else{
        outcome <- "dm_incident_1014"
        cols <- c("dm_incident_1014", "timedm_1014", i, j, "age_b","age_f","sex_b_rev","smk1_b","smk1_f","drk1_b","drk1_f","high_edu_b","paactive3_g_b","paactive3_g_f","high_fruveg","med_all7",
                  paste0(edc_traits8,"_detected"))
        
      }
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      colnames(phy_edc_dat_temp)[c(1,2,3,4)] <- c("CENSOR","TIME","EDC","MP4")
      
      ## 根据EDC浓度进行二分类 ##
      if(i %in% edc_traits8_log10){ # 检出率<50%, detected or not 分组
        edc_name <- gsub("_log10","",i)
        phy_edc_dat_temp$edc_group <- phy_edc_dat_temp[[paste0(edc_name,"_detected")]] # 0为EDC未检出组，1为EDC检出组
      }else{ # 检出率>50% or EDC index, 中位数分组
        edc_median <- median(phy_edc_dat_temp$EDC)
        phy_edc_dat_temp$edc_group <- ifelse(phy_edc_dat_temp$EDC <= edc_median, 0, 1) # 0为低EDC浓度组，1为高EDC浓度组
      }
      phy_edc_dat_temp$edc_group <- factor(phy_edc_dat_temp$edc_group) # EDC分组转换为因子
      # cox 分析 #
      {
        cox_fit_interaction <- coxph(Surv(TIME, CENSOR) ~ edc_group + age_b + sex_b_rev + smk1_b + drk1_b + high_edu_b + paactive3_g_b + high_fruveg, data = phy_edc_dat_temp)
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
        cox_results_interaction_temp$direction <- "2_2"
        row.names(cox_results_interaction_temp) <-NULL
        
        results_binary_exp <- bind_rows(results_binary_exp, cox_results_interaction_temp)
      }
      
    }
  }
}
results_binary_exp$exp.lci <- exp(results_binary_exp$coef - 1.96*results_binary_exp$se.coef.)
results_binary_exp$exp.uci <- exp(results_binary_exp$coef + 1.96*results_binary_exp$se.coef.)

results_binary_exp_short <- results_binary_exp[results_binary_exp$rowname %in% c("edc_group1","mp4_group1"),]
openxlsx::write.xlsx(results_binary_exp_short,"results/cox/interaction/cox_logistic_results_edc_mp4_binary_20260802.xlsx")
#### 二分类EDC和二分类GM分别对结局的关联 ####



# 读取stable中最终保留的17种交互组合（先运行“interaction_edc_mp4_outcome_for_table”）
keep_exp_med_out_17 <- read.table("results/cox/interaction/final_17_pairs_for_stable_20260802.txt", header = TRUE)
# 读取二分类EDC和二分类GM分别对结局的关联结果
results_binary_exp_short <- readxl::read_xlsx("results/cox/interaction/cox_logistic_results_edc_mp4_binary_20260802.xlsx")
table(results_binary_exp_short$outcome)

results_binary_exp_short$outcome2 <- ifelse(results_binary_exp_short$outcome == "cvd_incident_1421", "cvd_incident_1021", results_binary_exp_short$outcome)
results_binary_exp_short$outcome2 <- ifelse(results_binary_exp_short$outcome2 == "dm_incident_1014_no_self_report", "dm_incident_1014", results_binary_exp_short$outcome2)
results_binary_exp_short$keep_exp_med_out <- paste0(results_binary_exp_short$exposure,"|",results_binary_exp_short$mediator,"|",results_binary_exp_short$outcome2)
results_binary_exp_short <- results_binary_exp_short[results_binary_exp_short$keep_exp_med_out %in% keep_exp_med_out_17$x,]
results_binary_exp_short$exp.coef.dong <- exp(results_binary_exp_short$coef)
results_binary_exp_short$dif <- results_binary_exp_short$exp.coef.dong - results_binary_exp_short$exp.coef.
results_binary_exp_short$text <- sprintf("%.3f (%.3f, %.3f)", results_binary_exp_short$exp.coef.dong, results_binary_exp_short$exp.lci, results_binary_exp_short$exp.uci)


# 排序 #
results_binary_exp_short$sig_flag <- ifelse(results_binary_exp_short$Pr...z.. < 0.05, 1, 2)
results_binary_exp_short$exposure_group <- ifelse(results_binary_exp_short$exposure %in% c("edc_count2_edc14_f","edc_score_edc14_f"), 1,
                                                  ifelse(results_binary_exp_short$exposure %in% c("edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10), 2,
                                                         ifelse(results_binary_exp_short$exposure %in% c("edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10), 3,
                                                                ifelse(results_binary_exp_short$exposure %in% c("edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10), 4,
                                                                       ifelse(results_binary_exp_short$exposure %in% c("edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10), 5, 0)))))
results_binary_exp_short$exposure <- factor(results_binary_exp_short$exposure,
                                            levels = c("edc_count2_edc14_f","edc_score_edc14_f",
                                                       "edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10,
                                                       "edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10,
                                                       "edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10,
                                                       "edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10),
                                            labels = c("EDC Scoremedian (14 EDCs)", "EDC Scorequartile (14 EDCs)",
                                                       "EDC Scoremedian (PFAS)", "EDC Scorequartile (PFAS)", edc_traits2,
                                                       "EDC Scoremedian (PAEs)", "EDC Scorequartile (PAEs)", edc_traits3,
                                                       "EDC Scoremedian (antimicrobials)", "EDC Scorequartile (antimicrobials)", edc_traits5,
                                                       "EDC Scoremedian (bisphenols)", "EDC Scorequartile (bisphenols)", edc_traits4))
results_binary_exp_short$mediator <- trans_mp4_s_names(results_binary_exp_short$mediator)
results_binary_exp_short$outcome2 <- factor(results_binary_exp_short$outcome2, 
                                            levels = c("dm_incident_1014","ckd_incident_1014","cvd_incident_1021"),
                                            labels = c("Diabetes","CKD","CVD"))


# 提取EDC分组和结局关联
dat1 <- results_binary_exp_short[results_binary_exp_short$rowname == "edc_group1",]
dat1 <- dat1[,c("exposure","mediator","outcome2","text","Pr...z..")]
colnames(dat1)[c(4,5)] <- c("EDC (high vs. low) HR (95% CI)","P value")
# 提取GM分组和结局关联
dat2 <- results_binary_exp_short[results_binary_exp_short$rowname == "mp4_group1",]
dat2 <- dat2[,c("exposure","mediator","outcome2","text","Pr...z..")]
colnames(dat2)[c(4,5)] <- c("GM (high vs. low) OR/HR (95% CI)","P value")


#### 读取interaction stable主结果（保证binary结果顺序一致） ----
dat_all_stable <- readxl::read_xlsx(paste0("tables/interaction_results_20260802.xlsx"))  # 先运行“interaction_edc_mp4_outcome_for_table”
dat_all_stable <- dat_all_stable[,c("Analyte-based scores or individual analytes","Gut microbial species","Outcomes")]
#### 读取interaction stable主结果（保证binary结果顺序一致） ####

#### 读取clade name数据 ----
## 读取clade name ##
clade_name_mp4_s <- read.table("jiading/sourceDataTaxon/mpa4/JD.mp4.n4491_clade_name.txt", header = TRUE)
clade_name_mp4_s$species <- trans_mp4_s_names(clade_name_mp4_s$species)
clade_name_mp4_s <- clade_name_mp4_s[,c("phylum","species")]
table(clade_name_mp4_s$species)
table(clade_name_mp4_s$phylum)
#### 读取clade name数据 ####


dat_all_stable <- left_join(dat_all_stable, dat1, by=c("Analyte-based scores or individual analytes"="exposure","Gut microbial species"="mediator","Outcomes"="outcome2")) %>%
  left_join(dat2, by=c("Analyte-based scores or individual analytes"="exposure","Gut microbial species"="mediator","Outcomes"="outcome2")) %>%
  left_join(clade_name_mp4_s, by=c("Gut microbial species"="species"))
dat_all_stable$Phylum <- gsub("p__","",dat_all_stable$phylum)
colnames(dat_all_stable)
dat_all_stable <- dat_all_stable[,c("Analyte-based scores or individual analytes","Gut microbial species","Phylum","Outcomes",
                                    "EDC (high vs. low) HR (95% CI)","P value.x","GM (high vs. low) OR/HR (95% CI)","P value.y")]


openxlsx::write.xlsx(dat_all_stable, paste0("tables/interaction_results_binary_edc_mp4_outcome_20260802.xlsx"))
