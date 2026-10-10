library(data.table)
library(dplyr)
library(tidyr)
library(survey)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")

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
                "TCC","TCS",
                "BPA","BPS","BPF")
edc_traits2 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS")
edc_traits3 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP")
edc_traits3_q2 <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP") # 检出率>50%的PAE6
edc_traits3_q4 <- c("MEHP","MECPP","MEHHP","MEP") # 检出率>75%的PAE4
edc_traits4 <- c("BPA","BPS","BPF")
edc_traits4_q2 <- c("BPA") # 检出率>50%的BP1
edc_traits5 <- c("TCC","TCS")
edc_traits6 <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                 "TCC","TCS","BPA") # 检出率>50%
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
phy_out_cat <- c("dm_f","ckd_f","cvd_f")
phy_traits_cat <- c("dm_f","ckd_f","cvd_f","dm_b","ckd_b","cvd_b",
                    "ob_f","ob_b","abob_f","abob_b","ir_f","ir_b","dyslip_f","dyslip_b",
                    "mets_f","mets_b","nafld_f","nafld_b","hua_f","hua_b","hpt_f","hpt_b","as_imt_f","as_imt_b")
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


#### JD_2010 EDC 信息 ----
# 检测下限 #
jd_edc_low_limit_ngmL <- c(PFOA = 1, PFNA = 1, PFDA = 0.5, PFOS = 1, PFHxS = 0.2, 
                           MEHP = 1, MEHHP = 0.01, MEOHP = 0.02, MECPP = 0.1, MnBP = 0.5, MEP = 0.1, MBzP = 0.05, MCPP = 0.015, MiBP = 0.5, 
                           TCS = 0.2, TCC = 0.02,
                           BPA = 0.5, BPS = 0.1, BPF = 0.2)
dat_edc <- data.frame(edc_low_limit_ngmL = c(PFOA = 1, PFNA = 1, PFDA = 0.5, PFOS = 1, PFHxS = 0.2, 
                                             MEHP = 1, MEHHP = 0.01, MEOHP = 0.02, MECPP = 0.1, MnBP = 0.5, MEP = 0.1, MBzP = 0.05, MCPP = 0.015, MiBP = 0.5, 
                                             TCS = 0.2, TCC = 0.02,
                                             BPA = 0.5, BPS = 0.1, BPF = 0.2))
dat_edc$edc <- rownames(dat_edc)

# 检出率 #
for (i in edc_traits) {
  # i <- edc_traits[[1]]
  dat_edc[i,"dr"] <- sum(phy_edc_dat[[paste0(i,"_detected")]]) / nrow(phy_edc_dat)
}

# 检出样本的 median [q1, q3] #
for (i in edc_traits) {
  # i <- edc_traits[[1]]
  
  # 筛选出大于检测下限的数值，对他们进行四分位分类
  filtered_values <- phy_edc_dat[[i]][phy_edc_dat[[paste0(i,"_detected")]] > 0]
  # 四分位数转换
  q <- quantile(filtered_values, probs = c(0.25, 0.5, 0.75))
  
  dat_edc[i,"median"] <- q[2]
  dat_edc[i,"q1"] <- q[1]
  dat_edc[i,"q3"] <- q[3]
  
  dat_edc[i,"min"] <- min(filtered_values)
  dat_edc[i,"max"] <- max(filtered_values)
}

dat_edc$edc_low_limit_ngmL <- sprintf("%.2f", round(dat_edc$edc_low_limit_ngmL,2))
dat_edc$dr <- paste0(sprintf("%.2f", round(dat_edc$dr*100,2)),"%")
dat_edc$text <- paste0(sprintf("%.2f", round(dat_edc$median,2))," (",sprintf("%.2f", round(dat_edc$q1,2)),", ",sprintf("%.2f", round(dat_edc$q3,2)),")") 

dat_jd_edc_for_stable <- dat_edc[,c("edc","edc_low_limit_ngmL","dr","text")]

dat_jd_edc_for_stable$edc <- factor(dat_jd_edc_for_stable$edc, levels = edc_traits)
dat_jd_edc_for_stable <- dat_jd_edc_for_stable %>%
  arrange(edc)

colnames(dat_jd_edc_for_stable) <- c("EDCs or metabolites","Limit of quantification","Detection rate","Median (quartile 1-quartile 3)")
 
openxlsx::write.xlsx(dat_jd_edc_for_stable,"tables/(stable05)_edc_information_jd.xlsx")
#### JD_2010 EDC 信息 ####

#### NHANES EDC 信息 ----
# 读取NHANES信息
nhanes_dat <- read.csv("raw_data/nhanes_dat_2003-2018_for_analysis.csv")
# 13 EDCs (PFAS_5, PAEs_6, TC_1, BP_1); 2011-2014 (2c) #
{
  cols_edc <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                "TCS",
                "BPA")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat1 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]),]
  nhanes_dat1 <- nhanes_dat1[nhanes_dat1$release != "2015-2016",]
  table(nhanes_dat1$release)
}
# PAEs_6; 2003-2018 (8c) #
{
  cols_edc <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat5 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]),]
  table(nhanes_dat5$release)
}

# 检测下限 #
nhanes_edc_low_limit_ngmL <- c(PFOA = 0.1, PFNA = 0.1, PFDA = 0.1, PFOS = 0.1, PFHxS = 0.1, 
                               MEHP = 0.8, MEHHP = 0.4, MEOHP = 0.2, MECPP = 0.4, MnBP = 0.4, MEP = 1.2, MBzP = 0.3, MCPP = 0.4, MiBP = 0.8,
                               TCS = 1.7, TCC = 0.1,
                               BPA = 0.2, BPS = 0.1, BPF = 0.2) # 基于NHANES 2013-2014
dat_edc <- data.frame(edc_low_limit_ngmL = c(PFOA = 0.1, PFNA = 0.1, PFDA = 0.1, PFOS = 0.1, PFHxS = 0.1, 
                                             MEHP = 0.8, MEHHP = 0.4, MEOHP = 0.2, MECPP = 0.4, MnBP = 0.4, MEP = 1.2, MBzP = 0.3, MCPP = 0.4, MiBP = 0.8,
                                             TCS = 1.7, TCC = 0.1,
                                             BPA = 0.2, BPS = 0.1, BPF = 0.2)) # 基于NHANES 2013-2014
dat_edc$edc <- rownames(dat_edc)

# 检出率 #
edc_traits_1114 <- c("PFOA", "PFNA", "PFDA", "PFOS", "PFHxS",
                     "TCS", "BPA")
edc_traits_0318 <- c("MEHP", "MECPP", "MEHHP", "MEP", "MEOHP", "MiBP")
for (i in edc_traits_1114) {
  # i <- edc_traits_1114[[1]]
  
  nhanes_dat1[[paste0(i,"_detected")]] <- 1 - nhanes_dat1[[paste0(i,"_detected")]]
  dat_edc[i,"n"] <- length(nhanes_dat1[[i]])
  dat_edc[i,"dr"] <- sum(nhanes_dat1[[paste0(i,"_detected")]]) / length(nhanes_dat1[[i]])
  
  if(i %in% c("TCS", "BPA")){
    dat_edc[i,"sample_type"] <- "urine"
  }else{
    dat_edc[i,"sample_type"] <- "serum"
  }
  dat_edc[i,"cycle"] <- "2011-2014"
}
for (i in edc_traits_0318) {
  # i <- edc_traits_0318[[1]]
  
  nhanes_dat5[[paste0(i,"_detected")]] <- 1 - nhanes_dat5[[paste0(i,"_detected")]]
  dat_edc[i,"n"] <- length(nhanes_dat5[[i]])
  dat_edc[i,"dr"] <- sum(nhanes_dat5[[paste0(i,"_detected")]]) / length(nhanes_dat5[[i]])
  
  dat_edc[i,"sample_type"] <- "urine"
  dat_edc[i,"cycle"] <- "2003-2018"
}

# 检出样本的 median [q1, q3] #
for (i in edc_traits_1114) {
  # i <- edc_traits_1114[[1]]
  
  # 筛选出大于检测下限的数值，对他们进行四分位分类
  filtered_values <- nhanes_dat1[[i]][nhanes_dat1[[paste0(i,"_detected")]] == 1]
  # 四分位数转换
  q <- quantile(filtered_values, probs = c(0.25, 0.5, 0.75))
  
  dat_edc[i,"median"] <- q[2]
  dat_edc[i,"q1"] <- q[1]
  dat_edc[i,"q3"] <- q[3]
  
  dat_edc[i,"min"] <- min(filtered_values)
  dat_edc[i,"max"] <- max(filtered_values)
}
for (i in edc_traits_0318) {
  # i <- edc_traits_0318[[1]]
  
  # 筛选出大于检测下限的数值，对他们进行四分位分类
  filtered_values <- nhanes_dat5[[i]][nhanes_dat5[[paste0(i,"_detected")]] == 1]
  # 四分位数转换
  q <- quantile(filtered_values, probs = c(0.25, 0.5, 0.75))
  
  dat_edc[i,"median"] <- q[2]
  dat_edc[i,"q1"] <- q[1]
  dat_edc[i,"q3"] <- q[3]
  
  dat_edc[i,"min"] <- min(filtered_values)
  dat_edc[i,"max"] <- max(filtered_values)
}

dat_edc$edc_low_limit_ngmL <- sprintf("%.2f", round(dat_edc$edc_low_limit_ngmL,2))
dat_edc$dr <- paste0(sprintf("%.2f", round(dat_edc$dr*100,2)),"%")
dat_edc$text <- paste0(sprintf("%.2f", round(dat_edc$median,2))," (",sprintf("%.2f", round(dat_edc$q1,2)),", ",sprintf("%.2f", round(dat_edc$q3,2)),")") 

dat_nhanes_edc_for_stable <- dat_edc[,c("edc","n","edc_low_limit_ngmL","dr","text","sample_type","cycle")]

dat_nhanes_edc_for_stable$edc <- factor(dat_nhanes_edc_for_stable$edc, levels = edc_traits)
dat_nhanes_edc_for_stable <- dat_nhanes_edc_for_stable %>%
  arrange(edc)

colnames(dat_nhanes_edc_for_stable) <- c("EDCs or metabolites","Number of participants","Limit of quantification","Detection rate","Median (quartile 1-quartile 3)","Sample type","NHANES survey cycle")

openxlsx::write.xlsx(dat_nhanes_edc_for_stable,"tables/(stable11)_edc_information_nhanes.xlsx")
#### NHANES EDC 信息 ####
