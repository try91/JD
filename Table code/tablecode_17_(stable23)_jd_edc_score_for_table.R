library(data.table)
library(dplyr)
library(tidyr)
### EDC log10转换 ###
### 使用 MPA4 和 MPA3 (genus 和 species) ###

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

#### 变量整理 ----
# 菌群2014菌群 (分类和连续)
# 丰度>0.0001, 出现率>10%的微生物 (物种和属)
mp4_s_names <- read.table("jiading/sourceDataTaxon/mpa4/species_names_mp4_10%.txt")
mp4_s_names <- mp4_s_names[,1]
mp4_g_names <- read.table("jiading/sourceDataTaxon/mpa4/genus_names_mp4_10%.txt")
mp4_g_names <- mp4_g_names[,1]
mp3_s_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_10%.txt")
mp3_s_names <- mp3_s_names[,1]
mp3_g_names <- read.table("jiading/sourceDataTaxon/mpa3/genus_names_mp3_10%.txt")
mp3_g_names <- mp3_g_names[,1]
# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_zero <- paste0(mp4_g_names,"_zero") # 菌群MP4填补0值丰度 (属层面)

mp3_s_bin <- paste0(mp3_s_names,"_bin") # 菌群MP3出现与否的分类变量 (物种层面)
mp3_s_log10 <- paste0(mp3_s_names,"_log10") # 菌群MP3丰度的log10转换 (物种层面)
mp3_s_zero <- paste0(mp3_s_names,"_zero") # 菌群MP3填补0值丰度 (物种层面)

mp3_g_bin <- paste0(mp3_g_names,"_bin") # 菌群MP3出现与否的分类变量 (属层面)
mp3_g_log10 <- paste0(mp3_g_names,"_log10") # 菌群MP3丰度的log10转换 (属层面)
mp3_g_zero <- paste0(mp3_g_names,"_zero") # 菌群MP3填补0值丰度 (属层面)

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
phy_out_cat <- c("dm_f","ckd_f","cvd_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014新发case，横断面数据)
phy_traits_cat <- c("dm_f","dm_b","ckd_f","ckd_b","cvd_f","cvd_b",
                    "ob_f","ob_b","abob_f","abob_b","ir_f","ir_b","dyslip_f","dyslip_b",
                    "mets_f","mets_b","nafld_f","nafld_b","hua_f","hua_b","hpt_f","hpt_b","as_imt_f","as_imt_b")

# 2014、2010表型 (连续)
phy_traits_cont <- c("bmi_f","bmi_b","height_f","height_b","weight_f","weight_b","whr_f","whr_b","wc_f","wc_b","hc_f","hc_b",
                     "tg_f","tg_b","ldl_f","ldl_b","hdl_f","hdl_b","chol_f","chol_b","apoa_f","apoa_b","apob_f","apob_b","nonhdl_f","nonhdl_b",
                     "alt_f","alt_b","ast_f","ast_b","ggt_f","ggt_b","bia_f","bia_b",
                     "egfr_f","egfr_b","acr_f","acr_b","scr_f","scr_b","ua_f","ua_b",
                     "glu0_f","glu0_b","glu120_f","glu120_b","vhba1c_f","vhba1c_b",
                     "ins0_f","ins0_b","ins120_f","ins120_b","homair_f","homair_b","homab_f","homab_b",
                     "sbp_f","sbp_b","dbp_f","dbp_b","pr_f","pr_b",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "hgb_f","hgb_b","plt_f","plt_b","eos_f","lym_f","mon_f","neu_f",
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_f","wbc_b","crp_f")
                     
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

#### 读取EDC INDEX traits变量名 ----
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")

edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### 读取EDC INDEX traits变量名 ####


#### 读取嘉定人群信息 ----
phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))

# 总人群
phy_edc_dat1 <- phy_edc_temp_list[["phy_edc_temp"]] # 提取subgroup
phy_edc_dat1$sex <- phy_edc_dat1$sex_b

# # 菌群人群
# phy_edc_dat1 <- phy_edc_temp_list[["phy_edc_temp0_3"]] # 提取subgroup
# phy_edc_dat1$sex <- phy_edc_dat1$sex_b
#### 读取嘉定人群信息 ####

#### 读取EDC INDEX数据 ----
edc_count_b <- read.table("jiading/sourceDataEDCs/index/edc_count_2010_20260313.txt", header = TRUE)
edc_score_b <- read.table("jiading/sourceDataEDCs/index/edc_score_2010_20260313.txt", header = TRUE)

edc_count_f <- read.table("jiading/sourceDataEDCs/index/edc_count_2014_20260313.txt", header = TRUE)
edc_score_f <- read.table("jiading/sourceDataEDCs/index/edc_score_2014_20260313.txt", header = TRUE)
#### 读取EDC INDEX数据 ####

phy_edc_dat1 <- left_join(phy_edc_dat1,edc_count_b,by="ID") %>%
  left_join(edc_score_b,by="ID") %>%
  left_join(edc_count_f,by="ID") %>%
  left_join(edc_score_f,by="ID")


#### JD 连续变量 (mean(SD) or median(interquartile range)) ----
# 总人群 EDC SCORE #
{
  dat_temp <- phy_edc_dat1
  
  cols <- edc_index_b_keep
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列37个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- dat_temp[, ..keep_cols]
  # dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：定义统计函数（一次性计算所需指标）（所有连续变量的数量&Mean(SD)&Median） ==========
  stat_fun <- function(x) {
    x_clean <- na.omit(x)  # 剔除NA值
    tibble(
      n = length(x_clean),
      mean = round(mean(x_clean), 2),
      sd = round(sd(x_clean), 2),
      median = round(median(x_clean), 2),
      q1 = round(quantile(x_clean, 0.25), 2),
      q3 = round(quantile(x_clean, 0.75), 2),
      min = round(min(x_clean), 2),
      max = round(max(x_clean), 2)
    )
  }
  
  # ========== 步骤3：全局统计 ==========
  global_stats <- dat_temp %>%
    select(all_of(cols)) %>%
    summarise(across(everything(), stat_fun)) %>%
    pivot_longer(everything(), names_to = "连续变量名", values_to = "全局统计") %>%
    unnest(全局统计)
    # rename_with(~paste0("全局_", .), -连续变量名)
  
  # ========== 步骤4：分性别统计 ==========
  gender_stats <- dat_temp %>%
    group_by(!!sym(gender_var)) %>%
    select(all_of(cols), !!sym(gender_var)) %>%
    summarise(across(all_of(cols), stat_fun)) %>%
    pivot_longer(-!!sym(gender_var), names_to = "连续变量名", values_to = "性别统计") %>%
    unnest(性别统计) %>%
    pivot_wider(
      names_from = !!sym(gender_var),
      values_from = c(n, mean, sd, median, q1, q3, min, max),
      names_sep = "_"
    )
  
  # ========== 步骤5：整合全局统计 + 分性别统计 ==========
  final_stats_b <- global_stats %>%
    left_join(gender_stats, by = "连续变量名") %>%
    arrange(连续变量名)
  
  final_stats_b <- final_stats_b[,c("连续变量名", 
                                    "n","mean","sd","median","q1","q3","min","max", 
                                    "n_0","mean_0","sd_0","median_0","q1_0","q3_0","min_0","max_0", 
                                    "n_1","mean_1","sd_1","median_1","q1_1","q3_1","min_1","max_1")]
  final_stats_b$mean <- paste0(final_stats_b$mean," (",final_stats_b$sd,")")
  final_stats_b$mean_0 <- paste0(final_stats_b$mean_0," (",final_stats_b$sd_0,")")
  final_stats_b$mean_1 <- paste0(final_stats_b$mean_1," (",final_stats_b$sd_1,")")
  
  final_stats_b$median <- paste0(final_stats_b$median," (",final_stats_b$q1,"-",final_stats_b$q3,")")
  final_stats_b$median_0 <- paste0(final_stats_b$median_0," (",final_stats_b$q1_0,"-",final_stats_b$q3_0,")")
  final_stats_b$median_1 <- paste0(final_stats_b$median_1," (",final_stats_b$q1_1,"-",final_stats_b$q3_1,")")
  
  final_stats_b$min_max <- paste0("[",final_stats_b$min,", ",final_stats_b$max,"]")
  final_stats_b$min_max_0 <- paste0("[",final_stats_b$min_0,", ",final_stats_b$max_0,"]")
  final_stats_b$min_max_1 <- paste0("[",final_stats_b$min_1,", ",final_stats_b$max_1,"]")
  
  final_stats_b$连续变量名 <- factor(final_stats_b$连续变量名, levels = cols)
  final_stats_b1 <- final_stats_b %>%
    arrange(连续变量名)
  final_stats_b1 <- final_stats_b1[,c("连续变量名", "n","mean","median","min_max", 
                                      "n_0","mean_0","median_0","min_max_0", 
                                      "n_1","mean_1","median_1","min_max_1")]
}

# 菌群人群 EDC SCORE #
{
  dat_temp <- phy_edc_dat1[phy_edc_dat1$flag_mpa4_14 == 1,]
  
  cols <- edc_index_f_keep
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列37个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- dat_temp[, ..keep_cols]
  # dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：定义统计函数（一次性计算所需指标）（所有连续变量的数量&Mean(SD)&Median） ==========
  stat_fun <- function(x) {
    x_clean <- na.omit(x)  # 剔除NA值
    tibble(
      n = length(x_clean),
      mean = round(mean(x_clean), 2),
      sd = round(sd(x_clean), 2),
      median = round(median(x_clean), 2),
      q1 = round(quantile(x_clean, 0.25), 2),
      q3 = round(quantile(x_clean, 0.75), 2),
      min = round(min(x_clean), 2),
      max = round(max(x_clean), 2)
    )
  }
  
  # ========== 步骤3：全局统计 ==========
  global_stats <- dat_temp %>%
    select(all_of(cols)) %>%
    summarise(across(everything(), stat_fun)) %>%
    pivot_longer(everything(), names_to = "连续变量名", values_to = "全局统计") %>%
    unnest(全局统计)
  # rename_with(~paste0("全局_", .), -连续变量名)
  
  # ========== 步骤4：分性别统计 ==========
  gender_stats <- dat_temp %>%
    group_by(!!sym(gender_var)) %>%
    select(all_of(cols), !!sym(gender_var)) %>%
    summarise(across(all_of(cols), stat_fun)) %>%
    pivot_longer(-!!sym(gender_var), names_to = "连续变量名", values_to = "性别统计") %>%
    unnest(性别统计) %>%
    pivot_wider(
      names_from = !!sym(gender_var),
      values_from = c(n, mean, sd, median, q1, q3, min, max),
      names_sep = "_"
    )
  
  # ========== 步骤5：整合全局统计 + 分性别统计 ==========
  final_stats_b <- global_stats %>%
    left_join(gender_stats, by = "连续变量名") %>%
    arrange(连续变量名)
  
  final_stats_b <- final_stats_b[,c("连续变量名", 
                                    "n","mean","sd","median","q1","q3","min","max", 
                                    "n_0","mean_0","sd_0","median_0","q1_0","q3_0","min_0","max_0", 
                                    "n_1","mean_1","sd_1","median_1","q1_1","q3_1","min_1","max_1")]
  final_stats_b$mean <- paste0(final_stats_b$mean," (",final_stats_b$sd,")")
  final_stats_b$mean_0 <- paste0(final_stats_b$mean_0," (",final_stats_b$sd_0,")")
  final_stats_b$mean_1 <- paste0(final_stats_b$mean_1," (",final_stats_b$sd_1,")")
  
  final_stats_b$median <- paste0(final_stats_b$median," (",final_stats_b$q1,"-",final_stats_b$q3,")")
  final_stats_b$median_0 <- paste0(final_stats_b$median_0," (",final_stats_b$q1_0,"-",final_stats_b$q3_0,")")
  final_stats_b$median_1 <- paste0(final_stats_b$median_1," (",final_stats_b$q1_1,"-",final_stats_b$q3_1,")")
  
  final_stats_b$min_max <- paste0("[",final_stats_b$min,", ",final_stats_b$max,"]")
  final_stats_b$min_max_0 <- paste0("[",final_stats_b$min_0,", ",final_stats_b$max_0,"]")
  final_stats_b$min_max_1 <- paste0("[",final_stats_b$min_1,", ",final_stats_b$max_1,"]")
  
  final_stats_b$连续变量名 <- factor(final_stats_b$连续变量名, levels = cols)
  final_stats_b2 <- final_stats_b %>%
    arrange(连续变量名)
  final_stats_b2 <- final_stats_b2[,c("连续变量名", "n","mean","median","min_max", 
                                      "n_0","mean_0","median_0","min_max_0", 
                                      "n_1","mean_1","median_1","min_max_1")]
}

final_stats_b1 <- final_stats_b1[,c(1:5)]
colnames(final_stats_b1)[2:5] <- paste0(colnames(final_stats_b1)[2:5]," | JD_2010")
final_stats_b1$连续变量名 <- factor(final_stats_b1$连续变量名,
                               
                                    levels = c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b",
                                               "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),
                               
                                    labels = c("EDC Scoremedian (14 EDCs)","EDC Scoremedian (PFAS)","EDC Scoremedian (PAEs)","EDC Scoremedian (antimicrobials)","EDC Scoremedian (bisphenols)", 
                                               "EDC Scorequartile (14 EDCs)","EDC Scorequartile (PFAS)","EDC Scorequartile (PAEs)","EDC Scorequartile (antimicrobials)","EDC Scorequartile (bisphenols)"))
final_stats_b1 <- final_stats_b1 %>%
  arrange(连续变量名)

final_stats_b2 <- final_stats_b2[,c(1:5)]
colnames(final_stats_b2)[2:5] <- paste0(colnames(final_stats_b2)[2:5]," | JD_sub")
final_stats_b2$连续变量名 <- factor(final_stats_b2$连续变量名,
                               
                               levels = c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_tc_f","edc_count2_bp1_f",
                                          "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"),
                               
                               labels = c("EDC Scoremedian (14 EDCs)","EDC Scoremedian (PFAS)","EDC Scoremedian (PAEs)","EDC Scoremedian (antimicrobials)","EDC Scoremedian (bisphenols)", 
                                          "EDC Scorequartile (14 EDCs)","EDC Scorequartile (PFAS)","EDC Scorequartile (PAEs)","EDC Scorequartile (antimicrobials)","EDC Scorequartile (bisphenols)"))

final_stats_b_all <- left_join(final_stats_b1, final_stats_b2, by=c("连续变量名"))


colnames(final_stats_b_all) <- c("Analyte-based scores",
                                 "n","Mean (SD)","Median (quartile 1-quartile 3)","[Min., Max.]",
                                 "n","Mean (SD)","Median (quartile 1-quartile 3)","[Min., Max.]")
openxlsx::write.xlsx(final_stats_b_all, "tables/jd_edc_score_20260512.xlsx")
#### JD 连续变量 (mean(SD) or median(interquartile range)) ####
