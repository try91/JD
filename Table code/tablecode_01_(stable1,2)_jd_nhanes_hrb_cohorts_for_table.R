library(data.table)
library(dplyr)
library(tidyr)
library(survey)
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

# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
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
phy_out_cat <- c("dm_f","ckd_f","cvd_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014新发case，横断面数据)
phy_traits_cat <- c("dm_f","dm_b","ckd_f","ckd_b","cvd_f","cvd_b",
                    "ob_f","ob_b","abob_f","abob_b","ir_f","ir_b","dyslip_f","dyslip_b",
                    "mets_f","mets_b","nafld_f","nafld_b","hua_f","hua_b","hpt_f","hpt_b","as_imt_f","as_imt_b")

# 2014、2010表型 (连续)
phy_traits_cont <- c("bmi_f","bmi_b","height_f","height_b","weight_f","weight_b","whr_f","whr_b","wc_f","wc_b","hc_f","hc_b",
                     "tg_f","tg_b","ldl_f","ldl_b","hdl_f","hdl_b","chol_f","chol_b","apoa_f","apoa_b","apob_f","apob_b","nonhdl_f","nonhdl_b",
                     "alt_f","alt_b","ast_f","ast_b","ggt_f","ggt_b","bia_f","bia_b",
                     "egfr_f","egfr_b","scr_f","scr_b","ua_f","ua_b", #"acr_f","acr_b",
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

#### 读取嘉定人群信息 ----
phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))
# 总人群
phy_edc_dat1 <- phy_edc_temp_list[["phy_edc_temp"]] # 提取subgroup
phy_edc_dat1$sex <- phy_edc_dat1$sex_b
#### 读取嘉定人群信息 ####


#### JD 分类变量 (n(%)) ----
phy_traits_cat_b <- phy_traits_cat[grepl("_b",phy_traits_cat)]
phy_traits_cat_f <- phy_traits_cat[grepl("_f",phy_traits_cat)]

# 2010表型 总人群 #
{
  dat_temp <- phy_edc_dat1
  
  cols <- c("sex_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7", 
            "censorall_1021","censorall_1014","censorall_1421",
            "cvd_incident_1021","cvd_incident_1014","cvd_incident_1421",
            "ckd_incident_1014","dm_incident_1014",
            phy_traits_cat_b)
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列16个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- dat_temp[, ..keep_cols]
  dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：全局统计（所有分类变量的数量&占比） ==========
  global_summary <- dat_temp %>%
    # 宽表转长表：每个分类变量+取值为一行
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按变量+取值分组统计
    group_by(分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      p = round(n() / nrow(dat_temp) * 100, 1),
      .groups = "drop"
    )
  
  # ========== 步骤3：分性别统计（每个分类变量按性别统计数量&占比） ==========
  gender_summary <- dat_temp %>%
    # 宽表转长表
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按性别+变量+取值分组
    group_by(!!sym(gender_var), 分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      .groups = "drop_last"
    ) %>%
    # 计算该性别下的占比
    mutate(
      p = round(n / sum(n) * 100, 1)
    ) %>%
    ungroup() %>%
    # 转宽表：把不同性别的统计结果列出来（比如"男_数量""女_占比"）
    pivot_wider(
      names_from = !!gender_var,
      values_from = c(n, p),
      names_sep = "_",
      values_fill = 0  # 无数据的单元格填0
    )
  
  # ========== 步骤4：整合全局+分性别统计到一张表 ==========
  final_summary_b <- global_summary %>%
    left_join(gender_summary, by = c("分类变量名", "分类取值")) %>%
    # 按分类变量名排序
    arrange(分类变量名)
  
  final_summary_b <- final_summary_b[,c("分类变量名","分类取值", "n","p", "n_0","p_0", "n_1","p_1")]
  final_summary_b$n <- paste0(final_summary_b$n," (",final_summary_b$p,"%)")
  final_summary_b$n_0 <- paste0(final_summary_b$n_0," (",final_summary_b$p_0,"%)")
  final_summary_b$n_1 <- paste0(final_summary_b$n_1," (",final_summary_b$p_1,"%)")
  
  final_summary_b$分类变量名 <- factor(final_summary_b$分类变量名, levels = cols)
  final_summary_b1 <- final_summary_b %>%
    arrange(分类变量名, 分类取值)
  final_summary_b1 <- final_summary_b1[,c("分类变量名","分类取值", "n", "n_0", "n_1")]
}

# 2010表型 菌群人群 #
{
  dat_temp <- phy_edc_dat1[phy_edc_dat1$flag_mpa4_14 == 1,]
  
  cols <- c("sex_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg","med_all7",
            "censorall_1021","censorall_1014","censorall_1421",
            "cvd_incident_1021","cvd_incident_1014","cvd_incident_1421",
            "ckd_incident_1014","dm_incident_1014",
            phy_traits_cat_b)
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列16个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- dat_temp[, ..keep_cols]
  dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：全局统计（所有分类变量的数量&占比） ==========
  global_summary <- dat_temp %>%
    # 宽表转长表：每个分类变量+取值为一行
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按变量+取值分组统计
    group_by(分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      p = round(n() / nrow(dat_temp) * 100, 1),
      .groups = "drop"
    )
  
  # ========== 步骤3：分性别统计（每个分类变量按性别统计数量&占比） ==========
  gender_summary <- dat_temp %>%
    # 宽表转长表
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按性别+变量+取值分组
    group_by(!!sym(gender_var), 分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      .groups = "drop_last"
    ) %>%
    # 计算该性别下的占比
    mutate(
      p = round(n / sum(n) * 100, 1)
    ) %>%
    ungroup() %>%
    # 转宽表：把不同性别的统计结果列出来（比如"男_数量""女_占比"）
    pivot_wider(
      names_from = !!gender_var,
      values_from = c(n, p),
      names_sep = "_",
      values_fill = 0  # 无数据的单元格填0
    )
  
  # ========== 步骤4：整合全局+分性别统计到一张表 ==========
  final_summary_b <- global_summary %>%
    left_join(gender_summary, by = c("分类变量名", "分类取值")) %>%
    # 按分类变量名排序
    arrange(分类变量名)
  
  final_summary_b <- final_summary_b[,c("分类变量名","分类取值", "n","p", "n_0","p_0", "n_1","p_1")]
  final_summary_b$n <- paste0(final_summary_b$n," (",final_summary_b$p,"%)")
  final_summary_b$n_0 <- paste0(final_summary_b$n_0," (",final_summary_b$p_0,"%)")
  final_summary_b$n_1 <- paste0(final_summary_b$n_1," (",final_summary_b$p_1,"%)")
  
  final_summary_b$分类变量名 <- factor(final_summary_b$分类变量名, levels = cols)
  final_summary_b2 <- final_summary_b %>%
    arrange(分类变量名, 分类取值)
  final_summary_b2 <- final_summary_b2[,c("分类变量名","分类取值", "n", "n_0", "n_1")]
}

# 2014表型 菌群人群 #
{
  dat_temp <- phy_edc_dat1[phy_edc_dat1$flag_mpa4_14 == 1,]
  
  cols <- c("sex_b","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7",
            "censorall_1021","censorall_1014","censorall_1421",
            "cvd_incident_1021","cvd_incident_1014","cvd_incident_1421",
            "ckd_incident_1014","dm_incident_1014",
            phy_traits_cat_f)
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列16个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- dat_temp[, ..keep_cols]
  dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：全局统计（所有分类变量的数量&占比） ==========
  global_summary <- dat_temp %>%
    # 宽表转长表：每个分类变量+取值为一行
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按变量+取值分组统计
    group_by(分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      p = round(n() / nrow(dat_temp) * 100, 1),
      .groups = "drop"
    )
  
  # ========== 步骤3：分性别统计（每个分类变量按性别统计数量&占比） ==========
  gender_summary <- dat_temp %>%
    # 宽表转长表
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按性别+变量+取值分组
    group_by(!!sym(gender_var), 分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      .groups = "drop_last"
    ) %>%
    # 计算该性别下的占比
    mutate(
      p = round(n / sum(n) * 100, 1)
    ) %>%
    ungroup() %>%
    # 转宽表：把不同性别的统计结果列出来（比如"男_数量""女_占比"）
    pivot_wider(
      names_from = !!gender_var,
      values_from = c(n, p),
      names_sep = "_",
      values_fill = 0  # 无数据的单元格填0
    )
  
  # ========== 步骤4：整合全局+分性别统计到一张表 ==========
  final_summary_b <- global_summary %>%
    left_join(gender_summary, by = c("分类变量名", "分类取值")) %>%
    # 按分类变量名排序
    arrange(分类变量名)
  
  final_summary_b <- final_summary_b[,c("分类变量名","分类取值", "n","p", "n_0","p_0", "n_1","p_1")]
  final_summary_b$n <- paste0(final_summary_b$n," (",final_summary_b$p,"%)")
  final_summary_b$n_0 <- paste0(final_summary_b$n_0," (",final_summary_b$p_0,"%)")
  final_summary_b$n_1 <- paste0(final_summary_b$n_1," (",final_summary_b$p_1,"%)")
  
  final_summary_b$分类变量名 <- factor(final_summary_b$分类变量名, levels = cols)
  final_summary_b3 <- final_summary_b %>%
    arrange(分类变量名, 分类取值)
  final_summary_b3 <- final_summary_b3[,c("分类变量名","分类取值", "n", "n_0", "n_1")]
}

colnames(final_summary_b1)[3:5] <- paste0(colnames(final_summary_b1)[3:5]," | JD_2010")
final_summary_b1$分类变量名 <- gsub("_b|_f","",final_summary_b1$分类变量名)

colnames(final_summary_b2)[3:5] <- paste0(colnames(final_summary_b2)[3:5]," | JD_sub (2010)")
final_summary_b2$分类变量名 <- gsub("_b|_f","",final_summary_b2$分类变量名)

colnames(final_summary_b3)[3:5] <- paste0(colnames(final_summary_b3)[3:5]," | JD_sub (2014)")
final_summary_b3$分类变量名 <- gsub("_b|_f","",final_summary_b3$分类变量名)

final_summary_b_all <- left_join(final_summary_b1, final_summary_b2, by=c("分类变量名", "分类取值")) %>%
  left_join(final_summary_b3, by=c("分类变量名", "分类取值"))

openxlsx::write.xlsx(final_summary_b_all, "tables/jd_cohorts_cat_20260324.xlsx")
#### JD 分类变量 (n(%))  ####

#### JD 连续变量 (mean(SD) or median(interquartile range)) ----
phy_traits_cont_b <- phy_traits_cont[grepl("_b",phy_traits_cont)]
phy_traits_cont_f <- phy_traits_cont[grepl("_f",phy_traits_cont)]

# 2010表型 总人群 #
{
  dat_temp <- phy_edc_dat1
  
  cols <- c("age_b", phy_traits_cont_b,
            "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
            "crp_f",
            "eos_f","lym_f","mon_f","neu_f",
            "nlr_f","lmr_f","plr_f","sii_f","siri_f")
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
      q3 = round(quantile(x_clean, 0.75), 2)
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
      values_from = c(n, mean, sd, median, q1, q3),
      names_sep = "_"
    )
  
  # ========== 步骤5：整合全局统计 + 分性别统计 ==========
  final_stats_b <- global_stats %>%
    left_join(gender_stats, by = "连续变量名") %>%
    arrange(连续变量名)
  
  final_stats_b <- final_stats_b[,c("连续变量名", 
                                    "n","mean","sd","median","q1","q3", 
                                    "n_0","mean_0","sd_0","median_0","q1_0","q3_0", 
                                    "n_1","mean_1","sd_1","median_1","q1_1","q3_1")]
  final_stats_b$mean <- paste0(final_stats_b$mean," (",final_stats_b$sd,")")
  final_stats_b$mean_0 <- paste0(final_stats_b$mean_0," (",final_stats_b$sd_0,")")
  final_stats_b$mean_1 <- paste0(final_stats_b$mean_1," (",final_stats_b$sd_1,")")
  
  final_stats_b$median <- paste0(final_stats_b$median," (",final_stats_b$q1,"-",final_stats_b$q3,")")
  final_stats_b$median_0 <- paste0(final_stats_b$median_0," (",final_stats_b$q1_0,"-",final_stats_b$q3_0,")")
  final_stats_b$median_1 <- paste0(final_stats_b$median_1," (",final_stats_b$q1_1,"-",final_stats_b$q3_1,")")
  
  
  final_stats_b$连续变量名 <- factor(final_stats_b$连续变量名, levels = cols)
  final_stats_b1 <- final_stats_b %>%
    arrange(连续变量名)
  final_stats_b1 <- final_stats_b1[,c("连续变量名", "n","mean","median", "n_0","mean_0","median_0", "n_1","mean_1","median_1")]
}

# 2010表型 菌群人群 #
{
  dat_temp <- phy_edc_dat1[phy_edc_dat1$flag_mpa4_14 == 1,]
  
  cols <- c("age_b", phy_traits_cont_b,
            "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
            "crp_f",
            "eos_f","lym_f","mon_f","neu_f",
            "nlr_f","lmr_f","plr_f","sii_f","siri_f")
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
      q3 = round(quantile(x_clean, 0.75), 2)
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
      values_from = c(n, mean, sd, median, q1, q3),
      names_sep = "_"
    )
  
  # ========== 步骤5：整合全局统计 + 分性别统计 ==========
  final_stats_b <- global_stats %>%
    left_join(gender_stats, by = "连续变量名") %>%
    arrange(连续变量名)
  
  final_stats_b <- final_stats_b[,c("连续变量名", 
                                    "n","mean","sd","median","q1","q3", 
                                    "n_0","mean_0","sd_0","median_0","q1_0","q3_0", 
                                    "n_1","mean_1","sd_1","median_1","q1_1","q3_1")]
  final_stats_b$mean <- paste0(final_stats_b$mean," (",final_stats_b$sd,")")
  final_stats_b$mean_0 <- paste0(final_stats_b$mean_0," (",final_stats_b$sd_0,")")
  final_stats_b$mean_1 <- paste0(final_stats_b$mean_1," (",final_stats_b$sd_1,")")
  
  final_stats_b$median <- paste0(final_stats_b$median," (",final_stats_b$q1,"-",final_stats_b$q3,")")
  final_stats_b$median_0 <- paste0(final_stats_b$median_0," (",final_stats_b$q1_0,"-",final_stats_b$q3_0,")")
  final_stats_b$median_1 <- paste0(final_stats_b$median_1," (",final_stats_b$q1_1,"-",final_stats_b$q3_1,")")
  
  
  final_stats_b$连续变量名 <- factor(final_stats_b$连续变量名, levels = cols)
  final_stats_b2 <- final_stats_b %>%
    arrange(连续变量名)
  final_stats_b2 <- final_stats_b2[,c("连续变量名", "n","mean","median", "n_0","mean_0","median_0", "n_1","mean_1","median_1")]
}

# 2014表型 菌群人群 #
{
  dat_temp <- phy_edc_dat1[phy_edc_dat1$flag_mpa4_14 == 1,]
  
  cols <- c("age_f", phy_traits_cont_f)
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
      q3 = round(quantile(x_clean, 0.75), 2)
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
      values_from = c(n, mean, sd, median, q1, q3),
      names_sep = "_"
    )
  
  # ========== 步骤5：整合全局统计 + 分性别统计 ==========
  final_stats_b <- global_stats %>%
    left_join(gender_stats, by = "连续变量名") %>%
    arrange(连续变量名)
  
  final_stats_b <- final_stats_b[,c("连续变量名", 
                                    "n","mean","sd","median","q1","q3", 
                                    "n_0","mean_0","sd_0","median_0","q1_0","q3_0", 
                                    "n_1","mean_1","sd_1","median_1","q1_1","q3_1")]
  final_stats_b$mean <- paste0(final_stats_b$mean," (",final_stats_b$sd,")")
  final_stats_b$mean_0 <- paste0(final_stats_b$mean_0," (",final_stats_b$sd_0,")")
  final_stats_b$mean_1 <- paste0(final_stats_b$mean_1," (",final_stats_b$sd_1,")")
  
  final_stats_b$median <- paste0(final_stats_b$median," (",final_stats_b$q1,"-",final_stats_b$q3,")")
  final_stats_b$median_0 <- paste0(final_stats_b$median_0," (",final_stats_b$q1_0,"-",final_stats_b$q3_0,")")
  final_stats_b$median_1 <- paste0(final_stats_b$median_1," (",final_stats_b$q1_1,"-",final_stats_b$q3_1,")")
  
  
  final_stats_b$连续变量名 <- factor(final_stats_b$连续变量名, levels = cols)
  final_stats_b3 <- final_stats_b %>%
    arrange(连续变量名)
  final_stats_b3 <- final_stats_b3[,c("连续变量名", "n","mean","median", "n_0","mean_0","median_0", "n_1","mean_1","median_1")]
}

colnames(final_stats_b1)[2:10] <- paste0(colnames(final_stats_b1)[2:10]," | JD_2010")
final_stats_b1$连续变量名 <- gsub("_b|_f","",final_stats_b1$连续变量名)

colnames(final_stats_b2)[2:10] <- paste0(colnames(final_stats_b2)[2:10]," | JD_sub (2010)")
final_stats_b2$连续变量名 <- gsub("_b|_f","",final_stats_b2$连续变量名)

colnames(final_stats_b3)[2:10] <- paste0(colnames(final_stats_b3)[2:10]," | JD_sub (2014)")
final_stats_b3$连续变量名 <- gsub("_b|_f","",final_stats_b3$连续变量名)

final_stats_b_all <- left_join(final_stats_b1, final_stats_b2, by=c("连续变量名")) %>%
  left_join(final_stats_b3, by=c("连续变量名"))

openxlsx::write.xlsx(final_stats_b_all, "tables/jd_cohorts_cont_20260324.xlsx")
#### JD 连续变量 (mean(SD) or median(interquartile range)) ####



### 读取NHANES人群信息 ###
nhanes_dat <- read.csv("C:/TWang/DLiu/NHANES/nhanes_dat_2003-2018_for_analysis.csv")

#### Subgroup for analysis ----
# 13 EDCs (PFAS_5, PAEs_6, TC_1, BP_1); 2011-2014 (2c) #
{
  cols_edc <- c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP",
                "TCS",
                "BPA")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat1 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]), ]
  nhanes_dat1 <- nhanes_dat1[nhanes_dat1$release %in% c("2011-2012","2013-2014"),]
  
  table(nhanes_dat1$release)
}
# 重新计算权重
nhanes_dat1$WTMEC2YR_COMB <- nhanes_dat1$WTMEC2YR/2 # 两轮数据除以2

# PAEs_6; 2003-2018 (8c) #
{
  cols_edc <- c("MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP")
  
  # 筛选最终纳入分析的EDC (不为空值) 的样本
  nhanes_dat5 <- nhanes_dat[complete.cases(nhanes_dat[, cols_edc]), ]
  table(nhanes_dat5$release)
}
# 重新计算权重
nhanes_dat5$WTMEC2YR_COMB <- nhanes_dat5$WTMEC2YR/8 # 八轮数据除以8
#### Subgroup for analysis ####

#### 创建 2003-2018 NHANES调查设计对象 ----
nhanes_design_dat1 <- svydesign(
  id = ~SDMVPSU,          # 整群变量（PSU）
  strata = ~SDMVSTRA,     # 分层变量
  weights = ~WTMEC2YR_COMB,    # 合并调整后2年体检实验室权重
  data = nhanes_dat1,
  nest = TRUE             # NHANES的PSU嵌套在分层中，必须设为TRUE！
)

# 创建 2003-2018 NHANES调查设计对象
nhanes_design_dat5 <- svydesign(
  id = ~SDMVPSU,          # 整群变量（PSU）
  strata = ~SDMVSTRA,     # 分层变量
  weights = ~WTMEC2YR_COMB,    # 合并调整后2年体检实验室权重
  data = nhanes_dat5,
  nest = TRUE             # NHANES的PSU嵌套在分层中，必须设为TRUE！
)

# 创建 2003-2018 NHANES调查设计对象 (去除cvd = NA的case)
nhanes_design_dat6 <- svydesign(
  id = ~SDMVPSU,          # 整群变量（PSU）
  strata = ~SDMVSTRA,     # 分层变量
  weights = ~WTMEC2YR_COMB,    # 合并调整后2年体检实验室权重
  data = nhanes_dat5[!is.na(nhanes_dat5$cvd),],
  nest = TRUE             # NHANES的PSU嵌套在分层中，必须设为TRUE！
)

# 创建 2003-2018 NHANES调查设计对象 (去除ckd = NA的case)
nhanes_design_dat7 <- svydesign(
  id = ~SDMVPSU,          # 整群变量（PSU）
  strata = ~SDMVSTRA,     # 分层变量
  weights = ~WTMEC2YR_COMB,    # 合并调整后2年体检实验室权重
  data = nhanes_dat5[!is.na(nhanes_dat5$ckd),],
  nest = TRUE             # NHANES的PSU嵌套在分层中，必须设为TRUE！
)
#### 创建 2003-2018 NHANES调查设计对象 ####


# nhanes_dat_for_demographic <- nhanes_dat1 # 所有dat1的样本dat5中都包含
nhanes_dat_for_demographic <- nhanes_dat5
#### 非加权NHANES信息 (不包含strata, primary sampling units (PSUs), and sample weight) ----
# 统计分类变量基本情况
{
  nhanes_dat_for_demographic$sex <- nhanes_dat_for_demographic$sex_rev
  cols <- c("race_eth","sex_rev","high_edu","smk1","drk1","mvpa_g","kcal_q4","ucr_q5","release_cycle",
            "dm","ckd","cvd")
  keep_cols <- c("sex",cols)
  # ========== 步骤1：数据预处理（明确变量类型） ==========
  # 1. 定义性别变量名和分类变量名（自动识别，无需手动列16个变量）
  gender_var <- "sex"  # 你的性别变量名
  dat_temp <- nhanes_dat_for_demographic[, keep_cols]
  dat_temp[is.na(dat_temp)] <- 9999
  
  # ========== 步骤2：全局统计（所有分类变量的数量&占比） ==========
  global_summary <- dat_temp %>%
    # 宽表转长表：每个分类变量+取值为一行
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按变量+取值分组统计
    group_by(分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      p = round(n() / nrow(dat_temp) * 100, 1),
      .groups = "drop"
    )
  
  # ========== 步骤3：分性别统计（每个分类变量按性别统计数量&占比） ==========
  gender_summary <- dat_temp %>%
    # 宽表转长表
    pivot_longer(cols = all_of(cols), names_to = "分类变量名", values_to = "分类取值") %>%
    # 按性别+变量+取值分组
    group_by(!!sym(gender_var), 分类变量名, 分类取值) %>%
    summarise(
      n = n(),
      .groups = "drop_last"
    ) %>%
    # 计算该性别下的占比
    mutate(
      p = round(n / sum(n) * 100, 1)
    ) %>%
    ungroup() %>%
    # 转宽表：把不同性别的统计结果列出来（比如"男_数量""女_占比"）
    pivot_wider(
      names_from = !!gender_var,
      values_from = c(n, p),
      names_sep = "_",
      values_fill = 0  # 无数据的单元格填0
    )
  
  # ========== 步骤4：整合全局+分性别统计到一张表 ==========
  final_summary_nhanes <- global_summary %>%
    left_join(gender_summary, by = c("分类变量名", "分类取值")) %>%
    # 按分类变量名排序
    arrange(分类变量名)
  
  final_summary_nhanes <- final_summary_nhanes[,c("分类变量名","分类取值", "n","p", "n_0","p_0", "n_1","p_1")]
  final_summary_nhanes$n <- paste0(final_summary_nhanes$n," (",final_summary_nhanes$p,"%)")
  final_summary_nhanes$n_0 <- paste0(final_summary_nhanes$n_0," (",final_summary_nhanes$p_0,"%)")
  final_summary_nhanes$n_1 <- paste0(final_summary_nhanes$n_1," (",final_summary_nhanes$p_1,"%)")
  
  final_summary_nhanes$分类变量名 <- factor(final_summary_nhanes$分类变量名, levels = cols)
  final_summary_nhanes <- final_summary_nhanes %>%
    arrange(分类变量名, 分类取值)
  final_summary_nhanes <- final_summary_nhanes[,c("分类变量名","分类取值", "n", "n_0", "n_1")]
  
}
openxlsx::write.xlsx(final_summary_nhanes, "tables/nhanes_cat_NOsamplingweights_20260326.xlsx")


summary(nhanes_dat_for_demographic$age) # 20~80 years old (80:80 years of age and over)
sd(nhanes_dat_for_demographic$age)

summary(nhanes_dat$ucr_log10) # log10 transformed urinary creatinine concentrations (mg/dL)
sd(nhanes_dat_for_demographic$ucr_log10)

# table(nhanes_dat$race_eth, useNA = "always") # 0:Hispanic; 1:Non-Hispanic White; 2:Non-Hispanic Black; 3:Other
# table(nhanes_dat_for_demographic$sex_rev, useNA = "always") # 0:women, 1:men
# table(nhanes_dat_for_demographic$high_edu, useNA = "always") # 0:less than high school, 1:High school graduate or higher
# table(nhanes_dat_for_demographic$smk1, useNA = "always") # Do you now smoke cigarettes? 0:Not at All (Refused, Don't know, Missing), 1:Every day or Some days
# table(nhanes_dat_for_demographic$drk1, useNA = "always") # How often drink alcohol over past 12 mos? 0:Never (Don't know, Missing), 1:1 to 365 Range of Values
# table(nhanes_dat_for_demographic$mvpa_g, useNA = "always") # 0:Vigorous and Moderate recreational activities (less than 150 mins/per week), 1:Vigorous and Moderate recreational activities (greater or equal to 150 mins/per week)
# table(nhanes_dat_for_demographic$kcal_q4, useNA = "always") # 0:Q1, 1:Q2, 2:Q3, 3:Q4, 4:NA (女性男性分别计算)
# 
# table(nhanes_dat_for_demographic$dm, useNA = "always")
# table(nhanes_dat_for_demographic$ckd, useNA = "always")
# table(nhanes_dat_for_demographic$cvd, useNA = "always")
#### 非加权NHANES信息 (不包含strata, primary sampling units (PSUs), and sample weight) ####


# nhanes_design <- nhanes_design_dat1 # 所有dat1的样本dat5中都包含
nhanes_design <- nhanes_design_dat5 # 2003-2018 NHANES调查设计对象
# nhanes_design <- nhanes_design_dat6 # 2003-2018 NHANES调查设计对象 (去除cvd = NA的case)
# nhanes_design <- nhanes_design_dat7 # 2003-2018 NHANES调查设计对象 (去除ckd = NA的case)
#### 加权NHANES信息 (包含strata, primary sampling units (PSUs), and sample weight) ----
{
  # 用svytotal直接计算各分类的加权总人数
  # 1. 定义需要计算的分类变量列表（关键：和你数据中的变量名一致）
  var_list <- c(
    "race_eth",   # 种族/民族
    "sex_rev",    # 性别
    "high_edu",   # 高学历
    "smk1",       # 吸烟
    "drk1",       # 饮酒
    "mvpa_g",     # 中高强度运动
    "kcal_q4",    # 热量四分位数
    "ucr_q5",
    "release_cycle",
    "dm",         # 糖尿病
    "ckd",        # 慢性肾病
    "cvd"         # 心血管疾病
  )
  
  # 2. 初始化空列表，用于存储每个变量的结果
  result_list <- list()
  
  # 3. 循环计算每个变量的加权总人数
  for (var_name in var_list) {
    # 跳过不存在的变量（避免循环报错）
    if (!var_name %in% colnames(nhanes_design$variables)) {
      warning(paste("变量", var_name, "不存在，跳过"))
      next
    }
    
    # 构建公式（factor转换分类变量）
    formula_str <- paste("~factor(", var_name, ")")
    formula <- as.formula(formula_str)
    
    # 计算加权总人数（带标准误）
    svy_result <- svytotal(formula, design = nhanes_design)
    
    # 整理成数据框：分类标签 + 加权人数 + 标准误
    result_df <- data.frame(
      变量名 = var_name,
      分类标签 = names(coef(svy_result)),  # 分类（如1/2/男/女）
      加权总人数 = round(coef(svy_result)), # 四舍五入为整数
      标准误 = round(SE(svy_result))        # 标准误（保留整数）
    )
    
    # 将当前变量结果存入列表
    result_list[[var_name]] <- result_df
  }
  
  # 4. 合并所有变量的结果到一个数据框（核心输出）
  weighted_total_all <- bind_rows(result_list)
  
  # 5. 可选：美化分类标签（比如把1/2替换成中文，根据你的数据调整）
  weighted_total_all$分类标签 <- gsub(
    pattern = "factor\\(.*\\)(\\d+)",  # 匹配factor(变量名)1/2的格式
    replacement = "\\1",               # 提取数字部分
    x = weighted_total_all$分类标签
  )
}
openxlsx::write.xlsx(weighted_total_all, "tables/nhanes_cat_WITHsamplingweights_20260326.xlsx")


# ========== 年龄：加权均数、SE、95%CI ========== #
age_stats <- svymean(~age, design = nhanes_design, na.rm = TRUE)
age_var <- svyvar(~age, design = nhanes_design, na.rm = TRUE) # 计算加权方差
age_quantiles <- svyquantile(~age, design = nhanes_design, na.rm = TRUE,
                             quantiles = c(0.25, 0.5, 0.75), method = "constant")

# 提取具体数值
age_mean <- round(coef(age_stats)[1], 2)       # 加权均数（保留1位小数）
age_sd <- sqrt(coef(age_var)[1])               # 方差开平方=标准差
age_se <- round(SE(age_stats)[1], 1)           # 标准误
age_ci_lower <- round(confint(age_stats)[1,1], 1)  # 95%CI下限
age_ci_upper <- round(confint(age_stats)[1,2], 1)  # 95%CI上限

age_q1 <- round(age_quantiles$age[1], 1)       # 25%分位数
age_median <- round(age_quantiles$age[2], 1)   # 中位数
age_q3 <- round(age_quantiles$age[3], 1)       # 75%分位数

# ========== 尿肌酐：加权均数、SE、95%CI ========== #
ucr_stats <- svymean(~ucr_log10, design = nhanes_design, na.rm = TRUE)
ucr_var <- svyvar(~ucr_log10, design = nhanes_design, na.rm = TRUE) # 计算加权方差
ucr_quantiles <- svyquantile(~ucr_log10, design = nhanes_design, na.rm = TRUE,
                             quantiles = c(0.25, 0.5, 0.75), method = "constant")

# 提取具体数值
ucr_mean <- round(coef(ucr_stats)[1], 2)       # 加权均数（保留1位小数）
ucr_sd <- sqrt(coef(ucr_var)[1])               # 方差开平方=标准差
ucr_se <- round(SE(ucr_stats)[1], 1)           # 标准误
ucr_ci_lower <- round(confint(ucr_stats)[1,1], 1)  # 95%CI下限
ucr_ci_upper <- round(confint(ucr_stats)[1,2], 1)  # 95%CI上限

ucr_q1 <- round(ucr_quantiles$age[1], 1)       # 25%分位数
ucr_median <- round(ucr_quantiles$age[2], 1)   # 中位数
ucr_q3 <- round(ucr_quantiles$age[3], 1)       # 75%分位数
#### 加权NHANES信息 (包含strata, primary sampling units (PSUs), and sample weight) ####
