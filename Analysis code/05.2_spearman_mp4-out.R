library(data.table)
library(dplyr)
library(ppcor)
### EDC log10转换 ###
### 使用 MPA4 (genus 和 species) ###

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))

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

# 2021、2014新发表型 (分类)
phy_incident_cat <- c("cvd_incident_1021","cvd_incident_1014","cvd_incident_1421","ckd_incident_1014","dm_incident_1014")    # 新发 cvd, ckd, dm 去除基线 case (只做EDC对outcome，不做cvd_incident_1421)
phy_incident_time <- c("timecvd_1021","timecvd_1014","timecvd_1421","timeckd_1014","timedm_1014")
phy_censor_cat <- c("censorall_1021","censorall_1014")
phy_censor_time <- c("timeall_1021","timeall_1014")
# 2014、2010表型 (分类)
phy_out_cat <- c("cvd_f","ckd_f","dm_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014新发case，横断面数据)
phy_traits_cat <- c("as_imt_f","hpt_f","nafld_f","ob_f","abob_f","dyslip_f","hua_f","ir_f","mets_f")
# 2014、2010表型 (连续)
phy_traits_cont <- c("bmi_f","wc_f","hc_f","whr_f","height_f","weight_f",
                     "hdl_f","ldl_f","apoa_f","apob_f","chol_f","tg_f","nonhdl_f",
                     "alt_f","ast_f","ggt_f","scr_f","egfr_f","acr_f","ua_f","bia_f",
                     "glu0_f","glu120_f","vhba1c_f","ins0_f","ins120_f","homair_f","homab_f",
                     "sbp_f","dbp_f","pr_f",
                     "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                     "wbc_f","crp_f",
                     "plt_f","hgb_f","eos_f","lym_f","mon_f","neu_f",
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

#### 构建Spearman分析函数 ----
### 创建一个函数，用于计算连续变量之间的单变量Spearman相关性
spearman_test <- function(EXP, OUT, DAT) {
  # Spearman相关性测试
  results <- cor.test(DAT[[OUT]], DAT[[EXP]], method = "spearman")
  temp_spearman <- data.frame(
    exp_name = EXP,
    out_name = OUT,
    method = "spearman",
    estimate = results$estimate,
    p = results$p.value,
    sample = sample_name,
    n = length(na.omit(DAT[[OUT]]))
  )
  return(temp_spearman)
}
### 创建一个函数，用于计算连续变量之间的多变量Partial Spearman相关性 (Partial Spearman 中校正分类变量不能是因子)
partial_spearman_test <- function(EXP, OUT, COV, DAT) {
  DAT <- DAT[!is.na(DAT[[OUT]]) & !is.na(DAT[[EXP]]),]
  if(EXP == OUT){
    temp_spearman <- data.frame(
      exp_name = EXP,
      out_name = OUT,
      method = "partial spearman",
      estimate = NA,
      p = NA,
      sample = sample_name,
      n = length(na.omit(DAT[[OUT]]))
    )
    return(temp_spearman)
  }else{
    results <- pcor.test(DAT[[OUT]], DAT[[EXP]], DAT[,..COV], method = "spearman")
    temp_spearman <- data.frame(
      exp_name = EXP,
      out_name = OUT,
      method = "partial spearman",
      estimate = results$estimate,
      p = results$p.value,
      sample = sample_name,
      n = length(na.omit(DAT[[OUT]]))
    )
    return(temp_spearman)
  }
}
#### 构建Spearman分析函数 ####

#### MP4 || Outcomes ----
spearman_results_list <- list()
for (i in c("phy_edc_temp0_3")) {
  # i <- "phy_edc_temp0_3"

  phy_edc_dat <- phy_edc_temp_list[[i]]
  sample_name <- i
  
  #### Spearman (菌群和结局指标之间的关联) ----
  spearman_results <- data.frame()
  for (k in mp4_s_log10) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," |Spearman (exposure) ",which(c(mp4_s_log10) == k)," out of ",length(c(mp4_s_log10))))
    for (m in c(phy_traits_cont)) {
      
      ## 分析样本选取 (排除缺失的项) ##
      cols <- c(k, m)
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)

      result <- spearman_test(k, m, phy_edc_dat_temp)
      spearman_results <- rbind(spearman_results, result)
    }
  }
  spearman_results_list[[paste0("spearman_results-(",sample_name,")")]] <- spearman_results
  #### Spearman (菌群和结局指标之间的关联) ####
}
### 保存单个数据 ###
saveRDS(spearman_results_list, paste0("results/correlations/spearman/spearman_results_mp4-out_20260729.rds"))


partial_spearman_results_list <- list() # Partial Spearman 中校正分类变量不能是因子
for (i in c("phy_edc_temp0_3")) {
  # i <- "phy_edc_temp0_3"
  
  phy_edc_dat <- phy_edc_temp_list[[i]]
  sample_name <- i
  
  #### Partial Spearman (菌群和结局指标之间的关联，校正协变量) ----
  partial_spearman_results <- data.frame()
  for (k in mp4_s_log10) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," |Partial Spearman (exposure) ",which(c(mp4_s_log10) == k)," out of ",length(c(mp4_s_log10))))
    for (m in c(phy_traits_cont)) {
      
      ## 选取校正变量 ##
      # 菌群人群，菌群变量（校正用药）
      cov_traits_for_p_spearman <- c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7") # Partial Spearman 中校正分类变量不能是因子
      
      ## 分析样本选取 (排除缺失的项) ##
      cols <- c(k, m, cov_traits_for_p_spearman)
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      
      tryCatch(expr = {
        result <- partial_spearman_test(k, m, cov_traits_for_p_spearman, phy_edc_dat_temp)
        partial_spearman_results <- rbind(partial_spearman_results, result)
      },
      error = function(e) {
        message("错误信息: ", e$message)
        return(NA)  # 返回默认值
      })
      
    }
  }
  
  partial_spearman_results_list[[paste0("partial_spearman_results-(",sample_name,")")]] <- partial_spearman_results
  #### Partial Spearman (菌群和结局指标之间的关联，校正协变量) ####
}
### 保存单个数据 ###
saveRDS(partial_spearman_results_list, paste0("results/correlations/spearman/partial_spearman_results_mp4-out_20260729.rds"))
#### MP4 || Outcomes ####
