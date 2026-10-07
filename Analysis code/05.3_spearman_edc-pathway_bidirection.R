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
# # mp3中用于构建ma的菌的名称
# mp3_ma_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_ma.txt")
# mp3_ma_names <- mp3_ma_names$V1
# mp3_ma_names_log10 <- paste0(mp3_ma_names,"_log10")
# # microbial age (MA)
# # mean(phy_edc_dat$MA,na.rm = TRUE)

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

# edc index 变量名 #
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")

edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### 读取EDC INDEX数据 ####

#### 读取 Pathway变量名 ----
# 丰度>0.0001, 出现率>10%的pathway
pathway_mp3 <- read.table("jiading/sourceDataTaxon/mpa3/pathway_names_mp3_10%.txt",header = TRUE)
pathway_mp3_names_labels <- pathway_mp3[,1]
pathway_mp3_names <- pathway_mp3[,2]
# 转换后的pathway名称
pathway_mp3_log10 <- paste0(pathway_mp3_names,"_log10") # Pathway丰度的log10转换

# 完整386个pathway
pathway_mp3_clean <- fread("jiading/sourceDataTaxon/mpa3/WTG.humann3.pathway.clean_log10.csv")
#### 读取 Pathway变量名 ####


# 总体MP4变量名
mp4_traits <- c(mp4_s_log10)
# 总体分类变量名
cat_traits <- c(phy_incident_cat, phy_censor_cat, phy_out_cat, phy_traits_cat)
# 总体连续变量名
cont_traits <- c(phy_traits_cont)
# 总体EDC变量名
edc_traits <- c(edc_index_f_keep, edc_traits_log10)


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


#### EDC || EDC VS. Pathways ----
spearman_results_list <- list()
for (i in c("phy_edc_temp0_3")) {
  # i <- "phy_edc_temp0_3"
  
  phy_edc_dat <- phy_edc_temp_list[[i]]
  sample_name <- i
  phy_edc_dat <- left_join(phy_edc_dat,edc_count_b,by="ID") %>%
    left_join(edc_score_b,by="ID") %>%
    left_join(edc_count_f,by="ID") %>%
    left_join(edc_score_f,by="ID") %>%
    left_join(pathway_mp3_clean,by="id14_15")
  
  #### Spearman (菌群与Pathways之间的关联) ----
  spearman_results <- data.frame()
  for (k in edc_traits) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," |Spearman (exposure) ",which(c(edc_traits) == k)," out of ",length(c(edc_traits))))
    for (m in c(pathway_mp3_log10)) {
      
      ## 分析样本选取 (排除缺失的项) ##
      cols <- c(k, m)
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      
      result <- spearman_test(k, m, phy_edc_dat_temp)
      spearman_results <- rbind(spearman_results, result)
    }
  }
  spearman_results_list[[paste0("spearman_results-(",sample_name,")")]] <- spearman_results
  #### Spearman (菌群与Pathways之间的关联) ####
}
### 保存单个数据 ###
saveRDS(spearman_results_list, paste0("results/correlations/spearman/spearman_results_edc-pathway_20260313.rds"))


partial_spearman_results_list <- list() # Partial Spearman 中校正分类变量不能是因子
partial_spearman_results_list2 <- list() # Partial Spearman 中校正分类变量不能是因子
for (i in c("phy_edc_temp0_3")) {
  # i <- "phy_edc_temp0_3"
  
  phy_edc_dat <- phy_edc_temp_list[[i]]
  sample_name <- i
  phy_edc_dat <- left_join(phy_edc_dat,edc_count_b,by="ID") %>%
    left_join(edc_score_b,by="ID") %>%
    left_join(edc_count_f,by="ID") %>%
    left_join(edc_score_f,by="ID") %>%
    left_join(pathway_mp3_clean,by="id14_15")
  
  #### Partial Spearman(菌群与Pathways之间的关联，校正协变量) ----
  partial_spearman_results <- data.frame()
  partial_spearman_results2 <- data.frame()
  for (k in edc_traits) {
    print(paste0(round(Sys.time(),0)," || ",sample_name," |Partial Spearman (exposure) ",which(c(edc_traits) == k)," out of ",length(c(edc_traits))))
      
    for (m in c(pathway_mp3_log10)) {
      
      cov_traits_for_p_spearman <- c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7") # Partial Spearman 中校正分类变量不能是因子
      ## 分析样本选取 (排除缺失的项) ##
      cols <- c(k, m, cov_traits_for_p_spearman)
      phy_edc_dat_temp <- phy_edc_dat[,..cols]
      phy_edc_dat_temp <- na.omit(phy_edc_dat_temp)
      
      # EDC (EXP+校正变量) 与 Pathways (OUT)
      tryCatch(expr = {
        result <- partial_spearman_test(k, m, cov_traits_for_p_spearman, phy_edc_dat_temp)
        partial_spearman_results <- rbind(partial_spearman_results, result)
      },
      error = function(e) {
        message("错误信息: ", e$message)
        return(NA)  # 返回默认值
      })
      
      # # Pathways (EXP+校正变量) 与 EDC (OUT)
      # tryCatch(expr = {
      #   result2 <- partial_spearman_test(m, k, cov_traits_for_p_spearman, phy_edc_dat_temp)
      #   partial_spearman_results2 <- rbind(partial_spearman_results2, result2)
      # },
      # error = function(e) {
      #   message("错误信息: ", e$message)
      #   return(NA)  # 返回默认值
      # })
      
    }
  }
  
  partial_spearman_results_list[[paste0("partial_spearman_results-(",sample_name,")")]] <- partial_spearman_results
  # partial_spearman_results_list2[[paste0("partial_spearman_results-(",sample_name,")")]] <- partial_spearman_results2
  #### Partial Spearman(菌群与Pathways之间的关联，校正协变量) ####
}
### 保存单个数据 ###
saveRDS(partial_spearman_results_list, paste0("results/correlations/spearman/partial_spearman_results_edc-pathway_20260313.rds"))
# saveRDS(partial_spearman_results_list2, paste0("results/correlations/spearman/partial_spearman_results_pathway-edc_20260313.rds"))
#### EDC || EDC VS. Pathways ####
