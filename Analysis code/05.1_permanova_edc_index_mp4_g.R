library(data.table)
library(dplyr)
library(chemometrics)
library(vegan)
library(pcaPP)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

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


cov_traits <- c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")
cov_traits_cat <- c("sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")


# 汇总分类和连续变量
env_out_edc <- c("age_f", cov_traits_cat, edc_traits_log10, edc_index_f_keep)


#### 读取d.bray数据 ----
beta_diversity_list_mp4_g <- readRDS(paste0("results/indices/beta_diversity_mp4_genus.rds"))
#### 读取d.bray数据 ----


#### PERMANOVA (MP4, at the genus level) (分类表型和连续表型的genus水平的Bray-Curtis dissimilarity差异) ----
dat <- beta_diversity_list_mp4_g[[1]]
cov <- beta_diversity_list_mp4_g[[2]]
keep_col <- c("id14_15",env_out_edc)
cov <- cov[,keep_col]
# 分类校正变量必须转换为因子
for (i in cov_traits_cat) {
  cov[[i]] <- factor(cov[[i]]) 
}

# 筛选仅含一种唯一值的变量名（忽略NA）
temp <- cov[,cov_traits_cat]
one_value_vars <- names(temp)[sapply(temp, function(x) length(unique(na.omit(x))) == 1)]
# 从cov和变量名向量中删除只有一个值的变量名
env_out_edc <- env_out_edc[!env_out_edc %in% one_value_vars]
keep_col <- c("id14_15",env_out_edc)
cov <- cov[,keep_col]


permanova_results_all <- data.frame()
for (i in 1:length(env_out_edc)){
  # i <- 4
  
  # 筛选校正变量和主要分析变量
  if(env_out_edc[i] %in% cov_traits){
    col_temp <- c("id14_15",cov_traits)
    # 注意：主要分析元素要放在最后
    non_target <- col_temp[col_temp != env_out_edc[i]]       # 提取非目标元素
    col_temp <- c(non_target,env_out_edc[i])       # 合并
    cov_temp <- cov[,col_temp]
  }else{
    # 注意：主要分析元素要放在最后
    col_temp <- c("id14_15",cov_traits,env_out_edc[i])
    cov_temp <- cov[,col_temp]
  }
  
  # 删除主要分析变量的缺失值
  cov_temp <- cov_temp[!is.na(cov_temp[[env_out_edc[i]]]),]
  
  # 选取剩余的ID
  id_temp <- cov_temp[["id14_15"]]
  
  # 根据ID选取sub matrix
  dat_matrix_temp <- as.matrix(dat) # Bray-Curtis dissimilarity结果默认是dist对象，可以转换为矩阵查看
  dat_sub_matrix_temp <- dat_matrix_temp[id_temp,id_temp]
  
  # 把sub matrix转换回dist数据for analysis
  dat_temp <- as.dist(dat_sub_matrix_temp)
  
  # 删除协变量数据框ID
  cov_temp <- cov_temp[,-1]
  
  ### 多因素PERMANOVA (校正其他covariates) 注意：主要分析元素要放在最后 ###
  print(paste0(i," out of ",length(env_out_edc)," || ",env_out_edc[i]," (多因素，校正其他协变量) || ",Sys.time()))
  cov_order <- colnames(cov_temp) # 移动env_out_edc[i]元素到末尾 (已经在开头进行操作过了)
  f2 <- as.formula(paste0("dat_temp ~ ",paste(cov_order,collapse = " + "))) # 构建多因素分析公式
  set.seed(123)  # 固定随机种子
  adonis_result2 <- adonis2(f2, data = cov_temp, permutations = 999, parallel = 8)
  # print(adonis_result2)
  
  adonis_result2$n <- nrow(cov_temp)
  # 计算调整R方
  adonis_result2$adjusted_R2 <- RsquareAdj(adonis_result2$R2, adonis_result2$n, adonis_result2$Df)
  # 计算调整R方
  adonis_result2$cov <- row.names(adonis_result2)
  adonis_result2$trait <- env_out_edc[i]
  adonis_result2$adj <- "multivariable"
  adonis_result2$formula <- paste0("dat ~ ",paste(cov_order,collapse = " + "))
  permanova_results_all <- rbind(permanova_results_all,adonis_result2)
  ### 多因素PERMANOVA (校正其他covariates) 注意：主要分析元素要放在最后 ###
}
openxlsx::write.xlsx(permanova_results_all, paste0("results/permanova/permanova_edc_index-mp4_genus.xlsx"))
#### PERMANOVA (MP4, at the genus level) (分类表型和连续表型的genus水平的Bray-Curtis dissimilarity差异) ####
