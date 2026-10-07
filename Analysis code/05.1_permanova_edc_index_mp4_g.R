library(data.table)
library(dplyr)
library(chemometrics)
library(vegan)
library(pcaPP)
### EDC log10转换 ###
### 使用 MPA4 (genus 和 species) ###

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
# mp3_g_zero <- paste0(mp3_g_names,"_zero") # 菌群MP3填补0值丰度 (属层面)
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
                    "ob_f","ob_b","abob_f","abob_b","dyslip_f","dyslip_b","hua_f","hua_b","ir_f","ir_b","mets_f","mets_b",
                    
                    # "smk1_f","drk1_f","paactive3_g_f",
                    
                    "sitduration_f","sitduration_b","sleeptg_f",
                    "dm_treat_f","dm_treat_b","hpt_treat_f","hpt_treat_b","hpl_treat_f","hpl_treat_b",
                    "diet_score_g_f","high_fruveg_f","low_ssb_f","low_meat_f","high_fish_f")
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
                     "nlr_f","lmr_f","plr_f","sii_f","siri_f",
                     
                     "sleept_f","sittimet_f","sittimet_b","sum_met_f","sum_met_b",
                     "alco_f","alco_b","diet_score_f")
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

#### EDC INDEX traits ----
edc_index_b <- c("edc_count1_edc19_b","edc_count1_edc14_b","edc_count1_edc9_b","edc_count1_pfas_b","edc_count1_pae_b","edc_count1_pae6_b","edc_count1_pae4_b","edc_count1_bp_b","edc_count1_bp1_b","edc_count1_tc_b",
                 "edc_count2_edc19_b","edc_count2_edc14_b","edc_count2_edc9_b","edc_count2_pfas_b","edc_count2_pae_b","edc_count2_pae6_b","edc_count2_pae4_b","edc_count2_bp_b","edc_count2_bp1_b","edc_count2_tc_b",
                 "edc_count3_edc19_b","edc_count3_edc14_b","edc_count3_edc9_b","edc_count3_pfas_b","edc_count3_pae_b","edc_count3_pae6_b","edc_count3_pae4_b","edc_count3_bp_b","edc_count3_bp1_b","edc_count3_tc_b",
                 
                 "edc_score_edc19_b","edc_score_edc14_b","edc_score_edc9_b","edc_score_pfas_b","edc_score_pae_b","edc_score_pae6_b","edc_score_pae4_b","edc_score_bp_b","edc_score_bp1_b","edc_score_tc_b")

edc_index_f <- c("edc_count1_edc19_f","edc_count1_edc14_f","edc_count1_edc9_f","edc_count1_pfas_f","edc_count1_pae_f","edc_count1_pae6_f","edc_count1_pae4_f","edc_count1_bp_f","edc_count1_bp1_f","edc_count1_tc_f",
                 "edc_count2_edc19_f","edc_count2_edc14_f","edc_count2_edc9_f","edc_count2_pfas_f","edc_count2_pae_f","edc_count2_pae6_f","edc_count2_pae4_f","edc_count2_bp_f","edc_count2_bp1_f","edc_count2_tc_f",
                 "edc_count3_edc19_f","edc_count3_edc14_f","edc_count3_edc9_f","edc_count3_pfas_f","edc_count3_pae_f","edc_count3_pae6_f","edc_count3_pae4_f","edc_count3_bp_f","edc_count3_bp1_f","edc_count3_tc_f",
                 
                 "edc_score_edc19_f","edc_score_edc14_f","edc_score_edc9_f","edc_score_pfas_f","edc_score_pae_f","edc_score_pae6_f","edc_score_pae4_f","edc_score_bp_f","edc_score_bp1_f","edc_score_tc_f")


edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")

edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### EDC INDEX traits ####

cov_traits <- c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")
cov_traits_cat <- c("sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")


# 汇总分类和连续变量
env_out_edc <- c("age_f", cov_traits_cat, edc_traits_log10, edc_index_f_keep)


#### 读取d.bray数据 ----
beta_diversity_list_mp4_g <- readRDS(paste0("results/indices/beta_diversity_mp4_genus_20260313.rds"))
beta_diversity_list_mp4_s <- readRDS(paste0("results/indices/beta_diversity_mp4_species_20260313.rds"))
#### 读取d.bray数据 ----


#### PERMANOVA (MP4, at the genus level) (分类表型和连续表型的genus水平的Bray-Curtis dissimilarity差异) ----
dat <- beta_diversity_list_mp4_g[[1]]
# # Bray-Curtis dissimilarity结果默认是dist对象，可以转换为矩阵查看
# matrix_bray <- as.matrix(beta_diversity_list_mp4_g[[1]])
# # 查看前5x5矩阵
# print(matrix_bray[1:5, 1:5])
cov <- beta_diversity_list_mp4_g[[2]]
keep_col <- c("id14_15",env_out_edc)
cov <- cov[,..keep_col]
# 分类校正变量必须转换为因子
for (i in cov_traits_cat) {
  cov[[i]] <- factor(cov[[i]]) 
}
# ### 获取缺失列名及对应缺失数量 ###
# missing_stats <- colSums(is.na(cov))              # 计算每列缺失总数
# missing_cols <- names(missing_stats[missing_stats > 0])  # 筛选有缺失的列名
# result <- data.frame(Column = missing_cols, Missing = missing_stats[missing_cols])
# row.names(result) <- NULL
# print(result)


# 由于cat变量包含小于二分类的变量（如PFAS的detected），筛选仅含一种唯一值的变量名（忽略NA）
temp <- cov[,..cov_traits_cat]
one_value_vars <- names(temp)[sapply(temp, function(x) length(unique(na.omit(x))) == 1)]
# 从cov和变量名向量中删除只有一个值的变量名
env_out_edc <- env_out_edc[!env_out_edc %in% one_value_vars]
keep_col <- c("id14_15",env_out_edc)
cov <- cov[,..keep_col]


permanova_results_all <- data.frame()
for (i in 1:length(env_out_edc)){
  # i <- 4
  # i <- 20
  
  # 筛选校正变量和主要分析变量
  if(env_out_edc[i] %in% cov_traits){
    col_temp <- c("id14_15",cov_traits)
    # 注意：主要分析元素要放在最后
    non_target <- col_temp[col_temp != env_out_edc[i]]       # 提取非目标元素
    col_temp <- c(non_target,env_out_edc[i])       # 合并
    cov_temp <- cov[,..col_temp]
  }else{
    # 注意：主要分析元素要放在最后
    col_temp <- c("id14_15",cov_traits,env_out_edc[i])
    cov_temp <- cov[,..col_temp]
  }
  
  # 删除主要分析变量的缺失值
  cov_temp <- cov_temp[!is.na(cov_temp[[env_out_edc[i]]])]
  
  # 选取剩余的ID
  id_temp <- cov_temp[["id14_15"]]
  
  # 根据ID选取sub matrix
  dat_matrix_temp <- as.matrix(dat) # Bray-Curtis dissimilarity结果默认是dist对象，可以转换为矩阵查看
  dat_sub_matrix_temp <- dat_matrix_temp[id_temp,id_temp]
  
  # 把sub matrix转换回dist数据for analysis
  dat_temp <- as.dist(dat_sub_matrix_temp)
  
  # 删除协变量数据框ID
  cov_temp <- cov_temp[,-1]
  
  
  # ### 单因素PERMANOVA (不校正其他covariates) ###
  # print(paste0(i," out of ",length(env_out_edc)," || ",env_out_edc[i]," (单因素) || ",Sys.time()))
  # f1 <- as.formula(paste0("dat_temp ~ ",env_out_edc[i])) # 构建单因素分析公式
  # set.seed(123)  # 固定随机种子
  # adonis_result1 <- adonis2(f1, data = cov_temp, permutations = 999, parallel = 8)
  # # print(adonis_result1)
  # 
  # adonis_result1$n <- nrow(cov_temp)
  # # 计算调整R方
  # adonis_result1$adjusted_R2 <- RsquareAdj(adonis_result1$R2, adonis_result1$n, adonis_result1$Df)
  # # 计算调整R方
  # adonis_result1$cov <- row.names(adonis_result1)
  # adonis_result1$trait <- env_out_edc[i]
  # adonis_result1$adj <- "univariable"
  # adonis_result1$formula <- paste0("dat ~ ",env_out_edc[i])
  # permanova_results_all <- rbind(permanova_results_all,adonis_result1)
  # ### 单因素PERMANOVA (不校正其他covariates) ###
   
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
openxlsx::write.xlsx(permanova_results_all, paste0("results/permanova/permanova_edc_index-mp4_genus-20260313.xlsx"))
#### PERMANOVA (MP4, at the genus level) (分类表型和连续表型的genus水平的Bray-Curtis dissimilarity差异) ####
