library(data.table)
library(dplyr)
library(ggplot2)
library(corrplot)
library(linkET)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

phy_edc_temp_list <- readRDS(paste0("jiading/data_for_analysis/phy_edc_subset.rds"))
### 提取完整人群for EDC analysis ###
sample_name <- "phy_edc_temp"  # 提取subgroup的名称
phy_edc_dat <- phy_edc_temp_list[[sample_name]] # 提取subgroup

#### 配色 ----
PFAS_colors <- colorRampPalette(c("#f55d78","#FFFFFF"))(50)[c(1,8,15,22,29)]
PAE_colors <- colorRampPalette(c("#3498db","#FFFFFF"))(50)[c(1,6,11,16,21,26,31,36,41)]
BP_colors <- colorRampPalette(c("#f1c40f","#FFFFFF"))(50)[c(1,18,35)]
TC_colors <- colorRampPalette(c("#00b894","#FFFFFF"))(50)[c(1,21)]
my_colors <- c(PFAS_colors,PAE_colors,BP_colors,TC_colors)
#### 配色 ####

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
phy_out_cat <- c("cvd_f","ckd_f","dm_f") # cvd, ckd, dm 包括基线 case (2010基线case+2014或2021新发case，横断面数据)
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


#### EDC-EDC 相关性 (Partial Spearman) 填补后 ----
### 读取EDC-EDC关联性数据
p_spearman_results <- readRDS("results/correlations/spearman/partial_spearman_results_edc-edc_mp4_out_20260313.rds")
p_spearman_results <- p_spearman_results[[1]]
p_spearman_results <- p_spearman_results[p_spearman_results$exp_name %in% edc_traits_log10 & p_spearman_results$out_name %in% edc_traits_log10,]

p_spearman_results$exp_name <- factor(p_spearman_results$exp_name, 
                                      levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                        "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                        "TCC","TCS",
                                                        "BPA","BPS","BPF"),"_log10"),
                                      labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                 "TCC","TCS",
                                                 "BPA","BPS","BPF"))
p_spearman_results$out_name <- factor(p_spearman_results$out_name,
                                      levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                        "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                        "TCC","TCS",
                                                        "BPA","BPS","BPF"),"_log10"),
                                      labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                 "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                 "TCC","TCS",
                                                 "BPA","BPS","BPF"))


p_spearman_results <- p_spearman_results %>%
  arrange(exp_name,out_name)
# 以每个outcome表型为单位，校正exposure
p_spearman_results <- p_spearman_results %>%
  group_by(out_name) %>%  # 按outcome分组，校正exposure
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


df <- p_spearman_results
df <- df[df$exp_name != df$out_name,]
df <- df[,c("exp_name","out_name","estimate","p","p_adj_bh")]
# 核心：无向组合去重 #
df$key <- apply(df[, c("exp_name", "out_name")], 1, function(x) {
  paste(sort(x), collapse = "_")  # 排序后拼接成唯一键
})

# 保留每个唯一组合的第一行，删除重复的 (a,b) (b,a) #
df <- df[!duplicated(df$key), ]

# 删掉辅助列（可选） #
df$key <- NULL


colnames(df) <- c("Analyte 1","Analyte 2","Spearman's rho (two-sided)","P value","BH-adjusted P")
openxlsx::write.xlsx(df,"tables/edc_edc_correaltion_20260511.xlsx")
#### EDC-EDC 相关性 (Partial Spearman) 填补后 ####

#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ----
### 读取分类组间EDC差异数据
wilcoxon_results <- readRDS("results/correlations/wilcoxon/wilcoxon_results_edc_20260313.rds")
wilcoxon_results <- wilcoxon_results[[1]]
wilcoxon_results <- wilcoxon_results[wilcoxon_results$exp_name %in% edc_traits_log10 & wilcoxon_results$out_name %in% c("age_g_b","sex_b_rev","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg_f"),]
wilcoxon_results$exp_name <- factor(wilcoxon_results$exp_name,
                                    levels = paste0(c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                      "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                      "TCC","TCS",
                                                      "BPA","BPS","BPF"),"_log10"),
                                    labels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                               "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                               "TCC","TCS",
                                               "BPA","BPS","BPF"))
wilcoxon_results$out_name <- factor(wilcoxon_results$out_name, 
                                    levels = c("age_g_b","sex_b_rev","high_edu_b","smk1_b","drk1_b","paactive3_g_b","high_fruveg_f"), 
                                    labels = c("Older age (≥60 y)","Men","Higher educational attainment","Current smoking","Current drinking","Sufficient physical activity","Adequate fruits and vegetables intake"))
wilcoxon_results <- wilcoxon_results %>%
  arrange(out_name,exp_name)
# 以每个outcome表型为单位，校正exposure
wilcoxon_results <- wilcoxon_results %>%
  group_by(out_name) %>%  # 按outcome分组，校正exposure
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


df <- wilcoxon_results
# 当 Cliff's Delta > 0 时，表示第一组有更高的倾向获得更大的值。分析时第一组是0，第二组是1，此处统一方向加一个负号，表示Cliff's Delta > 0 时，表示第二组有更高的倾向获得更大的值。
df$cliff.delta <- -df$cliff.delta

df <- df[,c("out_name","exp_name","p","p_adj_bh","cliff.delta")]


colnames(df) <- c("Factor","Analyte","P value","BH-adjusted P","Effect size (Cliff's delta, two-sided Wilcoxon rank-sum test)")
openxlsx::write.xlsx(df,"tables/cov_edc_wilcoxon_20260511.xlsx")
#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ####
