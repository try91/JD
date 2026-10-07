library(data.table)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

#### 配色 ----
# 提取 RdBu 的 11 种颜色
rdbu_colors <- RColorBrewer::brewer.pal(11, "RdBu")
# 提取 BrBG 的 11 种颜色
brbg_colors <- RColorBrewer::brewer.pal(11, "BrBG")
# 设置11种sector颜色
sector_colors <- colorRampPalette(colors = c("#b71540", "#eb2f06", "#fa8231", "#fed330", "#26de81", "#45aaf2", "#cd84f1", "#7158e2"))(11)
# sector_colors <- c("#476066", sector_colors)
#### 配色 ####

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

#### 数据处理 (for circos heatmap & forest plot) ----
### 设置纳入图片的暴露和结局
## 暴露：EDC
exposure <- c("EDC_14","PFAS","PAE6","TC","BP1")
# 标准化名称
exposure_label <- c("EDCs (14)","PFAS (5)","PAEs (6)","Antimicrobials (2)","Bisphenols (1)")
exp_dat <- as.data.frame(exposure)

## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
outcome <- c("dm_incident_1014","ckd_incident_1014","cvd_incident_1021","cvd_incident_1014",
             "dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b",
             "bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
             "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
             "alt_b","ast_b","ggt_b","bia_b",
             "egfr_b","scr_b","ua_b", #"acr_b",
             "glu0_b","glu120_b","vhba1c_b",
             "ins0_b","ins120_b","homair_b","homab_b",
             "sbp_b","dbp_b","pr_b", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
             "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
             "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
             "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")
             # "diet_score_f","high_fruveg","high_fish_f","low_meat_f","low_ssb_f",
             # "smk1_b","drk1_b","paactive3_g_b","sum_met_b","sleept_f","sitduration_b","sittimet_b",
             # "high_edu_b","sex_b_rev","age_b")  # 重要，分析中(men=0,women=1)women为保护因素，但我希望图表中关系统一，别的结局都是正向的，因此这里使用sex_b_rev (women=0,men=1)
outcome_final <- outcome
out_dat <- as.data.frame(outcome_final)
# 标准化名称
outcome_label <- c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2021)","Incident CVD (2010-2014)",
                   "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT",
                   "BMI","Height","Weight","WHR","WC","HC",
                   "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                   "ALT","AST","GGT","Bile acid",
                   "eGFR","Serum creatinine","UA", #"UACR",
                   "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                   "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                   "SBP","DBP","PR", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
                   "FT3","FT4","TSH","TPOAb","TgAb",
                   "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                   "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")

# outcome_label <- c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2021)","Incident CVD (2010-2014)",
#                    "Diabetes (2010)","CKD (2010)","CVD (2010)","Obesity (2010)","Abdominal obesity (2010)","IR (2010)","Dyslipidemia (2010)","MetS (2010)","NAFLD (2010)","High UA (2010)","Hypertension (2010)","High CIMT (2010)",
#                    "BMI (2010)","Height (2010)","Weight (2010)","WHR (2010)","WC (2010)","HC (2010)",
#                    "TG (2010)","LDL-C (2010)","HDL-C (2010)","TC (2010)","ApoA-1 (2010)","ApoB (2010)","Non-HDL-C (2010)",
#                    "ALT (2010)","AST (2010)","GGT (2010)","Bile acid (2010)",
#                    "eGFR (2010)","Serum creatinine (2010)","UA (2010)", "UACR (2010)",
#                    "OGTT 0-h glucose (2010)","OGTT 2-h glucose (2010)","HbA1c (2010)",
#                    "OGTT 0-h insulin (2010)","OGTT 2-h insulin (2010)","HOMA-IR (2010)","HOMA-B (2010)",
#                    "SBP (2010)","DBP (2010)","PR (2010)", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
#                    "FT3 (2014)","FT4 (2014)","TSH (2014)","TPOAb (2014)","TgAb (2014)",
#                    "Hemoglobin (2010)","Platelet count (2010)","Eosinophil count (2014)","Lymphocyte count (2014)","Monocyte count (2014)","Neutrophil count (2014)",
#                    "NLR (2014)","LMR (2014)","PLR (2014)","SII (2014)","SIRI (2014)","WBC (2010)","Hs-CRP (2014)")
                   
## 构建有完整exposure和outcome的数据框 ##
out_dat1 <- out_dat
out_dat1$exposure <- exposure[1]
out_dat2 <- out_dat
out_dat2$exposure <- exposure[2]
out_dat3 <- out_dat
out_dat3$exposure <- exposure[3]
out_dat4 <- out_dat
out_dat4$exposure <- exposure[4]
out_dat5 <- out_dat
out_dat5$exposure <- exposure[5]
dat_exp_out_all <- rbind(out_dat1,out_dat2,out_dat3,out_dat4,out_dat5)
colnames(dat_exp_out_all) <- c("out", "exp")
## 构建有完整exposure和outcome的数据框 ##


# 设置人群名称
sample_name1 <- "phy_edc_temp" # 用总人群

# 读取数据 (总人群)
results_wqs_bin <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_bin_(",sample_name1,")_20260313.xlsx"))
# 读取数据 (总人群)
results_wqs_cont <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_cont_(",sample_name1,")_20260313.xlsx"))

results_wqs <- rbind(results_wqs_bin, results_wqs_cont) %>%
  filter(exp %in% exposure & out %in% outcome_final & cov == "wqs" & type == "Q2")

results_wqs_pos <- results_wqs[results_wqs$direction == "pos",] # 正向权重wqs分析结果
unique(results_wqs_pos$exp) # 缺BP1
unique(results_wqs_pos$out)
results_wqs_pos <- left_join(dat_exp_out_all,results_wqs_pos,by=c("exp", "out"))
results_wqs_pos$p <- ifelse(is.na(results_wqs_pos$p), 1, results_wqs_pos$p)

results_wqs_neg <- results_wqs[results_wqs$direction == "neg",] # 负向权重wqs分析结果
unique(results_wqs_neg$exp) # 缺BP1
unique(results_wqs_neg$out) # 缺whr_b
results_wqs_neg <- left_join(dat_exp_out_all,results_wqs_neg,by=c("exp", "out"))
results_wqs_neg$p <- ifelse(is.na(results_wqs_neg$p), 1, results_wqs_neg$p)

# FDR 校正
# 以每个exposure表型为单位，校正outcome
dat_wqs_pos <- results_wqs_pos %>%
  group_by(exp) %>%  # 按exposure分组，校正outcome
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()

dat_wqs_neg <- results_wqs_neg %>%
  group_by(exp) %>%  # 按exposure分组，校正outcome
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()


# EDC 权重数据整理(正)
dat_wqs_pos$numb_pos_weight <- rowSums(dat_wqs_pos[, 13:26] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
dat_wqs_pos$numb_neg_weight <- rowSums(dat_wqs_pos[, 13:26] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量

dat_wqs_pos$weight_threshold_pos <- 1 / dat_wqs_pos$numb_pos_weight
dat_wqs_pos$weight_threshold_neg <- 1 / dat_wqs_pos$numb_neg_weight

dat_wqs_pos_edc14 <- dat_wqs_pos[dat_wqs_pos$exp == "EDC_14",]
dat_wqs_pos_edc14_long <- tidyr::gather(dat_wqs_pos_edc14, edc, weight, 13:26, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！


# EDC 权重数据整理(负)
dat_wqs_neg$numb_pos_weight <- rowSums(dat_wqs_neg[, 13:26] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
dat_wqs_neg$numb_neg_weight <- rowSums(dat_wqs_neg[, 13:26] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量

dat_wqs_neg$weight_threshold_pos <- 1 / dat_wqs_neg$numb_pos_weight
dat_wqs_neg$weight_threshold_neg <- 1 / dat_wqs_neg$numb_neg_weight

dat_wqs_neg_edc14 <- dat_wqs_neg[dat_wqs_neg$exp == "EDC_14",]
dat_wqs_neg_edc14_long <- tidyr::gather(dat_wqs_neg_edc14, edc, weight, 13:26, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！


# Prepare data for plot
dat_wqs_pos_for_plot <- dat_wqs_pos
dat_wqs_pos_for_plot$exp <- factor(dat_wqs_pos_for_plot$exp, levels = exposure, labels = exposure_label)
dat_wqs_pos_for_plot$out <- factor(dat_wqs_pos_for_plot$out, levels = outcome_final, labels = outcome_label)
dat_wqs_pos_for_plot <- dat_wqs_pos_for_plot %>% arrange(out,exp)

dat_wqs_neg_for_plot <- dat_wqs_neg
dat_wqs_neg_for_plot$exp <- factor(dat_wqs_neg_for_plot$exp, levels = exposure, labels = exposure_label)
dat_wqs_neg_for_plot$out <- factor(dat_wqs_neg_for_plot$out, levels = outcome_final, labels = outcome_label)
dat_wqs_neg_for_plot <- dat_wqs_neg_for_plot %>% arrange(out,exp)
#### 数据处理 (for circos heatmap & forest plot) ####


############################################# 环形热图 (qgcomp) #############################################
# 删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]
dat_wqs_pos_for_plot1 <- dat_wqs_pos_for_plot[!dat_wqs_pos_for_plot$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2021)","Incident CVD (2010-2014)"),]
# 提取数据做EDC和表型关联性环状热图
dat_wqs_corr <- dat_wqs_pos_for_plot1[,c(2,1,9:12,27)]
# 提取数据做EDC权重环状热图
dat_wqs_weight <- dat_wqs_pos_for_plot1[dat_wqs_pos_for_plot1$exp == "EDCs (14)",c(2,1,13:26,27)]
dat_wqs_weight[dat_wqs_weight$p_adj_bh >= 0.05, c(3:16)] <- 0
# 正权重
direction <- "pos-weight"
# 提取 BrBG 的 11 种颜色
brbg_colors <- RColorBrewer::brewer.pal(11, "BrBG")
brbg_colors <- colorRampPalette(colors = brbg_colors[1:6])(11)
# # 查看颜色梯度
# scales::show_col(brbg_colors)

# # 删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]
# dat_wqs_neg_for_plot1 <- dat_wqs_neg_for_plot[!dat_wqs_neg_for_plot$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2021)","Incident CVD (2010-2014)"),]
# # 提取数据做EDC和表型关联性环状热图
# dat_wqs_corr <- dat_wqs_neg_for_plot1[,c(2,1,9:12,27)]
# # 提取数据做EDC权重环状热图
# dat_wqs_weight <- dat_wqs_neg_for_plot1[dat_wqs_neg_for_plot1$exp == "EDCs (14)",c(2,1,13:26,27)]
# dat_wqs_weight[dat_wqs_weight$p_adj_bh >= 0.05, c(3:16)] <- 0
# # 负权重
# direction <- "neg-weight"
# # 提取 BrBG 的 11 种颜色
# brbg_colors <- RColorBrewer::brewer.pal(11, "BrBG")
# brbg_colors <- colorRampPalette(colors = brbg_colors[11:6])(11)
# # # 查看颜色梯度
# # scales::show_col(brbg_colors)

#### circos heatmap (删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]) ----
## WQS 关联性
# 转换为wide格式data
dat_wqs_corr_3cols <- dat_wqs_corr
# dat_wqs_corr_3cols$z <- dat_wqs_corr_3cols$estimate
dat_wqs_corr_3cols <- dat_wqs_corr_3cols[,c("exp","out","z")]
# 由于有过大值和过小值，因此对原本数值再进行一次转换，减少数值间差异
dat_wqs_corr_3cols$z <- asinh(dat_wqs_corr_3cols$z) # 反双曲正弦变换（asinh）， 适用于正负值混合的数据，尤其适合处理尾部极端值，效果类似对数变换但对零值更平滑。
# 由于有过大值和过小值，因此对原本数值再进行一次log转换，减少数值间差异
dat1 <- tidyr::pivot_wider(dat_wqs_corr_3cols, names_from = exp, values_from = z) # 转换为行为“out”列为“exp”
dat1 <- data.frame(dat1, check.names = FALSE)  # 禁止自动修改列名
rownames(dat1) <- dat1$out
dat1 <- dat1[,-1]
dat1_mat <- as.matrix(dat1)

## WQS 关联性(P值)
dat_wqs_p_3cols <- dat_wqs_corr[,c("exp","out","p_adj_bh")]
dat1_p <- tidyr::pivot_wider(dat_wqs_p_3cols, names_from = exp, values_from = p_adj_bh) # 转换为行为“out”列为“exp”
dat1_p <- data.frame(dat1_p, check.names = FALSE)  # 禁止自动修改列名
rownames(dat1_p) <- dat1_p$out
# dat1_p <- as.character(dat1_p)
dat1_p <- dat1_p[,-1]
dat1_p_mat <- as.matrix(dat1_p)
# dat1_p_mat[dat1_p_mat < 0.05] <- "*"
# dat1_p_mat[dat1_p_mat != "*"] <- " "
dat1_mat[dat1_p_mat >= 0.05] <- NA

# 把weight数据转换为long data排序
dat_wqs_weight2 <- dat_wqs_weight[,c(2:16)]
dat_wqs_weight2_long <- tidyr::gather(dat_wqs_weight2, edc, weight, 2:15, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
dat_wqs_weight2_long$edc <- factor(dat_wqs_weight2_long$edc, 
                                   levels = edc_traits_log10,
                                   labels = edc_traits)
dat_wqs_weight2_long <- dat_wqs_weight2_long %>%
  arrange(edc)
# # 由于有过大值和过小值，因此对原本数值再进行一次转换，减少数值间差异
# dat_wqs_weight2_long$weight <- asinh(dat_wqs_weight2_long$weight) # 反双曲正弦变换（asinh）， 适用于正负值混合的数据，尤其适合处理尾部极端值，效果类似对数变换但对零值更平滑。
# # 由于有过大值和过小值，因此对原本数值再进行一次log转换，减少数值间差异
dat2 <- tidyr::pivot_wider(dat_wqs_weight2_long, names_from = edc, values_from = weight) # 转换为行为“out”列为“edc”
dat2 <- dat2[,-1]
dat2 <- data.frame(dat2, check.names = FALSE)  # 禁止自动修改列名
rownames(dat2) <- dat_wqs_weight2$out
dat2_mat <- as.matrix(dat2)

# 设置表型分组
phenotypes <- c("Disorder and disease",
                "Body measurement","Lipid and lipoprotein",
                "Liver function","Kidney function",
                "Glucose metabolism","Insulin metabolism","Blood pressure",
                "Thyroid function","Hematological trait","Inflammation")

split <- factor(c(rep(phenotypes[1],12),
                  rep(phenotypes[2],6),rep(phenotypes[3],7),
                  rep(phenotypes[4],4),rep(phenotypes[5],3),
                  rep(phenotypes[6],3),rep(phenotypes[7],4),
                  rep(phenotypes[8],3),rep(phenotypes[9],5),
                  rep(phenotypes[10],6),rep(phenotypes[11],7)),
                
                levels = phenotypes)


### 作图 ###
pdf(paste0("figures/supplementary_figures/circos_heatmap_edc_outcomes_",direction,"_20260714.pdf"), width = 9, height = 9)

circos.par(start.degree = 90, 
           points.overflow.warning = FALSE,
           # cell.padding = c(0.02, 0, 0.02, 0), # 控制每个热图单元格（cell）内部的边距（padding）
           track.margin = c(0.008, 0.005),
           gap.degree = c(rep(2,10),90)) # change to 1 to number of sector - 1

# 设置外圈颜色
corr_min <- min(dat1_mat,na.rm = TRUE)
corr_max <- max(dat1_mat,na.rm = TRUE)
# col_corr = colorRamp2(c(corr_min, 0, corr_max), c("#1082bd","white","#E63E56")) #外圈beta的取值范围
# col_corr = colorRamp2(c(corr_min, 0, corr_max), c("#369acc","#FFFFFF","#de324c")) #外圈beta的取值范围
# 创建颜色映射函数，均匀分配断点到 RdBu 颜色
col_corr <- colorRamp2(
  breaks = seq(corr_min, corr_max, length.out = 11), 
  colors = rev(rdbu_colors)
)

# 设置内圈颜色
weight_min <- min(dat2_mat)
weight_max <- max(dat2_mat)
# col_weight = colorRamp2(c(weight_min, 0, weight_max), c("#39A880","white","#EDAD3D")) #内圈beta的取值范围
# col_weight = colorRamp2(c(weight_min, 0, weight_max), c("#72b043","white","#f8cc1b")) #内圈beta的取值范围
# 创建颜色映射函数，均匀分配断点到 RdBu 颜色
col_weight <- colorRamp2(
  breaks = seq(weight_min, weight_max, length.out = 11), 
  colors = rev(brbg_colors)
)


# 1. 绘制第一圈热图
circos.heatmap(dat1_mat, cluster = FALSE, split = split, 
               col = col_corr, na.col = "#EEEAE7", 
               track.height = 0.18, 
               cell.border = "white", cell.lwd = 0.4,
               # bg.border = "grey", bg.lwd = 1, bg.lty = 2,
               rownames.cex = 0.9, rownames.side = "outside") # 外圈名称字体大小设置
# 2. 在第一圈热图第一个section左侧添加行名
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(dat1_mat)[5:1] # 编号反了，需要调整所以是5 to 1
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[1], n) + convert_x(-1, "mm"), # 左边界外1mm  # CELL_META$cell.xlim[2] (右侧添加)
                1:n - 0.5, 
                labels = cn,
                cex = 0.65, 
                adj = c(1, 0.5), # 右对齐、垂直居中
                facing = "downward",  # 强制垂直向下
                niceFacing = FALSE  # 禁用自动旋转
                )
  }
}, bg.border = NA)


# 3. 绘制第二圈热图
circos.heatmap(dat2_mat, cluster = FALSE, split = split, 
               col = col_weight, 
               track.height = 0.28,
               cell.border = "lightgrey", cell.lwd = 0.4
               # bg.border = "grey", bg.lwd = 1, bg.lty = 2
               )
# 4. 在第二圈热图第一个section左侧添加行名
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = colnames(dat2_mat)[14:1] # 编号反了，需要调整所以是14 to 1
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[1], n) + convert_x(-1, "mm"), # 左边界外0mm  # CELL_META$cell.xlim[2] (右侧添加)
                1:n - 0.5, 
                labels = cn,
                cex = 0.52, 
                adj = c(1, 0.5), # 右对齐、垂直居中
                facing = "downward",  # 强制垂直向下
                niceFacing = FALSE  # 禁用自动旋转
    )
  }
}, bg.border = NA)


# 5. 绘制第三圈sector分类
circos.track(ylim = c(0, 1), 
             track.height = 0.019,
             
             # panel.fun = function(x, y) {
             #   sector_index = CELL_META$sector.index
             #   xlim = CELL_META$xlim
             #   ylim = CELL_META$ylim
             #   circos.text(mean(xlim), mean(ylim), sector_index, cex = 0.7, col = "black",
             #               facing = "inside", 
             #               niceFacing = TRUE)
             # },
             
             bg.col = sector_colors, bg.border = "lightgrey", bg.lwd = 0.1)
# 6. 添加sector图例
legend(x = -0.6, y = 0.72, pch = 15, col = sector_colors, legend = phenotypes, 
       cex = 0.65,
       box.col = "white",
       ncol = 1, text.col = "black",
       title = " ", title.col = "black", title.adj = 0)


# 7. 创建外圈热图的图例 "ComplexHeatmap"包
min_neg_zscore <- round(min(dat1_mat, na.rm = TRUE),1)
max_neg_zscore <- round(max(dat1_mat[dat1_mat < 0], na.rm = TRUE),1)
min_pos_zscore <- round(min(dat1_mat[dat1_mat > 0], na.rm = TRUE),1)
max_pos_zscore <- round(max(dat1_mat, na.rm = TRUE),1)
legend_corr <- Legend(
  title = "Z", 
  col_fun = col_corr, 
  at = c(min_neg_zscore, max_neg_zscore, 0, min_pos_zscore, max_pos_zscore),  # 自定义刻度值
  title_gp = gpar(fontsize = 6),  # 标题字体大小
  labels_gp = gpar(fontsize = 6), # 刻度标签字体大小
  grid_width = unit(2.5, "mm"),     # 色块宽度
  legend_height = unit(5, "mm")     # 图例高度
)
# 8. 创建内圈热图的图例 "ComplexHeatmap"包
min_weight <- round(min(dat_wqs_weight2_long$weight),1)
max_weight <- round(max(dat_wqs_weight2_long$weight),1)
legend_weight <- Legend(
  title = "Weight", 
  col_fun = col_weight, 
  at = c(min_weight, round(min_weight/2,1), 0, round(max_weight/2,1), max_weight),  # 自定义刻度值
  title_gp = gpar(fontsize = 6),  # 标题字体大小
  labels_gp = gpar(fontsize = 6), # 刻度标签字体大小
  grid_width = unit(2.5, "mm"),     # 色块宽度
  legend_height = unit(5, "mm")     # 图例高度
)
# 9. 合并图例
combined_legend <- packLegend(legend_corr, legend_weight, 
                              direction = "horizontal")
# 10. 绘制图例（在circos.clear()之前调用）
grid.draw(combined_legend)


circos.clear()

dev.off()

#### circos heatmap (删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]) ####

# 逆变换恢复原始数据（反双曲正弦变换 asinh）
sinh(c(3.8,1.5,-1.6,-2.3))
sinh(c(2.3,1.7,-1.7,-3.5))

