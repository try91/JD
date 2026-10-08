library(data.table)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

#### 配色 ----
# 提取 RdBu 的 11 种颜色
rdbu_colors <- RColorBrewer::brewer.pal(11, "RdBu")
# 提取 BrBG 的 11 种颜色
brbg_colors <- RColorBrewer::brewer.pal(11, "BrBG")
brbg_colors <- colorRampPalette(colors = brbg_colors)(41) #内圈细分为41份 (1~0.5, 0, -0.5~-1)
# 构建权重图颜色
weight_colors <- colorRampPalette(colors = c("#82589F","white","#1289A7"))(41) #内圈细分为41份 (1~0.5, 0, -0.5~-1)
# 设置11种sector颜色
sector_colors <- colorRampPalette(colors = c("#b71540", "#eb2f06", "#fa8231", "#fed330", "#26de81", "#45aaf2", "#cd84f1", "#7158e2"))(11)
#### 配色 ####

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

#### 数据处理 (for circos heatmap & forest plot) ----
### 设置纳入图片的暴露和结局
## 暴露：EDC
exposure <- c("EDC_14","EDC_pos","EDC_neg","PFAS","PAE_6","TC","BP_1")
# 标准化名称
exposure_label <- c("EDCs (14)","EDCs (positive weight)","EDCs (negative weight)","PFAS (5)","PAEs (6)","Antimicrobials (2)","Bisphenols (1)")

## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
outcome <- c("dm_incident_1014","ckd_incident_1014","cvd_incident_1021","cvd_incident_1014",
             "dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b",
             "bmi_b","height_b","weight_b","whr_b","wc_b","hc_b",
             "tg_b","ldl_b","hdl_b","chol_b","apoa_b","apob_b","nonhdl_b",
             "alt_b","ast_b","ggt_b","bia_b",
             "egfr_b","scr_b","ua_b",
             "glu0_b","glu120_b","vhba1c_b",
             "ins0_b","ins120_b","homair_b","homab_b",
             "sbp_b","dbp_b","pr_b",
             "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
             "hgb_b","plt_b","eos_f","lym_f","mon_f","neu_f",
             "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_b","crp_f")
outcome_final <- outcome
# 标准化名称
outcome_label <- c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2021)","Incident CVD (2010-2014)",
                   "Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT",
                   "BMI","Height","Weight","WHR","WC","HC",
                   "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                   "ALT","AST","GGT","Bile acid",
                   "eGFR","Serum creatinine","UA",
                   "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                   "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                   "SBP","DBP","PR",
                   "FT3","FT4","TSH","TPOAb","TgAb",
                   "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                   "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")

# 设置人群名称
sample_name1 <- "phy_edc_temp" # 用总人群

# 读取数据 (总人群)
results_qg_all1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",sample_name1,").xlsx"))
# 读取数据 (总人群)
results_qg_edc1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(",sample_name1,").xlsx"))

results_qg <- rbind(results_qg_all1, results_qg_edc1) %>%
  filter(exp %in% exposure & out %in% outcome_final)
results_qg <- results_qg[!(results_qg$out %in% c("cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014") & results_qg$method == "qgcomp bin"),] # 删除"cvd_incident_1021","ckd_incident_1014","dm_incident_1014"的logistic结果
colnames(results_qg)[9:22] <- gsub("_log10","",colnames(results_qg)[9:22])

results_qg$z <- results_qg$estimate/results_qg$se

# FDR 校正
# 以每个exposure表型为单位，校正outcome
dat <- results_qg %>%
  group_by(exp) %>%  # 按exposure分组，校正outcome
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()

# EDC 权重数据整理
dat$numb_pos_weight <- rowSums(dat[, 9:22] > 0, na.rm = TRUE)# 计算每一行第9列到第22列中正数的数量
dat$numb_neg_weight <- rowSums(dat[, 9:22] < 0, na.rm = TRUE)# 计算每一行第9列到第22列中负数的数量

dat$weight_threshold_pos <- 1 / dat$numb_pos_weight
dat$weight_threshold_neg <- 1 / dat$numb_neg_weight

dat_edc <- dat[dat$exp == "EDC_14",]
dat_edc_19_long_qgcomp <- tidyr::gather(dat_edc, edc, weight, 9:22, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！

dat_qgcomp <- dat
dat_qgcomp$exp <- factor(dat_qgcomp$exp, levels = exposure, labels = exposure_label)
dat_qgcomp$out <- factor(dat_qgcomp$out, levels = outcome_final, labels = outcome_label)
dat_qgcomp <- dat_qgcomp %>% arrange(out,exp)
#### 数据处理 (for circos heatmap & forest plot) ####


############################################# 环形热图 (qgcomp) #############################################
#### circos heatmap (删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]) ----
# 提取数据做EDC和表型关联性环状热图
dat_qgcomp_corr <- dat_qgcomp[,c(1:8,23,24)]
dat_qgcomp_corr <- dat_qgcomp_corr[!dat_qgcomp_corr$out %in% c("Incident CVD (2010-2021)","Incident CVD (2010-2014)","Incident CKD (2010-2014)","Incident diabetes (2010-2014)") & !dat_qgcomp_corr$exp %in% c("EDCs (positive weight)","EDCs (negative weight)"),]
# 提取数据做EDC权重环状热图
dat_qgcomp_weight <- dat_qgcomp[dat_qgcomp$exp == "EDCs (14)",c(2,9:22,25:28)]
dat_qgcomp_weight <- dat_qgcomp_weight[!dat_qgcomp_weight$out %in% c("Incident CVD (2010-2021)","Incident CVD (2010-2014)","Incident CKD (2010-2014)","Incident diabetes (2010-2014)"),]

## QGCOMP 关联性
# 转换为wide格式data
dat_qgcomp_corr_3cols <- dat_qgcomp_corr
# dat_qgcomp_corr_3cols$z <- dat_qgcomp_corr_3cols$estimate
dat_qgcomp_corr_3cols <- dat_qgcomp_corr_3cols[,c("exp","out","z")]
# 由于有过大值和过小值，因此对原本数值再进行一次转换，减少数值间差异
dat_qgcomp_corr_3cols$z <- asinh(dat_qgcomp_corr_3cols$z) # 反双曲正弦变换（asinh）， 适用于正负值混合的数据，尤其适合处理尾部极端值，效果类似对数变换但对零值更平滑。
# 由于有过大值和过小值，因此对原本数值再进行一次log转换，减少数值间差异
dat1 <- tidyr::pivot_wider(dat_qgcomp_corr_3cols, names_from = exp, values_from = z) # 转换为行为“out”列为“exp”
dat1 <- data.frame(dat1, check.names = FALSE)  # 禁止自动修改列名
rownames(dat1) <- dat1$out
dat1 <- dat1[,-1]
dat1_mat <- as.matrix(dat1)

## QGCOMP 关联性(P值)
dat_qgcomp_p_3cols <- dat_qgcomp_corr[,c("exp","out","p_adj_bh")]
dat1_p <- tidyr::pivot_wider(dat_qgcomp_p_3cols, names_from = exp, values_from = p_adj_bh) # 转换为行为“out”列为“exp”
dat1_p <- data.frame(dat1_p, check.names = FALSE)  # 禁止自动修改列名
rownames(dat1_p) <- dat1_p$out
# dat1_p <- as.character(dat1_p)
dat1_p <- dat1_p[,-1]
dat1_p_mat <- as.matrix(dat1_p)
# dat1_p_mat[dat1_p_mat < 0.05] <- "*"
# dat1_p_mat[dat1_p_mat != "*"] <- " "
dat1_mat[dat1_p_mat >= 0.05] <- NA
  
# 把weight数据转换为long data排序
dat_qgcomp_weight2 <- dat_qgcomp_weight[,c(1:15)]
dat_qgcomp_weight2_long <- tidyr::gather(dat_qgcomp_weight2, edc, weight, 2:15, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
dat_qgcomp_weight2_long$edc <- factor(dat_qgcomp_weight2_long$edc, levels = c("PFOS","PFOA","PFNA","PFDA","PFHxS",
                                                                              "MEHP","MECPP","MEHHP","MEP","MEOHP","MiBP","MnBP","MCPP","MBzP",
                                                                              "TCC","TCS",
                                                                              "BPA","BPS","BPF"))
dat_qgcomp_weight2_long <- dat_qgcomp_weight2_long %>%
  arrange(edc)
# # 由于有过大值和过小值，因此对原本数值再进行一次转换，减少数值间差异
# dat_qgcomp_weight2_long$weight <- asinh(dat_qgcomp_weight2_long$weight) # 反双曲正弦变换（asinh）， 适用于正负值混合的数据，尤其适合处理尾部极端值，效果类似对数变换但对零值更平滑。
# # 由于有过大值和过小值，因此对原本数值再进行一次log转换，减少数值间差异
dat2 <- tidyr::pivot_wider(dat_qgcomp_weight2_long, names_from = edc, values_from = weight) # 转换为行为“out”列为“edc”
dat2 <- dat2[,-1]
dat2 <- data.frame(dat2, check.names = FALSE)  # 禁止自动修改列名
rownames(dat2) <- dat_qgcomp_weight2$out
dat2_mat <- as.matrix(dat2)

# 设置表型分组
phenotypes <- c("Disorder and disease",
                "Body measurement","Lipid and lipoprotein",
                "Liver function","Kidney function",
                "Glucose metabolism","Insulin metabolism","Blood pressure",
                "Thyroid function","Hematological trait","Inflammation")
                
phenotypes[1]
split <- factor(c(rep(phenotypes[1],12),
                  rep(phenotypes[2],6),rep(phenotypes[3],7),
                  rep(phenotypes[4],4),rep(phenotypes[5],3),
                  rep(phenotypes[6],3),rep(phenotypes[7],4),
                  rep(phenotypes[8],3),rep(phenotypes[9],5),
                  rep(phenotypes[10],6),rep(phenotypes[11],7)),
                levels = phenotypes)


### 作图 ###
pdf("figures/main_figures/(fig3a)_circos_heatmap_edc_outcomes.pdf", width = 9, height = 9)

circos.par(start.degree = 90, 
           points.overflow.warning = FALSE,
           track.margin = c(0.008, 0.005),
           gap.degree = c(rep(2,10),90)) # change to 1 to number of sector - 1

# 设置外圈颜色
corr_min <- min(dat1_mat,na.rm = TRUE)
corr_max <- max(dat1_mat,na.rm = TRUE)
# 创建颜色映射函数，均匀分配断点到 RdBu 颜色
col_corr <- colorRamp2(
  breaks = seq(corr_min, corr_max, length.out = 11), 
  colors = rev(rdbu_colors)
)

# 设置内圈颜色
weight_min <- min(dat2_mat) # -0.846
weight_max <- max(dat2_mat) # 0.770
# 创建颜色映射函数，均匀分配断点到 RdBu 颜色
col_weight <- colorRamp2(
  breaks = seq(weight_min, weight_max, length.out = 33), 
  colors = rev(brbg_colors[6:38])  # 0.77~-0.85
)


# 1. 绘制第一圈热图
circos.heatmap(dat1_mat, cluster = FALSE, split = split, 
               col = col_corr, na.col = "#EEEAE7", 
               track.height = 0.19, 
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
               track.height = 0.30,
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
min_weight <- round(min(dat_qgcomp_weight2_long$weight),1)
max_weight <- round(max(dat_qgcomp_weight2_long$weight),1)
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


# # 逆变换恢复原始数据（反双曲正弦变换 asinh）
# sinh(c(3.7,1.5,-1.5,-3.6))  

#### circos heatmap (删除 4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]) ####


############################################# 森林图 (qgcomp) #############################################
#### 森林图+权重图 (4 incidence outcome ["dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021"]) ----
dat_forest <- dat_qgcomp[dat_qgcomp$out %in% c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"),]
dat_forest$hr <- exp(dat_forest$estimate)
dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
dat_forest$out <- factor(dat_forest$out, levels = c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)"))
dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 

dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
dat_forest2 <- dat_forest[dat_forest$out == "Incident CVD (2010-2021)",]
dat_forest3 <- dat_forest[dat_forest$out == "Incident CKD (2010-2014)",]
dat_forest4 <- dat_forest[dat_forest$out == "Incident diabetes (2010-2014)",]

### 森林图
# Incident CVD (2010-2014)
{
  f_forest1 <- ggplot(data=dat_forest1, aes(x=estimate, y=exp)) +
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
    geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 5,
              position = position_dodge(width = 0.7)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.3)
    ) +
    
    scale_y_discrete(expand = c(0.1, 0.1)) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(
      title = "Incident CVD (2010-2014)",
      x = "HR (95% CI)") +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = "EDC groups",  # 修改图例标题
      values = c("#f1c40f", "#00b894", "#3498db", "#f55d78", "#5758BB", "#EE5A24", "#476066"),
      guide = guide_legend(reverse = TRUE)  # 逆序排列
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.position = "none")
}
# Incident CVD (2010-2021)
{
  f_forest2 <- ggplot(data=dat_forest2, aes(x=estimate, y=exp)) +
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
    geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 5,
              position = position_dodge(width = 0.7)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.3)
    ) +
    
    scale_y_discrete(expand = c(0.1, 0.1)) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(
      title = "Incident CVD (2010-2021)",
      x = "HR (95% CI)") +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = "EDC groups",  # 修改图例标题
      values = c("#f1c40f", "#00b894", "#3498db", "#f55d78", "#5758BB", "#EE5A24", "#476066"),
      guide = guide_legend(reverse = TRUE)  # 逆序排列
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.position = "none")
}
# Incident CKD (2010-2014)
{
  f_forest3 <- ggplot(data=dat_forest3, aes(x=estimate, y=exp)) +
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
    geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 5,
              position = position_dodge(width = 0.7)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.3)
    ) +
    
    scale_y_discrete(expand = c(0.1, 0.1)) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(
      title = "Incident CKD (2010-2014)",
      x = "HR (95% CI)") +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = "EDC groups",  # 修改图例标题
      values = c("#f1c40f", "#00b894", "#3498db", "#f55d78", "#5758BB", "#EE5A24", "#476066"),
      guide = guide_legend(reverse = TRUE)  # 逆序排列
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.position = "none")
}
# Incident DM (2010-2014)
{
  f_forest4 <- ggplot(data=dat_forest4, aes(x=estimate, y=exp)) +
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
    geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 5,
              position = position_dodge(width = 0.7)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.3)
    ) +
    
    scale_y_discrete(expand = c(0.1, 0.1)) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(
      title = "Incident diabetes (2010-2014)",
      x = "HR (95% CI)") +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = "EDC groups",  # 修改图例标题
      values = c("#f1c40f", "#00b894", "#3498db", "#f55d78", "#5758BB", "#EE5A24", "#476066"),
      guide = guide_legend(reverse = TRUE)  # 逆序排列
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.position = "none")
}
### 森林图


dat_weight <- dat_forest[dat_forest$exp == "EDCs (14)",]
dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 9:22, na.rm=TRUE, factor_key=TRUE)
dat_weight_long$weight_abs <- abs(dat_weight_long$weight)

dat_weight_long1 <- dat_weight_long[dat_weight_long$out == "Incident CVD (2010-2014)",]%>%
  arrange(rev(out),-weight_abs)
dat_weight_long1$edc <- factor(dat_weight_long1$edc, levels = rev(dat_weight_long1$edc))

dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "Incident CVD (2010-2021)",]%>%
  arrange(rev(out),-weight_abs)
dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))

dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "Incident CKD (2010-2014)",]%>%
  arrange(rev(out),-weight_abs)
dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))

dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Incident diabetes (2010-2014)",]%>%
  arrange(rev(out),-weight_abs)
dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))

### 权重图
# Incident CVD (2010-2014)
{
  f_weight1 <- ggplot(dat_weight_long1, aes(y = edc, x = weight, fill = weight)) +
    # 绘制分簇柱状图
    geom_bar(stat = "identity",
             width = 0.9,
             color = "black",
             linewidth = 0.2) +
    
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    
    scale_x_continuous(limits = c(-0.4, 0.4),
                       breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    
    # 可视化美化
    labs(
      title = " ",
      x = "Weight", y = "Outcomes") +
    
    scale_fill_gradientn(
      # colours = rev(weight_colors[c(13:33)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      colours = rev(weight_colors[c(7:41)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      limits = c(min(dat_weight_long$weight), max(dat_weight_long$weight)),
      name = "Weight"
    ) +
    
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.text = element_text(size = 13, colour = "black"),   # 调整legend文本大小
          legend.title = element_text(size = 15, colour = "black"),
          legend.position = "right"  # 图例位置（可选）
    )
}
# Incident CVD (2010-2021)
{
  f_weight2 <- ggplot(dat_weight_long2, aes(y = edc, x = weight, fill = weight)) +
    # 绘制分簇柱状图
    geom_bar(stat = "identity",
             width = 0.9,
             color = "black",
             linewidth = 0.2) +
    
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    
    scale_x_continuous(limits = c(-0.4, 0.4),
                       breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    
    # 可视化美化
    labs(
      title = " ",
      x = "Weight", y = "Outcomes") +
    
    scale_fill_gradientn(
      # colours = rev(weight_colors[c(13:33)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      colours = rev(weight_colors[c(7:41)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      limits = c(min(dat_weight_long$weight), max(dat_weight_long$weight)),
      name = "Weight"
    ) +
    
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.text = element_text(size = 13, colour = "black"),   # 调整legend文本大小
          legend.title = element_text(size = 15, colour = "black"),
          legend.position = "right"  # 图例位置（可选）
    )
}
# Incident CKD (2010-2014)
{
  f_weight3 <- ggplot(dat_weight_long3, aes(y = edc, x = weight, fill = weight)) +
    # 绘制分簇柱状图
    geom_bar(stat = "identity",
             width = 0.9,
             color = "black",
             linewidth = 0.2) +
    
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    
    scale_x_continuous(limits = c(-0.4, 0.4),
                       breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    
    # 可视化美化
    labs(
      title = " ",
      x = "Weight", y = "Outcomes") +
    
    scale_fill_gradientn(
      # colours = rev(weight_colors[c(13:33)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      colours = rev(weight_colors[c(7:41)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      limits = c(min(dat_weight_long$weight), max(dat_weight_long$weight)),
      name = "Weight"
    ) +
    
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.text = element_text(size = 13, colour = "black"),   # 调整legend文本大小
          legend.title = element_text(size = 15, colour = "black"),
          legend.position = "right"  # 图例位置（可选）
    )
}
# Incident DM (2010-2014)
{
  f_weight4 <- ggplot(dat_weight_long4, aes(y = edc, x = weight, fill = weight)) +
    # 绘制分簇柱状图
    geom_bar(stat = "identity",
             width = 0.9,
             color = "black",
             linewidth = 0.2) +
    
    geom_vline(xintercept = 0,
               linetype = "solid",
               linewidth = 0.5) +
    
    scale_x_continuous(limits = c(-0.7, 0.35),
                       breaks = c(-0.7, -0.35, 0, 0.35)) +
    
    # 可视化美化
    labs(
      title = " ",
      x = "Weight", y = "Outcomes") +
    
    scale_fill_gradientn(
      # colours = rev(weight_colors[c(13:33)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      colours = rev(weight_colors[c(7:41)]),  # 因为最小值到最大值对于-0.58~0.39，因此也从41种颜色中选取对应范围的颜色
      limits = c(min(dat_weight_long$weight), max(dat_weight_long$weight)),
      name = "Weight"
    ) +
    
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 14, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 14, colour = "black"),
          
          legend.text = element_text(size = 13, colour = "black"),   # 调整legend文本大小
          legend.title = element_text(size = 15, colour = "black"),
          legend.position = "right"  # 图例位置（可选）
    )
}


f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest1, f_forest2, nrow = 4, rel_heights = c(1,1,1,1))
f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight1, f_weight2, nrow = 4, rel_heights = c(1,1,1,1))
f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.2,0.4,1.2))
ggsave(f_fw, filename=paste0("figures/main_figures/(fig3b)_forest_weight_4incidence.pdf"), width = 10, height = 14, limitsize = FALSE)
#### 森林图+权重图 (4 incidence outcome ["cvd_incident_1021","cvd_incident_1014","ckd_incident_1014","dm_incident_1014"]) ####
