library(data.table)
library(dplyr)
library(ggplot2)
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
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp4_s_names_short <- mp4_s_names[!grepl("_GGB",mp4_s_names)] # 排除未分类的菌属（GGB）, 未分类菌种（SGB）先保留
mp4_g_names_short <- mp4_g_names[!grepl("_GGB",mp4_g_names)] # 排除未分类的菌属（GGB）
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp3_s_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_10%.txt")
mp3_s_names <- mp3_s_names[,1]
mp3_g_names <- read.table("jiading/sourceDataTaxon/mpa3/genus_names_mp3_10%.txt")
mp3_g_names <- mp3_g_names[,1]
# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)
mp4_g_zero <- paste0(mp4_g_names,"_zero") # 菌群MP4填补0值丰度 (属层面)

mp3_s_bin <- paste0(mp3_s_names,"_bin") # 菌群MP3出现与否的分类变量 (物种层面)
mp3_s_log10 <- paste0(mp3_s_names,"_log10") # 菌群MP3丰度的log10转换 (物种层面)
mp3_s_zero <- paste0(mp3_s_names,"_zero") # 菌群MP3填补0值丰度 (物种层面)

mp3_g_bin <- paste0(mp3_g_names,"_bin") # 菌群MP3出现与否的分类变量 (属层面)
mp3_g_log10 <- paste0(mp3_g_names,"_log10") # 菌群MP3丰度的log10转换 (属层面)
mp3_g_zero <- paste0(mp3_g_names,"_zero") # 菌群MP3填补0值丰度 (属层面)
# mp3中用于构建ma的菌的名称
mp3_ma_names <- read.table("jiading/sourceDataTaxon/mpa3/species_names_mp3_ma.txt")
mp3_ma_names <- mp3_ma_names$V1
mp3_ma_names_log10 <- paste0(mp3_ma_names,"_log10")
# microbial age (MA)
# mean(phy_edc$MA,na.rm = TRUE)

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
                    "ob_f","ob_b","abob_f","abob_b","dyslip_f","dyslip_b","hua_f","hua_b","ir_f","ir_b","mets_f","mets_b")
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

# 十类药物 (使用人数>20, 包括Statins)
med_cat10 <- c("med_dm1_f","med_dm2_f","med_dm3_f","med_dm4_f",
               "med_hbp1_f","med_hbp2_f","med_hbp3_6_f","med_hbp4_f","med_hbp5_f",
               "med_lip1_f")
med_cat10_label <- c("Sulfonylureas","Biguanides","Thiazolidinediones (TZDs)","Alpha-glucosidase inhibitors (AGIs)",
                     "Angiotensin II receptor blockers (ARBs)","Angiotensin-converting enzyme\ninhibitors (ACEIs)","Beta blockers","Calcium antagonists (CCBs)","Diuretics",
                     "Statins")


cov_traits <- c("age_f","sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")
cov_traits_cat <- c("sex_b_rev","smk1_f","drk1_f","high_edu_b","paactive3_g_f","high_fruveg","med_all7")

edc_traits_cat <- c(paste0(edc_traits,"_detected"),
                    paste0(edc_traits,"_median"),
                    paste0(edc_traits,"_quantile"),
                    "num_edc1",
                    "num_edc2",
                    "num_edc3")

# 总体分类变量名
cat_traits <- c(cov_traits_cat)
# 总体连续变量名
cont_traits <- c("age_f", phy_traits_cont, 
                 "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f",
                 "edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                 edc_traits_log10)


#### 读取Permanova分析结果 ----
permanova_results_med_mp3_g <- readxl::read_xlsx("results/permanova/medication10_permanova_beta_diversity_mp3_genus_20251017.xlsx")
permanova_results_med_mp4_g <- readxl::read_xlsx("results/permanova/medication10_permanova_beta_diversity_mp4_genus_20260313.xlsx")
permanova_results_med_mp4_s <- readxl::read_xlsx("results/permanova/medication10_permanova_beta_diversity_mp4_species_20260313.xlsx")

permanova_results_env_out_edc_mp4_s <- readxl::read_xlsx("results/permanova/permanova_edc_index-mp4_species-20260313.xlsx")
dat1 <- permanova_results_env_out_edc_mp4_s
dat1 <- dat1[dat1$trait %in% c(cat_traits,cont_traits),]
permanova_results_env_out_edc_mp4_g <- readxl::read_xlsx("results/permanova/permanova_edc_index-mp4_genus-20260313.xlsx")
dat2 <- permanova_results_env_out_edc_mp4_g
dat2 <- dat2[dat2$trait %in% c(cat_traits,cont_traits),]
#### 读取Permanova分析结果 ####

#### 筛选药物 ----
permanova_results_med_mp3_g <- permanova_results_med_mp3_g[permanova_results_med_mp3_g$cov == permanova_results_med_mp3_g$trait,]
permanova_results_med_mp3_g <- permanova_results_med_mp3_g[permanova_results_med_mp3_g$adj == "multivariable",]
permanova_results_med_mp3_g$p_adj_bh <- p.adjust(permanova_results_med_mp3_g$`Pr(>F)`, method = "BH", n = length(permanova_results_med_mp3_g$`Pr(>F)`))
permanova_results_med_mp3_g$sig_flag <- ifelse(permanova_results_med_mp3_g$p_adj_bh < 0.05, "TRUE", "FALSE")

permanova_results_med_mp4_g <- permanova_results_med_mp4_g[permanova_results_med_mp4_g$cov == permanova_results_med_mp4_g$trait,]
permanova_results_med_mp4_g <- permanova_results_med_mp4_g[permanova_results_med_mp4_g$adj == "multivariable",]
permanova_results_med_mp4_g$p_adj_bh <- p.adjust(permanova_results_med_mp4_g$`Pr(>F)`, method = "BH", n = length(permanova_results_med_mp4_g$`Pr(>F)`))
permanova_results_med_mp4_g$sig_flag <- ifelse(permanova_results_med_mp4_g$p_adj_bh < 0.05, "TRUE", "FALSE")

permanova_results_med_mp4_s <- permanova_results_med_mp4_s[permanova_results_med_mp4_s$cov == permanova_results_med_mp4_s$trait,]
permanova_results_med_mp4_s <- permanova_results_med_mp4_s[permanova_results_med_mp4_s$adj == "multivariable",]
permanova_results_med_mp4_s$p_adj_bh <- p.adjust(permanova_results_med_mp4_s$`Pr(>F)`, method = "BH", n = length(permanova_results_med_mp4_s$`Pr(>F)`))
permanova_results_med_mp4_s$sig_flag <- ifelse(permanova_results_med_mp4_s$p_adj_bh < 0.05, "TRUE", "FALSE")


permanova_results_med_mp3_g$trait <- factor(permanova_results_med_mp3_g$trait, levels = med_cat10, labels = med_cat10_label)
permanova_results_med_mp3_g <- permanova_results_med_mp3_g %>%
  arrange(-adjusted_R2)
permanova_results_med_mp3_g$trait <- factor(permanova_results_med_mp3_g$trait, levels = permanova_results_med_mp3_g$trait)

permanova_results_med_mp4_g$trait <- factor(permanova_results_med_mp4_g$trait, levels = med_cat10, labels = med_cat10_label)
permanova_results_med_mp4_g <- permanova_results_med_mp4_g %>%
  arrange(-adjusted_R2)
permanova_results_med_mp4_g$trait <- factor(permanova_results_med_mp4_g$trait, levels = permanova_results_med_mp4_g$trait)

permanova_results_med_mp4_s$trait <- factor(permanova_results_med_mp4_s$trait, levels = med_cat10, labels = med_cat10_label)
permanova_results_med_mp4_s <- permanova_results_med_mp4_s %>%
  arrange(-adjusted_R2)
permanova_results_med_mp4_s$trait <- factor(permanova_results_med_mp4_s$trait, levels = permanova_results_med_mp4_s$trait)
#### 筛选药物 ####

#### 整合Permanova分析结果 (MP4-Species Genus) ----
## MP4 Species ##
dat1 <- dat1[dat1$cov == dat1$trait,]
dat1_multi <- dat1[dat1$adj == "multivariable" & dat1$trait %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f",
                                                                   edc_traits6_log10),]
dat1_multi$p_adj_bh <- p.adjust(dat1_multi$`Pr(>F)`, method = "BH", n = length(dat1_multi$`Pr(>F)`))
dat1_multi$sig_flag <- ifelse(dat1_multi$p_adj_bh < 0.05, "TRUE", "FALSE")
dat1_multi$edc_cat <- ifelse(dat1_multi$trait %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"), "EDC Score", 
                             ifelse(dat1_multi$trait %in% edc_traits2_log10, "PFAS", 
                                    ifelse(dat1_multi$trait %in% edc_traits3_log10, "PAE", 
                                           ifelse(dat1_multi$trait %in% edc_traits5_log10, "TC", "BP"))))
dat1_multi$edc_cat <- factor(dat1_multi$edc_cat, levels = c("EDC Score","PFAS","PAE","TC","BP"))

dat1_multi <- dat1_multi %>%
  arrange(edc_cat, -adjusted_R2)
dat1_multi$trait <- factor(dat1_multi$trait, levels = dat1_multi$trait, labels = gsub("_log10","",dat1_multi$trait))


## MP4 Genus ##
dat2 <- dat2[dat2$cov == dat2$trait,]
dat2_multi <- dat2[dat2$adj == "multivariable" & dat2$trait %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f",
                                                                   edc_traits6_log10),]
dat2_multi$p_adj_bh <- p.adjust(dat2_multi$`Pr(>F)`, method = "BH", n = length(dat2_multi$`Pr(>F)`))
dat2_multi$sig_flag <- ifelse(dat2_multi$p_adj_bh < 0.05, "TRUE", "FALSE")
dat2_multi$edc_cat <- ifelse(dat2_multi$trait %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"), "EDC Score", 
                             ifelse(dat2_multi$trait %in% edc_traits2_log10, "PFAS", 
                                    ifelse(dat2_multi$trait %in% edc_traits3_log10, "PAE", 
                                           ifelse(dat2_multi$trait %in% edc_traits5_log10, "TC", "BP"))))
dat2_multi$edc_cat <- factor(dat2_multi$edc_cat, levels = c("EDC Score","PFAS","PAE","TC","BP"))

dat2_multi <- dat2_multi %>%
  arrange(edc_cat, -adjusted_R2)
dat2_multi$trait <- factor(dat2_multi$trait, levels = dat2_multi$trait, labels = gsub("_log10","",dat2_multi$trait))
#### 整合Permanova分析结果 (MP4-Species Genus) ####


#### 作图柱状图 ----
### 构建用药、EDC菌属差异作图函数 ###
barchart_function1 <- function(DAT, TITLE){
  
  # 调整adjusted_R2为负数的情况
  DAT$adjusted_R2 <- ifelse(DAT$adjusted_R2 < 0, 0, DAT$adjusted_R2)
  
  # 显著添加星号
  DAT <- DAT %>%
    mutate(text = case_when(
      is.na(p_adj_bh) ~ paste(" "), 
      !is.na(p_adj_bh) & p_adj_bh < 5e-2 ~ paste("*")))
  
  barchart <- ggplot(DAT, aes(x = trait, y = adjusted_R2)) +
    geom_hline(yintercept = 0,
               linetype = "solid",
               linewidth = 0.2) +
    
    geom_col(aes(fill = sig_flag), width = 0.6) +  # 关键修改：根据sig_flag分组填充颜色
    scale_fill_manual(
      # values = c("TRUE" = "#845ec2", "FALSE" = "#00c9a7"),
      values = c("TRUE" = "#3c6382", "FALSE" = "#7f7f7f"),
      # 直接在scale_color_manual中设置标签
      labels = c("TRUE" = "BH-adjusted P<0.05",
                 "FALSE" = "Non-significant"),
      # 通过breaks参数明确指定图例顺序
      breaks = c("TRUE", "FALSE")
      ) +
    
    coord_flip() +  # 关键：翻转坐标轴
    labs(x = "Traits", y = "Explained variance (adjusted R2)") +
    ggtitle(paste0(TITLE)) +
    
    theme_classic() +
    theme(
      # plot.margin = margin(5, 5, 5, 5, "mm"),
      # panel.border = element_rect(color = "black", fill = NA, size = 0.5), # 面板区域边框
      
      plot.title = element_text(size = 14, colour = "black"),  # 关键修改：设置标题大小和加粗
      axis.ticks = element_line(size = 0.5, colour = "black"),  # 调整 X 轴刻度线, 线宽(默认0.5)
      
      axis.ticks.length = unit(1, "mm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      axis.title.x = element_text(size = 12, colour = "black", margin = margin(t = 3, unit = "mm")),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(size = 12, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 12, colour = "black"),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 8, colour = "black"), # 调整legend文本大小
      # 图例键大小
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(4, "mm")
      # legend.position = "none"  # 这行代码删除了图例
    ) +
    scale_x_discrete(limits = rev(levels(DAT$trait)))  # 反转类别顺序
  
  return(barchart)
}


# 用药菌属差异 #
# f1_1 <- barchart_function1(permanova_results_med_mp3_g,"PERMANOVA Drug (genus level-MP3)")
f1_2 <- barchart_function1(permanova_results_med_mp4_g,"PERMANOVA Drug (genus level")
# f1_3 <- barchart_function1(permanova_results_med_mp4_s,"PERMANOVA Drug (species level)")
# f1 <- cowplot::plot_grid(f1_1, f1_2, f1_3, 
#                          ncol = 3,
#                          rel_widths = c(1, 1, 1))
# ggsave(paste0("C:/TWang/DLiu/EDC_Micro/figures/main_figures/barchart_permanova_drug_all_(phy_edc_temp0_3)_20251216.pdf"), 
#        f1, width = 90, height = 18, limitsize = FALSE)
ggsave(paste0("C:/TWang/DLiu/EDC_Micro/figures/main_figures/barchart_permanova_drug_mp4_g_(phy_edc_temp0_3)_20260313.pdf"), 
       f1_2, width = 10, height = 6, limitsize = FALSE)

# EDC菌属差异 #
f2 <- barchart_function1(dat2_multi,"")
ggsave(paste0("C:/TWang/DLiu/EDC_Micro/figures/main_figures/barchart_permanova_edc_mp4_g_(phy_edc_temp0_3)_20260313.pdf"), 
       f2, width = 6, height = 8, limitsize = FALSE)
#### 作图柱状图 ####
