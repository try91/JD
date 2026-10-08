library(data.table)
library(dplyr)
library(ggplot2)
library(corrplot)

setwd("your_file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

#### 配色 ----
PFAS_colors <- colorRampPalette(c("#f55d78","#FFFFFF"))(50)[c(1,8,15,22,29)]
PAE_colors <- colorRampPalette(c("#3498db","#FFFFFF"))(50)[c(1,6,11,16,21,26,31,36,41)]
BP_colors <- colorRampPalette(c("#f1c40f","#FFFFFF"))(50)[c(1,18,35)]
TC_colors <- colorRampPalette(c("#00b894","#FFFFFF"))(50)[c(1,21)]
my_colors <- c(PFAS_colors,PAE_colors,BP_colors,TC_colors)
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


############################################# EDC相关性热图 #############################################
#### EDC-EDC 相关性热图 (Partial Spearman) 填补后 ----
### 读取EDC-EDC关联性数据
p_spearman_results <- readRDS("results/correlations/spearman/partial_spearman_results_edc-edc_mp4_out.rds")
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

dat <- p_spearman_results
data_min <- min(dat$estimate, na.rm = TRUE)
data_max <- max(dat$estimate, na.rm = TRUE)


# EDC和EDC关联heatmap
heatmap_function1 <- function(DAT, TITLE){
  
  filtered_dat <- DAT
  
  # 色块显示内容
  filtered_dat <- filtered_dat %>%
    mutate(text = case_when(
      is.na(p_adj_bh) ~ paste(" "),
      !is.na(p_adj_bh) & p_adj_bh < 5e-2 ~ paste("*")))
  
  heatmap <- ggplot(filtered_dat, aes(exp_name,out_name)) +
    
    geom_tile(aes(fill = estimate), colour = "darkgrey", linewidth = 0.1) +
    
    # scale_fill_gradient2(low = "#0984e3", mid = "white", high = "#d63031", na.value = "#EDEDED") +
    scale_fill_gradientn(
      colours = c("#0984e3", "white", "#d63031"), # 定义颜色序列
      values = scales::rescale(c(min(filtered_dat$estimate,na.rm = TRUE), 0, max(filtered_dat$estimate,na.rm = TRUE))), # 关键：将颜色精确对应到数据点
      na.value = "#EDEDED"
      # guide = guide_colorbar(
      #   nbin = 400,  # 增加分段数使过渡更平滑
      #   ticks.at = c(min(filtered_dat$estimate,na.rm = TRUE), max(filtered_dat$estimate,na.rm = TRUE))  # 强制在极值处显示刻度
      # )
    ) +
    
    geom_text(aes(label=text),col ="black",size = 14, vjust = 0.85) +
    
    labs(fill = paste0("Spearman's rho")) +   # 修改 legend 内容
    # ggtitle(paste0(TITLE)) +
    
    theme_minimal() + # 不要背景
    theme(
      plot.title = element_text(size = 25),  # 设置标题大小和加粗
      axis.title.x=element_blank(), # 去掉 x轴title
      axis.title.y=element_blank(), # 去掉 y轴title
      axis.ticks=element_blank(), # 去掉刻度线
      panel.grid =element_blank(),  # 删去网格线
      axis.text.x = element_text(angle = 45, hjust = 0, size = 21, color ="black"), # 调整x轴文字
      axis.text.y = element_text(size = 21, color ="black"), #调整y轴文字
      
      # 添加以下 legend 设置
      legend.title = element_text(size = 15),  # 设置legend标题大小
      legend.text = element_text(size = 15),   # 设置legend文本大小
      
      # 设置color panel的长宽
      legend.key.width = unit(6, "mm"),   # 设置颜色条宽度
      legend.key.height = unit(6, "mm")   # 设置颜色条高度
    ) +
    
    scale_x_discrete(position = "top") + # 将 X 轴放置在最上面
    scale_y_discrete(limits = rev)  # 添加以下行以实现Y轴逆序
  
  return(heatmap)
}

# EDC和EDC关联 (作图)
f_heatmap1 <- heatmap_function1(dat, "")
ggsave(f_heatmap1, filename=paste0("figures/main_figures/(fig2c)_heatmap_19edc-edc.pdf"), width = 15, height = 13, limitsize = FALSE)
#### EDC-EDC 相关性热图 (Partial Spearman) 填补后 ####

#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ----
### 读取分类组间EDC差异数据
wilcoxon_results <- readRDS("results/correlations/wilcoxon/wilcoxon_results_edc.rds")
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
                                    labels = c("Age>=60 (2010)","Men (2010)","Higher educational attainment (2010)","Current smoking (2010)","Current drinking (2010)","Sufficient physical activity (2010)","Adequate fruit and vegetable (2014)"))
wilcoxon_results <- wilcoxon_results %>%
  arrange(exp_name,out_name)
# 以每个outcome表型为单位，校正exposure
wilcoxon_results <- wilcoxon_results %>%
  group_by(out_name) %>%  # 按outcome分组，校正exposure
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()

dat <- wilcoxon_results


# EDC和COV关联heatmap
heatmap_function2 <- function(DAT, TITLE){
  
  filtered_dat <- DAT
  
  # 色块显示内容
  filtered_dat <- filtered_dat %>%
    mutate(text = case_when(
      is.na(p_adj_bh) ~ paste(" "),
      !is.na(p_adj_bh) & p_adj_bh < 5e-2 & p_adj_bh >= 1e-10 ~ paste("1"),
      !is.na(p_adj_bh) & p_adj_bh < 1e-10 & p_adj_bh >= 1e-20 ~ paste("2"),
      !is.na(p_adj_bh) & p_adj_bh < 1e-20 & p_adj_bh >= 1e-50 ~ paste("3"),
      !is.na(p_adj_bh) & p_adj_bh < 1e-50 ~ paste("4")))
  
  filtered_dat <- filtered_dat %>%
    mutate(text = case_when(
      p_adj_bh < 5e-2 & estimate > 0 ~ paste0("+",text),
      p_adj_bh < 5e-2 & estimate < 0 ~ paste0("-",text)))
  
  heatmap <- ggplot(filtered_dat, aes(exp_name,out_name)) +
    
    geom_tile(aes(fill = estimate), colour = "black", linewidth = 0.05) +
    
    scale_fill_gradientn(
      colours = c("#1dd1a1", "white", "#ff9f43"), # 定义颜色序列
      values = scales::rescale(c(min(filtered_dat$estimate), 0, max(filtered_dat$estimate))), # 关键：将颜色精确对应到数据点
      na.value = "#EDEDED"
    ) +
    
    geom_text(aes(label=text),col ="black",size = 5, vjust = 0.85) +
    # ggtitle(paste0(TITLE)) +
    
    theme_minimal() + # 不要背景
    theme(plot.title = element_text(size = 25),  # 设置标题大小和加粗
          axis.title.x=element_blank(), # 去掉 x轴title
          axis.title.y=element_blank(), # 去掉 y轴title
          axis.ticks=element_blank(), # 去掉刻度线
          panel.grid =element_blank(),  # 删去网格线
          axis.text.x = element_text(angle = 90, hjust = 0, size = 12), # 调整x轴文字
          axis.text.y = element_text(size = 13)) + #调整y轴文字
    labs(fill = paste0("Estimate")) +   # 修改 legend 内容
    scale_x_discrete(position = "top") + # 将 X 轴放置在最上面
    scale_y_discrete(limits = rev)  # 添加以下行以实现Y轴逆序
  
  return(heatmap)
}

# EDC和COV关联 (作图)
# 当 Cliff's Delta > 0 时，表示第一组有更高的倾向获得更大的值。分析时第一组是0，第二组是1，此处统一方向加一个负号，表示Cliff's Delta > 0 时，表示第二组有更高的倾向获得更大的值。
dat$cliff.delta <- -dat$cliff.delta
colnames(dat)[9] <- "estimate"

f_heatmap2 <- heatmap_function2(dat, "")
ggsave(f_heatmap2, filename=paste0("figures/main_figures/(fig2c)_heatmap_19edc-cov.pdf"), width = 15, height = 5, limitsize = FALSE)
#### EDC-outcome 相关性热图 (Wilcoxon) 填补后 ----
