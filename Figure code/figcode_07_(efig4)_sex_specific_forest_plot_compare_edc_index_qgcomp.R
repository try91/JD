library(data.table)
library(dplyr)
library(ggplot2)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

#### 配色 ----
PFAS_colors <- colorRampPalette(c("#f55d78","#FFFFFF"))(50)[c(1,8,15,22,29)]
PAE_colors <- colorRampPalette(c("#3498db","#FFFFFF"))(50)[c(1,6,11,16,21,26,31,36,41)]
BP_colors <- colorRampPalette(c("#f1c40f","#FFFFFF"))(50)[c(1,18,35)]
TC_colors <- colorRampPalette(c("#00b894","#FFFFFF"))(50)[c(1,21)]
my_colors <- c(PFAS_colors,PAE_colors,BP_colors,TC_colors)

# 创建25色渐变色标尺
my_palette <- colorRampPalette(colors = c("#cf6a87", "#f19066", "#f5cd79", "#33d9b2", "#63cdda", "#34ace0", "#786fa6"))(25)
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

#### 结局变量名汇总 ----
## 结局 横断面研究结果
select_out_cat <- c("dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b")
# 标准化名称
select_out_cat_labels <- c("Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT")
#### 结局变量名汇总 ####


# # 男性 (总人群, n=3930)
# sample_name <- "phy_edc_temp0_1"  # 男性 (总人群, n=3930)
# # 女性 (总人群, n=6394)
# sample_name <- "phy_edc_temp0_2"  # 女性 (总人群, n=6394)
for(sample_name in c("phy_edc_temp0_1","phy_edc_temp0_2")){
  #### 数据处理 ----
  ## 读取COX & Logistic分析结果
  cox_results_edc_incident_all <- readxl::read_xlsx("results/cox/sex_specific_cox_results_edc_incident.xlsx")
  cox_results_edc_index_incident_all <- readxl::read_xlsx("results/cox/sex_specific_cox_results_edc_index_incident.xlsx")
  
  logistic_results_edc_incident_all <- readxl::read_xlsx("results/glm/sex_specific_logistic_results_edc_incident.xlsx")
  logistic_results_edc_index_incident_all <- readxl::read_xlsx("results/glm/sex_specific_logistic_results_edc_index_incident.xlsx")
  
  ## 读取Qgcomp分析结果
  qg_results_q2 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(",sample_name,").xlsx"))
  qg_results_q2_pn <- readxl::read_xlsx(paste0("results/correlations/qgcomp/qgcomp_results_(q2)_(edc+-)_(",sample_name,").xlsx"))
  
  
  ## 整合分析结果
  cox_results <- rbind(cox_results_edc_incident_all,cox_results_edc_index_incident_all)
  cox_results <- cox_results[cox_results$sample == sample_name,]
  cox_results <- cox_results[cox_results$adjust == "adj",]
  cox_results <- cox_results[cox_results$outcome %in% phy_incident_cat,]
  cox_results <- cox_results[cox_results$rowname %in% edc_index_b_keep,]
  cox_results <- cox_results[,c(1,3:13)]
  colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
  out_cox <- unique(cox_results$outcome) # 提取结局变量
  
  logistic_results <- rbind(logistic_results_edc_incident_all,logistic_results_edc_index_incident_all)
  logistic_results <- logistic_results[logistic_results$sample == sample_name,]
  logistic_results <- logistic_results[logistic_results$adjust == "adj",]
  logistic_results <- logistic_results[logistic_results$outcome %in% select_out_cat,]
  logistic_results <- logistic_results[logistic_results$rowname %in% edc_index_b_keep,]
  colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
  out_logistic <- unique(logistic_results$outcome) # 提取结局变量
  
  qg_results_q2 <- rbind(qg_results_q2,qg_results_q2_pn)
  qg_results_q2 <- qg_results_q2[qg_results_q2$sample == sample_name,]
  qg_results_q2 <- qg_results_q2[,c(1:8)]
  qg_results_q2$z <- qg_results_q2$estimate / qg_results_q2$se 
  colnames(qg_results_q2)[1:2] <- c("rowname","outcome")
  qg_results_q2$outcome <- ifelse((qg_results_q2$outcome %in% c(phy_incident_cat,phy_censor_cat)) & (qg_results_q2$method == "qgcomp bin"), paste0(qg_results_q2$outcome,"_logistic"), qg_results_q2$outcome)
  qg_results_q2 <- qg_results_q2[qg_results_q2$outcome %in% c(phy_incident_cat,select_out_cat),]
  qg_results_q2 <- qg_results_q2[qg_results_q2$rowname %in% c("EDC_14", "PFAS", "PAE_6", "BP_1", "TC"),]
  
  
  dat_result1 <- rbind(cox_results,logistic_results)
  dat_result1_1 <- dat_result1[dat_result1$rowname %in% c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_1 <- dat_result1_1 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_2 <- dat_result1[dat_result1$rowname %in% c("edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),]
  # 以每个OUT表型为单位进行校正
  dat_result1_2 <- dat_result1_2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  # 以每个OUT表型为单位进行校正
  qg_results_q2 <- qg_results_q2 %>%
    group_by(outcome) %>%  # 按outcome分组
    mutate(p_adj_bh = p.adjust(p, method = "BH"))  # 对每个分组的P值进行FDR校正
  
  dat_result1_1$method_type <- "edc_count_index"
  dat_result1_2$method_type <- "edc_score"
  qg_results_q2$method_type <- "qgcomp"
  
  
  dat_result_all <- bind_rows(dat_result1_1,dat_result1_2,qg_results_q2)
  
  dat_result_all$exposure_type <- ifelse(dat_result_all$rowname %in% c("EDC_14","edc_count2_edc14_b","edc_score_edc14_b"), "EDCs (14)",
                                         ifelse(dat_result_all$rowname %in% c("PFAS","edc_count2_pfas_b","edc_score_pfas_b"), "PFAS (5)",
                                                ifelse(dat_result_all$rowname %in% c("PAE_6","edc_count2_pae6_b","edc_score_pae6_b"), "PAEs (6)",
                                                       ifelse(dat_result_all$rowname %in% c("TC","edc_count2_tc_b","edc_score_tc_b"), "Antimicrobials (2)", "Bisphenols (1)"))))
  
  
  table(dat_result1_1$rowname,dat_result1_1$outcome)
  table(dat_result1_2$rowname,dat_result1_2$outcome)
  table(qg_results_q2$rowname,qg_results_q2$outcome)
  
  dat_result_all$rowname <- factor(dat_result_all$rowname, levels = c("EDC_14", "PFAS", "PAE_6", "TC", "BP_1",
                                                                      "edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_tc_b","edc_count2_bp1_b",
                                                                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_tc_b","edc_score_bp1_b"),
                                   labels = c("EDCs (14)", "PFAS (5)", "PAEs (6)", "Antimicrobials (2)", "Bisphenols (1)",
                                              "EDC Scoremedian (14 EDCs)", "EDC Scoremedian (PFAS)", "EDC Scoremedian (PAEs)", "EDC Scoremedian (antimicrobials)", "EDC Scoremedian (bisphenols)",    # EDC count index
                                              "EDC Scorequartile (14 EDCs)", "EDC Scorequartile (PFAS)", "EDC Scorequartile (PAEs)", "EDC Scorequartile (antimicrobials)", "EDC Scorequartile (bisphenols)"))
  
  dat_result_all$outcome <- factor(dat_result_all$outcome,
                                   levels = rev(c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021",select_out_cat)),
                                   labels = rev(c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)",select_out_cat_labels)))
  
  dat_result_all$method_type <- factor(dat_result_all$method_type, levels = rev(c("qgcomp", "edc_count_index", "edc_score")))
  
  
  # 计算95%CI
  dat_result_all$z.lci <- dat_result_all$z - 1.96
  dat_result_all$z.uci <- dat_result_all$z + 1.96
  
  dat_result_all$lci <- dat_result_all$estimate - 1.96*dat_result_all$se
  dat_result_all$uci <- dat_result_all$estimate + 1.96*dat_result_all$se
  #### 数据处理 ####
  
  
  ############################################# 森林图 #############################################
  #### forest plot (对比edc index和qgcomp发现对4个主要outcome和12个次要outcome的结果) ----
  # 构建绘制森林图的函数 (根据每一个污染物index分开，5列for附图)
  forest_function1 <- function(DAT,X_TITLE){
    
    f <- ggplot(data=DAT, aes(x=z, y=outcome, col=method_type)) +
      geom_errorbar(aes(xmin=z.lci, xmax=z.uci, col=method_type), width=0, cex=0.8, position = position_dodge(0.55)) +
      geom_point(aes(col=method_type), cex = 2.4, position = position_dodge(0.55)) +
      geom_text(aes(label = text, group = method_type),
                hjust = -0.5,
                size = 3.6,
                position = position_dodge(width = 0.65)) +
      
      scale_x_continuous(
        labels = function(x) round(x, 1)  # 显示为x，保留1位小数
      ) +
      
      geom_vline(xintercept = 0,
                 linetype = "dashed",
                 linewidth = 0.5) +
      
      # 将颜色、图例标题、整合到 scale_colour_manual
      scale_colour_manual(
        name = "",  # legend title
        
        # 设置颜色
        values = c("qgcomp" = "#476066",
                   "edc_count_index" = "#d62c2c",
                   "edc_score" = "#591129"),
        # 直接在scale_color_manual中设置标签
        labels = c("qgcomp" = "Mixture effectqg-comp",  
                   "edc_count_index" = "EDC Scoremedian",
                   "edc_score" = "EDC Scorequartile"),
        # 通过breaks参数明确指定图例顺序
        breaks = c("qgcomp",
                   "edc_count_index",
                   "edc_score")
      ) +
      
      labs(x="Z (95% CI)") +
      
      # 纵向排列面板
      facet_wrap(~exposure_type, ncol = 5, scales = "free_x") +   # x轴自由缩放
      
      theme_classic() +
      
      theme(
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        
        strip.text = element_text(size = 15, colour = "black"), # facet_wrap面板字体大小
        
        axis.line = element_line(size = 0.5, colour = "black"),
        strip.background = element_rect(size = 0.5),
        
        panel.spacing.x = unit(8, "mm"),  # 横向间距
        panel.spacing.y = unit(8, "mm"), # 纵向间距
        
        axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
        axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
        
        axis.text.x = element_text(size = 15, colour = "black"), # 调整x轴文字
        axis.text.y = element_text(size = 15, colour = "black"),
        
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 15, colour = "black")   # 调整legend文本大小
      )
    
    return(f)
  }
  
  # 绘制森林图比较构建的EDC index和qgcomp结果
  dat_result_all$exposure_type <- factor(dat_result_all$exposure_type, levels = c("EDCs (14)", "PFAS (5)", "PAEs (6)", "Antimicrobials (2)", "Bisphenols (1)"))
  dat_result_all$text <- paste0(sprintf("%.3f", dat_result_all$z)," (",sprintf("%.3f", dat_result_all$z.lci),", ",sprintf("%.3f", dat_result_all$z.uci),")") 
  
  f_forest1 <- forest_function1(dat_result_all)
  
  if(sample_name == "phy_edc_temp0_1"){
    ggsave(f_forest1, filename=paste0("figures/supplementary_figures/(efig4a)_forest_compare_edc_qgcomp_outcome_sex_specific_(",sample_name,").pdf"), width = 13, height = 12, limitsize = FALSE)
  }
  
  if(sample_name == "phy_edc_temp0_2"){
    ggsave(f_forest1, filename=paste0("figures/supplementary_figures/(efig4b)_forest_compare_edc_qgcomp_outcome_sex_specific_(",sample_name,").pdf"), width = 13, height = 12, limitsize = FALSE)
  }
  #### forest plot (对比edc index和qgcomp发现对4个主要outcome和12个次要outcome的结果) ####
}

