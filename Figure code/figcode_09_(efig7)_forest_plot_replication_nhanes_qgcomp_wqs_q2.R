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


############################# QGCOMP #############################
### QGCOMP: nhanes_dat1 (13EDC) (未纳入NHANES权重) ###
{
  forest_plot_weight1 <- list()
  for (i in c("nhanes_dat1_1114")) {
    sample_name1 <- i # (13 EDC)
    #### 数据处理 (qgcomp results for circos heatmap & forest plot) ----
    ### 设置纳入图片的暴露和结局
    ## 暴露：EDC
    exposure <- c("EDC_13","EDC_pos","EDC_neg","PFAS","PAE_6","TC_1","BP_1")
    # 标准化名称
    exposure_label <- c("EDCs (13)","EDCs (positive weight)","EDCs (negative weight)","PFAS (5)","PAEs (6)","Antimicrobials (1)","Bisphenols (1)")
    
    ## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
    outcome <- c("dm","ckd","cvd")
    outcome_final <- outcome
    # 标准化名称
    outcome_label <- c("Diabetes","CKD","CVD")
    
    # 读取数据 (总人群)
    results_qg_all1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name1,")_NOsamplingweights.xlsx"))
    # 读取数据 (总人群)
    results_qg_edc1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(edc+-)_(",sample_name1,")_NOsamplingweights.xlsx"))
    
    results_qg <- rbind(results_qg_all1, results_qg_edc1) %>%
      filter(exp %in% exposure & out %in% outcome_final)
    colnames(results_qg)[9:21] <- gsub("_log10","",colnames(results_qg)[9:21])
    
    results_qg$z <- results_qg$estimate/results_qg$se
    
    # FDR 校正
    # 以每个exposure表型为单位，校正outcome
    dat <- results_qg %>%
      group_by(exp) %>%  # 按exposure分组，校正outcome
      mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
      ungroup()
    
    # EDC 权重数据整理
    dat$numb_pos_weight <- rowSums(dat[, 9:21] > 0, na.rm = TRUE)# 计算每一行第9列到第22列中正数的数量
    dat$numb_neg_weight <- rowSums(dat[, 9:21] < 0, na.rm = TRUE)# 计算每一行第9列到第22列中负数的数量
    
    dat$weight_threshold_pos <- 1 / dat$numb_pos_weight
    dat$weight_threshold_neg <- 1 / dat$numb_neg_weight
    
    dat_edc <- dat[dat$exp == "EDC_13",]
    dat_edc_19_long_qgcomp <- tidyr::gather(dat_edc, edc, weight, 9:21, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
    
    dat_qgcomp <- dat
    dat_qgcomp$exp <- factor(dat_qgcomp$exp, levels = exposure, labels = exposure_label)
    dat_qgcomp$out <- factor(dat_qgcomp$out, levels = outcome_final, labels = outcome_label)
    dat_qgcomp <- dat_qgcomp %>% arrange(out,exp)
    #### 数据处理 (qgcomp results for circos heatmap & forest plot) ####
    ############################################# 森林图 (qgcomp) #############################################
    #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
    dat_forest <- dat_qgcomp[dat_qgcomp$out %in% c("Diabetes","CKD","CVD"),]
    dat_forest$hr <- exp(dat_forest$estimate)
    dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
    dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
    dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
    dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
    dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
    
    # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
    dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
    dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
    dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
    
    ### 森林图
    # CVD
    {
      f_forest2 <- ggplot(data=dat_forest2, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("CVD"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
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
    # CKD
    {
      f_forest3 <- ggplot(data=dat_forest3, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("CKD"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
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
    # DM
    {
      f_forest4 <- ggplot(data=dat_forest4, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("Diabetes"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
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
    
    
    dat_weight <- dat_forest[dat_forest$exp == "EDCs (13)",]
    dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 9:21, na.rm=TRUE, factor_key=TRUE)
    dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
    
    dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
    
    dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
    
    dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
    
    ### 权重图
    # CVD
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    # CKD
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    # DM
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    
    
    f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
    f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
    f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1))
    # ggsave(f_fw, filename=paste0("figures/supplementary_figures/forest_qgcomp_(",sample_name1,")_20260214.pdf"), width = 10, height = 10, limitsize = FALSE)
    
    forest_plot_weight1[[sample_name1]] <- f_fw
    #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
  }
  forest_plot_weight1_all <- cowplot::plot_grid(forest_plot_weight1[[1]])
  ggsave(forest_plot_weight1_all, filename=paste0("figures/supplementary_figures/(efig7a)_forest_qgcomp_nhanes_13edc_NOsamplingweights.pdf"), width = 9, height = 12, limitsize = FALSE)
}

### QGCOMP: nhanes_dat5 (PAE_6) (未纳入NHANES权重) ###
{
  forest_plot_weight2 <- list()
  for (i in c("nhanes_dat5_0318")) {
    sample_name1 <- i # (PAE_6)
    #### 数据处理 (qgcomp results for circos heatmap & forest plot) ----
    ### 设置纳入图片的暴露和结局
    ## 暴露：EDC
    exposure <- c("PAE_6","EDC_pos","EDC_neg")
    # 标准化名称
    exposure_label <- c("PAEs (6)","EDCs (positive weight)","EDCs (negative weight)")
    
    ## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
    outcome <- c("dm","ckd","cvd")
    outcome_final <- outcome
    # 标准化名称
    outcome_label <- c("Diabetes","CKD","CVD")
    
    # 读取数据 (总人群)
    results_qg_all1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(",sample_name1,")_NOsamplingweights.xlsx"))
    # 读取数据 (总人群)
    results_qg_edc1 <- readxl::read_xlsx(paste0("results/correlations/qgcomp/nhanes/qgcomp_results_(q2)_(edc+-)_(",sample_name1,")_NOsamplingweights.xlsx"))
    
    results_qg <- rbind(results_qg_all1, results_qg_edc1) %>%
      filter(exp %in% exposure & out %in% outcome_final)
    colnames(results_qg)[9:14] <- gsub("_log10","",colnames(results_qg)[9:14])
    
    results_qg$z <- results_qg$estimate/results_qg$se
    
    # FDR 校正
    # 以每个exposure表型为单位，校正outcome
    dat <- results_qg %>%
      group_by(exp) %>%  # 按exposure分组，校正outcome
      mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
      ungroup()
    
    # EDC 权重数据整理
    dat$numb_pos_weight <- rowSums(dat[, 9:14] > 0, na.rm = TRUE)# 计算每一行第9列到第22列中正数的数量
    dat$numb_neg_weight <- rowSums(dat[, 9:14] < 0, na.rm = TRUE)# 计算每一行第9列到第22列中负数的数量
    
    dat$weight_threshold_pos <- 1 / dat$numb_pos_weight
    dat$weight_threshold_neg <- 1 / dat$numb_neg_weight
    
    dat_edc <- dat[dat$exp == "EDC_13",]
    dat_edc_19_long_qgcomp <- tidyr::gather(dat_edc, edc, weight, 9:14, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
    
    dat_qgcomp <- dat
    dat_qgcomp$exp <- factor(dat_qgcomp$exp, levels = exposure, labels = exposure_label)
    dat_qgcomp$out <- factor(dat_qgcomp$out, levels = outcome_final, labels = outcome_label)
    dat_qgcomp <- dat_qgcomp %>% arrange(out,exp)
    #### 数据处理 (qgcomp results for circos heatmap & forest plot) ####
    ############################################# 森林图 (qgcomp) #############################################
    #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
    dat_forest <- dat_qgcomp[dat_qgcomp$out %in% c("Diabetes","CKD","CVD"),]
    dat_forest$hr <- exp(dat_forest$estimate)
    dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
    dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
    dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
    dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
    dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
    
    # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
    dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
    dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
    dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
    
    ### 森林图
    # CVD
    {
      f_forest2 <- ggplot(data=dat_forest2, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("CVD"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
        # 将颜色、图例标题、逆序整合到 scale_colour_manual
        scale_colour_manual(
          name = "EDC groups",  # 修改图例标题
          values = c("#5758BB", "#EE5A24", "#3498db"),
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
    # CKD
    {
      f_forest3 <- ggplot(data=dat_forest3, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("CKD"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
        # 将颜色、图例标题、逆序整合到 scale_colour_manual
        scale_colour_manual(
          name = "EDC groups",  # 修改图例标题
          values = c("#5758BB", "#EE5A24", "#3498db"),
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
    # DM
    {
      f_forest4 <- ggplot(data=dat_forest4, aes(x=estimate, y=exp)) +
        geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, col=exp), width=0, cex=1, position = position_dodge(0.7)) +
        geom_point(aes(col=exp), cex = 3, position = position_dodge(0.7)) +
        geom_text(aes(label = text),
                  hjust = -0.5,
                  size = 5,
                  position = position_dodge(width = 0.7)) +
        
        scale_x_continuous(
          labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
        ) +
        
        scale_y_discrete(expand = c(0.1, 0.1)) +
        
        geom_vline(xintercept = 0,
                   linetype = "dashed",
                   linewidth = 0.5) +
        
        labs(
          title = paste0("Diabetes"," (",sample_name1,")"),
          x = "OR (95% CI)") +
        
        # 将颜色、图例标题、逆序整合到 scale_colour_manual
        scale_colour_manual(
          name = "EDC groups",  # 修改图例标题
          values = c("#5758BB", "#EE5A24", "#3498db"),
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
    
    
    dat_weight <- dat_forest[dat_forest$exp == "PAEs (6)",]
    dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 9:14, na.rm=TRUE, factor_key=TRUE)
    dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
    
    dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
    
    dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
    
    dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
      arrange(rev(out),-weight_abs)
    dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
    
    ### 权重图
    # CVD
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    # CKD
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    # DM
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
        
        # 可视化美化
        labs(
          title = " ",
          x = "Weight", y = "Outcomes") +
        
        scale_fill_gradientn(
          # colours = rev(brbg_colors[c(4:8)]),
          colours = rev(weight_colors),
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
    
    
    f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
    f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
    f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1))
    # ggsave(f_fw, filename=paste0("figures/supplementary_figures/forest_qgcomp_(",sample_name1,")_20260212.pdf"), width = 10, height = 7, limitsize = FALSE)
    
    forest_plot_weight2[[sample_name1]] <- f_fw
    #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
  }
  forest_plot_weight2_all <- cowplot::plot_grid(forest_plot_weight2[[1]])
  ggsave(forest_plot_weight2_all, filename=paste0("figures/supplementary_figures/(efig7c)_forest_qgcomp_nhanes_6pae_NOsamplingweights.pdf"), width = 9, height = 7, limitsize = FALSE)
}


############################# WQS #############################
### WQS: nhanes_dat1 (13EDC) (未纳入NHANES权重) ###
for (i in c("nhanes_dat1_1114")) {
  # i <- "nhanes_dat1_1114"
  # 设置人群名称
  sample_name1 <- i
  # 读取数据
  results_wqs_bin <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_bin_(",sample_name1,")_NOsamplingweights.xlsx"))
  
  #### 数据处理 (WQS results for circos heatmap & forest plot) ----
  ### 设置纳入图片的暴露和结局
  ## 暴露：EDC
  exposure <- c("EDC_13","PFAS","PAE_6","TC_1","BP_1")
  # 标准化名称
  exposure_label <- c("EDCs (13)","PFAS (5)","PAEs (6)","Antimicrobials (1)","Bisphenols (1)")
  
  ## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
  outcome <- c("dm","ckd","cvd")
  outcome_final <- outcome
  out_dat <- as.data.frame(outcome_final)
  # 标准化名称
  outcome_label <- c("Diabetes","CKD","CVD")
  
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
  
  
  results_wqs <- results_wqs_bin %>%
    filter(exp %in% exposure & out %in% outcome_final & cov == "wqs" & type == "Q2")
  
  results_wqs_pos <- results_wqs[results_wqs$direction == "pos",] # 正向权重wqs分析结果
  unique(results_wqs_pos$exp) # 缺TC_1 BP_1
  unique(results_wqs_pos$out)
  results_wqs_pos <- left_join(dat_exp_out_all,results_wqs_pos,by=c("exp", "out"))
  results_wqs_pos$p <- ifelse(is.na(results_wqs_pos$p), 1, results_wqs_pos$p)
  
  results_wqs_neg <- results_wqs[results_wqs$direction == "neg",] # 负向权重wqs分析结果
  unique(results_wqs_neg$exp) # 缺TC_1 BP_1
  unique(results_wqs_neg$out)
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
  dat_wqs_pos$numb_pos_weight <- rowSums(dat_wqs_pos[, 13:25] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
  dat_wqs_pos$numb_neg_weight <- rowSums(dat_wqs_pos[, 13:25] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量
  
  dat_wqs_pos$weight_threshold_pos <- 1 / dat_wqs_pos$numb_pos_weight
  dat_wqs_pos$weight_threshold_neg <- 1 / dat_wqs_pos$numb_neg_weight
  
  dat_wqs_pos_edc13 <- dat_wqs_pos[dat_wqs_pos$exp == "EDC_13",]
  dat_wqs_pos_edc13_long <- tidyr::gather(dat_wqs_pos_edc13, edc, weight, 13:25, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
  dat_wqs_pos_edc13_long$edc <- factor(dat_wqs_pos_edc13_long$edc, levels = )
  
  # EDC 权重数据整理(负)
  dat_wqs_neg$numb_pos_weight <- rowSums(dat_wqs_neg[, 13:25] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
  dat_wqs_neg$numb_neg_weight <- rowSums(dat_wqs_neg[, 13:25] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量
  
  dat_wqs_neg$weight_threshold_pos <- 1 / dat_wqs_neg$numb_pos_weight
  dat_wqs_neg$weight_threshold_neg <- 1 / dat_wqs_neg$numb_neg_weight
  
  dat_wqs_neg_edc13 <- dat_wqs_neg[dat_wqs_neg$exp == "EDC_13",]
  dat_wqs_neg_edc13_long <- tidyr::gather(dat_wqs_neg_edc13, edc, weight, 13:25, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
  
  
  # Prepare data for plot
  dat_wqs_pos_for_plot <- dat_wqs_pos
  dat_wqs_pos_for_plot$exp <- factor(dat_wqs_pos_for_plot$exp, levels = exposure, labels = exposure_label)
  dat_wqs_pos_for_plot$out <- factor(dat_wqs_pos_for_plot$out, levels = outcome_final, labels = outcome_label)
  dat_wqs_pos_for_plot <- dat_wqs_pos_for_plot %>% arrange(out,exp)
  
  dat_wqs_neg_for_plot <- dat_wqs_neg
  dat_wqs_neg_for_plot$exp <- factor(dat_wqs_neg_for_plot$exp, levels = exposure, labels = exposure_label)
  dat_wqs_neg_for_plot$out <- factor(dat_wqs_neg_for_plot$out, levels = outcome_final, labels = outcome_label)
  dat_wqs_neg_for_plot <- dat_wqs_neg_for_plot %>% arrange(out,exp)
  #### 数据处理 (WQS results for circos heatmap & forest plot) ####
  
  
  for (j in c("positive","negative")) {
    if(j == "positive"){
      # 正向模型 #
      dat_forest <- dat_wqs_pos_for_plot
      dat_forest <- dat_forest[dat_forest$exp %in% c("EDCs (13)","PFAS (5)","PAEs (6)"),]
      
      forest_color <- c("#3498db", "#f55d78", "#476066")
      color_wqs_weight <- "#82589F"
      export <- paste0("figures/supplementary_figures/(efig7b)_forest_wqs_pos_weight_3_nhanes_(",sample_name1,")_NOsamplingweights.pdf")
      
      ############################################# 森林图 (WQS) #############################################
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
      dat_forest$hr <- exp(dat_forest$estimate)
      dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
      dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
      dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
      dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
      dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
      
      # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
      dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
      dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
      dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
      
      ### 森林图
      # CVD
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
            title = "CVD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # CKD
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
            title = "CKD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # DM
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
            title = "Diabetes",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      
      
      dat_weight <- dat_forest[dat_forest$exp == "EDCs (13)",]
      dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 13:25, na.rm=TRUE, factor_key=TRUE)
      dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
      
      dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
      
      dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
      
      dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
      
      ### 权重图
      # CVD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # CKD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # DM
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      
      f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
      f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
      f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1.2))
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
      ggsave(f_fw, filename=export, width = 6, height = 12, limitsize = FALSE)
      
    }else if(j == "negative"){
      # 负向模型 #
      dat_forest <- dat_wqs_neg_for_plot
      dat_forest <- dat_forest[dat_forest$exp %in% c("EDCs (13)","PFAS (5)","PAEs (6)"),]
      
      forest_color <- c("#3498db", "#f55d78", "#476066")
      color_wqs_weight <- "#1289A7"
      export <- paste0("figures/supplementary_figures/(efig7b)_forest_wqs_neg_weight_3_nhanes_(",sample_name1,")_NOsamplingweights.pdf")
      
      ############################################# 森林图 (WQS) #############################################
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
      dat_forest$hr <- exp(dat_forest$estimate)
      dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
      dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
      dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
      dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
      dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
      
      # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
      dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
      dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
      dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
      
      ### 森林图
      # CVD
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
            title = "CVD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # CKD
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
            title = "CKD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # DM
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
            title = "Diabetes",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      
      
      dat_weight <- dat_forest[dat_forest$exp == "EDCs (13)",]
      dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 13:25, na.rm=TRUE, factor_key=TRUE)
      dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
      
      dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
      
      dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
      
      dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
      
      ### 权重图
      # CVD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # CKD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # DM
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      
      f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
      f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
      f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1.2))
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
      ggsave(f_fw, filename=export, width = 6, height = 12, limitsize = FALSE)
      
    }
  }
}

### WQS: nhanes_dat5 (PAE_6) (未纳入NHANES权重) ###
for (i in c("nhanes_dat5_0318")) {
  # 设置人群名称
  # sample_name1 <- "nhanes_dat5_0318"
  sample_name1 <- i
  # 读取数据
  results_wqs_bin <- readxl::read_xlsx(paste0("results/correlations/wqs/wqs_results_bin_(",sample_name1,")_NOsamplingweights.xlsx"))
  
  #### 数据处理 (WQS results for circos heatmap & forest plot) ----
  ### 设置纳入图片的暴露和结局
  ## 暴露：EDC
  exposure <- c("PAE_6")
  # 标准化名称
  exposure_label <- c("PAEs (6)")
  
  ## 结局：12类biomarker（用总人群（用10年指标，没有的用14年的补充））
  outcome <- c("dm","ckd","cvd")
  outcome_final <- outcome
  out_dat <- as.data.frame(outcome_final)
  # 标准化名称
  outcome_label <- c("Diabetes","CKD","CVD")
  
  
  results_wqs <- results_wqs_bin %>%
    filter(exp %in% exposure & out %in% outcome_final & cov == "wqs" & type == "Q2")
  
  results_wqs_pos <- results_wqs[results_wqs$direction == "pos",] # 正向权重wqs分析结果
  unique(results_wqs_pos$exp) # 缺TC_1 BP_1
  unique(results_wqs_pos$out)
  results_wqs_pos$p <- ifelse(is.na(results_wqs_pos$p), 1, results_wqs_pos$p)
  
  results_wqs_neg <- results_wqs[results_wqs$direction == "neg",] # 负向权重wqs分析结果
  unique(results_wqs_neg$exp) # 缺TC_1 BP_1
  unique(results_wqs_neg$out)
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
  dat_wqs_pos$numb_pos_weight <- rowSums(dat_wqs_pos[, 13:18] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
  dat_wqs_pos$numb_neg_weight <- rowSums(dat_wqs_pos[, 13:8] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量
  
  dat_wqs_pos$weight_threshold_pos <- 1 / dat_wqs_pos$numb_pos_weight
  dat_wqs_pos$weight_threshold_neg <- 1 / dat_wqs_pos$numb_neg_weight
  
  dat_wqs_pos_edc13 <- dat_wqs_pos[dat_wqs_pos$exp == "PAE_6",]
  dat_wqs_pos_edc13_long <- tidyr::gather(dat_wqs_pos_edc13, edc, weight, 13:18, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
  
  
  # EDC 权重数据整理(负)
  dat_wqs_neg$numb_pos_weight <- rowSums(dat_wqs_neg[, 13:18] > 0, na.rm = TRUE)# 计算每一行第13列到第26列中正数的数量
  dat_wqs_neg$numb_neg_weight <- rowSums(dat_wqs_neg[, 13:18] < 0, na.rm = TRUE)# 计算每一行第13列到第26列中负数的数量
  
  dat_wqs_neg$weight_threshold_pos <- 1 / dat_wqs_neg$numb_pos_weight
  dat_wqs_neg$weight_threshold_neg <- 1 / dat_wqs_neg$numb_neg_weight
  
  dat_wqs_neg_edc13 <- dat_wqs_neg[dat_wqs_neg$exp == "PAE_6",]
  dat_wqs_neg_edc13_long <- tidyr::gather(dat_wqs_neg_edc13, edc, weight, 13:18, na.rm=TRUE, factor_key=TRUE)  #注意根据数据格式进行修改，容易忘记！！！！！
  
  
  # Prepare data for plot
  dat_wqs_pos_for_plot <- dat_wqs_pos
  dat_wqs_pos_for_plot$exp <- factor(dat_wqs_pos_for_plot$exp, levels = exposure, labels = exposure_label)
  dat_wqs_pos_for_plot$out <- factor(dat_wqs_pos_for_plot$out, levels = outcome_final, labels = outcome_label)
  dat_wqs_pos_for_plot <- dat_wqs_pos_for_plot %>% arrange(out,exp)
  
  dat_wqs_neg_for_plot <- dat_wqs_neg
  dat_wqs_neg_for_plot$exp <- factor(dat_wqs_neg_for_plot$exp, levels = exposure, labels = exposure_label)
  dat_wqs_neg_for_plot$out <- factor(dat_wqs_neg_for_plot$out, levels = outcome_final, labels = outcome_label)
  dat_wqs_neg_for_plot <- dat_wqs_neg_for_plot %>% arrange(out,exp)
  #### 数据处理 (WQS results for circos heatmap & forest plot) ####
  
  for (j in c("positive","negative")) {
    if(j == "positive"){
      # 正向模型 #
      dat_forest <- dat_wqs_pos_for_plot
      dat_forest <- dat_forest[dat_forest$exp %in% c("PAEs (6)"),]
      
      forest_color <- c("#3498db")
      color_wqs_weight <- "#82589F"
      export <- paste0("figures/supplementary_figures/(efig7d)_forest_wqs_pos_weight_3_nhanes_(",sample_name1,")_NOsamplingweights.pdf")
      
      ############################################# 森林图 (WQS) #############################################
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
      dat_forest$hr <- exp(dat_forest$estimate)
      dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
      dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
      dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
      dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
      dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
      
      # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
      dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
      dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
      dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
      
      ### 森林图
      # CVD
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
            title = "CVD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # CKD
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
            title = "CKD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # DM
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
            title = "Diabetes",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      
      
      dat_weight <- dat_forest[dat_forest$exp == "PAEs (6)",]
      dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 13:18, na.rm=TRUE, factor_key=TRUE)
      dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
      
      dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
      
      dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
      
      dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
      
      ### 权重图
      # CVD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # CKD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # DM
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      
      f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
      f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
      f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1.2))
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
      ggsave(f_fw, filename=export, width = 6, height = 7, limitsize = FALSE)
      
    }else if(j == "negative"){
      # 负向模型 #
      dat_forest <- dat_wqs_neg_for_plot
      dat_forest <- dat_forest[dat_forest$exp %in% c("PAEs (6)"),]
      
      forest_color <- c("#3498db")
      color_wqs_weight <- "#1289A7"
      export <- paste0("figures/supplementary_figures/(efig7d)_forest_wqs_neg_weight_3_nhanes_(",sample_name1,")_NOsamplingweights.pdf")
      
      ############################################# 森林图 (WQS) #############################################
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ----
      dat_forest$hr <- exp(dat_forest$estimate)
      dat_forest$lci <- exp(dat_forest$estimate - 1.96*dat_forest$se)
      dat_forest$uci <- exp(dat_forest$estimate + 1.96*dat_forest$se)
      dat_forest$out <- factor(dat_forest$out, levels = c("Diabetes","CKD","CVD"))
      dat_forest$exp <- factor(dat_forest$exp, levels = rev(exposure_label))
      dat_forest$text <- paste0(sprintf("%.3f", dat_forest$hr)," (",sprintf("%.3f", dat_forest$lci),", ",sprintf("%.3f", dat_forest$uci),")") 
      
      # dat_forest1 <- dat_forest[dat_forest$out == "Incident CVD (2010-2014)",]
      dat_forest2 <- dat_forest[dat_forest$out == "CVD",]
      dat_forest3 <- dat_forest[dat_forest$out == "CKD",]
      dat_forest4 <- dat_forest[dat_forest$out == "Diabetes",]
      
      ### 森林图
      # CVD
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
            title = "CVD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # CKD
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
            title = "CKD",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      # DM
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
            title = "Diabetes",
            subtitle = sample_name1,
            x = "OR (95% CI)") +
          
          # 将颜色、图例标题、逆序整合到 scale_colour_manual
          scale_colour_manual(
            name = "EDC groups",  # 修改图例标题
            values = forest_color,
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
      
      
      dat_weight <- dat_forest[dat_forest$exp == "PAEs (6)",]
      dat_weight_long <- tidyr::gather(dat_weight, edc, weight, 13:18, na.rm=TRUE, factor_key=TRUE)
      dat_weight_long$weight_abs <- abs(dat_weight_long$weight)
      
      dat_weight_long2 <- dat_weight_long[dat_weight_long$out == "CVD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long2$edc <- factor(dat_weight_long2$edc, levels = rev(dat_weight_long2$edc))
      
      dat_weight_long3 <- dat_weight_long[dat_weight_long$out == "CKD",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long3$edc <- factor(dat_weight_long3$edc, levels = rev(dat_weight_long3$edc))
      
      dat_weight_long4 <- dat_weight_long[dat_weight_long$out == "Diabetes",]%>%
        arrange(rev(out),-weight_abs)
      dat_weight_long4$edc <- factor(dat_weight_long4$edc, levels = rev(dat_weight_long4$edc))
      
      ### 权重图
      # CVD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # CKD
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      # DM
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
          
          # 可视化美化
          labs(
            title = " ",
            x = "Weight", y = "Outcomes") +
          
          # 设置连续颜色渐变
          scale_fill_gradient2(
            low = "white",            # 低点 对应白色
            high = color_wqs_weight,  # 高点 对应X色
            midpoint = 0,      # 中间点（默认是 0，可不写）
            limits = c(0, max(dat_weight_long$weight)), # 确保颜色范围覆盖数据
            name = "Weight"    # 图例标题
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
      
      f_forest <- cowplot::plot_grid(f_forest4, f_forest3, f_forest2, nrow = 3, rel_heights = c(1,1,1))
      f_weight <- cowplot::plot_grid(f_weight4, f_weight3, f_weight2, nrow = 3, rel_heights = c(1,1,1))
      f_fw <- cowplot::plot_grid(f_forest, NULL, f_weight, ncol = 3, rel_widths = c(1.1,0.5,1.2))
      #### 森林图+权重图 (3 outcome ["dm","ckd","cvd"]) ####
      ggsave(f_fw, filename=export, width = 6, height = 7, limitsize = FALSE)
      
    }
  }
}
