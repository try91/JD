library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ppcor)
library(survival)
library(rms)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/analyte_measurements_dat_20261006.txt", header = TRUE)
micro_dat <- read.table("raw_data/gut_microbial_composition_function_pathway_profiles_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp"  # 提取subgroup的名称

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

#### 配色 ----
my_colors <- c(rep("#f55d78",5),rep("#3498db",9),rep("#00b894",2),rep("#f1c40f",3))

my_colors_list <- list()
for (i in 1:19) {
  my_colors_list[[edc_traits_log10[i]]] <- my_colors[i]
}
#### 配色 ####


############################################# forest图 #############################################
#### 读取COX Logistic分析结果 ----
### 读取汇总COX结果数据 ###
cox_results_all <- readxl::read_xlsx("results/cox/cox_results_edc_incident.xlsx")
cox_results_all_temp <- cox_results_all[cox_results_all$rowname %in% c(edc_traits_log10) & cox_results_all$adjust == "adj" & cox_results_all$sample == "phy_edc_temp",]
cox_results_all_temp <- cox_results_all_temp[cox_results_all_temp$outcome %in% c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021"),]
cox_results_all_temp$outcome <- factor(cox_results_all_temp$outcome, 
                                       levels = rev(c("dm_incident_1014","ckd_incident_1014","cvd_incident_1014","cvd_incident_1021")),
                                       labels = rev(c("Incident diabetes (2010-2014)","Incident CKD (2010-2014)","Incident CVD (2010-2014)","Incident CVD (2010-2021)")))
cox_results_all_temp$exposure <- factor(cox_results_all_temp$exposure, levels = edc_traits)
cox_results_all_temp <- cox_results_all_temp %>%
  arrange(exposure,outcome)

cox_results_all_temp$hr <- exp(cox_results_all_temp$coef)
cox_results_all_temp$lci <- exp(cox_results_all_temp$coef - 1.96*cox_results_all_temp$se.coef.)
cox_results_all_temp$uci <- exp(cox_results_all_temp$coef + 1.96*cox_results_all_temp$se.coef.)
cox_results_all_temp$text <- paste0(sprintf("%.3f", cox_results_all_temp$hr)," (",sprintf("%.3f", cox_results_all_temp$lci),", ",sprintf("%.3f", cox_results_all_temp$uci),")")

### 读取汇总Logistic结果数据 ###
logistic_results_all <- readxl::read_xlsx("results/glm/logistic_results_edc_incident.xlsx")
logistic_results_all_temp <- logistic_results_all[logistic_results_all$rowname %in% c(edc_traits_log10) & logistic_results_all$adjust == "adj" & logistic_results_all$sample == "phy_edc_temp",]
logistic_results_all_temp <- logistic_results_all_temp[logistic_results_all_temp$outcome %in% c("dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b"),]
logistic_results_all_temp$outcome <- factor(logistic_results_all_temp$outcome, 
                                            levels = rev(c("dm_b","ckd_b","cvd_b","ob_b","abob_b","ir_b","dyslip_b","mets_b","nafld_b","hua_b","hpt_b","as_imt_b")),
                                            labels = rev(c("Diabetes","CKD","CVD","Obesity","Abdominal obesity","IR","Dyslipidemia","MetS","NAFLD","High UA","Hypertension","High CIMT")))
logistic_results_all_temp$exposure <- factor(logistic_results_all_temp$exposure, levels = edc_traits)
logistic_results_all_temp <- logistic_results_all_temp %>%
  arrange(exposure,outcome)

logistic_results_all_temp$or <- exp(logistic_results_all_temp$coef)
logistic_results_all_temp$lci <- exp(logistic_results_all_temp$coef - 1.96*logistic_results_all_temp$se.coef.)
logistic_results_all_temp$uci <- exp(logistic_results_all_temp$coef + 1.96*logistic_results_all_temp$se.coef.)
logistic_results_all_temp$text <- paste0(sprintf("%.3f", logistic_results_all_temp$or)," (",sprintf("%.3f", logistic_results_all_temp$lci),", ",sprintf("%.3f", logistic_results_all_temp$uci),")") 
#### 读取COX Logistic分析结果 ----

#### Forest plot (19种EDC (log10)) ----
# 构建绘制 log10 暴露森林图的函数 (根据每一个污染物独立坐标轴分开)
forest_plot1 <- function(DAT,X_TITLE){
  n_colors <- length(unique(DAT$outcome))
  j <- length(unique(DAT$exposure))  # 计算类别数
  
  f <- ggplot(data=DAT, aes(x=coef, y=outcome, col=outcome)) +
    geom_errorbar(aes(xmin=coef-1.96*se.coef., xmax=coef+1.96*se.coef., col=outcome), width=0, cex=1.5, position = position_dodge(0.7)) +
    geom_point(aes(col=outcome), cex = 5, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 7.5,
              position = position_dodge(width = 0.75)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      ) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(x=X_TITLE) +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = "Outcome",  # 修改图例标题
      values = rev(scales::hue_pal()(n_colors)), # 默认颜色
      guide = guide_legend(reverse = TRUE)  # 逆序排列
    ) +
    
    # 纵向排列面板
    facet_wrap(~exposure, ncol = 5, scales = "free_x") + # x轴自由缩放
    
    theme_classic() +
    
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),  # 四周各 1cm 边距
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 25, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          strip.text = element_text(size = 22, colour = "black"), # facet_wrap面板字体大小
          
          panel.border = element_blank(),
          
          panel.spacing.x = unit(8, "mm"),  # 横向间距
          panel.spacing.y = unit(8, "mm"), # 纵向间距
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 22, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 22, colour = "black"),
          
          legend.text = element_text(size = 22, colour = "black", margin = margin(l = 10, t = 5, b = 5)),   # 调整legend文本大小
          legend.title = element_text(size = 25, colour = "black"),
          legend.position = "none")
  
  return(f)
}
# 构建绘制 log10 暴露森林图的函数 (单个EDC for RCS)
forest_plot5 <- function(DAT,EXP,X_TITLE){
  n_colors <- length(unique(DAT$outcome))
  DAT_temp <-  DAT[DAT[["exposure"]] ==  EXP,]
  
  f <- ggplot(data=DAT_temp, aes(x=coef, y=outcome, col=outcome)) +
    geom_errorbar(aes(xmin=coef-1.96*se.coef., xmax=coef+1.96*se.coef., col=outcome), width=0, cex=1.5, position = position_dodge(0.7)) +
    geom_point(aes(col=outcome), cex = 5, position = position_dodge(0.7)) +
    geom_text(aes(label = text),
              hjust = -0.5,
              size = 7.5,
              position = position_dodge(width = 0.75)) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f", exp(x))  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
    ) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(title = sub("_log10","",EXP), x=X_TITLE) +
    
    # 将颜色、图例标题、逆序整合到 scale_colour_manual
    scale_colour_manual(
      name = " ",  # 修改图例标题
      values = c(
        "Incident diabetes (2010-2014)" = "#195e31",       
        "Incident CKD (2010-2014)" = "#01449b",
        "Incident CVD (2010-2014)" = "#ea0c22",       
        "Incident CVD (2010-2021)" = "#831a1f"
      ),
      breaks = c("Incident diabetes (2010-2014)", "Incident CKD (2010-2014)", "Incident CVD (2010-2014)", "Incident CVD (2010-2021)")
    ) +
    
    theme_classic() +
    
    theme(plot.title = element_text(size = 18),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 18, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
          axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
          
          axis.text.x = element_text(size = 18, colour = "black"), # 调整x轴文字，字体加粗
          axis.text.y = element_text(size = 18, colour = "black"),
          
          legend.text = element_text(size = 18, colour = "black"),   # 调整legend文本大小
          legend.title = element_text(size = 18, colour = "black"),
          
          legend.position = "none")
  
  return(f)
}


### Logistic 作图 (考虑19 EDC和所有二分类结果横向展示，放在extended fig) ###
logistic_log10_plot1 <- forest_plot1(logistic_results_all_temp,"OR (95% CI)")
ggsave(logistic_log10_plot1, filename=paste0("figures/supplementary_figures/(efig6)_logistic_edc1_(phy_edc_temp).pdf"), width = 18, height = 24, limitsize = FALSE)

### COX 作图 (所有EDC单独分开，组合RCS) ###
forest_plot_list <- list()
for (i in edc_traits) {
  forest_plot_list[[i]] <- forest_plot5(cox_results_all_temp, i, "HR (95% CI)")
}
#### Forest plot (19种EDC (log10)) ####


############################################# RCS + forest图 #############################################
#### 读取绘制RCS图 + forest plot ----
# 构建绘制 RCS 图的函数
rcs_plot <- function(DAT,EXP,OUT,COLOR){
  nonlinear <- anova(DAT) # 检验自变量与因变量的非线性关系
  nonlinear <- data.frame(nonlinear)
  nonlinear$P <- ifelse(nonlinear$P < 0.001, "<0.001", paste0("=", sprintf("%.3f", round(nonlinear$P,3))))
  
  # 预测风险比（HR）并绘制RCS曲线
  HR <- Predict(DAT, EDC, fun = exp, type = "predictions", ref.zero = TRUE, np = 200)
  
  # 创建图形
  f <- ggplot() +
    
    # HR 图层
    geom_line(data=HR, aes(EDC, yhat), linetype="solid", size=1, alpha=0.7, colour=COLOR) +
    geom_ribbon(data=HR, aes(EDC, ymin=lower, ymax=upper), alpha=0.1, fill=COLOR) +
    
    geom_hline(yintercept=0, linetype = "dashed", colour = "black", linewidth = 1) +
    
    scale_x_continuous(
      labels = function(x) sprintf("%.1f",x)  # 显示为 x，保留1位小数 (保留小数点后最后一位的0)
    ) +
    scale_y_continuous(
      labels = function(y) sprintf("%.1f",exp(y))  # 显示为 exp(y)，保留1位小数 (保留小数点后最后一位的0)
    ) +
    
    # 添左上角左对齐文本
    annotate("text",
             x = -Inf, y = Inf,
             label = OUT,
             hjust = -0.05,  # 0左对齐 (数字越小越靠右)
             vjust = 1.35,   # 0上对齐 (数字越小越靠上)
             size = 5,       # 字体大小
             colour = "black") +
    
    # 添加右上角右对齐文本
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("  P",nonlinear[1,3]),
             hjust = 1.15,   # 0右对齐
             vjust = 2,      # 0上对齐
             size = 5,       # 字体大小
             colour = "black") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("P-nonlinear",nonlinear[2,3]),
             hjust = 1.07,   # 0右对齐
             vjust = 4,      # 0上对齐
             size = 5,       # 字体大小
             colour = "black") +
    
    labs(x = paste0("Log10(",sub("_log10","",EXP),") (ng/mL)"), y = "HR") +
    
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1), # 面板区域边框
          
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          
          axis.text.y = element_text(color= "black", size = 18),
          axis.text.x = element_text(color= "black", size = 18),
          
          axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 0, b = 0, l = 5)),
          axis.title.x = element_text(size = 18, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          
          plot.title = element_blank()
    )
  
  return(f)
}


### 读取汇总RCS结果数据 ###
rcs_results_list <- readRDS("results/cox/rcs_results.rds") 

rcs_plot_all <- list()
rcs_plot_knot3 <- list()
rcs_plot_knot4 <- list()
# 创建单个EDC的单个结局图
for (i in c("phy_edc_temp")) {
  # i <- "phy_edc_temp"
  
  sample_name <- i  # 提取subgroup的名称
  
  for (j in 1:length(rcs_results_list)) {
    # j<-1
    
    rcs_fit <- rcs_results_list[[j]]
    list_name <- names(rcs_results_list[j])
    
    sample <- sub("\\|.*", "", list_name)
    exp <- sub(".*?\\|(.*?)\\|.*", "\\1", list_name)
    out <- sub(".*?\\|.*?\\|(.*?)\\|.*", "\\1", list_name)
    knot <- sub(".*\\|", "", list_name)
    
    if(!exp %in% edc_traits_log10){
      next
    }
    
    if(out == "cvd_incident_1014"){
      out <- "Incident CVD\n(2010-2014)"
    }else if(out == "cvd_incident_1021"){
      out <- "Incident CVD\n(2010-2021)"
    }else if(out == "ckd_incident_1014"){
      out <- "Incident CKD\n(2010-2014)"
    }else if(out == "dm_incident_1014"){
      out <- "Incident diabetes\n(2010-2014)"
    }
    
    
    plot_temp <- rcs_plot(rcs_fit,exp,out,my_colors_list[[exp]])
    
    if(knot == 3){
      rcs_plot_knot3[[list_name]] <- plot_temp
    }else{
      rcs_plot_knot4[[list_name]] <- plot_temp
    }
    
    rcs_plot_all[[list_name]] <- plot_temp
  }
}


# 创建单个EDC的3个结局图
rcs_plot_knot3_incident <- list()
rcs_plot_knot4_incident <- list()
for (i in 1:19) {
  rcs_plot_knot3_incident[[edc_traits[i]]] <- cowplot::plot_grid(forest_plot_list[[i]], NULL,
                                                                 rcs_plot_knot3[[(i-1)*4+1]], NULL,
                                                                 rcs_plot_knot3[[(i-1)*4+2]], NULL,
                                                                 rcs_plot_knot3[[(i-1)*4+3]], NULL,
                                                                 rcs_plot_knot3[[(i-1)*4+4]], ncol = 9, rel_widths = c(1.5, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1))

  rcs_plot_knot4_incident[[edc_traits[i]]] <- cowplot::plot_grid(forest_plot_list[[i]], NULL,
                                                                 rcs_plot_knot4[[(i-1)*4+1]], NULL,
                                                                 rcs_plot_knot4[[(i-1)*4+2]], NULL,
                                                                 rcs_plot_knot4[[(i-1)*4+3]], NULL,
                                                                 rcs_plot_knot4[[(i-1)*4+4]], ncol = 9, rel_widths = c(1.5, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1))
}


# 创建19个EDC的4个结局汇总图
rcs_19edc_knot3 <- cowplot::plot_grid(plotlist = rcs_plot_knot3_incident[1:19],  # 提取前19个子图
                                      nrow = 19,                                 # 指定19行
                                      rel_heights = rep(1, 19),                  # 每行高度权重相同
                                      align = "v"                                # 垂直对齐（可选）
)

rcs_19edc_knot4 <- cowplot::plot_grid(plotlist = rcs_plot_knot4_incident[1:19],  # 提取前19个子图
                                      nrow = 19,                                 # 指定19行
                                      rel_heights = rep(1, 19),                  # 每行高度权重相同
                                      align = "v"                                # 垂直对齐（可选）
)

ggsave(rcs_19edc_knot4, filename=paste0("figures/supplementary_figures/(efig5)_rcs_knot4_(phy_edc_temp)1.pdf"), width = 22, height = 65, limitsize = FALSE)
#### 读取绘制RCS图 + forest plot ####
