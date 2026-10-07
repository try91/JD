library(data.table)
library(dplyr)
library(ggplot2)

setwd("C:/TWang/DLiu/EDC_Micro/") # Windows路径

#### 构建函数把species_name转换为可以作图的标准化名称 ----
trans_mp4_s_names <- function(STRING){
  
  string_temp <- gsub("s__","",STRING)
  string_temp <- gsub("_log10","",string_temp)
  # # 删除所有"_" (一般情况下第一个"_"分隔属和种)
  # string_temp <- gsub("_"," ",string_temp)
  # 删除"_" (一般情况下第一个"_"分隔属和种)
  string_temp <- sub("_"," ",string_temp) # 使用 sub()替换第一个匹配项。
  # 删除特殊位置的"_"
  string_temp <- gsub("_bacterium"," bacterium",string_temp)
  string_temp <- gsub("bacterium_","bacterium ",string_temp)
  string_temp <- gsub("_unclassified"," unclassified",string_temp)
  string_temp <- gsub("unclassified_","unclassified ",string_temp)
  string_temp <- gsub("copri_","copri ",string_temp)
  string_temp <- gsub("Cibionibacter_","Cibionibacter ",string_temp)
  string_temp <- gsub("Marseille_","Marseille ",string_temp)
  string_temp <- gsub("oral_taxon_","oral taxon ",string_temp)
  # # 添加特殊位置的"_"
  # string_temp <- gsub(" Family","_Family",string_temp)
  # string_temp <- gsub(" phylum","_phylum",string_temp)
  # string_temp <- gsub("Candidatus ","Candidatus_",string_temp)
  # 处理sp
  string_temp <- gsub("sp_","sp. ",string_temp)
  string_temp <- gsub("_sp."," sp.",string_temp)
  return(string_temp)
}
#### 构建函数把species_name转换为可以作图的标准化名称 ####

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

# 转换后的菌的名称
mp4_s_bin <- paste0(mp4_s_names,"_bin") # 菌群MP4出现与否的分类变量 (物种层面)
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)
mp4_s_zero <- paste0(mp4_s_names,"_zero") # 菌群MP4填补0值丰度 (物种层面)

mp4_g_bin <- paste0(mp4_g_names,"_bin") # 菌群MP4出现与否的分类变量 (属层面)
mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)
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
                 "TCC","TCS",
                 "BPA") # 检出率>50%
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

#### 读取+处理interaction分析结果 ----
# "相乘交互作用"显著(Multiplicative scale P<0.05), 同时"相加交互作用"显著(以RERI显著为标准)的组
results_interaction_sig1 <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig_20260728.xlsx")
# "相乘交互作用"显著的组(Multiplicative scale P<0.05)
results_interaction_sig2 <- readxl::read_xlsx("results/cox/interaction/multi_interaction_sig_20260728.xlsx")
results_interaction_sig2 <- results_interaction_sig2[!results_interaction_sig2$keep_exp_med_out %in% results_interaction_sig1$keep_exp_med_out,] # 排除相乘和相加交互同时显著的结果
# "相加交互作用"显著的组(以RERI显著为标准)
results_interaction_sig3 <- readxl::read_xlsx("results/cox/interaction/add_interaction_sig_20260728.xlsx")
results_interaction_sig3 <- results_interaction_sig3[!results_interaction_sig3$keep_exp_med_out %in% results_interaction_sig1$keep_exp_med_out,] # 排除相乘和相加交互同时显著的结果

results_interaction_sig1 <- results_interaction_sig1[results_interaction_sig1$mediator %in% mp4_s_log10_short,]


results_interaction_sig1_dm <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "DM",]
results_interaction_sig1_dm$hr_95ci <- sprintf("%.3f (%.3f, %.3f)", results_interaction_sig1_dm$exp.coef., results_interaction_sig1_dm$exp.lci, results_interaction_sig1_dm$exp.uci)
results_interaction_sig1_dm$text <- paste0(results_interaction_sig1_dm$rowname,": ",results_interaction_sig1_dm$hr_95ci)
results_interaction_sig1_dm$text <- ifelse(results_interaction_sig1_dm$text == "SI: NA (NA, NA)", "", results_interaction_sig1_dm$text)
results_interaction_sig1_dm$p_text <- ifelse(results_interaction_sig1_dm$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_dm$Pr...z..,3)))
results_interaction_sig1_dm <- results_interaction_sig1_dm[results_interaction_sig1_dm$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_dm <- results_interaction_sig1_dm[!(results_interaction_sig1_dm$rowname %in% c("EDC high") & results_interaction_sig1_dm$group %in% c("interaction_recode")),]

results_interaction_sig1_ckd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CKD",]
results_interaction_sig1_ckd$hr_95ci <- sprintf("%.3f (%.3f, %.3f)", results_interaction_sig1_ckd$exp.coef., results_interaction_sig1_ckd$exp.lci, results_interaction_sig1_ckd$exp.uci)
results_interaction_sig1_ckd$text <- paste0(results_interaction_sig1_ckd$rowname,": ",results_interaction_sig1_ckd$hr_95ci)
results_interaction_sig1_ckd$text <- ifelse(results_interaction_sig1_ckd$text == "SI: NA (NA, NA)", "", results_interaction_sig1_ckd$text)
results_interaction_sig1_ckd$p_text <- ifelse(results_interaction_sig1_ckd$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_ckd$Pr...z..,3)))
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[!(results_interaction_sig1_ckd$rowname %in% c("EDC high") & results_interaction_sig1_ckd$group %in% c("interaction_recode")),]

results_interaction_sig1_cvd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CVD",]
results_interaction_sig1_cvd$hr_95ci <-sprintf("%.3f (%.3f, %.3f)", results_interaction_sig1_cvd$exp.coef., results_interaction_sig1_cvd$exp.lci, results_interaction_sig1_cvd$exp.uci)
results_interaction_sig1_cvd$text <- paste0(results_interaction_sig1_cvd$rowname,": ",results_interaction_sig1_cvd$hr_95ci)
results_interaction_sig1_cvd$text <- ifelse(results_interaction_sig1_cvd$text == "SI: NA (NA, NA)", NA, results_interaction_sig1_cvd$text)
results_interaction_sig1_cvd$p_text <- ifelse(results_interaction_sig1_cvd$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_cvd$Pr...z..,3)))
results_interaction_sig1_cvd <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_cvd <- results_interaction_sig1_cvd[!(results_interaction_sig1_cvd$rowname %in% c("EDC high") & results_interaction_sig1_cvd$group %in% c("interaction_recode")),]



results_interaction_sig1_dm$point_color <- ifelse(results_interaction_sig1_dm$rowname == "EDC high" & results_interaction_sig1_dm$group == "mp4_low", "cat1", 
                                                  ifelse(results_interaction_sig1_dm$rowname == "EDC high" & results_interaction_sig1_dm$group == "mp4_high", "cat2", 
                                                         ifelse(results_interaction_sig1_dm$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_dm$point_color <- factor(results_interaction_sig1_dm$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_dm <- results_interaction_sig1_dm %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_dm1 <- results_interaction_sig1_dm[results_interaction_sig1_dm$direction == "1_3",]
results_interaction_sig1_dm2 <- results_interaction_sig1_dm[results_interaction_sig1_dm$direction == "1_4", "text"]
colnames(results_interaction_sig1_dm2) <- "additive_text"
results_interaction_sig1_dm3 <- cbind(results_interaction_sig1_dm1,results_interaction_sig1_dm2)
results_interaction_sig1_dm3$hr_95ci <- paste0("HR: ",results_interaction_sig1_dm3$hr_95ci)


results_interaction_sig1_ckd$point_color <- ifelse(results_interaction_sig1_ckd$rowname == "EDC high" & results_interaction_sig1_ckd$group == "mp4_low", "cat1", 
                                                   ifelse(results_interaction_sig1_ckd$rowname == "EDC high" & results_interaction_sig1_ckd$group == "mp4_high", "cat2", 
                                                          ifelse(results_interaction_sig1_ckd$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_ckd$point_color <- factor(results_interaction_sig1_ckd$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_ckd <- results_interaction_sig1_ckd %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_ckd1 <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$direction == "1_3",] 
results_interaction_sig1_ckd2 <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$direction == "1_4", "text"]
colnames(results_interaction_sig1_ckd2) <- "additive_text"
results_interaction_sig1_ckd3 <- cbind(results_interaction_sig1_ckd1,results_interaction_sig1_ckd2)
results_interaction_sig1_ckd3$hr_95ci <- paste0("HR: ",results_interaction_sig1_ckd3$hr_95ci)


results_interaction_sig1_cvd$point_color <- ifelse(results_interaction_sig1_cvd$rowname == "EDC high" & results_interaction_sig1_cvd$group == "mp4_low", "cat1", 
                                                   ifelse(results_interaction_sig1_cvd$rowname == "EDC high" & results_interaction_sig1_cvd$group == "mp4_high", "cat2", 
                                                          ifelse(results_interaction_sig1_cvd$rowname == "EDC high:GM high", "cat3", "cat4")))
results_interaction_sig1_cvd$point_color <- factor(results_interaction_sig1_cvd$point_color, levels = c("cat1","cat2","cat3","cat4"))
results_interaction_sig1_cvd <- results_interaction_sig1_cvd %>%
  arrange(exposure,mediator,point_color)
results_interaction_sig1_cvd1 <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$direction == "1_3",] 
results_interaction_sig1_cvd2 <- results_interaction_sig1_cvd[results_interaction_sig1_cvd$direction == "1_4", "text"]
colnames(results_interaction_sig1_cvd2) <- "additive_text"
results_interaction_sig1_cvd3 <- cbind(results_interaction_sig1_cvd1,results_interaction_sig1_cvd2)
results_interaction_sig1_cvd3$hr_95ci <- paste0("HR: ",results_interaction_sig1_cvd3$hr_95ci)
#### 读取+处理interaction分析结果 ####

############################################# 森林图 (根据"both_interaction_dm_ckd_20251203.pdf" 挑选4个代表) #############################################
#### forest plot (挑选4个代表, 高低GM组中EDC对结局的影响 + multiplicative & additive interaction) ----
# 构建绘制森林图的函数
forest_function <- function(DAT,SUBTITLE,TITLE){
  
  DAT$point_color <- factor(DAT$point_color, levels = rev(c("cat2","cat1","cat3","cat4")))
  
  DAT$coef <- log(DAT$exp.coef.)
  DAT$lci <- log(DAT$exp.lci)
  DAT$uci <- log(DAT$exp.uci)
  
  f <- ggplot(data=DAT, aes(x=coef, y=exposure, col=point_color)) +
    geom_errorbar(aes(xmin=lci, xmax=uci, col=point_color), width=0, cex=0.8, position = position_dodge(width = 0.75)) +
    geom_point(aes(col=point_color), cex = 2.4, position = position_dodge(width = 0.75)) +
    
    geom_text(aes(label = hr_95ci, col=point_color),
              hjust = 1.5,
              size = 4,
              position = position_dodge(width = 0.75)) +
    
    geom_text(aes(label = p_text, col=point_color),
              hjust = -1,
              size = 4,
              position = position_dodge(width = 0.75)) +
    
    geom_text(aes(label = additive_text, col=point_color),
              hjust = -2.5,
              size = 4,
              position = position_dodge(width = 0.75)) +
    
    
    # 关键：expand 控制轴线与数据两端空白
    scale_x_continuous(
      # labels = function(x) round(exp(x), 1),  # 显示为 exp(x)，保留1位小数 (不保留小数点后最后一位的0)
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = expansion(mult = c(0.05, 0.05))  # x 轴左右留白 5%
      ) + 
    
    scale_y_discrete(
      expand = expansion(mult = c(0.1, 0.1))  # y 轴上下留白 10%
      ) +
  
  
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    # 将颜色、图例标题、整合到 scale_colour_manual
    scale_colour_manual(
      name = "",  # legend title
      
      # 设置颜色
      # values = c("cat1" = "#f1c40f",
      #            "cat2" = "#d63031",
      #            "cat3" = "#5f27cd",
      #            "cat4" = "#48dbfb"),
      values = c("cat2" = "#5e2a2a",
                 "cat1" = "#b5be26",
                 "cat3" = "#5f27cd",
                 "cat4" = "#f1c40f"),
      # 直接在scale_color_manual中设置标签
      labels = c("cat2" = "High abundance",
                 "cat1" = "Low abundance",
                 "cat3" = "Multiplicative scale",
                 "cat4" = "Multiplicative scale"),
      # 通过breaks参数明确指定图例顺序
      breaks = c("cat2",
                 "cat1",
                 "cat3",
                 "cat4")
    ) +
    
    labs(title = TITLE,
         subtitle = SUBTITLE,
         x="HR (95%CI)") +
    
    theme_classic() +
    theme(
      # plot.margin = margin(1, 1, 1, 1, "cm"),  # 四周各 1cm 边距
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 14),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 14, margin = margin(t = 5, r = 0, b = 0, l = 0)),
      
      panel.spacing.x = unit(8, "mm"),  # 横向间距
      panel.spacing.y = unit(8, "mm"), # 纵向间距
      
      axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      
      axis.text.x = element_text(size = 12, colour = "black"), # 调整x轴文字
      # axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.y = element_blank(),
      
      # legend.position = "bottom",
      # legend.title = element_text(size = 15, colour = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 10, colour = "black")   # 调整legend文本大小
    )
  
  return(f)
}

# 绘制森林图比较构建的EDC index和qgcomp结果
results_interaction_sig1_dm3$mediator <- trans_mp4_s_names(results_interaction_sig1_dm3$mediator)
results_interaction_sig1_dm4 <- results_interaction_sig1_dm3[results_interaction_sig1_dm3$exposure %in% c("edc_score_edc14_f") & results_interaction_sig1_dm3$mediator %in% c("Prevotella copri clade_A") |
                                                               results_interaction_sig1_dm3$exposure %in% c("edc_score_edc14_f") & results_interaction_sig1_dm3$mediator %in% c("Flavonifractor plautii"),]
dm_plot_list <- list()
for (i in c("Flavonifractor plautii","Prevotella copri clade_A")) {
  dat_for_plot <- results_interaction_sig1_dm4[results_interaction_sig1_dm4$mediator == i,]
  f_forest1 <- forest_function(dat_for_plot, i, paste0("Risk of diabetes associated with high EDC Scorequartile (",dat_for_plot$exposure[[2]],")"))
  dm_plot_list[[i]] <- f_forest1
}


results_interaction_sig1_ckd3$mediator <- trans_mp4_s_names(results_interaction_sig1_ckd3$mediator)
results_interaction_sig1_ckd4 <- results_interaction_sig1_ckd3[results_interaction_sig1_ckd3$exposure %in% c("edc_score_edc14_f") & results_interaction_sig1_ckd3$mediator %in% c("Anaerobutyricum hallii") |
                                                                 results_interaction_sig1_ckd3$exposure %in% c("MiBP_log10") & results_interaction_sig1_ckd3$mediator %in% c("Granulicatella SGB8255"),]
ckd_plot_list <- list()
for (i in c("Anaerobutyricum hallii","Granulicatella SGB8255")) {
  dat_for_plot <- results_interaction_sig1_ckd4[results_interaction_sig1_ckd4$mediator == i,]
  f_forest2 <- forest_function(dat_for_plot, i, paste0("Risk of CKD associated with high EDC Scorequartile (",dat_for_plot$exposure[[2]],")"))
  ckd_plot_list[[i]] <- f_forest2
}


results_interaction_sig1_cvd3$mediator <- trans_mp4_s_names(results_interaction_sig1_cvd3$mediator)
results_interaction_sig1_cvd4 <- results_interaction_sig1_cvd3[results_interaction_sig1_cvd3$exposure %in% c("edc_score_pae6_f") & results_interaction_sig1_cvd3$mediator %in% c("Dysosmobacter sp. BX15") |
                                                                 results_interaction_sig1_cvd3$exposure %in% c("edc_score_pae6_f") & results_interaction_sig1_cvd3$mediator %in% c("Clostridium SGB4751"),]
cvd_plot_list <- list()
for (i in c("Dysosmobacter sp. BX15","Clostridium SGB4751")) {
  dat_for_plot <- results_interaction_sig1_cvd4[results_interaction_sig1_cvd4$mediator == i,]
  f_forest3 <- forest_function(dat_for_plot, i, paste0("Risk of CVD associated with high EDC Scorequartile (",dat_for_plot$exposure[[2]],")"))
  cvd_plot_list[[i]] <- f_forest3
}


f_all <- cowplot::plot_grid(plotlist = c(dm_plot_list, ckd_plot_list, cvd_plot_list), 
                            nrow = 6, rel_heights = c(1, 1, 1, 1, 1, 1),
                            align = "hv")
ggsave(paste0("figures/main_figures/forest_interaction_dm_ckd_cvd_20260702.pdf"),
       f_all, width = 3.8, height = 12, limitsize = FALSE)
#### forest plot (挑选6个代表, 高低GM组中EDC对结局的影响 + multiplicative & additive interaction) ####
