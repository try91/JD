library(data.table)
library(dplyr)
library(tidyr)
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


#### 读取 (JD HRB) Replication结果 ----
## JD结果 ##
{
  results_jd <- 
    results_jd <- 
    
  results_jd <- read.csv("results/replication/mp4_s_for_replication/disease_interation_related_mp4_species_20260728.csv")
  
  dm_jd <- results_jd[,c("species","estimate_dm","se_dm","p_dm")]
  colnames(dm_jd) <- c("species","estimate_dm_jd","se_dm_jd","p_dm_jd")
  dm_jd$estimate_dm_jd <- as.numeric(dm_jd$estimate_dm_jd)
  dm_jd$p_dm_jd <- as.numeric(dm_jd$p_dm_jd)
  nrow(dm_jd[dm_jd$p_dm_jd < 0.05,]) # 35 dm-related sig species in JD (no self-report cases)
  
  ckd_jd <- results_jd[,c("species","estimate_ckd","se_ckd","p_ckd")]
  colnames(ckd_jd) <- c("species","estimate_ckd_jd","se_ckd_jd","p_ckd_jd")
  ckd_jd$estimate_ckd_jd <- as.numeric(ckd_jd$estimate_ckd_jd)
  ckd_jd$p_ckd_jd <- as.numeric(ckd_jd$p_ckd_jd)
  nrow(ckd_jd[ckd_jd$p_ckd_jd < 0.05,]) # 84 ckd-related sig species in JD
  
  cvd_jd <- results_jd[,c("species","estimate_cvd","se_cvd","p_cvd")]
  colnames(cvd_jd) <- c("species","estimate_cvd_jd","se_cvd_jd","p_cvd_jd")
  cvd_jd$estimate_cvd_jd <- as.numeric(cvd_jd$estimate_cvd_jd)
  cvd_jd$p_cvd_jd <- as.numeric(cvd_jd$p_cvd_jd)
  nrow(cvd_jd[cvd_jd$p_cvd_jd < 0.05,]) # 26 cvd-related sig species in JD
}

## HEB结果 ##
{
  dm_hrb <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "T2D")
  dm_hrb <- dm_hrb[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(dm_hrb) <- c("species","estimate_dm_hrb","se_dm_hrb","p_dm_hrb","p_fdr_dm_hrb")
  dm_hrb <- dm_hrb %>%
    mutate(p_adj_bh_dm_hrb = p.adjust(p_dm_hrb, method = "BH"))
  
  ckd_hrb <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "CKD")
  ckd_hrb <- ckd_hrb[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(ckd_hrb) <- c("species","estimate_ckd_hrb","se_ckd_hrb","p_ckd_hrb","p_fdr_ckd_hrb")
  ckd_hrb <- ckd_hrb %>%
    mutate(p_adj_bh_ckd_hrb = p.adjust(p_ckd_hrb, method = "BH"))
}

## JieZ ACVD结果 ##
{
  cvd_gd <- readxl::read_xlsx("results/replication/T2D_CKD_ASCVD_validation_com_20260726_n359.xlsx", sheet = "CVD") # 之前由分组变量设置了“Control”和“ASCVD”默认ASCVD作为ref，导致所有结果反向，已更正
  cvd_gd <- cvd_gd[,c("...1","logit.coeff","logit.coeff.SE","logit.coeff.pval","logit.fdr")]
  colnames(cvd_gd) <- c("species","estimate_cvd_gd","se_cvd_gd","p_cvd_gd","p_fdr_cvd_gd")
  cvd_gd <- cvd_gd %>%
    mutate(p_adj_bh_cvd_gd = p.adjust(p_cvd_gd, method = "BH"))
}
#### 读取 (JD HRB) Replication结果 ####

#### (JD HRB) Replication结果处理 ----
# DM结果处理 #
{
  results_dm <- left_join(dm_jd, dm_hrb, by = "species") %>%
    filter(p_dm_jd < 0.05) %>%
    arrange(desc(estimate_dm_jd))
    
  # 筛选JD和HRB方向一致的结果
  results_dm_same_direction <- results_dm[results_dm$estimate_dm_jd * results_dm$estimate_dm_hrb > 0,]
  # results_dm_same_direction_hrb_sig <- results_dm_same_direction[results_dm_same_direction$p_dm_hrb < 0.05,]
  results_dm_same_direction_hrb_sig_fdr <- results_dm_same_direction[results_dm_same_direction$p_fdr_dm_hrb < 0.2,]
  
  results_dm_estimate <- results_dm[,c("species","estimate_dm_jd","estimate_dm_hrb")]
  colnames(results_dm_estimate) <- c("species","JD","HRB")
  
  results_dm_se <- results_dm[,c("species","se_dm_jd","se_dm_hrb")]
  colnames(results_dm_se) <- c("species","JD","HRB")
  
  results_dm_p <- results_dm[,c("species","p_dm_jd","p_dm_hrb")]
  colnames(results_dm_p) <- c("species","JD","HRB")
  
  # 转换成long data
  results_dm_estimate_long <- results_dm_estimate %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "estimate"           # 新列，存放原列中的数值
    )
  results_dm_se_long <- results_dm_se %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "se"           # 新列，存放原列中的数值
    )
  results_dm_p_long <- results_dm_p %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "p"           # 新列，存放原列中的数值
    )
  results_dm_long <- left_join(results_dm_estimate_long, results_dm_se_long, by=c("species","study")) %>%
    left_join(results_dm_p_long, by=c("species","study"))
}

# CKD结果处理 #
{
  results_ckd <- left_join(ckd_jd, ckd_hrb, by = "species") %>%
    filter(p_ckd_jd < 0.05) %>%
    arrange(desc(estimate_ckd_jd))
    
  # 筛选JD和HRB方向一致的结果
  results_ckd_same_direction <- results_ckd[results_ckd$estimate_ckd_jd * results_ckd$estimate_ckd_hrb > 0,]
  # results_ckd_same_direction_hrb_sig <- results_ckd_same_direction[results_ckd_same_direction$p_ckd_hrb < 0.05,]
  results_ckd_same_direction_hrb_sig_fdr <- results_ckd_same_direction[results_ckd_same_direction$p_fdr_ckd_hrb < 0.2,]
  
  results_ckd_estimate <- results_ckd[,c("species","estimate_ckd_jd","estimate_ckd_hrb")]
  colnames(results_ckd_estimate) <- c("species","JD","HRB")
  
  results_ckd_se <- results_ckd[,c("species","se_ckd_jd","se_ckd_hrb")]
  colnames(results_ckd_se) <- c("species","JD","HRB")
  
  results_ckd_p <- results_ckd[,c("species","p_ckd_jd","p_ckd_hrb")]
  colnames(results_ckd_p) <- c("species","JD","HRB")
  # 转换成long data
  results_ckd_estimate_long <- results_ckd_estimate %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "estimate"           # 新列，存放原列中的数值
    )
  results_ckd_se_long <- results_ckd_se %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "se"           # 新列，存放原列中的数值
    )
  results_ckd_p_long <- results_ckd_p %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "p"           # 新列，存放原列中的数值
    )
  results_ckd_long <- left_join(results_ckd_estimate_long, results_ckd_se_long, by=c("species","study")) %>%
    left_join(results_ckd_p_long, by=c("species","study"))
}

# ASCVD结果处理 #
{
  results_acvd <- left_join(cvd_jd, cvd_gd, by = "species") %>%
    filter(p_cvd_jd < 0.05) %>%
    arrange(desc(estimate_cvd_jd))
    
  # 筛选JD和GD方向一致的结果
  results_acvd_same_direction <- results_acvd[results_acvd$estimate_cvd_jd * results_acvd$estimate_cvd_gd > 0,]
  # results_acvd_same_direction_acvd_sig <- results_acvd_same_direction[results_acvd_same_direction$p_cvd_gd < 0.05,]
  results_acvd_same_direction_acvd_sig_fdr <- results_acvd_same_direction[results_acvd_same_direction$p_fdr_cvd_gd < 0.2,]
  
  results_acvd_estimate <- results_acvd[,c("species","estimate_cvd_jd","estimate_cvd_gd")]
  colnames(results_acvd_estimate) <- c("species","JD","GD")
  
  results_acvd_se <- results_acvd[,c("species","se_cvd_jd","se_cvd_gd")]
  colnames(results_acvd_se) <- c("species","JD","GD")
  
  results_acvd_p <- results_acvd[,c("species","p_cvd_jd","p_cvd_gd")]
  colnames(results_acvd_p) <- c("species","JD","GD")
  # 转换成long data
  results_acvd_estimate_long <- results_acvd_estimate %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "estimate"           # 新列，存放原列中的数值
    )
  results_acvd_se_long <- results_acvd_se %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "se"           # 新列，存放原列中的数值
    )
  results_acvd_p_long <- results_acvd_p %>%
    pivot_longer(
      cols = -species,                   # 保留xx列不动
      names_to = "study",              # 新列，存放原列名（如“数学_期中”）
      values_to = "p"           # 新列，存放原列中的数值
    )
  results_acvd_long <- left_join(results_acvd_estimate_long, results_acvd_se_long, by=c("species","study")) %>%
    left_join(results_acvd_p_long, by=c("species","study"))
}
#### (JD HRB) Replication结果处理 ----

#### 保留验证队列方向一致的结果 ----
# for fig 6a #
results_dm_same_direction_for_fig6a <- results_dm_same_direction
results_ckd_same_direction_for_fig6a <- results_ckd_same_direction
results_acvd_same_direction_for_fig6a <- results_acvd_same_direction

results_dm_same_direction_for_fig6a$species <- paste0(results_dm_same_direction_for_fig6a$species,"_log10")
results_ckd_same_direction_for_fig6a$species <- paste0(results_ckd_same_direction_for_fig6a$species,"_log10")
results_acvd_same_direction_for_fig6a$species <- paste0(results_acvd_same_direction_for_fig6a$species,"_log10")

openxlsx::write.xlsx(results_dm_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_dm_20260728.xlsx")
openxlsx::write.xlsx(results_ckd_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_ckd_20260728.xlsx")
openxlsx::write.xlsx(results_acvd_same_direction_for_fig6a, "results/replication/replication_same_direction/replication_same_direction_acvd_20260728.xlsx")

# 验证队列中方向一致且FDR-P<0.2的结果
results_dm_same_direction_hrb_sig <- results_dm_same_direction_hrb_sig_fdr
results_ckd_same_direction_hrb_sig <- results_ckd_same_direction_hrb_sig_fdr
results_acvd_same_direction_acvd_sig <- results_acvd_same_direction_acvd_sig_fdr

results_dm_same_direction_hrb_sig$species <- paste0(results_dm_same_direction_hrb_sig$species,"_log10")
results_ckd_same_direction_hrb_sig$species <- paste0(results_ckd_same_direction_hrb_sig$species,"_log10")
results_acvd_same_direction_acvd_sig$species <- paste0(results_acvd_same_direction_acvd_sig$species,"_log10")

openxlsx::write.xlsx(results_dm_same_direction_hrb_sig, "results/replication/replication_same_direction/replication_same_direction_sig_dm_20260728.xlsx")
openxlsx::write.xlsx(results_ckd_same_direction_hrb_sig, "results/replication/replication_same_direction/replication_same_direction_sig_ckd_20260728.xlsx")
openxlsx::write.xlsx(results_acvd_same_direction_acvd_sig, "results/replication/replication_same_direction/replication_same_direction_sig_acvd_20260728.xlsx")
#### 保留验证队列方向一致的结果 ####

#### 作图Forest (JD HRB) ----
# Logistic分析结果(MP4-DM CKD)作图 #
# 构建绘制森林图的函数
forest_function <- function(DAT, TITLE){
  
  f <- ggplot(data=DAT, aes(x=estimate, y=species)) +
    
    geom_errorbar(aes(xmin=estimate-1.96*se, xmax=estimate+1.96*se, 
                      group = study,              # 关键1：按study分组
                      color = errorbar_color),    # 关键2：误差线颜色映射到errorbar_color变量
                  width=0, cex=0.9, position = position_dodge(0.75)) +
    
    geom_point(aes(group = study,            # 关键1：按study分组
                   color = point_color),     # 关键2：点颜色映射到point_color变量
               size = 3.8,
                   position = position_dodge(0.75)) +
    
    geom_text(aes(label = text, group = study),
              hjust = -0.5,
              size = 6.6,
              position = position_dodge(width = 0.9)) +
    
    # 将颜色、图例标题、整合到 scale_colour_manual
    scale_colour_manual(
      name = "Cohort",  # legend title
      
      # 设置颜色
      values = c("jd" = "#8e44ad",
                 "rep" = "#2980b9"),
      # 直接在scale_color_manual中设置标签
      labels = c("jd" = "JD_sub", # "Discovery: JD_sub"
                 "rep" = "Validation"), # "Validation: (a) HRB_sub1 (diabetes), (b) HRB_sub2 (CKD), (c) ACVD_2017 (CVD)"
      # 通过breaks参数明确指定图例顺序
      breaks = c("jd", "rep")
    ) +
    
    scale_x_continuous(
      # labels = function(x) round(exp(x), 1),  # 显示为 exp(x)，保留1位小数 (不保留小数点后最后一位的0)
      labels = function(x) sprintf("%.1f", exp(x)),  # 显示为 exp(x)，保留1位小数 (保留小数点后最后一位的0)
      expand = c(0.1, 0.1)
    ) +
    
    geom_vline(xintercept = 0,
               linetype = "dashed",
               linewidth = 0.5) +
    
    labs(title = TITLE,
         x="OR/HR (95%CI)") +
    
    theme_classic() +
    
    theme(
      # plot.margin = margin(1, 1, 1, 1, "cm"),  # 四周各 1cm 边距
      plot.title = element_text(size = 25),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 25, margin = margin(t = 5, r = 0, b = 0, l = 0)),
      
      panel.spacing.x = unit(8, "mm"),  # 横向间距
      panel.spacing.y = unit(8, "mm"), # 纵向间距
      
      axis.ticks = element_line(size = 0.8, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.15, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      
      axis.text.x = element_text(size = 25, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 23, colour = "black"),
      
      # legend.position = "bottom",
      # legend.title = element_text(size = 15, colour = "black"),
      legend.key.height = unit(10, "mm"),
      legend.key.width = unit(10, "mm"),
      legend.title = element_blank(),
      legend.text = element_text(size = 25, colour = "black")   # 调整legend文本大小
    )
  
  return(f)
}

# DM数据整理for plot #
{
  dat_dm_for_forest <- results_dm_long[results_dm_long$species %in% results_dm_same_direction_hrb_sig_fdr$species,]
  dat_dm_for_forest$hror <- exp(dat_dm_for_forest$estimate)
  dat_dm_for_forest$lci <- exp(dat_dm_for_forest$estimate - 1.96*dat_dm_for_forest$se)
  dat_dm_for_forest$uci <- exp(dat_dm_for_forest$estimate + 1.96*dat_dm_for_forest$se)
  
  ### 如果p<0.05需要确保UCI<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_dm_for_forest[dat_dm_for_forest$p < 0.05 & (dat_dm_for_forest$uci >= 0.9995 & dat_dm_for_forest$uci < 1),]
    dat_dm_for_forest$uci <- ifelse(dat_dm_for_forest$p < 0.05 & (dat_dm_for_forest$uci >= 0.9995 & dat_dm_for_forest$uci < 1), 0.999, dat_dm_for_forest$uci)
  }
  
  dat_dm_for_forest$text <- paste0(sprintf("%.3f", dat_dm_for_forest$hror)," (",sprintf("%.3f", dat_dm_for_forest$lci),", ",sprintf("%.3f", dat_dm_for_forest$uci),")") 
  
  dat_dm_for_forest$point_color <- ifelse(dat_dm_for_forest$study == "JD", "jd", "rep")
  dat_dm_for_forest$errorbar_color <- ifelse(dat_dm_for_forest$study == "JD", "jd", "rep")
  
  # dat_dm_for_forest$species <- ifelse(dat_dm_for_forest$species %in% gsub("_log10","",results_dm_same_direction_hrb_sig$species), paste0(dat_dm_for_forest$species,"*"), dat_dm_for_forest$species)
}
dat_dm_for_forest$species <- factor(dat_dm_for_forest$species, 
                                    levels = rev(unique(dat_dm_for_forest$species)), 
                                    labels = rev(trans_mp4_s_names(unique(dat_dm_for_forest$species))))
f1 <- forest_function(dat_dm_for_forest, "")

# CKD数据整理for plot #
{
  dat_ckd_for_forest <- results_ckd_long[results_ckd_long$species %in% results_ckd_same_direction_hrb_sig_fdr$species,]
  dat_ckd_for_forest$hror <- exp(dat_ckd_for_forest$estimate)
  dat_ckd_for_forest$lci <- exp(dat_ckd_for_forest$estimate - 1.96*dat_ckd_for_forest$se)
  dat_ckd_for_forest$uci <- exp(dat_ckd_for_forest$estimate + 1.96*dat_ckd_for_forest$se)
  
  ### 如果p<0.05需要确保UCI<1，因为四舍五入原因可能存在0.9995 变成 1.000，因此筛选出来赋值为0.999 ###
  {
    test <- dat_ckd_for_forest[dat_ckd_for_forest$p < 0.05 & (dat_ckd_for_forest$uci >= 0.9995 & dat_ckd_for_forest$uci < 1),]
    dat_ckd_for_forest$uci <- ifelse(dat_ckd_for_forest$p < 0.05 & (dat_ckd_for_forest$uci >= 0.9995 & dat_ckd_for_forest$uci < 1), 0.999, dat_ckd_for_forest$uci)
  }
  
  dat_ckd_for_forest$text <- paste0(sprintf("%.3f", dat_ckd_for_forest$hror)," (",sprintf("%.3f", dat_ckd_for_forest$lci),", ",sprintf("%.3f", dat_ckd_for_forest$uci),")") 
  
  dat_ckd_for_forest$point_color <- ifelse(dat_ckd_for_forest$study == "JD", "jd", "rep")
  dat_ckd_for_forest$errorbar_color <- ifelse(dat_ckd_for_forest$study == "JD", "jd", "rep")
  
  # dat_ckd_for_forest$species <- ifelse(dat_ckd_for_forest$species %in% gsub("_log10","",results_ckd_same_direction_hrb_sig$species), paste0(dat_ckd_for_forest$species,"*"), dat_ckd_for_forest$species)
}
dat_ckd_for_forest$species <- factor(dat_ckd_for_forest$species, 
                                     levels = rev(unique(dat_ckd_for_forest$species)),
                                     labels = rev(trans_mp4_s_names(unique(dat_ckd_for_forest$species))))
f2 <- forest_function(dat_ckd_for_forest, "")

# CVD数据整理for plot #
{
  dat_acvd_for_forest <- results_acvd_long[results_acvd_long$species %in% results_acvd_same_direction_acvd_sig_fdr$species,]
  dat_acvd_for_forest$hror <- exp(dat_acvd_for_forest$estimate)
  dat_acvd_for_forest$lci <- exp(dat_acvd_for_forest$estimate - 1.96*dat_acvd_for_forest$se)
  dat_acvd_for_forest$uci <- exp(dat_acvd_for_forest$estimate + 1.96*dat_acvd_for_forest$se)
  dat_acvd_for_forest$text <- paste0(sprintf("%.3f", dat_acvd_for_forest$hror)," (",sprintf("%.3f", dat_acvd_for_forest$lci),", ",sprintf("%.3f", dat_acvd_for_forest$uci),")") 
  
  dat_acvd_for_forest$point_color <- ifelse(dat_acvd_for_forest$study == "JD", "jd", "rep")
  dat_acvd_for_forest$errorbar_color <- ifelse(dat_acvd_for_forest$study == "JD", "jd", "rep")
  
  # dat_acvd_for_forest$species <- ifelse(dat_acvd_for_forest$species %in% gsub("_log10","",results_acvd_same_direction_acvd_sig$species), paste0(dat_acvd_for_forest$species,"*"), dat_acvd_for_forest$species)
}
dat_acvd_for_forest$species <- factor(dat_acvd_for_forest$species, 
                                    levels = rev(unique(dat_acvd_for_forest$species)), 
                                    labels = rev(trans_mp4_s_names(unique(dat_acvd_for_forest$species))))
f3 <- forest_function(dat_acvd_for_forest, "")


f1_2 <- cowplot::plot_grid(f1, NULL,
                           nrow = 2,
                           rel_heights = c(1, 0.25),
                           align = "v")
f3_2 <- cowplot::plot_grid(f3, NULL,
                           nrow = 2,
                           rel_heights = c(1, 2.6),
                           align = "v")

f <- cowplot::plot_grid(f1_2, f2, f3_2,
                        ncol = 3,
                        rel_widths = c(1, 1, 0.8))

ggsave(paste0("C:/TWang/DLiu/EDC_Micro/figures/main_figures/forest_plot_replication_dm-ckd-cvd_20260728.pdf"),
       f, width = 34, height = 20, limitsize = FALSE)
#### 作图Forest (JD HRB) ####


#### 补充信息 (不同疾病相关菌交集) ----
results_dm_temp <- results_dm %>% select(species)
results_ckd_temp <- results_ckd %>% select(species)
results_acvd_temp <- results_acvd %>% select(species)

dat_dm_ckd <- inner_join(results_dm_temp, results_ckd_temp, by="species")
dat_dm_ckd$significant <- "Diabetes and CKD"
dat_dm_acvd <- inner_join(results_dm_temp, results_acvd_temp, by="species")
dat_dm_acvd$significant <- "Diabetes and CVD"
dat_ckd_acvd <- inner_join(results_ckd_temp, results_acvd_temp, by="species")
dat_ckd_acvd$significant <- "CKD and CVD"

dat_dm_ckd_acvd <- inner_join(results_dm_temp, results_ckd_temp, by="species") %>%
  inner_join(results_acvd_temp, by="species")

dat_all <- rbind(dat_dm_ckd, dat_dm_acvd, dat_ckd_acvd) %>%
  arrange(significant, species)
dat_all <- dat_all[,c("species","significant")]
# openxlsx::write.xlsx(dat_all, "C:/TWang/DLiu/EDC_Micro/figures/main_figures/fig4c_supp_information_20260617.xlsx")
#### 补充信息 (不同疾病相关菌交集) ####


#### 添加菌的门信息 (作图数据处理) ----
{
  dat_dm_for_forest <- dat_dm_for_forest[dat_dm_for_forest$study == "JD",]
  dat_ckd_for_forest <- dat_ckd_for_forest[dat_ckd_for_forest$study == "JD",]
  dat_acvd_for_forest <- dat_acvd_for_forest[dat_acvd_for_forest$study == "JD",]
  
  #### 读取clade name数据 ----
  ## 读取clade name ##
  clade_name_mp4_s <- read.table("jiading/sourceDataTaxon/mpa4/JD.mp4.n4491_clade_name.txt", header = TRUE)
  clade_name_mp4_s$species <- trans_mp4_s_names(clade_name_mp4_s$species)
  clade_name_mp4_s <- clade_name_mp4_s[,c("phylum","species")]
  table(clade_name_mp4_s$species)
  table(clade_name_mp4_s$phylum)
  #### 读取clade name数据 ####
  
  dat_dm_for_phy_heatmap <- data.frame(species = unique(dat_dm_for_forest$species)) 
  dat_dm_for_phy_heatmap <- left_join(dat_dm_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_dm_for_phy_heatmap$y_axis <- 1
  
  dat_ckd_for_phy_heatmap <- data.frame(species = unique(dat_ckd_for_forest$species)) 
  dat_ckd_for_phy_heatmap <- left_join(dat_ckd_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_ckd_for_phy_heatmap$y_axis <- 1
  
  dat_acvd_for_phy_heatmap <- data.frame(species = unique(dat_acvd_for_forest$species)) 
  dat_acvd_for_phy_heatmap <- left_join(dat_acvd_for_phy_heatmap, clade_name_mp4_s, by=c("species"))
  dat_acvd_for_phy_heatmap$y_axis <- 1
  
  
  dat_dm_for_phy_heatmap$species <- factor(dat_dm_for_phy_heatmap$species,
                                           levels = c(dat_dm_for_phy_heatmap$species))
  dat_ckd_for_phy_heatmap$species <- factor(dat_ckd_for_phy_heatmap$species,
                                            levels = c(dat_ckd_for_phy_heatmap$species))
  dat_acvd_for_phy_heatmap$species <- factor(dat_acvd_for_phy_heatmap$species,
                                             levels = c(dat_acvd_for_phy_heatmap$species))
  
  phylum_colors <- c(
    "p__Actinobacteria"  = "#cf6a87",
    "p__Bacteroidetes"   = "#f19066",
    "p__Candidatus_Saccharibacteria" = "#f5cd79",
    "p__Firmicutes"      = "#2ecc71",
    "p__Proteobacteria"  = "#34ace0"
  )
}
# 添加菌的门信息 (作图)
f_heatmap_phylum_dm <- ggplot(dat_dm_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_dm, filename=paste0("figures/main_figures/forest_mp4_out_phylum_dm_20260728.pdf"), width = 14, height = 3.8, limitsize = FALSE)


f_heatmap_phylum_ckd1 <- ggplot(dat_ckd_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_ckd1, filename=paste0("figures/main_figures/forest_mp4_out_phylum_ckd_20260728.pdf"), width = 14, height = 3.8, limitsize = FALSE)


f_heatmap_phylum_cvd <- ggplot(dat_acvd_for_phy_heatmap, aes(x = species, y = y_axis)) +
  # 热图方块：fill 根据 phylum 自动上色
  geom_tile(aes(fill = phylum), 
            color="#EDEDED",  # 添加边框颜色
            linewidth = 0.05) +
  
  # 手动分配5种颜色
  scale_fill_manual(values = phylum_colors) +
  
  # 标签设置
  labs(
    x = "Species", 
    y = "",           # 隐藏Y轴标题（只有1行无需标题）
    fill = "Phylum"   # 图例标题
  ) +
  
  # 主题美化
  theme_minimal() + # 不要背景
  theme(
    plot.title = element_text(size = 25),  # 设置标题大小和加粗
    axis.title.x=element_blank(), # 去掉 x轴title
    axis.title.y=element_blank(), # 去掉 y轴title
    axis.ticks=element_blank(), # 去掉刻度线
    panel.grid =element_blank(),  # 删去网格线
    axis.text.x = element_text(angle = 300, hjust = 0, size = 12, color = "black"), # 调整x轴文字
    axis.text.y = element_text(size = 13, color = "black") #调整y轴文字
  )
ggsave(f_heatmap_phylum_cvd, filename=paste0("figures/main_figures/forest_mp4_out_phylum_cvd_20260728.pdf"), width = 14, height = 3.8, limitsize = FALSE)
#### 添加菌的门信息 (作图数据处理) ####
