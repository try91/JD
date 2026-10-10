library(data.table)
library(dplyr)
library(ggplot2)
library(corrplot)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

phenotype_dat <- read.table("raw_data/clinical_phenotypes_dat_20261006.txt", header = TRUE)
edc_dat <- read.table("raw_data/EDC_analytes_dat_20261006.txt", header = TRUE)
# 统计每人血液中detected的EDC和大于中位数的人数 #
{
  # 设置检测上下限
  limits <- list(
    low = c(MEHP = 1, MEHHP = 0.01, MEOHP = 0.02, MECPP = 0.1, MnBP = 0.5, MEP = 0.1, MBzP = 0.05, MCPP = 0.015, MiBP = 0.5, 
            TCS = 0.2, TCC = 0.02,
            BPA = 0.5, BPS = 0.1, BPF = 0.2, 
            PFOA = 1, PFNA = 1, PFDA = 0.5, PFOS = 1, PFHxS = 0.2),
    high = c(MEHP = 250, MEHHP = 2.5, MEOHP = 5, MECPP = 25, MnBP = 125, MEP = 25, MBzP = 25, MCPP = 30, MiBP = 125,
             TCS = 50, TCC = 5,
             BPA = 125, BPS = 25, BPF = 50, 
             PFOA = 250, PFNA = 250, PFDA = 125, PFOS = 250, PFHxS = 50)
  )
  
  # 将下限添加到数据框
  for (edc_name in names(limits$low)) {
    edc_dat[[paste0(edc_name, "_llimit")]] <- limits$low[[edc_name]]
  }
  
  # 通过自定义quantile_function函数，划分未检出和四分位数
  quantile_function <- function(dat, edc){
    if ((sum(dat[[paste0(edc, "_detected")]]) / nrow(dat)) > 0.75) {  # detection rate > 75%, 对所有进行分组
      # 四分位数转换
      q <- quantile(dat[[edc]], probs = c(0.25, 0.5, 0.75))
      # 定义Q1, Q2, Q3
      Q1 <- q[1]
      Q2 <- q[2]
      Q3 <- q[3]
      # 将每个数值分类到对应的四分位数
      dat[[paste0(edc,"_quantile")]] <- ifelse(dat[[edc]] <= Q1, 1, 
                                               ifelse(dat[[edc]] <= Q2, 2, 
                                                      ifelse(dat[[edc]] <= Q3, 3, 4)))
    }else{  # detection rate < 75%, 只对检出数据进行分组
      # 筛选出大于最低限度的数值，对他们进行四分位分类
      low_limit <- dat[[paste0(edc, "_llimit")]][1]
      filtered_values <- dat[[edc]][dat[[edc]] > low_limit]
      q <- quantile(filtered_values, probs = c(0.25, 0.5, 0.75))
      # 定义Q1, Q2, Q3
      Q1 <- q[1]
      Q2 <- q[2]
      Q3 <- q[3]
      # 将每个数值分类到对应的四分位数
      dat[[paste0(edc,"_quantile")]] <- ifelse(dat[[edc]] <= low_limit, 0, 
                                               ifelse(dat[[edc]] <= Q1, 1, 
                                                      ifelse(dat[[edc]] <= Q2, 2, 
                                                             ifelse(dat[[edc]] <= Q3, 3, 4))))
    }
    
    return(dat)
  }
  for (i in edc_traits) {
    edc_dat <- quantile_function(edc_dat,i)
  }
  
  # 通过自定义median_function函数，划分未检出和中位数
  median_function <- function(dat, edc){
    if ((sum(dat[[paste0(edc, "_detected")]]) / nrow(dat)) > 0.5) {  # detection rate > 50%, 对所有进行分组
      # 取中位数
      Median <- median(dat[[edc]])
      # 将每个数值分类到对应的中位数
      dat[[paste0(edc,"_median")]] <- ifelse(dat[[edc]] <= Median, 1, 2)
    }else{  # detection rate < 50%, 只对检出数据进行分组
      # 筛选出大于最低限度的数值，对他们取中位数
      low_limit <- dat[[paste0(edc, "_llimit")]][1]
      filtered_values <- dat[[edc]][dat[[edc]] > low_limit]
      Median <- median(filtered_values)
      # 将每个数值分类到对应的中位数
      dat[[paste0(edc,"_median")]] <- ifelse(dat[[edc]] <= low_limit, 0, 
                                             ifelse(dat[[edc]] <= Median, 1, 2))
    }
    
    return(dat)
  }
  for (i in edc_traits) {
    edc_dat <- median_function(edc_dat,i)
  }
  
  
  # 统计每人血液中detected的EDC大于中位数和位于Q4的人数
  for (i in 1:19) {
    # i <- 1
    colnames(edc_dat)
    EDC <- c(edc_traits[i],paste0(edc_traits[i],"_detected"))
    temp <- edc_dat[, EDC]
    temp2 <- temp[temp[[paste0(edc_traits[i],"_detected")]] == 1,]
    
    q <- quantile(temp2[[1]], probs = c(0.25, 0.5, 0.75))
    # 定义Q1, Q2, Q3
    Q1 <- q[1]
    Q2 <- q[2]
    Q3 <- q[3]
    
    edc_dat[[paste0(edc_traits[i],"_Q2")]] <- Q2
    edc_dat[[paste0(edc_traits[i],"_Q2_flag")]] <- ifelse(edc_dat[[paste0(edc_traits[i])]] > edc_dat[[paste0(edc_traits[i],"_Q2")]], 1, 0)
    edc_dat[[paste0(edc_traits[i],"_Q3")]] <- Q3
    edc_dat[[paste0(edc_traits[i],"_Q3_flag")]] <- ifelse(edc_dat[[paste0(edc_traits[i])]] > edc_dat[[paste0(edc_traits[i],"_Q3")]], 1, 0)
  }
  
  # 血液中可检测的EDC数量
  columns_containing_de <- grep("_detected", names(edc_dat), value = TRUE) # 19 EDC
  edc_dat$num_edc1 <- rowSums(edc_dat[, columns_containing_de], na.rm = TRUE) # 19 EDC
  
  # 血液污染物浓度>中位数的EDC数量
  columns_containing_q2 <- grep("_Q2_flag", names(edc_dat), value = TRUE) # 19 EDC
  edc_dat$num_edc2 <- rowSums(edc_dat[, columns_containing_q2], na.rm = TRUE) # 19 EDC
  
  # 血液污染物浓度位于Q4的EDC数量
  columns_containing_q3 <- grep("_Q3_flag", names(edc_dat), value = TRUE) # 19 EDC
  edc_dat$num_edc3 <- rowSums(edc_dat[, columns_containing_q3], na.rm = TRUE) # 19 EDC
}
micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

phy_edc_dat <- left_join(phenotype_dat, edc_dat, by = "ID") %>%
  left_join(micro_dat, by = "ID")
sample_name <- "phy_edc_temp"  # 提取subgroup的名称

#### 配色 ----
PFAS_colors <- colorRampPalette(c("#f55d78","#FFFFFF"))(50)[c(1,8,15,22,29)]
PAE_colors <- colorRampPalette(c("#3498db","#FFFFFF"))(50)[c(1,6,11,16,21,26,31,36,41)]
TC_colors <- colorRampPalette(c("#00b894","#FFFFFF"))(50)[c(1,21)]
BP_colors <- colorRampPalette(c("#f1c40f","#FFFFFF"))(50)[c(1,18,35)]
my_colors <- c(PFAS_colors,PAE_colors,TC_colors,BP_colors)
#### 配色 ####

#### 变量整理 ----
# 菌群2014菌群 (分类和连续)
# 丰度>0.0001, 出现率>10%的微生物 (物种和属)
mp4_s_names <- colnames(micro_dat)[3:361]
mp4_g_names <- colnames(micro_dat)[721:912]
# 排除未分类的菌属（GGB）和菌种（SGB） #
mp4_s_names_short <- mp4_s_names[!grepl("_GGB",mp4_s_names)] # 排除未分类的菌属（GGB）, 未分类菌种（SGB）先保留
mp4_g_names_short <- mp4_g_names[!grepl("_GGB",mp4_g_names)] # 排除未分类的菌属（GGB）
# 排除未分类的菌属（GGB）和菌种（SGB） #

# 转换后的菌的名称
mp4_s_log10 <- paste0(mp4_s_names,"_log10") # 菌群MP4丰度的log10转换 (物种层面)
mp4_s_log10_short <- paste0(mp4_s_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 物种层面)

mp4_g_log10 <- paste0(mp4_g_names,"_log10") # 菌群MP4丰度的log10转换 (属层面)
mp4_g_log10_short <- paste0(mp4_g_names_short,"_log10") # 菌群MP4丰度的log10转换 (有鉴定菌属, 属层面)


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


#### 数据处理 ----
# 分别挑选EDC浓度和是否被检出两个数据框
cols_edc_traits <- c(edc_traits,paste0(edc_traits,"_log10"),paste0(edc_traits,"_detected"),"sex_b_rev","age_b","smk1_b","drk1_b","high_edu_b","paactive3_g_b","high_fruveg")
phy_edc_temp <- phy_edc_dat[,cols_edc_traits]

cols_edc_detected <- c(paste0(edc_traits,"_detected"),"sex_b_rev")
phy_edc_temp_detected <- phy_edc_dat[,cols_edc_detected]

for (i in edc_traits){
  # 根据EDC_detected列中的数值来删除对应EDC和EDC_log10中的数值。具体来说，如果EDC_detected中的单元格是1，则对应EDC的数值保持不变；如果是0，则将对应EDC的数值改为NA。
  phy_edc_temp[[i]] <- ifelse(phy_edc_temp[[paste0(i,"_detected")]] == 0, NA, phy_edc_temp[[i]])
  phy_edc_temp[[paste0(i,"_log10")]] <- ifelse(phy_edc_temp[[paste0(i,"_detected")]] == 0, NA, phy_edc_temp[[paste0(i,"_log10")]])
}

# 统计第1到第19列的非NA值个数 (detected EDC数量)
phy_edc_temp$num_edc1 <- rowSums(!is.na(phy_edc_temp[, 1:19]))
table(phy_edc_temp$num_edc1)
#### 数据处理 ####


############################################# EDC人群分布特征图 #############################################
#### 小提琴图 ----
### 创建一个函数，用于画EDC小提琴图
violin_plot <- function(DAT,X_TEXT,COLOR){
  plot_edc <- ggplot(DAT, aes(x=edc, y=concentration, fill=edc)) + 
    geom_violin(trim=FALSE, size = 0.4) +
    geom_boxplot(
      width=0.12,                # 设置箱体宽度
      outlier.shape = 19,        # 实心圆点
      outlier.size = 0.6,        # 设置异常点大小
      fatten = 1.2,              # 中位线相对粗细比例
      fill="white") +
    
    geom_text(  # 添加检出率标签
      data = X_TEXT,
      aes(x = edc, y = MAX_CONCENTRATION, label = paste0("DR=",DR,"%")),
      hjust = -0.1,
      # vjust = -0.7,               # 标签上移避免重叠
      size = 4,                     # 字体大小
      inherit.aes = FALSE           # 忽略主图映射的 fill
    ) +
    
    scale_y_continuous(
      labels = scales::math_format(10^.x)  # 标签显示为 10^y  # 修正科学记数法标签
    ) +
    
    scale_fill_manual(values = COLOR) +  # 手动指定颜色
    
    # labs(title="Plot of concentration by EDC",x="EDC", y = "Concentration") +
    labs(x="EDC", y = "Serum concentration (ng/ml)") +
    
    theme_classic() +
    theme(
      plot.margin = margin(2, 2, 2, 2, "mm"), # t,r,b,l
      
      axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      
      axis.title.x = element_blank(), 
      axis.title.y = element_text(size = 15, colour = "black"),
      
      axis.text.x = element_text(size = 15, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 15, colour = "black"),
      
      # legend.title = element_blank(),
      # legend.text = element_text(size = 8, colour = "black"), # 调整legend文本大小
      legend.position = "none"
    )
  
  
  return(plot_edc)
}

## 总人群 ##
phy_edc_temp_long <- tidyr::gather(phy_edc_temp, edc, concentration, PFOS_log10:BPF_log10, na.rm=TRUE, factor_key=TRUE)
phy_edc_temp_long$edc <- gsub("_log10","",phy_edc_temp_long$edc)
phy_edc_temp_long$edc <- factor(phy_edc_temp_long$edc,levels = edc_traits)

# 计算每一种EDC的detected rate、max concentration为label定位
phy_edc_temp_long_distinct <- phy_edc_temp_long %>%
  group_by(edc) %>%
  mutate(
    N = n(),
    DR = round(n()/nrow(phy_edc_temp),4)*100,
    MAX_CONCENTRATION = max(concentration)) %>%
  ungroup() %>% 
  distinct(edc,N,DR,MAX_CONCENTRATION)

pfas_for_plot <- phy_edc_temp_long[phy_edc_temp_long$edc %in% edc_traits2,c("edc","concentration")]
pae_for_plot <- phy_edc_temp_long[phy_edc_temp_long$edc %in% edc_traits3,c("edc","concentration")]
bp_for_plot <- phy_edc_temp_long[phy_edc_temp_long$edc %in% edc_traits4,c("edc","concentration")]
tc_for_plot <- phy_edc_temp_long[phy_edc_temp_long$edc %in% edc_traits5,c("edc","concentration")]

pfas_for_plot_label <- phy_edc_temp_long_distinct[phy_edc_temp_long_distinct$edc %in% edc_traits2,]
pae_for_plot_label <- phy_edc_temp_long_distinct[phy_edc_temp_long_distinct$edc %in% edc_traits3,]
bp_for_plot_label <- phy_edc_temp_long_distinct[phy_edc_temp_long_distinct$edc %in% edc_traits4,]
tc_for_plot_label <- phy_edc_temp_long_distinct[phy_edc_temp_long_distinct$edc %in% edc_traits5,]
## 总人群 ##

f1 <- violin_plot(pfas_for_plot,pfas_for_plot_label,PFAS_colors)
f2 <- violin_plot(pae_for_plot,pae_for_plot_label,PAE_colors)
f3 <- violin_plot(tc_for_plot,tc_for_plot_label,TC_colors)
f4 <- violin_plot(bp_for_plot,bp_for_plot_label,BP_colors)

f_edc1 <- cowplot::plot_grid(f1, f2, ncol = 2, rel_widths = c(1, 1.8))
f_edc2 <- cowplot::plot_grid(f3, f4, ncol = 3, rel_widths = c(1.6, 2, 4.8))
f_edc3 <- cowplot::plot_grid(f_edc1, f_edc2, nrow = 2, rel_heights = c(1, 1))
ggsave(f_edc3, filename=paste0("figures/main_figures/(fig2a)_violin_plot_4_edc_group.pdf"), width = 17, height = 8, limitsize = FALSE)
#### 小提琴图 ####


#### 堆叠柱状图1 (统计每人血液中detected EDC种类。同时显示在每一个数量的人中，每一种EDC的比例) ----
# EDC计数堆叠柱状图标准：检出率>50%，直接根据中位数分为2组（如果检出率<50%，则在检出的样本中根据中位数分为2组），检出率>75%，直接根据四分位数分为4组（如果检出率<75%，则在检出的样本根据四分位数分为4组） #
### 创建一个函数，用于画堆叠柱状图
stack_plot <- function(DAT,X_TEXT,X_TITLE){
  f <- ggplot(DAT, aes(
    x = factor(num_edc1, levels = 4:18),  
    y = total_num, fill=edc))+
    # geom_col(position = 'stack', color="black", size=0.1) +  # 添加黑色边框
    geom_col(position = 'stack') +
    geom_text(  # 添加总数标签
      data = X_TEXT,
      aes(x = x_position, y = n, label = n),
      vjust = -0.7,                 # 标签上移避免重叠
      size = 5,                     # 字体大小
      inherit.aes = FALSE           # 忽略主图映射的 fill
    ) +
    scale_x_discrete(breaks = 4:18) +      # 确保显示所有刻度
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +      # Y 轴从下部不扩展，而上部扩展5%
    labs(x = X_TITLE, y = "No. of participants") +
    scale_fill_manual(values = my_colors) +  # 手动指定颜色
    theme_classic() +
    theme(
      axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      axis.title.x = element_text(size = 20, colour = "black"),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.x = element_text(size = 20, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 20, colour = "black"),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 10, colour = "black"), # 调整legend文本大小
      # 图例键大小
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(4, "mm")
    )
  
  return(f)
}

## 查找问题是否有遗漏的人，因为没有检测出EDC ##
phy_edc_temp_detected$num_edc1 <- rowSums(phy_edc_temp_detected[,1:19] != 0)
table(phy_edc_temp_detected$num_edc1)

## 总人群 ##
dat_for_plot <- phy_edc_temp_detected %>%
  group_by(num_edc1) %>%
  mutate(n = n(),
         total_num_edc = sum(num_edc1),
         total_num_pfas1 = sum(PFOS_detected),
         total_num_pfas2 = sum(PFOA_detected),
         total_num_pfas3 = sum(PFNA_detected),
         total_num_pfas4 = sum(PFDA_detected),
         total_num_pfas5 = sum(PFHxS_detected),
         
         total_num_pae1 = sum(MEHP_detected),
         total_num_pae2 = sum(MECPP_detected),
         total_num_pae3 = sum(MEHHP_detected),
         total_num_pae4 = sum(MEP_detected),
         total_num_pae5 = sum(MEOHP_detected),
         total_num_pae6 = sum(MiBP_detected),
         total_num_pae7 = sum(MnBP_detected),
         total_num_pae8 = sum(MCPP_detected),
         total_num_pae9 = sum(MBzP_detected),
         
         total_num_tc1 = sum(TCC_detected),
         total_num_tc2 = sum(TCS_detected),
         
         total_num_bp1 = sum(BPA_detected),
         total_num_bp2 = sum(BPS_detected),
         total_num_bp3 = sum(BPF_detected),
         
         
         participants_pfas1_percent = total_num_pfas1/total_num_edc *n,
         participants_pfas2_percent = total_num_pfas2/total_num_edc *n,
         participants_pfas3_percent = total_num_pfas3/total_num_edc *n,
         participants_pfas4_percent = total_num_pfas4/total_num_edc *n,
         participants_pfas5_percent = total_num_pfas5/total_num_edc *n,
         
         participants_pae1_percent = total_num_pae1/total_num_edc *n,
         participants_pae2_percent = total_num_pae2/total_num_edc *n,
         participants_pae3_percent = total_num_pae3/total_num_edc *n,
         participants_pae4_percent = total_num_pae4/total_num_edc *n,
         participants_pae5_percent = total_num_pae5/total_num_edc *n,
         participants_pae6_percent = total_num_pae6/total_num_edc *n,
         participants_pae7_percent = total_num_pae7/total_num_edc *n,
         participants_pae8_percent = total_num_pae8/total_num_edc *n,
         participants_pae9_percent = total_num_pae9/total_num_edc *n,
         
         participants_tc1_percent = total_num_tc1/total_num_edc *n,
         participants_tc2_percent = total_num_tc2/total_num_edc *n,
         
         participants_bp1_percent = total_num_bp1/total_num_edc *n,
         participants_bp2_percent = total_num_bp2/total_num_edc *n,
         participants_bp3_percent = total_num_bp3/total_num_edc *n) %>%
  distinct(num_edc1,.keep_all = TRUE) %>%
  ungroup

colnames(dat_for_plot)[43:61] <- edc_traits

dat_for_plot_long <- tidyr::gather(dat_for_plot, edc, total_num, PFOS:BPF, na.rm=TRUE, factor_key=TRUE)
# 统计每个edc数量组的总数，标注在图中
total_labels <- dat_for_plot_long %>% distinct(num_edc1,n)
total_labels$x_position <-  total_labels$num_edc1 - 3 # 由于原本图x轴从4开始，因此num_edc1=4其实是1，通过-3进行调整对齐
## 总人群 ##

f <- stack_plot(dat_for_plot_long,total_labels,"No. of detected EDCs")
ggsave(f, filename=paste0("figures/main_figures/(fig2b)_stack_num_edc1.pdf"), width = 8, height = 6, limitsize = FALSE)
#### 堆叠柱状图1 (统计每人血液中detected EDC种类。同时显示在每一个数量的人中，每一种EDC的比例) ####

#### 堆叠柱状图2 (统计每人血液中浓度>median EDC种类。同时显示在每一个数量的人中，每一种EDC的比例) ----
# EDC计数堆叠柱状图标准：检出率>50%，直接根据中位数分为2组（如果检出率<50%，则在检出的样本中根据中位数分为2组），检出率>75%，直接根据四分位数分为4组（如果检出率<75%，则在检出的样本根据四分位数分为4组） #
### 创建一个函数，用于画堆叠柱状图
stack_plot2 <- function(DAT,X_TEXT,X_TITLE){
  f <- ggplot(DAT, aes(
    x = factor(num_edc2),  # 将数值转为因子并强制定义所有分类
    y = total_num, fill=edc))+
    geom_col(position = 'stack') +
    geom_text(  # 添加总数标签
      data = X_TEXT,
      aes(x = x_position, y = n, label = n),
      vjust = -0.7,                 # 标签上移避免重叠
      size = 5,                     # 字体大小
      inherit.aes = FALSE           # 忽略主图映射的 fill
    ) +
    scale_x_discrete(breaks = 0:16) +  # 确保显示所有刻度
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +      # Y 轴从下部不扩展，而上部扩展5%
    labs(x = X_TITLE, y = "No. of participants") +
    scale_fill_manual(values = my_colors) +  # 手动指定颜色
    theme_classic() +
    theme(
      axis.ticks = element_line(size = 0.5, colour = "black"),# 调整 X 轴刻度线, 线宽(默认0.5)
      axis.ticks.length = unit(0.1, "cm"),  # 可选：调整刻度线长度(需配合 axis.ticks.length), 长度(默认0.2cm)
      axis.title.x = element_text(size = 20, colour = "black"),
      axis.title.y = element_text(size = 20, colour = "black"),
      axis.text.x = element_text(size = 20, colour = "black"), # 调整x轴文字
      axis.text.y = element_text(size = 20, colour = "black"),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 10, colour = "black"), # 调整legend文本大小
      # 图例键大小
      legend.key.height = unit(4, "mm"),
      legend.key.width = unit(4, "mm")
    )
  
  return(f)
}


cols_edc_median <- c(paste0(edc_traits,"_Q2_flag"),"sex_b_rev")
phy_edc_temp_median <- phy_edc_dat[,cols_edc_median]
phy_edc_temp_median$num_edc2 <- rowSums(phy_edc_temp_median[,1:19] != 0)
table(phy_edc_temp_median$num_edc2)

## 总人群 ##
dat_for_plot <- phy_edc_temp_median %>%
  group_by(num_edc2) %>%
  mutate(n = n(),
         total_num_edc = sum(num_edc2),
         total_num_pfas1 = sum(PFOS_Q2_flag),
         total_num_pfas2 = sum(PFOA_Q2_flag),
         total_num_pfas3 = sum(PFNA_Q2_flag),
         total_num_pfas4 = sum(PFDA_Q2_flag),
         total_num_pfas5 = sum(PFHxS_Q2_flag),
         
         total_num_pae1 = sum(MEHP_Q2_flag),
         total_num_pae2 = sum(MECPP_Q2_flag),
         total_num_pae3 = sum(MEHHP_Q2_flag),
         total_num_pae4 = sum(MEP_Q2_flag),
         total_num_pae5 = sum(MEOHP_Q2_flag),
         total_num_pae6 = sum(MiBP_Q2_flag),
         total_num_pae7 = sum(MnBP_Q2_flag),
         total_num_pae8 = sum(MCPP_Q2_flag),
         total_num_pae9 = sum(MBzP_Q2_flag),
         
         total_num_tc1 = sum(TCC_Q2_flag),
         total_num_tc2 = sum(TCS_Q2_flag),
         
         total_num_bp1 = sum(BPA_Q2_flag),
         total_num_bp2 = sum(BPS_Q2_flag),
         total_num_bp3 = sum(BPF_Q2_flag),
         
         
         participants_pfas1_percent = total_num_pfas1/total_num_edc *n,
         participants_pfas2_percent = total_num_pfas2/total_num_edc *n,
         participants_pfas3_percent = total_num_pfas3/total_num_edc *n,
         participants_pfas4_percent = total_num_pfas4/total_num_edc *n,
         participants_pfas5_percent = total_num_pfas5/total_num_edc *n,
         
         participants_pae1_percent = total_num_pae1/total_num_edc *n,
         participants_pae2_percent = total_num_pae2/total_num_edc *n,
         participants_pae3_percent = total_num_pae3/total_num_edc *n,
         participants_pae4_percent = total_num_pae4/total_num_edc *n,
         participants_pae5_percent = total_num_pae5/total_num_edc *n,
         participants_pae6_percent = total_num_pae6/total_num_edc *n,
         participants_pae7_percent = total_num_pae7/total_num_edc *n,
         participants_pae8_percent = total_num_pae8/total_num_edc *n,
         participants_pae9_percent = total_num_pae9/total_num_edc *n,
         
         participants_tc1_percent = total_num_tc1/total_num_edc *n,
         participants_tc2_percent = total_num_tc2/total_num_edc *n,
         
         participants_bp1_percent = total_num_bp1/total_num_edc *n,
         participants_bp2_percent = total_num_bp2/total_num_edc *n,
         participants_bp3_percent = total_num_bp3/total_num_edc *n) %>%
  distinct(num_edc2,.keep_all = TRUE)

colnames(dat_for_plot)[43:61] <- edc_traits

dat_for_plot_long <- tidyr::gather(dat_for_plot, edc, total_num, PFOS:BPF, na.rm=FALSE, factor_key=TRUE)
dat_for_plot_long$total_num <- ifelse(dat_for_plot_long$num_edc2 == 0 & dat_for_plot_long$edc == "PFOS", 91, 
                                      ifelse(dat_for_plot_long$num_edc2 == 0, 0, dat_for_plot_long$total_num))
# 统计每个edc数量组的总数，标注在图中
total_labels <- dat_for_plot_long %>% distinct(num_edc2,n)
total_labels$x_position <-  total_labels$num_edc2 + 1 # 由于原本图x轴从0开始，因此num_edc1=0其实是1，通过+1进行调整对齐
## 总人群 ##

f <- stack_plot2(dat_for_plot_long,total_labels,"No. of EDCs with a concentration above the median level")
ggsave(f, filename=paste0("figures/main_figures/(fig2b)_stack_num_edc2.pdf"), width = 8, height = 6, limitsize = FALSE)
#### 堆叠柱状图2 (统计每人血液中浓度>median EDC种类。同时显示在每一个数量的人中，每一种EDC的比例) ####
