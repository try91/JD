library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

#### 配色 ----
my_palette <- colorRampPalette(colors = c("#cf6a87", "#f19066", "#2ecc71", "#34ace0"))(4)
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


#### 处理菌群数据 ----
## 读取clade name ##
clade_name_mp4_s <- read.table("raw_data/JD.mp4.n4491_clade_name.txt", header = TRUE)
clade_name_mp4_g <- clade_name_mp4_s %>%
  distinct(kingdom, phylum, class, order, family, genus)

clade_name_mp4_s <- clade_name_mp4_s[clade_name_mp4_s$species %in% mp4_s_names,]
row.names(clade_name_mp4_s) <- NULL
clade_name_mp4_s$var_sort <- as.numeric(row.names(clade_name_mp4_s))

clade_name_mp4_g <- clade_name_mp4_g[clade_name_mp4_g$genus %in% mp4_g_names,]
row.names(clade_name_mp4_g) <- NULL
clade_name_mp4_g$var_sort <- as.numeric(row.names(clade_name_mp4_g))


## 读取所有赛选后门、纲、目、科信息 (丰度>0.0001检出率>10%) ##
mp4_p_names <- read.table("raw_data/phylum_names_mp4_10%.txt")
mp4_p_names <- mp4_p_names[,1]
mp4_c_names <- read.table("raw_data/class_names_mp4_10%.txt")
mp4_c_names <- mp4_c_names[,1]
mp4_o_names <- read.table("raw_data/order_names_mp4_10%.txt")
mp4_o_names <- mp4_o_names[,1]
mp4_f_names <- read.table("raw_data/family_names_mp4_10%.txt")
mp4_f_names <- mp4_f_names[,1]
#### 处理菌群数据 ####

#### 处理菌群-outcome数据 ----
## 读取MP4-Outcomes COX & logistic分析结果 ##
cox_results <- readxl::read_xlsx("results/cox/cox_results_mp4_incident.xlsx")
cox_results <- cox_results[cox_results$adjust == "adj" & cox_results$outcome == "cvd_incident_1421",]
cox_results <- cox_results[,c(1,3:13)]
colnames(cox_results)[c(1:4)] <- c("estimate","se","z","p")
out_cox <- unique(cox_results$outcome) # 提取结局变量
exp_cox <- unique(cox_results$exposure) # 提取暴露变量

logistic_results <- readxl::read_xlsx("results/glm/logistic_results_mp4_incident.xlsx")
logistic_results <- logistic_results[logistic_results$adjust == "adj" & logistic_results$outcome %in% c("ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"),]
logistic_results <- logistic_results[logistic_results$adjust == "adj",]
colnames(logistic_results)[c(1:4)] <- c("estimate","se","z","p")
out_logistic <- unique(logistic_results$outcome) # 提取结局变量
exp_logistic <- unique(logistic_results$exposure) # 提取暴露变量

cox_logistic_results <- rbind(cox_results,logistic_results)
cox_logistic_results <- cox_logistic_results[cox_logistic_results$rowname %in% c(mp4_s_log10),]
# 以每个OUT表型为单位进行校正（以outcome为组，校正每个菌）
cox_logistic_results <- cox_logistic_results %>%
  group_by(outcome) %>%  # 按outcome分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
cox_logistic_short <- cox_logistic_results[,c("rowname","outcome","estimate","se","p","p_adj_bh")]
colnames(cox_logistic_short) <- c("exp_name","out_name","estimate_cox_logistic","se_cox_logistic","p","p_adj_bh")

cox_logistic_short$out_name <- factor(cox_logistic_short$out_name, levels = c("cvd_incident_1421","ckd_incident_1014_no_self_report","dm_incident_1014_no_self_report"))
cox_logistic_short <- cox_logistic_short %>%
  arrange(exp_name,out_name)

# effect size 转换为 HR/OR
cox_logistic_short$exp_estimate <- exp(cox_logistic_short$estimate_cox_logistic)

## 分别提取三个结局结果 ##
cox_logistic_short_cvd <- cox_logistic_short[cox_logistic_short$out_name == "cvd_incident_1421",]
cox_logistic_short_ckd <- cox_logistic_short[cox_logistic_short$out_name == "ckd_incident_1014_no_self_report",]
cox_logistic_short_dm <- cox_logistic_short[cox_logistic_short$out_name == "dm_incident_1014_no_self_report",]

gm_sig_cvd <-unique(cox_logistic_short_cvd[cox_logistic_short_cvd$p < 0.05,]) # 26个CVD显著相关的菌种(新结果)
gm_sig_ckd <-unique(cox_logistic_short_ckd[cox_logistic_short_ckd$p < 0.05,]) # 84个CKD显著相关的菌种(新结果)
gm_sig_dm <-unique(cox_logistic_short_dm[cox_logistic_short_dm$p < 0.05,]) # 35个DM显著相关的菌种(新结果)
#### 处理菌群-outcome数据 ####


############################################# 曼哈顿图 #############################################
#### 筛选和结局相关的菌 ----
# 每组数量如果小于等于100则全保留，如果大于100则保留前100以及剩余随机保留10%数值
dat_3_outcome_gm <- cox_logistic_short[cox_logistic_short$out_name %in% c("dm_incident_1014_no_self_report","ckd_incident_1014_no_self_report","cvd_incident_1421"),]
dat_3_outcome_gm$OUTCOME <- ifelse(dat_3_outcome_gm$out_name %in% c("cvd_incident_1421"), "CVD",
                                   ifelse(dat_3_outcome_gm$out_name %in% c("ckd_incident_1014_no_self_report"), "CKD", "DM"))
dat_3_outcome_gm$exp_name <- gsub("_log10","",dat_3_outcome_gm$exp_name)
dat_3_outcome_gm$p_trans <- -log10(dat_3_outcome_gm$p)

# 匹配上每个species的排列序号，方便后期作图
clade_name_mp4_s_sort <- clade_name_mp4_s[,c("species","var_sort")]
dat_sort <- left_join(dat_3_outcome_gm, clade_name_mp4_s_sort, by=c("exp_name" = "species"))

dat_sort <- dat_sort %>%
  arrange(var_sort)
dat_sort <- dat_sort[dat_sort$p < 0.05,]
dat_sort_uniq <- dat_sort %>%
  distinct(exp_name,var_sort)
dat_sort_uniq$sort_new <- as.numeric(row.names(dat_sort_uniq))

# 匹配上每个species的排列序号，方便后期作图
dat_sort <- inner_join(dat_3_outcome_gm, dat_sort_uniq, by="exp_name")
#### 筛选和结局相关的菌 ####


#### Stack plot (x轴sort 1:127，y轴 phylum:species) ----
clade_name_mp4_s_for_plot <- inner_join(clade_name_mp4_s, dat_sort_uniq, by=c("species" = "exp_name"))
clade_name_mp4_s_for_plot <- clade_name_mp4_s_for_plot[,c(1:7,10)]
colnames(clade_name_mp4_s_for_plot)[8] <- "var_sort"

## 每一个clade层次设置颜色 ##
{
  # 构建函数
  color_palette_function <- function(DAT){
    
    # 为每一个二级分类设置颜色 #
    dat_temp <- DAT
    
    color_temp <- vector()
    for (i in c(1:nrow(dat_temp))) {
      
      if(dat_temp$n[i] == 1){
        my_palette_temp <- dat_temp$clade_name_color[i]
      }else{
        if(i == nrow(dat_temp)){
          my_palette_temp <- colorRampPalette(colors = c(dat_temp$clade_name_color[i], "#a29bfe"))(dat_temp$n[i]+1)
          my_palette_temp <- my_palette_temp[-c(length(my_palette_temp))]
        }else{
          my_palette_temp <- colorRampPalette(colors = c(dat_temp$clade_name_color[i], dat_temp$clade_name_color[i+1]))(dat_temp$n[i]+1) # 选取上层颜色和下一个上层颜色产生渐变色
          my_palette_temp <- my_palette_temp[-c(length(my_palette_temp))]  # 删除下一个上层颜色，包颜色不会重叠
        }
      }
      
      color_temp <- c(color_temp,my_palette_temp)
    }
    
    return(color_temp)
  }
  
  
  # 为每一个phylum设置颜色 #
  unique_p <- clade_name_mp4_s_for_plot %>%  # 1:5
    distinct(class,.keep_all = TRUE) %>%
    group_by(phylum) %>%
    mutate(n = n()) %>%
    distinct(phylum,n) %>%
    ungroup()
  unique_p$color <- my_palette
  clade_name_color_p <- unique_p[,c("phylum","color","n")]
  colnames(clade_name_color_p) <- c("clade_name","clade_name_color","n")
  
  
  # 为每一个class设置颜色 #
  unique_c <- clade_name_mp4_s_for_plot %>%  # 1:23
    distinct(order,.keep_all = TRUE) %>%
    group_by(class) %>%
    mutate(n = n()) %>%
    distinct(class,n) %>%
    ungroup()
  color_all <- color_palette_function(clade_name_color_p)
  # # 查看颜色梯度
  # scales::show_col(color_all)
  unique_c$color <- color_all
  clade_name_color_c <- unique_c[,c("class","color","n")]
  colnames(clade_name_color_c) <- c("clade_name","clade_name_color","n")
  
  
  # 为每一个order设置颜色 #
  unique_o <- clade_name_mp4_s_for_plot %>%  # 1:29
    distinct(family,.keep_all = TRUE) %>%
    group_by(order) %>%
    mutate(n = n()) %>%
    distinct(order,n) %>%
    ungroup()
  color_all <- color_palette_function(clade_name_color_c)
  unique_o$color <- color_all
  clade_name_color_o <- unique_o[,c("order","color","n")]
  colnames(clade_name_color_o) <- c("clade_name","clade_name_color","n")
  
  
  # 为每一个family设置颜色 #
  unique_f <- clade_name_mp4_s_for_plot %>%  # 1:40
    distinct(genus,.keep_all = TRUE) %>%
    group_by(family) %>%
    mutate(n = n()) %>%
    distinct(family,n) %>%
    ungroup()
  color_all <- color_palette_function(clade_name_color_o)
  unique_f$color <- color_all
  clade_name_color_f <- unique_f[,c("family","color","n")]
  colnames(clade_name_color_f) <- c("clade_name","clade_name_color","n")
  
  
  # 为每一个genus设置颜色 #
  unique_g <- clade_name_mp4_s_for_plot %>%  # 1:87
    distinct(species,.keep_all = TRUE) %>%
    group_by(genus) %>%
    mutate(n = n()) %>%
    distinct(genus,n) %>%
    ungroup()
  color_all <- color_palette_function(clade_name_color_f)
  unique_g$color <- color_all
  clade_name_color_g <- unique_g[,c("genus","color","n")]
  colnames(clade_name_color_g) <- c("clade_name","clade_name_color","n")
  
  
  # 为每一个species设置颜色 #
  unique_s <- clade_name_mp4_s_for_plot %>%  # 1:127
    group_by(species) %>%
    mutate(n = n()) %>%
    distinct(species,n) %>%
    ungroup()
  color_all <- color_palette_function(clade_name_color_g)
  unique_s$color <- color_all
  clade_name_color_s <- unique_s[,c("species","color","n")]
  colnames(clade_name_color_s) <- c("clade_name","clade_name_color","n")
  
  clade_name_color_matche_all <- rbind(clade_name_color_p,clade_name_color_c,clade_name_color_o,
                                       clade_name_color_f,clade_name_color_g,clade_name_color_s)
}


# 整合每个层级clade name中species的数量
{
  clade_name_for_stack_k <- clade_name_mp4_s_for_plot%>%
    group_by(kingdom) %>%
    mutate(n=n()) %>%
    distinct(kingdom,n,.keep_all = TRUE) %>%
    dplyr::select(kingdom,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_p <- clade_name_mp4_s_for_plot %>%
    group_by(phylum) %>%
    mutate(n=n()) %>%
    distinct(phylum,n,.keep_all = TRUE) %>%
    dplyr::select(phylum,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_c <- clade_name_mp4_s_for_plot %>%
    group_by(class) %>%
    mutate(n=n()) %>%
    distinct(class,n,.keep_all = TRUE) %>%
    dplyr::select(class,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_o <- clade_name_mp4_s_for_plot %>%
    group_by(order) %>%
    mutate(n=n()) %>%
    distinct(order,n,.keep_all = TRUE) %>%
    dplyr::select(order,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_f <- clade_name_mp4_s_for_plot %>%
    group_by(family) %>%
    mutate(n=n()) %>%
    distinct(family,n,.keep_all = TRUE) %>%
    dplyr::select(family,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_g <- clade_name_mp4_s_for_plot %>%
    group_by(genus) %>%
    mutate(n=n()) %>%
    distinct(genus,n,.keep_all = TRUE) %>%
    dplyr::select(genus,var_sort,n) %>%
    ungroup()
  
  clade_name_for_stack_s <- clade_name_mp4_s_for_plot %>%
    group_by(species) %>%
    mutate(n=n()) %>%
    distinct(species,n,.keep_all = TRUE) %>%
    dplyr::select(species,var_sort,n) %>%
    ungroup()
  
  
  colnames(clade_name_for_stack_k) <- c("clade_name","var_sort","count")
  clade_name_for_stack_k$clade_cat <- "kingdom"
  colnames(clade_name_for_stack_p) <- c("clade_name","var_sort","count")
  clade_name_for_stack_p$clade_cat <- "phylum"
  colnames(clade_name_for_stack_c) <- c("clade_name","var_sort","count")
  clade_name_for_stack_c$clade_cat <- "class"
  colnames(clade_name_for_stack_o) <- c("clade_name","var_sort","count")
  clade_name_for_stack_o$clade_cat <- "order"
  colnames(clade_name_for_stack_f) <- c("clade_name","var_sort","count")
  clade_name_for_stack_f$clade_cat <- "family"
  colnames(clade_name_for_stack_g) <- c("clade_name","var_sort","count")
  clade_name_for_stack_g$clade_cat <- "genus"
  colnames(clade_name_for_stack_s) <- c("clade_name","var_sort","count")
  clade_name_for_stack_s$clade_cat <- "species"
  
  clade_name_for_stack_all <- rbind(clade_name_for_stack_p,clade_name_for_stack_c,clade_name_for_stack_o,clade_name_for_stack_f,clade_name_for_stack_g,clade_name_for_stack_s)
  clade_name_color_matche_all <- clade_name_color_matche_all[,c("clade_name","clade_name_color")]
  clade_name_for_stack_all_color <- left_join(clade_name_for_stack_all,clade_name_color_matche_all,by="clade_name")
}


# 选取在图中显示的clade_name
clade_name_for_stack_all_color$clade_name_text <- ifelse(clade_name_for_stack_all_color$count >= 4, clade_name_for_stack_all_color$clade_name, "")
clade_name_for_stack_all_color$clade_name_text <- gsub("(p__|c__|o__|f__|g__|s__)\\s*", "", clade_name_for_stack_all_color$clade_name_text)
# Negativicutes, Lachnospiraceae_unclassified 和 Clostridiaceae_unclassified 与边上的在图中重叠，删除一个
clade_name_for_stack_all_color$clade_name_text <- ifelse(clade_name_for_stack_all_color$clade_name_text %in% c("Negativicutes","Lachnospiraceae_unclassified","Clostridiaceae_unclassified"), "", clade_name_for_stack_all_color$clade_name_text)

clade_name_for_stack_all_color$clade_cat <- factor(clade_name_for_stack_all_color$clade_cat,
                                                   levels = rev(c("phylum","class","order","family","genus","species")),
                                                   labels = rev(c("Phylum","Class","Order","Family","Genus","Species")))
# 逆序排列clade_name
clade_name_for_stack_all_color <- clade_name_for_stack_all_color %>%
  arrange(desc(clade_cat),-var_sort)
clade_name_for_stack_all_color$clade_name <- factor(clade_name_for_stack_all_color$clade_name, levels = clade_name_for_stack_all_color$clade_name)


# 构建stack plot函数
stack_function <- function(DAT){
  # 创建用于文本标签的数据框
  Label_DATA <- DAT %>%
    group_by(clade_cat) %>% # 按类别分组计算
    arrange(clade_name) %>% # 确保顺序与堆叠顺序一致（重要！）
    mutate(
      cumulative_count = cumsum(count),                       # 计算累积高度（每个堆叠块的上边界）
      mid_position = 127 - (cumulative_count - (count / 2))   # 计算每个堆叠块的中心位置 (实际图中text为右向左排列，所以需要用总长度127减去左边的位置)（重要！）
    ) %>%
    ungroup()
  
  f_temp <- ggplot(DAT, aes(
    x = count,  
    y = clade_cat, fill=clade_name))+
    geom_col(position = 'stack', color="#4b4b4b",  # 添加边框颜色
             size=0.1,
             width = 0.85) +    # 这个参数依然控制整体的组间间隔
    scale_fill_manual(values = DAT$clade_name_color) + # 应用自定义颜色
  
    # 添加文本标签
    geom_text(
      data = Label_DATA, # 指定新的标签数据
      aes(
        x = mid_position, # 标签放在每个堆叠块的中间
        y = clade_cat,
        label = clade_name_text
      ),
      size = 4.5,
      inherit.aes = FALSE # 忽略主图映射
    ) +
    
    theme_minimal() + # 不要背景
    
    theme(
      plot.margin = margin(10, 5, 5, 5), # t,r,b,l
      
      axis.ticks = element_blank(), # 去掉刻度线
      panel.grid = element_blank(),  # 删去网格线
      
      axis.title = element_blank(),
      
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 15, colour = "black"),
      
      legend.position = "none"
    ) +
    # 关键：expand 控制轴线与数据两端空白
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.03))) +  # x 轴左右留白 1%
    scale_y_discrete(expand = expansion(mult = c(0, 0)))    # y 轴上下留白 1%
  
  return(f_temp)
}
f_stack <- stack_function(clade_name_for_stack_all_color)
#### Stack plot (x轴sort 1:127，y轴 phylum:species) ####


#### Manhattan plot (x轴菌种，y轴相对丰度) ----
### 构建函数把species_name转换为可以作图的标准化名称 ###
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
### 构建函数把species_name转换为可以作图的标准化名称 ###

# 设置点的颜色分类
species_color <- clade_name_for_stack_all_color[,c("clade_name","clade_name_color")]
dat_sort_for_plot <- left_join(dat_sort, species_color, by = c("exp_name" = "clade_name"))
dat_sort_for_plot$exp_name_for_plot <- trans_mp4_s_names(dat_sort_for_plot$exp_name)

dat_sort_for_plot$clade_name_color <- ifelse(dat_sort$p >= 0.05, "non-sig", dat_sort_for_plot$clade_name_color)
dat_sort_for_plot$clade_name_color <- ifelse(dat_sort$p < 0.05 & dat_sort$estimate_cox_logistic > 0, "+", dat_sort_for_plot$clade_name_color)
dat_sort_for_plot$clade_name_color <- ifelse(dat_sort$p < 0.05 & dat_sort$estimate_cox_logistic < 0, "-", dat_sort_for_plot$clade_name_color)

dat_sort_for_plot$point_shape <- ifelse(dat_sort_for_plot$OUTCOME == "CVD", "cvd", "")
dat_sort_for_plot$point_shape <- ifelse(dat_sort_for_plot$OUTCOME == "CKD", "ckd", dat_sort_for_plot$point_shape)
dat_sort_for_plot$point_shape <- ifelse(dat_sort_for_plot$OUTCOME == "DM", "dm", dat_sort_for_plot$point_shape)

dat_sort_for_plot$point_size <- ifelse(dat_sort_for_plot$p >= 0.05, "cat1", "")
dat_sort_for_plot$point_size <- ifelse(dat_sort_for_plot$p < 0.05 & dat_sort_for_plot$p_adj_bh >= 0.05, "cat2", dat_sort_for_plot$point_size)
dat_sort_for_plot$point_size <- ifelse(dat_sort_for_plot$p_adj_bh < 0.05, "cat3", dat_sort_for_plot$point_size)

# 筛选各个疾病top5菌
top_gm <- dat_sort_for_plot %>%
  group_by(OUTCOME) %>%
  arrange(p) %>%
  filter(exp_name %in% mp4_s_names_short) %>%
  slice(1:5) %>%
  ungroup() %>%
  distinct(OUTCOME, exp_name)
top_gm$top_gm_flag <- 1
dat_sort_for_plot <- left_join(dat_sort_for_plot,top_gm,by=c("exp_name","OUTCOME"))
dat_sort_for_plot$text <- ifelse(dat_sort_for_plot$top_gm_flag == 1, dat_sort_for_plot$exp_name_for_plot, "")
dat_sort_for_plot$text <- ifelse(dat_sort_for_plot$p_adj_bh < 0.05, dat_sort_for_plot$exp_name_for_plot, dat_sort_for_plot$text)


# 构建曼哈顿图函数
manhattan_function <- function(DAT){
  f <- ggplot(DAT, aes(x = sort_new, y = p_trans, color = clade_name_color, shape = point_shape, size = point_size)) +
    geom_vline(xintercept = c(1:127), linetype = "solid", linewidth = 0.1, color = "#d7dbdc") +
    
    geom_point(alpha = 0.8) +  
    
    # 添加效应值正负标签
    geom_text(aes(label = text), color = "black", vjust = -2,
              size = 5) +

    # 设置颜色映射
    scale_color_manual(
      name = "",
      values = c(
        "+" = "#c44569",       
        "-" = "#546de5",
        "non-sig" = "#7f7f7f"
      ),
      labels = c(
        "+" = "Positive association",              
        "-" = "Negative association",
        "non-sig" = "Non-significant"
      ),
      breaks = c("+", "-", "non-sig")
    ) +
    
    # 设置形状映射
    scale_shape_manual(
      name = "",
      values = c(
        "dm" = 17,       # 三角形
        "ckd" = 15,      # 方形
        "cvd" = 16       # 圆形
      ),
      labels = c(
        "dm" = "Diabetes-related gut microbiome",              
        "ckd" = "CKD-related gut microbiome",
        "cvd" = "CVD-related gut microbiome"
      ),
      breaks = c("dm", "ckd", "cvd")
    ) +
    
    # 设置大小映射
    scale_size_manual(
      name = "",
      values = c(
        "cat1" = 2,      # 小
        "cat2" = 3,      # 中
        "cat3" = 5       # 大
      ),
      labels = c(
        "cat1" = "P>=0.05",              
        "cat2" = "P<0.05 and BH-adjusted P>=0.05",
        "cat3" = "BH-adjusted P<0.05"
      ),
      breaks = c("cat3", "cat2", "cat1")
    ) +
    
    # 添加阈值线
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.6, color = "red") +
    
    labs(x = "Species", y = "-log10(p-value)") +
    
    # 关键：expand 控制轴线与数据两端空白
    scale_x_continuous(
      breaks = unique(DAT$sort_new),  # 设置刻度位置为sort_new的值
      labels = unique(DAT$exp_name_for_plot),  # 设置刻度标签为对应的exp_name
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.005, 0.2))) +   # y 轴上下留白 1%
  
    # 设置主题
    theme_classic() +
    
    theme(
      plot.margin = margin(10, 40, 3, 25), # t,r,b,l
      
      panel.border = element_rect(color = "black", fill = NA, size = 0.5), # 面板区域边框
      
      axis.ticks.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.5),
      
      axis.text.y = element_text(color= "black", size = 15),
      axis.text.x = element_text(color= "black", angle = 300, hjust = 0, size = 10),
      
      axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 10)),
      axis.title.x = element_blank(),
      
      legend.position = "none",
      
      plot.title = element_blank()
    ) 
  
  return(f)
}

f_manhattan1 <- manhattan_function(dat_sort_for_plot)
#### Manhattan plot (x轴菌种，y轴相对丰度) ####

f1 <- cowplot::plot_grid(f_manhattan1, f_stack, 
                         nrow = 2,
                         rel_heights = c(10, 1.8))
ggsave(f1, filename="figures/main_figures/(fig5a)_manhattan_heatmap_clade.pdf", width = 25, height = 12, limitsize = FALSE)
