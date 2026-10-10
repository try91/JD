library(data.table)
library(dplyr)
library(ggplot2)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

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

#### 限制条件1: 读取replication验证方向一致(且FDR-P<0.2)的结果 ----
# 验证队列中方向一致,且FDR-P<0.2的结果
results_dm_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_dm.xlsx")
results_dm_replication <- results_dm_replication %>%
  distinct(species)
results_ckd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_ckd.xlsx")
results_ckd_replication <- results_ckd_replication %>%
  distinct(species)
results_acvd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_acvd.xlsx")
results_acvd_replication <- results_acvd_replication %>%
  distinct(species)
#### 限制条件1: 读取replication验证方向一致(且FDR-P<0.2)的结果 ####

#### 限制条件️2: 读取sensitivity interaction分析验证后的结果 ----
dat_validate_same_direction_sig <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig_validated_resutls.xlsx")
#### 限制条件2: 读取sensitivity interaction分析验证后的结果 ####


#### 读取+处理interaction分析结果 ----
# "相乘交互作用"显著(Multiplicative scale P<0.05), 同时"相加交互作用"显著(以RERI显著为标准)的组
results_interaction_sig1 <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig.xlsx")
results_interaction_sig1$outcome <- gsub("_no_self_report","",results_interaction_sig1$outcome)
results_interaction_sig1$keep_exp_med_out <- gsub("_no_self_report","",results_interaction_sig1$keep_exp_med_out)
unique(results_interaction_sig1$keep_exp_med_out)

# 限制条件1: 只保留验证结果与主结果菌种对于结局方向一致(且FDR-P<0.2)的菌种结果 #
results_interaction_sig1 <- results_interaction_sig1[(results_interaction_sig1$OUTCOME == "DM" & results_interaction_sig1$mediator %in% results_dm_replication$species) |
                                                       (results_interaction_sig1$OUTCOME == "CKD" & results_interaction_sig1$mediator %in% results_ckd_replication$species) |
                                                       (results_interaction_sig1$OUTCOME == "CVD" & results_interaction_sig1$mediator %in% results_acvd_replication$species),]
unique(results_interaction_sig1$keep_exp_med_out)
# 限制条件2: 只保留敏感性分析与主结果一致的菌种 #
results_interaction_sig1 <- results_interaction_sig1[results_interaction_sig1$keep_exp_med_out %in% dat_validate_same_direction_sig$keep_exp_med_out,] # 只选取验证后的结果
unique(results_interaction_sig1$keep_exp_med_out)
unique(results_interaction_sig1$mediator) # 13个菌,对应17组EDC-菌-Outcome组合
# 保留有明确菌属信息的菌种结果 #
results_interaction_sig1 <- results_interaction_sig1[results_interaction_sig1$mediator %in% mp4_s_log10_short,]
unique(results_interaction_sig1$keep_exp_med_out)


results_interaction_sig1_dm <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "DM",]
results_interaction_sig1_dm$hr_95ci <- sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_dm$exp.coef., results_interaction_sig1_dm$exp.lci, results_interaction_sig1_dm$exp.uci)
results_interaction_sig1_dm$text <- paste0(results_interaction_sig1_dm$rowname,": ",results_interaction_sig1_dm$hr_95ci)
results_interaction_sig1_dm$text <- ifelse(results_interaction_sig1_dm$text == "SI: NA (NA, NA)", "", results_interaction_sig1_dm$text)
results_interaction_sig1_dm$p_text <- ifelse(results_interaction_sig1_dm$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_dm$Pr...z..,3)))
results_interaction_sig1_dm <- results_interaction_sig1_dm[results_interaction_sig1_dm$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_dm <- results_interaction_sig1_dm[!(results_interaction_sig1_dm$rowname %in% c("EDC high") & results_interaction_sig1_dm$group %in% c("interaction_recode")),]

results_interaction_sig1_ckd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CKD",]
results_interaction_sig1_ckd$hr_95ci <- sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_ckd$exp.coef., results_interaction_sig1_ckd$exp.lci, results_interaction_sig1_ckd$exp.uci)
results_interaction_sig1_ckd$text <- paste0(results_interaction_sig1_ckd$rowname,": ",results_interaction_sig1_ckd$hr_95ci)
results_interaction_sig1_ckd$text <- ifelse(results_interaction_sig1_ckd$text == "SI: NA (NA, NA)", "", results_interaction_sig1_ckd$text)
results_interaction_sig1_ckd$p_text <- ifelse(results_interaction_sig1_ckd$Pr...z.. < 0.001, "P for interaction<0.001", paste0("P for interaction=",round(results_interaction_sig1_ckd$Pr...z..,3)))
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[results_interaction_sig1_ckd$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","RERI","AP","SI"),]
results_interaction_sig1_ckd <- results_interaction_sig1_ckd[!(results_interaction_sig1_ckd$rowname %in% c("EDC high") & results_interaction_sig1_ckd$group %in% c("interaction_recode")),]

results_interaction_sig1_cvd <- results_interaction_sig1[results_interaction_sig1$OUTCOME == "CVD",]
results_interaction_sig1_cvd$hr_95ci <-sprintf("%.2f (%.2f, %.2f)", results_interaction_sig1_cvd$exp.coef., results_interaction_sig1_cvd$exp.lci, results_interaction_sig1_cvd$exp.uci)
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

############################################# Table for main figure #############################################
#### Both multiplicative and additive interaction for main Fig 6a ----
interaction_both <- results_interaction_sig1

# interaction类型
interaction_both$interaction_type <- "both"

# MP4的高低
interaction_both$mp4_level <- ifelse(interaction_both$group == "mp4_high", "high", 
                                            ifelse(interaction_both$group == "mp4_low","low", NA))

interaction_both <- interaction_both[!is.na(interaction_both$mp4_level),]
interaction_both <- interaction_both[interaction_both$exp.coef. > 1 & interaction_both$sig_flag == "sig",]

interaction_both <- interaction_both[,c("exposure","mediator","OUTCOME","mp4_level")]
interaction_both$exposure_group <- ifelse(interaction_both$exposure %in% c("edc_count2_edc14_f","edc_score_edc14_f"), 1,
                                          ifelse(interaction_both$exposure %in% c("edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10), 2,
                                                 ifelse(interaction_both$exposure %in% c("edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10), 3,
                                                        ifelse(interaction_both$exposure %in% c("edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10), 4,
                                                               ifelse(interaction_both$exposure %in% c("edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10), 5, 0)))))
interaction_both$exposure <- factor(interaction_both$exposure,
                                    levels = c("edc_count2_edc14_f","edc_score_edc14_f",
                                               "edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10,
                                               "edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10,
                                               "edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10,
                                               "edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10),
                                    labels = c("EDC Score (14)", "EDC Score (14)",
                                               "EDC Score (PFAS)", "EDC Score (PFAS)", edc_traits2,
                                               "EDC Score (PAEs)", "EDC Score (PAEs)", edc_traits3,
                                               "EDC Score (TC)", "EDC Score (TC)", edc_traits5,
                                               "EDC Score (BP)", "EDC Score (BP)", edc_traits4))
interaction_both$mediator <- trans_mp4_s_names(interaction_both$mediator)
interaction_both$mp4_level <- factor(interaction_both$mp4_level, 
                                     levels = c("high","low"),
                                     labels = c("High abundance","Low abundance"))
interaction_both$OUTCOME <- factor(interaction_both$OUTCOME, 
                                   levels = c("DM","CKD","CVD"),
                                   labels = c("Diabetes","CKD","CVD"))
interaction_both <- interaction_both %>%
  arrange(exposure_group,mp4_level,exposure,mediator,OUTCOME)

# 找出：除了第一列以外，其他所有列的数值都完全一样的重复行，给这些行标记 same，否则为空。
interaction_both <- interaction_both %>%
  group_by(across(-1)) %>%
  mutate(
    same_flag = ifelse(n() >= 2, "same", "")
  ) %>%
  ungroup() %>%
  # 生成连续不跳号的配对ID
  group_by(across(-1)) %>%
  mutate(
    pair_id = ifelse(
      same_flag == "same",
      as.character(cur_group_id()),
      ""
    )
  ) %>%
  ungroup()

# 删除所有完全相同的重复行（只保留第一次出现的行）
interaction_both <- interaction_both %>%
  distinct(.keep_all = TRUE)

openxlsx::write.xlsx(interaction_both, "figures/main_figures/(fig6a)_both_interaction.xlsx")
#### Both multiplicative and additive interaction for main Fig 6a ####
