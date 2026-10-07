library(data.table)
library(dplyr)

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

#### 读取EDC INDEX 变量名 ----
edc_index_b_keep <- c("edc_count2_edc14_b","edc_count2_pfas_b","edc_count2_pae6_b","edc_count2_bp1_b","edc_count2_tc_b",
                      "edc_score_edc14_b","edc_score_pfas_b","edc_score_pae6_b","edc_score_bp1_b","edc_score_tc_b")
edc_index_f_keep <- c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f")
#### 读取EDC INDEX数据 ####

#### 读取interaction分析结果 ----
# 根据外来物质有害理论，我们只保留EDC高对于结局有危害效应的结果 #
results_both_interaction_sig_keep <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig_20260728.xlsx")
unique_both_interaction <- unique(results_both_interaction_sig_keep$keep_exp_med_out)
results_multiplicative_interaction_sig_keep <- readxl::read_xlsx("results/cox/interaction/multi_interaction_sig_20260728.xlsx")
unique_multiplicative_interaction <- unique(results_multiplicative_interaction_sig_keep$keep_exp_med_out)
results_additive_interaction_sig_keep <- readxl::read_xlsx("results/cox/interaction/add_interaction_sig_20260728.xlsx")
unique_additive_interaction <- unique(results_additive_interaction_sig_keep$keep_exp_med_out)


results_multiplicative_interaction_sig_keep_only <- results_multiplicative_interaction_sig_keep[!results_multiplicative_interaction_sig_keep$keep_exp_med_out %in% unique_both_interaction,]
unique_multiplicative_interaction <- unique(results_multiplicative_interaction_sig_keep_only$keep_exp_med_out)
results_additive_interaction_sig_keep_only <- results_additive_interaction_sig_keep[!results_additive_interaction_sig_keep$keep_exp_med_out %in% unique_both_interaction,]
unique_additive_interaction <- unique(results_additive_interaction_sig_keep_only$keep_exp_med_out)

all_results_interaction <- rbind(results_both_interaction_sig_keep,results_multiplicative_interaction_sig_keep_only,results_additive_interaction_sig_keep_only)
unique_all_interaction <- unique(all_results_interaction$keep_exp_med_out)

all_results_interaction$interaction_type <- ifelse(all_results_interaction$keep_exp_med_out %in% unique_both_interaction, "Both interaction",
                                                   ifelse(all_results_interaction$keep_exp_med_out %in% unique_multiplicative_interaction, "Multiplicative interaction only",
                                                          ifelse(all_results_interaction$keep_exp_med_out %in% unique_additive_interaction, "Additive interaction only","")))
all_results_interaction$outcome <- gsub("_no_self_report","",all_results_interaction$outcome)
all_results_interaction$keep_exp_med_out <- gsub("_no_self_report","",all_results_interaction$keep_exp_med_out)
unique(all_results_interaction$keep_exp_med_out)
#### 读取interaction分析结果 ####


#### 限制条件1: 读取replication验证方向一致,且FDR-P<0.2的结果 ----
# results_dm_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_dm_20260408.xlsx")
# results_dm_replication <- results_dm_replication %>%
#   distinct(species)
# 
# results_ckd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_ckd_20260408.xlsx")
# results_ckd_replication <- results_ckd_replication %>%
#   distinct(species)
# 
# results_acvd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_acvd_20260424.xlsx")
# results_acvd_replication <- results_acvd_replication %>%
#   distinct(species)

# 验证队列中方向一致,且FDR-P<0.2的结果
results_dm_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_dm_20260728.xlsx")
results_dm_replication <- results_dm_replication %>%
  distinct(species)

results_ckd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_ckd_20260728.xlsx")
results_ckd_replication <- results_ckd_replication %>%
  distinct(species)

results_acvd_replication <- readxl::read_xlsx("results/replication/replication_same_direction/replication_same_direction_sig_acvd_20260728.xlsx")
results_acvd_replication <- results_acvd_replication %>%
  distinct(species)

# 只保留验证结果与主结果菌种对于结局方向一致,且FDR-P<0.2的菌种结果 #
all_results_interaction <- all_results_interaction[(all_results_interaction$OUTCOME == "DM" & all_results_interaction$mediator %in% results_dm_replication$species) |
                                                     (all_results_interaction$OUTCOME == "CKD" & all_results_interaction$mediator %in% results_ckd_replication$species) |
                                                     (all_results_interaction$OUTCOME == "CVD" & all_results_interaction$mediator %in% results_acvd_replication$species),]
#### 限制条件1: 读取replication验证方向一致,且FDR-P<0.2的结果 ####

#### 限制条件️2: 读取sensitivity interaction分析验证后的结果 ----
dat_validate_same_direction_sig <- readxl::read_xlsx("results/cox/interaction/both_interaction_sig_validated_resutls_20260728.xlsx")

# 只保留敏感性分析与主结果一致的菌种 #
all_results_interaction <- all_results_interaction[all_results_interaction$keep_exp_med_out %in% dat_validate_same_direction_sig$keep_exp_med_out,] # 只选取验证后的结果
#### 限制条件2: 读取sensitivity interaction分析验证后的结果 ####


#### interaction结果数据处理 ----
all_results_interaction$text <- sprintf("%.3f (%.3f, %.3f)", all_results_interaction$exp.coef., all_results_interaction$exp.lci, all_results_interaction$exp.uci)

all_results_interaction <- all_results_interaction[all_results_interaction$rowname %in% c("EDC high","EDC high:GM high","EDC high:GM low","Multiplicative scale","RERI","AP","SI"),]
all_results_interaction <- all_results_interaction[!(all_results_interaction$rowname %in% c("EDC high") & all_results_interaction$group %in% c("interaction_recode")),]

all_results_interaction_short <- all_results_interaction[,c("rowname","keep_exp_med_out","group","interaction_type", "exposure","mediator","outcome","n","text","Pr...z..")]
unique(all_results_interaction_short$keep_exp_med_out)


# 排序 #
all_results_interaction_short$sig_flag <- ifelse(all_results_interaction_short$Pr...z.. < 0.05, 1, 2)
all_results_interaction_short$exposure_group <- ifelse(all_results_interaction_short$exposure %in% c("edc_count2_edc14_f","edc_score_edc14_f"), 1,
                                                       ifelse(all_results_interaction_short$exposure %in% c("edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10), 2,
                                                              ifelse(all_results_interaction_short$exposure %in% c("edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10), 3,
                                                                     ifelse(all_results_interaction_short$exposure %in% c("edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10), 4,
                                                                            ifelse(all_results_interaction_short$exposure %in% c("edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10), 5, 0)))))
all_results_interaction_short$exposure <- factor(all_results_interaction_short$exposure,
                                    levels = c("edc_count2_edc14_f","edc_score_edc14_f",
                                               "edc_count2_pfas_f","edc_score_pfas_f",edc_traits2_log10,
                                               "edc_count2_pae6_f","edc_score_pae6_f",edc_traits3_log10,
                                               "edc_count2_tc_f","edc_score_tc_f",edc_traits5_log10,
                                               "edc_count2_bp1_f","edc_score_bp1_f",edc_traits4_log10),
                                    labels = c("EDC Scoremedian (14 EDCs)", "EDC Scorequartile (14 EDCs)",
                                               "EDC Scoremedian (PFAS)", "EDC Scorequartile (PFAS)", edc_traits2,
                                               "EDC Scoremedian (PAEs)", "EDC Scorequartile (PAEs)", edc_traits3,
                                               "EDC Scoremedian (antimicrobials)", "EDC Scorequartile (antimicrobials)", edc_traits5,
                                               "EDC Scoremedian (bisphenols)", "EDC Scorequartile (bisphenols)", edc_traits4))
all_results_interaction_short$mediator <- trans_mp4_s_names(all_results_interaction_short$mediator)
all_results_interaction_short$outcome <- factor(all_results_interaction_short$outcome, 
                                                levels = c("dm_incident_1014","ckd_incident_1014","cvd_incident_1021"),
                                                labels = c("Diabetes","CKD","CVD"))

# 挑选高菌群组EDC效应 #
dat1 <- all_results_interaction_short[all_results_interaction_short$group == "mp4_high",c("keep_exp_med_out","interaction_type", "exposure","mediator","outcome","n","text","Pr...z..", "exposure_group","sig_flag")]
dat1 <- dat1 %>%
  arrange(exposure_group,sig_flag,exposure,mediator,outcome)

dat1 <- dat1[,c("keep_exp_med_out","interaction_type", "exposure","mediator","outcome","n","text","Pr...z..")]
colnames(dat1)[2:8] <- c("Interaction types","Analyte-based scores or individual analytes","Gut microbial species","Outcomes","N-high_mp4","HR (95% CI)-high_mp4","P value")
# 挑选低菌群组EDC效应 #
dat2 <- all_results_interaction_short[all_results_interaction_short$group == "mp4_low",c("keep_exp_med_out","n","text","Pr...z..")]
colnames(dat2)[2:4] <- c("N-low_mp4","HR (95% CI)-low_mp4","P value")

# 挑选multiplicative interaction效应 #
dat3 <- all_results_interaction_short[all_results_interaction_short$rowname %in% c("EDC high:GM high","EDC high:GM low"),c("keep_exp_med_out","rowname","text","Pr...z..")]
dat3$rowname <- ifelse(dat3$rowname == "EDC high:GM high", "Low EDC and low species abundance", 
                       ifelse(dat3$rowname == "EDC high:GM low", "Low EDC and high species abundance", ""))
colnames(dat3)[2:4] <- c("Reference groups for interaction","Multiplicative scale (95% CI)","P for interaction")
# dat3_1 <- all_results_interaction_short[all_results_interaction_short$rowname %in% c("Multiplicative scale"),c("keep_exp_med_out","text","Pr...z..")]
# colnames(dat3_2)[2:3] <- c("HR (95% CI)-multiplicative","P for interaction")
# dat3_2 <- left_join(dat3, dat3_1, by="keep_exp_med_out")

# 挑选RERI效应 #
dat4 <- all_results_interaction_short[all_results_interaction_short$rowname %in% c("RERI"),c("keep_exp_med_out","text")]
colnames(dat4)[2] <- c("RERI (95% CI)")
# 挑选AP效应 #
dat5 <- all_results_interaction_short[all_results_interaction_short$rowname %in% c("AP"),c("keep_exp_med_out","text")]
colnames(dat5)[2] <- c("AP (95% CI)")
# 挑选SI效应 #
dat6 <- all_results_interaction_short[all_results_interaction_short$rowname %in% c("SI"),c("keep_exp_med_out","text")]
dat6$text <- ifelse(dat6$text == "NA (NA, NA)", "/", dat6$text)
colnames(dat6)[2] <- c("SI (95% CI)")

# 合并各项结果
dat_all <- left_join(dat1, dat2, by="keep_exp_med_out") %>%
  left_join(dat3, by="keep_exp_med_out") %>%
  left_join(dat4, by="keep_exp_med_out") %>%
  left_join(dat5, by="keep_exp_med_out") %>%
  left_join(dat6, by="keep_exp_med_out")

write.table(dat_all$keep_exp_med_out, "results/cox/interaction/final_17_pairs_for_stable_20260802.txt", row.names = FALSE)

dat_all[,c(1:2,6,9)] <- NULL
#### interaction结果数据处理 ####


openxlsx::write.xlsx(dat_all,paste0("tables/interaction_results_20260802.xlsx"))
