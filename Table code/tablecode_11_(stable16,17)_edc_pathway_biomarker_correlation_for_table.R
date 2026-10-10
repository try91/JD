library(data.table)
library(dplyr)

setwd("file_path") # File path includes "raw_data", "results", "figures", and "tables" folders

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
phy_out_cat <- c("dm_f","ckd_f","cvd_f")
phy_traits_cat <- c("dm_f","ckd_f","cvd_f","dm_b","ckd_b","cvd_b",
                    "ob_f","ob_b","abob_f","abob_b","ir_f","ir_b","dyslip_f","dyslip_b",
                    "mets_f","mets_b","nafld_f","nafld_b","hua_f","hua_b","hpt_f","hpt_b","as_imt_f","as_imt_b")
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

#### Pathway变量名 ----
# 读取Pathway名称
pathway_mp3 <- read.table("raw_data/pathway_names_mp3_10%.txt",header = TRUE)
pathway_mp3_names_labels <- pathway_mp3[,1]
pathway_mp3_names <- pathway_mp3[,2]
pathway_mp3_log10 <- paste0(pathway_mp3_names,"_log10") # Pathway丰度的log10转换

# 读取Pathway数据
pathway_mp3_anno <- fread("raw_data/Pathway.annotation_for_edc_micro.csv",header = TRUE)
pathway_mp3_anno <- pathway_mp3_anno[pathway_mp3_anno$V1 %in% pathway_mp3_names_labels,]
pathway_mp3_anno$L2 <- ifelse(pathway_mp3_anno$L2 == "Superpathways" & !is.na(pathway_mp3_anno$L3), pathway_mp3_anno$L3, pathway_mp3_anno$L2)
pathway_mp3_anno$L2 <- ifelse(pathway_mp3_anno$L2 == "Other Biosynthesis" & !is.na(pathway_mp3_anno$L3), pathway_mp3_anno$L3, pathway_mp3_anno$L2)
pathway_mp3_anno$L2 <- ifelse(pathway_mp3_anno$L2 == "Degradation/Utilization/Assimilation - Other" & !is.na(pathway_mp3_anno$L3), pathway_mp3_anno$L3, pathway_mp3_anno$L2)
pathway_mp3_anno$L2 <- ifelse(pathway_mp3_anno$L2 == "nylon-6 oligomer degradation", "Nylon-6 Oligomer Degradation", pathway_mp3_anno$L2)
pathway_mp3_anno$L2 <- ifelse(pathway_mp3_anno$L2 == "8-Amino-7-oxononanoate Biosynthesis", "8-Amino-7-Oxononanoate Biosynthesis", pathway_mp3_anno$L2)
# 筛选Degradation Pathway
pathway_mp3_anno_deg <- pathway_mp3_anno %>%
  filter(if_any(L2, ~ grepl("Degradation|degradation|oxidation", .))) %>%
  arrange(L2,V1)
# 筛选Biosynthesis Pathway
pathway_mp3_anno_bio <- pathway_mp3_anno %>%
  filter(if_any(L2, ~ grepl("Biosynthesis|biosynthesis", .))) %>%
  arrange(L2,V1)
# 筛选Other Pathway
pathway_mp3_anno_other <- pathway_mp3_anno[!(pathway_mp3_anno$V1 %in% pathway_mp3_anno_deg$V1 | pathway_mp3_anno$V1 %in% pathway_mp3_anno_bio$V1),] %>%
  arrange(L2,V1)

# Pathway重排列
pathway_mp3_anno <- rbind(pathway_mp3_anno_deg,pathway_mp3_anno_bio,pathway_mp3_anno_other) %>%
  distinct(V1,.keep_all = TRUE)
pathway_l2_sort <- unique(pathway_mp3_anno$L2)

pathway_mp3_anno_short <- pathway_mp3_anno[,c("V1","L2")] 
#### Pathway变量名 ####


#### 读取EDC-Pathway相关性分析结果 ----
spearman_edc_results <- readRDS("results/correlations/spearman/partial_spearman_results_edc-pathway.rds")
spearman_edc_results <- spearman_edc_results[[1]]
spearman_edc_results <- spearman_edc_results[spearman_edc_results$exp_name %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"),]

# 以每个OUT表型为单位进行校正
spearman_edc_results <- spearman_edc_results %>%
  group_by(out_name) %>%  # 按outcome分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup() %>%
  group_by(out_name) %>%  # 按outcome分组
  mutate(p_adj_bf = p.adjust(p, method = "bonferroni")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
spearman_edc_results <- spearman_edc_results[spearman_edc_results$out_name %in% pathway_mp3_log10,]

spearman_edc_results$exp_name_label <- factor(spearman_edc_results$exp_name, 
                                              levels = c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"), 
                                              labels = c("EDC Scorequartile (14 EDCs)", 
                                                         "EDC Scorequartile (PFAS)", 
                                                         "EDC Scorequartile (PAEs)", 
                                                         "EDC Scorequartile (antimicrobials)", 
                                                         "EDC Scorequartile (bisphenols)"))
spearman_edc_results$pathway_name_label <- factor(spearman_edc_results$out_name, 
                                                  levels = pathway_mp3_log10, 
                                                  labels = pathway_mp3_names_labels)
spearman_edc_results <- left_join(spearman_edc_results, pathway_mp3_anno_short, by=c("pathway_name_label" = "V1"))
#### 读取EDC-Pathway相关性分析结果 ####
unique(spearman_edc_results$pathway_name_label)

#### 挑选与EDC Score (14 EDCs) 相关的 ----
pair_edc_mp4_s <- spearman_edc_results[,c("exp_name","out_name","estimate","p","p_adj_bh","exp_name_label","pathway_name_label","L2")]
pair_edc_mp4_s <- pair_edc_mp4_s[pair_edc_mp4_s$exp_name %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f"),]

# 保留与14 EDC score显著相关菌
pair_edc_mp4_s_mp4_14edcscore <- pair_edc_mp4_s[pair_edc_mp4_s$exp_name == "edc_score_edc14_f",]
pair_edc_mp4_s_mp4_14edcscore_sig <- pair_edc_mp4_s_mp4_14edcscore[pair_edc_mp4_s_mp4_14edcscore$p_adj_bh < 0.05,]

pair_edc_mp4_s_mp4_14edcscore_sig$L2 <- factor(pair_edc_mp4_s_mp4_14edcscore_sig$L2,
                                               levels = pathway_l2_sort)
pair_edc_mp4_s_mp4_14edcscore_sig <- pair_edc_mp4_s_mp4_14edcscore_sig %>%
  arrange(L2, -estimate)

pathway_v1_sort <- unique(pair_edc_mp4_s_mp4_14edcscore_sig$pathway_name_label)
pathway_v1_l2_sort <- pair_edc_mp4_s_mp4_14edcscore_sig[,c("pathway_name_label","L2")]
pathway_l2_sort <- as.character(unique(pair_edc_mp4_s_mp4_14edcscore_sig$L2)) 

# 挑选所有EDC score中 14 EDC score显著的结果
pair_edc_mp4_s4 <- pair_edc_mp4_s[pair_edc_mp4_s$pathway_name_label %in% pathway_v1_sort,]
pair_edc_mp4_s4$exp_name_label <- factor(pair_edc_mp4_s4$exp_name_label, 
                                         levels = c("EDC Scorequartile (14 EDCs)", 
                                                    "EDC Scorequartile (PFAS)", 
                                                    "EDC Scorequartile (PAEs)", 
                                                    "EDC Scorequartile (antimicrobials)", 
                                                    "EDC Scorequartile (bisphenols)"))
pair_edc_mp4_s4$pathway_name_label <- factor(pair_edc_mp4_s4$pathway_name_label, 
                                             levels = pathway_v1_sort)
unique(pair_edc_mp4_s4$pathway_name_label)
#### 挑选与EDC Score (14 EDCs) 相关的 ####


#### 读取Pathway-Biomarker相关性分析结果 ----
## 结局：10类biomarker（用菌群人群，用14年指标）
outcome <- c("bmi_f","height_f","weight_f","whr_f","wc_f","hc_f",
             "tg_f","ldl_f","hdl_f","chol_f","apoa_f","apob_f","nonhdl_f",
             "alt_f","ast_f","ggt_f","bia_f",
             "egfr_f","scr_f","ua_f", #"acr_f",
             "glu0_f","glu120_f","vhba1c_f",
             "ins0_f","ins120_f","homair_f","homab_f",
             "sbp_f","dbp_f","pr_f", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
             "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
             "hgb_f","plt_f","eos_f","lym_f","mon_f","neu_f",
             "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_f","crp_f")

# 标准化名称
outcome_label <- c("BMI","Height","Weight","WHR","WC","HC",
                   "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                   "ALT","AST","GGT","Bile acid",
                   "eGFR","Serum creatinine","UA", #"UACR",
                   "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                   "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                   "SBP","DBP","PR", # 2026.02.13 备注: 脉率(Pulse Rate，PR)不是脉压差 #
                   "FT3","FT4","TSH","TPOAb","TgAb",
                   "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                   "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")


spearman_out_results <- readRDS("results/correlations/spearman/partial_spearman_results_pathway-out.rds")
spearman_out_results <- spearman_out_results[[1]]
spearman_out_results <- spearman_out_results[spearman_out_results$out_name %in% c("bmi_f","height_f","weight_f","whr_f","wc_f","hc_f",
                                                                                  "tg_f","ldl_f","hdl_f","chol_f","apoa_f","apob_f","nonhdl_f",
                                                                                  "alt_f","ast_f","ggt_f","bia_f",
                                                                                  "egfr_f","scr_f","ua_f",
                                                                                  "glu0_f","glu120_f","vhba1c_f",
                                                                                  "ins0_f","ins120_f","homair_f","homab_f",
                                                                                  "sbp_f","dbp_f","pr_f",
                                                                                  "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
                                                                                  "hgb_f","plt_f","eos_f","lym_f","mon_f","neu_f",
                                                                                  "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_f","crp_f"),]

# 以每个EXP表型为单位进行校正
spearman_out_results <- spearman_out_results %>%
  group_by(exp_name) %>%  # 按exp分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup() %>%
  group_by(exp_name) %>%  # 按exp分组
  mutate(p_adj_bf = p.adjust(p, method = "bonferroni")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
spearman_out_results <- spearman_out_results[spearman_out_results$exp_name %in% pathway_mp3_log10,]
spearman_out_results <- spearman_out_results[spearman_out_results$out_name %in% outcome,]

spearman_out_results$pathway_name_label <- factor(spearman_out_results$exp_name, 
                                                  levels = pathway_mp3_log10, 
                                                  labels = pathway_mp3_names_labels)
spearman_out_results$out_name_label <- factor(spearman_out_results$out_name, levels = outcome, labels = outcome_label)
spearman_out_results <- left_join(spearman_out_results, pathway_mp3_anno_short, by=c("pathway_name_label" = "V1"))
#### 读取Pathway-Biomarker相关性分析结果 ####

#### Pathway和Biomarker关联筛选 ----
spearman_out_results <- spearman_out_results[spearman_out_results$pathway_name_label %in% pathway_v1_sort,]
spearman_out_results <- spearman_out_results[,c("exp_name","out_name","estimate","p","p_adj_bh","pathway_name_label","out_name_label","L2")]

spearman_out_results <- spearman_out_results[spearman_out_results$exp_name %in% pair_edc_mp4_s4$out_name,]

spearman_out_results$pathway_name_label <- factor(spearman_out_results$pathway_name_label, 
                                                  levels = pathway_v1_sort)

spearman_out_results$out_name_label <- factor(spearman_out_results$out_name_label,
                                              levels = outcome_label)
#### Pathway和Biomarker关联筛选 ----


dat_edc_pathway_for_stable <- pair_edc_mp4_s4[,c("exp_name_label","pathway_name_label","L2","estimate","p","p_adj_bh")]
dat_edc_pathway_for_stable$L2 <- factor(dat_edc_pathway_for_stable$L2, levels = pathway_l2_sort)
dat_edc_pathway_for_stable <- dat_edc_pathway_for_stable %>%
  arrange(exp_name_label, L2, pathway_name_label)
# unique(dat_edc_pathway_for_stable$exp_name_label)
# unique(dat_edc_pathway_for_stable$pathway_name_label)
# unique(dat_edc_pathway_for_stable$L2)


dat_pathway_biomarker_for_stable <- spearman_out_results[,c("pathway_name_label","L2","out_name_label","estimate","p","p_adj_bh")]
dat_pathway_biomarker_for_stable$L2 <- factor(dat_pathway_biomarker_for_stable$L2, levels = pathway_l2_sort)
dat_pathway_biomarker_for_stable <- dat_pathway_biomarker_for_stable %>%
  arrange(L2, pathway_name_label, out_name_label)
# unique(dat_pathway_biomarker_for_stable$pathway_name_label)
# unique(dat_pathway_biomarker_for_stable$L2)
# unique(dat_pathway_biomarker_for_stable$out_name_label)


colnames(dat_edc_pathway_for_stable) <- c("Analyte-based scores","Gut microbial pathways","Pathway categories","Spearman's rho (two-sided)","P value","BH-adjusted P")
colnames(dat_pathway_biomarker_for_stable) <- c("Gut microbial pathways","Pathway categories","Clinical biomarkers","Spearman's rho (two-sided)","P value","BH-adjusted P")

openxlsx::write.xlsx(dat_edc_pathway_for_stable,"tables/(stable16)_pathway_edc_correaltion.xlsx")
openxlsx::write.xlsx(dat_pathway_biomarker_for_stable,"tables/(stable17)_pathway_biomarker_correaltion.xlsx")
