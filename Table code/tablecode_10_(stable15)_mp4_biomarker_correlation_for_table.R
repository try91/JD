library(data.table)
library(dplyr)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

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


#### 筛选 EDC-菌 pairs ----
## 读取EDC-菌 Spearman分析结果 ##
{
  p_spearman_results_mp4 <- readRDS("results/correlations/spearman/partial_spearman_results_edc-edc_mp4_out.rds")
  p_spearman_results_edc_mp4 <- p_spearman_results_mp4[[2]]
  p_spearman_results_edc_mp4 <- p_spearman_results_edc_mp4[p_spearman_results_edc_mp4$sample == "phy_edc_temp0_3",]
  # 选取两个EDC index和19种单个EDC以及359种菌种
  p_spearman_results_edc_mp4 <- p_spearman_results_edc_mp4[p_spearman_results_edc_mp4$exp_name %in% c("edc_count2_edc14_f","edc_count2_pfas_f","edc_count2_pae6_f","edc_count2_bp1_f","edc_count2_tc_f",
                                                                                                      "edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_bp1_f","edc_score_tc_f",
                                                                                                      edc_traits_log10) &
                                                             p_spearman_results_edc_mp4$out_name %in% c(mp4_s_log10),] # 359 species
  # 以每个EXP表型为单位进行校正（以EDC为组，校正每个菌）
  p_spearman_results_edc_mp4 <- p_spearman_results_edc_mp4 %>%
    group_by(exp_name) %>%  # 按exposure分组
    mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
    ungroup()
  p_spearman_short <- p_spearman_results_edc_mp4[,c("exp_name","out_name","estimate","p","p_adj_bh")]
}
# spearman 的EDC-菌 pairs
pair_edc_mp4_s <- p_spearman_short
table(pair_edc_mp4_s$exp_name)
#### 筛选 EDC-菌 pairs ####

#### 筛选 菌-Outcome pairs ----
## 结局：10类biomarker（用菌群人群，用14年指标）
outcome <- c("bmi_f","height_f","weight_f","whr_f","wc_f","hc_f",
             "tg_f","ldl_f","hdl_f","chol_f","apoa_f","apob_f","nonhdl_f",
             "alt_f","ast_f","ggt_f","bia_f",
             "egfr_f","scr_f","ua_f",
             "glu0_f","glu120_f","vhba1c_f",
             "ins0_f","ins120_f","homair_f","homab_f",
             "sbp_f","dbp_f","pr_f",
             "ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f",
             "hgb_f","plt_f","eos_f","lym_f","mon_f","neu_f",
             "nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_f","crp_f")

# 标准化名称
outcome_label <- c("BMI","Height","Weight","WHR","WC","HC",
                   "TG","LDL-C","HDL-C","TC","ApoA-1","ApoB","Non-HDL-C",
                   "ALT","AST","GGT","Bile acid",
                   "eGFR","Serum creatinine","UA",
                   "OGTT 0-h glucose","OGTT 2-h glucose","HbA1c",
                   "OGTT 0-h insulin","OGTT 2-h insulin","HOMA-IR","HOMA-B",
                   "SBP","DBP","PR",
                   "FT3","FT4","TSH","TPOAb","TgAb",
                   "Hemoglobin","Platelet count","Eosinophil count","Lymphocyte count","Monocyte count","Neutrophil count",
                   "NLR","LMR","PLR","SII","SIRI","WBC","Hs-CRP")


## 读取EDC-菌 Spearman分析结果 ##
p_spearman_results_mp4_out <- readRDS("results/correlations/spearman/partial_spearman_results_mp4-out.rds")
p_spearman_results_mp4_out <- p_spearman_results_mp4_out[[1]]
p_spearman_results_mp4_out <- p_spearman_results_mp4_out[p_spearman_results_mp4_out$sample == "phy_edc_temp0_3",]
p_spearman_results_mp4_out <- p_spearman_results_mp4_out[p_spearman_results_mp4_out$out_name %in% outcome,]

# 以每个EXP表型为单位进行校正（以EDC为组，校正每个菌）
p_spearman_results_mp4_out <- p_spearman_results_mp4_out %>%
  group_by(exp_name) %>%  # 按exposure分组
  mutate(p_adj_bh = p.adjust(p, method = "BH")) %>% # 对每个分组的P值进行FDR校正
  ungroup()
p_spearman_results_mp4_out_short <- p_spearman_results_mp4_out[,c("exp_name","out_name","estimate","p","p_adj_bh")]
#### 筛选 菌-Outcome pairs ####

# 挑选与EDC Score quartile (14 EDCs) 显著相关的菌
colnames(pair_edc_mp4_s)[c(3,4,5)] <- c("estimate","p","p_adj_bh")
pair_edc_mp4_s <- pair_edc_mp4_s[pair_edc_mp4_s$exp_name %in% c("edc_score_edc14_f"),]
pair_edc_mp4_s_sig <- pair_edc_mp4_s[pair_edc_mp4_s$p_adj_bh < 0.05,] 

unique(pair_edc_mp4_s_sig[pair_edc_mp4_s_sig$out_name %in% mp4_s_log10_short,]$out_name) # 62 有明确属名称 for fig4b
unique(pair_edc_mp4_s_sig$out_name) # 75 包括无明确属名称 for stable

p_spearman_results_mp4_out_short_keep <- p_spearman_results_mp4_out_short[p_spearman_results_mp4_out_short$exp_name %in% unique(pair_edc_mp4_s_sig$out_name),]
table(p_spearman_results_mp4_out_short_keep$out_name)



dat_for_stable <- p_spearman_results_mp4_out_short_keep
dat_for_stable$exp_name <- factor(dat_for_stable$exp_name, 
                                  levels = mp4_s_log10,
                                  labels = trans_mp4_s_names(mp4_s_log10))
dat_for_stable$out_name <- factor(dat_for_stable$out_name, 
                                  levels = outcome,
                                  labels = outcome_label)

colnames(dat_for_stable) <- c("Gut microbial species","Clinical biomarkers","Spearman's rho (two-sided)","P value","BH-adjusted P")
openxlsx::write.xlsx(dat_for_stable,"tables/(stable15)_mp4_biomarker_correaltion.xlsx")
