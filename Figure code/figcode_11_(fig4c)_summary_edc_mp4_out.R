library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("C:/TWang/DLiu/EDC_Micro/submission") # File path includes "raw_data", "results", "figures", and "tables" folders

micro_dat <- read.table("raw_data/microbial_composition_pathway_dat_20261006.txt", header = TRUE)

#### 配色 ----
PFAS_colors <- colorRampPalette(c("#f55d78","#FFFFFF"))(50)[c(1,8,15,22,29)]
PAE_colors <- colorRampPalette(c("#3498db","#FFFFFF"))(50)[c(1,6,11,16,21,26,31,36,41)]
BP_colors <- colorRampPalette(c("#f1c40f","#FFFFFF"))(50)[c(1,18,35)]
TC_colors <- colorRampPalette(c("#00b894","#FFFFFF"))(50)[c(1,21)]
my_colors <- c(PFAS_colors,PAE_colors,BP_colors,TC_colors)
# 创建25色渐变色标尺
my_palette <- colorRampPalette(colors = c("#cf6a87", "#f19066", "#f5cd79", "#2ecc71", "#34ace0"))(5)
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
pair_edc_mp4_s_sig <- pair_edc_mp4_s[pair_edc_mp4_s$p_adj_bh < 0.05,] # (fdr_p显著)
table(pair_edc_mp4_s_sig$exp_name)
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

#### 读取clade name数据 ----
## 读取clade name ##
clade_name_mp4_s <- read.table("raw_data/JD.mp4.n4491_clade_name.txt", header = TRUE)
#### 读取clade name数据 ####


############################################# 柱状图 #############################################
#### barchart (EDC vs. MP4 + MP4 vs. Outcome) ----
# 作图只使用EDC score和检出率>50%单个EDC (14)
colnames(pair_edc_mp4_s)[c(3,4,5)] <- c("estimate","p","p_adj_bh")
pair_edc_mp4_s <- pair_edc_mp4_s[pair_edc_mp4_s$exp_name %in% c("edc_score_edc14_f","edc_score_pfas_f","edc_score_pae6_f","edc_score_tc_f","edc_score_bp1_f",
                                                                edc_traits6_log10),]

# 所有EDC和菌(门)关联 (作图数据处理)
{
  # 保留与14 EDC score显著相关菌
  pair_edc_mp4_s_mp4_14edcscore <- pair_edc_mp4_s[pair_edc_mp4_s$exp_name == "edc_score_edc14_f",]
  pair_edc_mp4_s_mp4_14edcscore_sig <- pair_edc_mp4_s_mp4_14edcscore[pair_edc_mp4_s_mp4_14edcscore$p_adj_bh < 0.05,]
  pair_edc_mp4_s_mp4_14edcscore_sig$out_name_new <- gsub("_log10","",pair_edc_mp4_s_mp4_14edcscore_sig$out_name)
  pair_edc_mp4_s_mp4_14edcscore_sig <- left_join(pair_edc_mp4_s_mp4_14edcscore_sig, clade_name_mp4_s, by=c("out_name_new"="species"))
}
dat_edc_mp4_for_fig <- pair_edc_mp4_s_mp4_14edcscore_sig[,c("phylum","out_name","exp_name","estimate")]
colnames(dat_edc_mp4_for_fig) <- c("phylum","species","edc","estimate1")

# 菌(门)和Outcome关联 (作图数据处理)
{
  p_spearman_results_mp4_out_short <- p_spearman_results_mp4_out_short[p_spearman_results_mp4_out_short$exp_name %in% pair_edc_mp4_s_mp4_14edcscore_sig$out_name,]
  p_spearman_results_mp4_out_short <- p_spearman_results_mp4_out_short[p_spearman_results_mp4_out_short$p_adj_bh < 0.05,]
  p_spearman_results_mp4_out_short$exp_name_new <- gsub("_log10","",p_spearman_results_mp4_out_short$exp_name)
  p_spearman_results_mp4_out_short <- left_join(p_spearman_results_mp4_out_short, clade_name_mp4_s, by=c("exp_name_new"="species"))
}
dat_mp4_out_for_fig <- p_spearman_results_mp4_out_short[,c("phylum","exp_name","out_name","estimate")]
colnames(dat_mp4_out_for_fig) <- c("phylum","species","biomarker","estimate2")





dat_for_fig <- left_join(dat_mp4_out_for_fig, dat_edc_mp4_for_fig, by=c("phylum","species"))
unique(dat_for_fig$phylum)

dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("bmi_f","height_f","weight_f","whr_f","wc_f","hc_f"), "Body measurement", "")
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("tg_f","ldl_f","hdl_f","chol_f","apoa_f","apob_f","nonhdl_f"), "Lipid and lipoprotein", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("alt_f","ast_f","ggt_f","bia_f"), "Liver function", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("egfr_f","scr_f","ua_f"), "Kidney function", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("glu0_f","glu120_f","vhba1c_f"), "Glucose metabolism", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("ins0_f","ins120_f","homair_f","homab_f"), "Insulin metabolism", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("sbp_f","dbp_f","pr_f"), "Blood pressure", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("ft3_f","ft4_f","tsh_f","tpoab_f","tgab_f"), "Thyroid function", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("hgb_f","plt_f","eos_f","lym_f","mon_f","neu_f"), "Hematological trait", dat_for_fig$biomarker_group)
dat_for_fig$biomarker_group <- ifelse(dat_for_fig$biomarker %in% c("nlr_f","lmr_f","plr_f","sii_f","siri_f","wbc_f","crp_f"), "Inflammation", dat_for_fig$biomarker_group)

dat_for_fig$direction_group <- ifelse(dat_for_fig$estimate1 > 0 & dat_for_fig$estimate2 > 0, 1,
                                      ifelse(dat_for_fig$estimate1 < 0 & dat_for_fig$estimate2 < 0, 2, 3)) 

dat_for_fig <- dat_for_fig %>%
  group_by(phylum, biomarker_group, direction_group) %>%
  mutate(n=n()) %>%
  distinct(phylum, biomarker_group, direction_group, .keep_all = TRUE) %>%
  ungroup()

all_biomarker_cat1 <- data.frame(phylum = c(rep("p__Actinobacteria",10), rep("p__Firmicutes",10), rep("p__Bacteroidetes",10), rep("p__Proteobacteria",10)),
                                 biomarker_group = rep(c("Body measurement","Lipid and lipoprotein",
                                                         "Liver function","Kidney function",
                                                         "Glucose metabolism","Insulin metabolism","Blood pressure",
                                                         "Thyroid function","Hematological trait","Inflammation"),4),
                                 direction_group = 1)
all_biomarker_cat2 <- data.frame(phylum = c(rep("p__Actinobacteria",10), rep("p__Firmicutes",10), rep("p__Bacteroidetes",10), rep("p__Proteobacteria",10)),
                                 biomarker_group = rep(c("Body measurement","Lipid and lipoprotein",
                                                         "Liver function","Kidney function",
                                                         "Glucose metabolism","Insulin metabolism","Blood pressure",
                                                         "Thyroid function","Hematological trait","Inflammation"),4),
                                 direction_group = 2)

all_biomarker_cat <- rbind(all_biomarker_cat1, all_biomarker_cat2)

dat_for_fig <- left_join(all_biomarker_cat, dat_for_fig, by=c("phylum","biomarker_group","direction_group"))
dat_for_fig$n <- ifelse(is.na(dat_for_fig$n), 0, dat_for_fig$n)

dat_for_fig$phylum <- factor(dat_for_fig$phylum,
                             levels = rev(c("p__Firmicutes","p__Bacteroidetes","p__Actinobacteria","p__Proteobacteria")),
                             labels = rev(c("Firmicutes","Bacteroidetes","Actinobacteria","Proteobacteria")))
dat_for_fig$biomarker_group <- factor(dat_for_fig$biomarker_group,
                                      levels = rev(c("Body measurement","Lipid and lipoprotein",
                                                     "Liver function","Kidney function",
                                                     "Glucose metabolism","Insulin metabolism","Blood pressure",
                                                     "Thyroid function","Hematological trait","Inflammation")))

dat_for_fig <- dat_for_fig %>%
  group_by(phylum, direction_group) %>%
  mutate(n_total=sum(n),
         percent=n/n_total) %>%
  ungroup() %>%
  arrange(phylum, direction_group, biomarker_group)
dat_for_fig$percent <- ifelse(is.na(dat_for_fig$percent), 0, dat_for_fig$percent)
dat_for_fig$percent <- ifelse(dat_for_fig$direction_group == 2, -dat_for_fig$percent, dat_for_fig$percent)
unique(dat_for_fig$n_total)


# 设置11种sector颜色
sector_colors <- colorRampPalette(colors = c("#b71540", "#eb2f06", "#fa8231", "#fed330", "#26de81", "#45aaf2", "#cd84f1", "#7158e2"))(11)

f <- ggplot(dat_for_fig, aes(x =  percent, y =phylum, fill = biomarker_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  
  # X=0垂直参考线
  geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
  
  # 自定义biomarker_group颜色，把values替换成你自己的色号
  scale_fill_manual(
    values = rev(sector_colors[2:11])
  ) +
  
  # x轴放到顶部，关闭底部坐标轴
  scale_x_continuous(
    expand = c(0.05, 0.05),
    limits = c(-0.4, 0.21),   # 坐标轴范围
    breaks = c(-0.4, -0.2, 0, 0.2),
    labels = c("40%", "20%", "0", "20%")
  ) +
  
  # 图例逆序，不改变柱子绘制顺序
  guides(fill = guide_legend(reverse = TRUE)) +
  
  labs(x="Percent (%)", y="Phylum", fill="Biomarker category") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 12),
    
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    
    panel.grid = element_blank()
  )

ggsave(f, filename=paste0("figures/main_figures/(fig4c)_barchart_summary_edc_mp4_out.pdf"), width = 7, height = 7, limitsize = FALSE)

