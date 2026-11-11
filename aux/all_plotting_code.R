height <- read.csv("data/h2_results/meta/all_height_rintres_foldenr.csv")
bmi <- read.csv("data/h2_results/meta/all_bmi_rintres_foldenr.csv")
whr <- read.csv("data/h2_results/meta/all_whr_rintres_foldenr.csv")
whradj <- read.csv("data/h2_results/meta/all_WHRadjBMI_rintres_foldenr.csv")
hgb <- read.csv("data/h2_results/meta/all_median_HGB_rintres_foldenr.csv")
t2d <- read.csv("data/h2_results/meta/all_t2d_finalcontrols_foldenr.csv")
hba1c <- read.csv("data/h2_results/meta/all_hba1c_noDM_rintres_v8_foldenr.csv")

height_sepcod <- read.csv("data/h2_results/meta/all_height_rintres_sepcod_foldenr.csv")
bmi_sepcod <- read.csv("data/h2_results/meta/all_bmi_rintres_sepcod_foldenr.csv")
height_sepRNA <- read.csv("data/h2_results/meta/all_height_rintres_sepRNA_foldenr.csv")
hba1c_sepRNA <- read.csv("data/h2_results/meta/all_hba1c_noDM_rintres_v8_sepRNA_foldenr.csv")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpattern)

############################################################################################################
# FIG1: total h2WGS (per ancestry, per trait) + rare vs common h2WGS
############################################################################################################

peranc_total <- ggplot(dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF.bin=="total"), 
       aes(fill=Ancestry,y=h2WGS, x=Trait ) ) +
  geom_bar(position=position_dodge(), width=.7, stat="identity") +
  scale_fill_manual(values = c("#2596be","#7ad239","#a124bd", "#1c3d50"), labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_x_discrete() +
  scale_y_continuous(limits=c(-0.05,0.8), n.breaks=6, oob = scales::squish) +
  geom_errorbar(aes(ymin=h2WGS-1.96*SE, ymax=h2WGS+1.96*SE),
                width=0, linewidth=0.5,position=position_dodge(.7)) +
  labs(
    title = paste("META-analysed total SNP heritability per trait and MAF range"),
    x="",
    y=bquote(paste({h^2}[WGS])) ) +
  theme_linedraw(base_size=18) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        panel.grid.major.y = element_line(colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(hjust=0.5),
        legend.position = "bottom")  
peranc_total

plot_df_rc_peranc_forscatter <- dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF_range=="rare")
plot_df_rc_peranc_forscatter$h2_total <- dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF_range=="total")$h2
colnames(plot_df_rc_peranc_forscatter)[4:6] <- c("h2_rare", "ci_rare", "p_rare")
plot_df_rc_peranc_forscatter$h2_common <- dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF_range=="common")$h2
plot_df_rc_peranc_forscatter$ci_common<- dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF_range=="common")$ci
plot_df_rc_peranc_forscatter$p_common <- dplyr::filter(plot_df_rc_peranc, plot_df_rc_peranc$MAF_range=="common")$p


anc_cols <- c("#2596be","#7ad239","#a124bd", "#1c3d50")

height_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="Height"),
       aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y=bquote(paste({h^2}[R])), x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="left") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank())
bmi_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="BMI"),
                            aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
whr_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="WHR"),
                         aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
whradj_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="WHRadjBMI"),
                         aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x=bquote(paste({h^2}[C]))) +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
hgb_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="HGB"),
                         aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
t2d_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="T2D"),
                         aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())
hba1c_rc_scatter <- ggplot(dplyr::filter(plot_df_rc_peranc_forscatter, plot_df_rc_peranc_forscatter$Trait=="HbA1c"),
                         aes(col=Ancestry,x= h2_common, y=h2_rare ) ) +
  geom_abline(slope=1, col='grey50', linetype=2, linewidth=0.5) +
  geom_errorbar(aes(y=h2_rare, xmin=h2_common-ci_common, xmax=h2_common+ci_common), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_errorbar(aes(x=h2_common, ymin=h2_rare-ci_rare, ymax=h2_rare+ci_rare), 
                width=0, linewidth=0.4, col=anc_cols) +
  geom_point(size=2, pch=16) +
  labs(y="", x="") +
  scale_color_manual(values = anc_cols, labels=c("AFR", "AMR", "EUR", "META"), name="Genetic similarity group") +
  scale_y_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish, position="right") +
  scale_x_continuous(labels=abs, limits=c(0,0.5), n.breaks=2, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  guides(color = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_blank(), axis.title=element_blank(), axis.text.y=element_blank())

ggarrange(peranc_total, ggarrange(plotlist=c(height_rc_scatter,bmi_rc_scatter,whr_rc_scatter, whradj_rc_scatter, hgb_rc_scatter, t2d_rc_scatter, hba1c_rc_scatter), 
                                  widths=c(1, 0.55,0.55,0.55,0.55,0.55,0.55), 
                                  align='h',ncol=7),
          nrow=2, common.legend = T, legend="bottom", heights=c(1, 0.25))

ggsave("plots/Fig1.png", width=11, height=11)


############################################################################################################
# FIG2A: STACKED HERITABILITY (all bins)
############################################################################################################

# make DF of meta results
annot <- c('1.coding rare','2.coding comm', '3.RNA rare', '4.RNA comm', "5.proximal UTR rare", 
           "6.proximal UTR comm", "7.proximal non-UTR rare", "8.proximal non-UTR comm", 
           "9.intronic rare", "10.intronic common", "11.distal regulatory rare", "12.distal regulatory comm",
           "13.un-annotated rare", "14.un-annotated common")
annotation <- c(rep(factor(annot, levels=annot), 7))
ann <- c('coding','coding', 'ncRNA', 'ncRNA', "proximal UTR", 
           "proximal UTR", "proximal non-UTR", "proximal non-UTR", 
           "intronic", "intronic", "distal regulatory", "distal regulatory",
           "un-annotated", "un-annotated")
a <- c(rep(factor(ann, levels=c('coding', 'ncRNA', "proximal UTR", "proximal non-UTR", 
                                "intronic", "distal regulatory","un-annotated")), 7))

meta_stack_h2 <- data.frame(full_annotation=annotation, annotation = a,
                    h2 = c(height$META_h2[1:14], bmi$META_h2[1:14], whr$META_h2[1:14], whradj$META_h2[1:14], hgb$META_h2[1:14], t2d$META_h2[1:14], hba1c$META_h2[1:14]), 
                    se = c(height$META_h2_SE[1:14], bmi$META_h2_SE[1:14], whr$META_h2_SE[1:14], whradj$META_h2_SE[1:14], hgb$META_h2[1:14], t2d$META_h2_SE[1:14], hba1c$META_h2_SE[1:14]), 
                    pval=c(height$META_h2_p[1:14], bmi$META_h2_p[1:14], whr$META_h2_p[1:14], whradj$META_h2_p[1:14], hgb$META_h2_p[1:14], t2d$META_h2_p[1:14], hba1c$META_h2_p[1:14]),
                    trait=factor(c(rep("1.height", 14), rep("2.BMI", 14), rep("3.WHR", 14), rep("4.WHRadjBMI", 14), rep("5.HGB", 14), rep("6.T2D", 14),  rep("7.HbA1c", 14)) )  )
colnames(meta_stack_h2)[1] <- "full_annotation"
rarecomm <- rep(c('1.rare', '2.comm'), 7)
meta_stack_h2$pattern <-  rep(rarecomm, 7)

# add ymin,ymax,ymid for error bars
meta_stack_pos <- meta_stack_h2 %>%
  group_by(trait) %>%
  arrange(trait, desc(full_annotation)) %>%   # ensure consistent stacking order
  mutate(
    ymin = cumsum(h2) - 1.96*se,
    ymax = cumsum(h2) + 1.96*se) %>%
  ungroup()

#meta_stack_pos_sggd <- dplyr::filter(meta_stack_pos, meta_stack_pos$trait != "1.height" | meta_stack_pos$trait != "3.WHR" |meta_stack_pos$trait != "4.WHRadjBMI" |meta_stack_pos$trait != "6.HGB" |meta_stack_pos$trait != "7.HbA1c" )

# make all values <0 small +ve
meta_stack_h2$h2[meta_stack$h2 < 0] <- 0.001
meta_stack_h2$se[meta_stack$se < 0] <- 0.001

hjust <- c(0.5,0.5,0.5,0.8,0.5, 0.5)
h2_meta <- ggplot(meta_stack_h2, aes(fill=annotation, y=h2, x=trait, pattern=pattern)) + 
  geom_bar_pattern(position="stack", stat="identity", colour="grey20", linewidth=0.2,
                    pattern_spacing = 0.02, pattern_frequency = 2, pattern_angle = 45, pattern_fill='white', pattern_color='white', pattern_size=0.01, pattern_alpha=0.5) + 
  scale_fill_manual(
    values = c(
      "coding" = "#fe6200",
      "ncRNA" = "#ffc72d",
      "proximal UTR" = "#2596be",
      "proximal non-UTR" = "#1c3d50",
      "intronic" = "#785ef0",
      "distal regulatory" = "#6ebd33",
      "un-annotated" = "#3d691d"
    ), name = "Functional annotation"
  ) +
  scale_pattern_manual(
    values = c("stripe", "none"),
    name = "MAF range", labels=c("rare", "common")
  ) +
  guides(lines=NULL, fill = guide_legend(
    override.aes = list(pattern = rep('none',7) ))
    ) +
  #geom_errorbar(
    #aes(ymin = ymin, ymax = ymax), 
      #col=rep(rev(c("green4","green4","darkgoldenrod3","darkgoldenrod3","black", "black", 
      #"brown4","brown4", "purple4", "purple4", "deepskyblue4",  "deepskyblue4", 
      #"black", "black")), 5),
    #width = 0, linewidth=0.3, position = position_dodge(width = 0.8)) +
  labs(
    y = (bquote({h^2}[WGS])) ,
    x = "Trait", 
    title = bquote(paste("AoU META-analysed cumulative ", {h^2}[WGS] ," per functional category")) ) +
  scale_x_discrete(
    labels=c("Height", "BMI", "WHR", "WHRadj", "HGB", "T2D", "HbA1c")) + 
  theme_linedraw(base_size = 18) + 
  ylim(0,0.7) +
  #theme(axis.text.x = element_text(angle = 45, hjust=0.5, vjust=0.8) ) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
                                            plot.background = element_blank(), panel.background = element_blank(),
                                            legend.background = element_blank(), 
                                            axis.ticks.x = element_blank()
        ) 
h2_meta

ggsave('plots/stacked_meta_h2_newvep_2.png', width=14, height=7)


############################################################################################################
# FIG2B: ENRICHMENT (all bins)
############################################################################################################

# prep
annot <- c('01.coding rare','02.coding common', '03.ncRNA rare', '04.ncRNA common', "05.proximal UTR rare", 
           "06.proximal UTR common", "07.proximal non-UTR rare", "08.proximal non-UTR common", 
           "09.intronic rare", "10.intronic common", "11.distal regulatory rare", "12.distal regulatory common",
           "13.un-annotated rare", "14.un-annotated common")
annotation <- c(rep(factor(annot, levels=annot), 7))

ann <- c('coding','coding', 'ncRNA', 'ncRNA', "proximal\nUTR",
         "proximal\nUTR", "proximal\nnon-UTR", "proximal\nnon-UTR",
         "intronic", "intronic", "distal\nregulatory", "distal\nregulatory",
         "non-\nannotated", "un-\nannotated")


meta_stack <- data.frame(Annotation=annotation,
                         enr = c(height$META_enr[1:14], bmi$META_enr[1:14], whr$META_enr[1:14], whradj$META_enr[1:14], hgb$META_enr[1:14], t2d$META_enr[1:14], hba1c$META_enr[1:14]),
                         enr_se = c(height$META_enr_SE[1:14], bmi$META_enr_SE[1:14], whr$META_enr_SE[1:14], whradj$META_enr_SE[1:14], hgb$META_enr_SE[1:14], t2d$META_enr_SE[1:14], hba1c$META_enr_SE[1:14]),
                         enr_pval = c(height$META_enr_p[1:14], bmi$META_enr_p[1:14], whr$META_enr_p[1:14], whradj$META_enr_p[1:14], hgb$META_enr_p[1:14], t2d$META_enr_p[1:14], hba1c$META_enr_p[1:14]),
                         enr_Q_P = c(height$META_enr_Q_p[1:14], bmi$META_enr_Q_p[1:14], whr$META_enr_Q_p[1:14], whradj$META_enr_Q_p[1:14], hgb$META_enr_Q_p[1:14], t2d$META_enr_Q_p[1:14], hba1c$META_enr_Q_p[1:14]),
                         h2 = c(height$META_h2[1:14], bmi$META_h2[1:14], whr$META_h2[1:14], whradj$META_h2[1:14], hgb$META_h2[1:14], t2d$META_h2[1:14], hba1c$META_h2[1:14]),
                         h2_se = c(height$META_h2_SE[1:14], bmi$META_h2_SE[1:14], whr$META_h2_SE[1:14], whradj$META_h2_SE[1:14], hgb$META_h2_SE[1:14], t2d$META_h2_SE[1:14], hba1c$META_h2_SE[1:14]),
                         h2_pval = c(height$META_h2_p[1:14], bmi$META_h2_p[1:14], whr$META_h2_p[1:14], whradj$META_h2_p[1:14], hgb$META_h2_p[1:14], t2d$META_h2_p[1:14], hba1c$META_h2_p[1:14]),
                         Trait = factor(c(rep("Height", 14), rep("BMI", 14), rep("WHR", 14), rep("WHRadjBMI", 14), rep("HGB", 14), rep("T2D", 14), rep("HbA1c", 14)), 
                                      levels=c("Height", "BMI", "WHR", "WHRadjBMI", "HGB", "T2D", "HbA1c"))
                         )


# add proper annotation labels
ann <- c('coding','coding', 'ncRNA', 'ncRNA', "proximal\nUTR",
         "proximal\nUTR", "proximal\nnon-UTR", "proximal\nnon-UTR",
         "intronic", "intronic", "distal\nregulatory", "distal\nregulatory",
         "un-\nannotated", "un-\nannotated")
meta_stack$annsmall <- c(rep(factor(ann, levels=c('coding', 'ncRNA', "proximal\nUTR", "proximal\nnon-UTR",
                                                  "intronic", "distal\nregulatory","un-\nannotated")), 7))
# graphical labels
rarecomm <- rep(c('1.rare', '2.comm'), 7)
meta_stack$pattern <-  rep(rarecomm, 7)

# CI calculation
# move SE to log scale
meta_stack$log_SE <- meta_stack$enr_se / meta_stack$enr
meta_stack$enr_low_ci <- exp( log(meta_stack$enr) - 1.96*meta_stack$log_SE )
meta_stack$enr_high_ci <- exp( log(meta_stack$enr) + 1.96*meta_stack$log_SE )

###### PYRAMID PLOT
rare_enr <- dplyr::filter(meta_stack, meta_stack$pattern == "1.rare")
comm_enr <- dplyr::filter(meta_stack, meta_stack$pattern == "2.comm")

annot_cols <- c("#fe6200","#ffc72d","#2596be","#1c3d50", "#785ef0","#6ebd33","#3d691d")
trait_cols <- c("#fe6200",   "#1c3d50","#2596be","#88ccee",  "#3d691d","#6ebd33","#bde99c", "grey60")

rareplot_depl <- ggplot(rare_enr, aes(x = enr, fill = annsmall, y=Trait)) + 
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  labs(
    subtitle = "Rare variants (1e-4≤MAF<0.01)",
    x="fold enrichment",
    y=NULL) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(labels=abs, limits=c(0, 60), n.breaks=6, oob = scales::squish) +
  ggplot2::scale_fill_manual(values = c(annot_cols), name="Functional annotation" ) +
  geom_errorbar(aes(xmin=enr-1.96*enr_se, xmax=enr+1.96*enr_se, alpha=0.5), 
                width=0, linewidth=0.4,
                position=position_dodge(.7, reverse=T)) +
  theme_linedraw(base_size=22) +
  geom_vline(xintercept = 1, col='black', linewidth=0.3) +
  
  guides(fill = guide_legend(nrow=1, ncol=8, byrow=F) ) + 
  
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(linetype=5, colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(), 
        legend.text = element_text(size = 22), legend.title = element_text(size = 24),
        legend.position = "bottom", 
        axis.text.y = element_text(hjust=0.5,  margin = margin(r = 14.5)))


commplot_depl <- ggplot(comm_enr, aes(x = -enr, fill = annsmall, y=Trait)) + 
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  labs(
    subtitle = "Common variants (0.01≤MAF<0.5)",
    x="fold enrichment",
    y=NULL) +
  scale_y_discrete(limits = rev, labels = NULL) +
  scale_x_continuous(limits=c(-60, 0), labels=abs, n.breaks=6, oob = scales::squish) +
  ggplot2::scale_fill_manual(values = c(annot_cols), name="Functional annotation") +
  geom_errorbar(aes(xmin=-(enr-1.96*enr_se), xmax=-(enr+1.96*enr_se)), 
                width=0, linewidth=0.4,
                position=position_dodge(.7, reverse=T)) +
  geom_vline(xintercept = -1, col='black', linewidth=0.3) +
  theme_linedraw(base_size=22) +
  guides(fill = guide_legend(nrow=1, ncol=8, byrow=F) ) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_line(linetype=5, colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),legend.box="vertical",
        legend.text = element_text(size = 22), legend.title = element_text(size = 24),
        axis.ticks.y = element_blank(),
        legend.position = "bottom", 
        axis.text.y = element_text(hjust=0.5, margin = margin(r = 14.5)))


ggarrange(commplot_depl, rareplot_depl, ncol=2,nrow=1, common.legend=T, legend="bottom", widths=c(1, 1.18))

ggsave('plots/commrare.png', width=16, height=14)


############################################################################################################
# FIG3A, B: sepcod + sepRNA pyramids
############################################################################################################
anc_cols <- c("#2596be","#7ad239","#a124bd", "#1c3d50")


#################################################### sepcod #################################################
height_sepcod <- read.csv("data/h2_results/meta/all_height_rintres_sepcod_gem_metafix.csv")
bmi_sepcod <- read.csv("data/h2_results/meta/all_bmi_rintres_sepcod_gem_metafix.csv")
whr_sepcod <- read.csv("data/h2_results/meta/all_whr_rintres_sepcod_gem_metafix.csv")
whradj_sepcod <- read.csv("data/h2_results/meta/all_WHRadjBMI_rintres_sepcod_gem_metafix.csv")
hgb_sepcod <- read.csv("data/h2_results/meta/all_median_hgb_rintres_sepcod_gem_metafix.csv")
t2d_sepcod <- read.csv("data/h2_results/meta/all_t2d_finalcontrols_sepcod_gem_metafix.csv")
hba1c_sepcod <- read.csv("data/h2_results/meta/all_hba1c_noDM_rintres_v8_sepcod_gem_metafix.csv")

annots_sepcod <- factor(height_sepcod$annotation )
meta_stack_sepcod <- data.frame(Annotation= rep(annots_sepcod[1:8], 7),
                        h2 = c(height_sepcod$META_h2[1:8], bmi_sepcod$META_h2[1:8], whr_sepcod$META_h2[1:8], whradj_sepcod$META_h2[1:8], 
                                hgb_sepcod$META_h2[1:8], t2d_sepcod$META_h2[1:8], hba1c_sepcod$META_h2[1:8]),
                        h2_se = c(height_sepcod$META_h2_SE[1:8], bmi_sepcod$META_h2_SE[1:8], whr_sepcod$META_h2_SE[1:8], whradj_sepcod$META_h2_SE[1:8], 
                                  hgb_sepcod$META_h2_SE[1:8], t2d_sepcod$META_h2_SE[1:8], hba1c_sepcod$META_h2_SE[1:8]),
                        h2_pval = c(height_sepcod$META_h2_p[1:8], bmi_sepcod$META_h2_p[1:8], whr_sepcod$META_h2_p[1:8], whradj_sepcod$META_h2_p[1:8], 
                                  hgb_sepcod$META_h2_p[1:8], t2d_sepcod$META_h2_p[1:8], hba1c_sepcod$META_h2_p[1:8]),
                        enr = c(height_sepcod$META_enr[1:8], bmi_sepcod$META_enr[1:8], whr_sepcod$META_enr[1:8], whradj_sepcod$META_enr[1:8], 
                                  hgb_sepcod$META_enr[1:8], t2d_sepcod$META_enr[1:8], hba1c_sepcod$META_enr[1:8]),
                        enr_se = c(height_sepcod$META_enr_SE[1:8], bmi_sepcod$META_enr_SE[1:8], whr_sepcod$META_enr_SE[1:8], whradj_sepcod$META_enr_SE[1:8], 
                                  hgb_sepcod$META_enr_SE[1:8], t2d_sepcod$META_enr_SE[1:8], hba1c_sepcod$META_enr_SE[1:8]),
                        pval = c(height_sepcod$META_enr_p[1:8],
                                  bmi_sepcod$META_enr_p[1:8], whr_sepcod$META_enr_p[1:8], whradj_sepcod$META_enr_p[1:8], hgb_sepcod$META_enr_p[1:8], t2d_sepcod$META_enr_p[1:8], hba1c_sepcod$META_enr_p[1:8]),
                        Q_P = c(height_sepcod$META_enr_Q_p[1:8], bmi_sepcod$META_enr_Q_p[1:8], whr_sepcod$META_enr_Q_p[1:8], whradj_sepcod$META_enr_Q_p[1:8], 
                                hgb_sepcod$META_enr_Q_p[1:8], t2d_sepcod$META_enr_Q_p[1:8], hba1c_sepcod$META_enr_Q_p[1:8]),
                        Trait = factor(c(rep("Height", 8), rep("BMI", 8), rep("WHR", 8), rep("WHRadjBMI", 8), 
                                         rep("HGB", 8), rep("T2D", 8), rep("HbA1c", 8)), 
                                        levels=c("Height", "BMI", "WHR", "WHRadjBMI", "HGB", "T2D", "HbA1c"))
                      )

meta_stack_sepcod$annsmall <- c(rep(factor(c('rare\npLoF','common\npLoF', 'rare\nsplicing', 'common\nsplicing', "rare\nmissense", "common\nmissense", "rare\nsynonymous", "common\nsynonymous"), 
                                           levels=c('rare\npLoF','common\npLoF', 'rare\nsplicing', 'common\nsplicing', "rare\nmissense", "common\nmissense", "rare\nsynonymous", "common\nsynonymous")), 7))


all_sepcod_pyramid <- ggplot(dplyr::filter(meta_stack_sepcod, meta_stack_sepcod$ancestry=="META"), 
                             aes(x = log10(enr), fill = Trait, y=annsmall)) +
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  labs(
    title = "Meta-analysis across genetic similarity groups: coding",
    x="log10(fold enrichment)",
    y=NULL) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits=c(0, 3.5), n.breaks=6, oob = scales::squish) +
  scale_fill_manual(values = c("#68c0f2", "#0072B2", "#08bf8d","#096e52", "#f0d50a", "darkorange", "darkorange4"), name="Genetic similarity group" ) +
  geom_errorbar(aes(xmin=(enr-1.96*enr_se), xmax=(enr+1.96*enr_se)),
                width=0, linewidth=0.3,
                position=position_dodge(.7, reverse=T)) + theme_linedraw(base_size=18) +
  geom_vline(xintercept = 0, col='black', linewidth=0.3) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(linetype=5, colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2))

#################################################### sepRNA #################################################
height_sepRNA <- read.csv("data/h2_results/meta/all_height_rintres_sepRNA_gem_metafix.csv")
bmi_sepRNA <- read.csv("data/h2_results/meta/all_bmi_rintres_sepRNA_gem_metafix.csv")
whr_sepRNA <- read.csv("data/h2_results/meta/all_whr_rintres_sepRNA_gem_metafix.csv")
whradj_sepRNA <- read.csv("data/h2_results/meta/all_WHRadjBMI_rintres_sepRNA_gem_metafix.csv")
hgb_sepRNA <- read.csv("data/h2_results/meta/all_median_hgb_rintres_sepRNA_gem_metafix.csv")
hba1c_sepRNA <- read.csv("data/h2_results/meta/all_hba1c_noDM_rintres_v8_sepRNA_gem_metafix.csv")

annots_sepRNA <- factor( height_sepRNA$annotation )
meta_stack_sepRNA <- data.frame(Annotation= rep(annots_sepRNA[13:16], 6),
                                enr = c(height_sepRNA$META_enr[13:16], 
                                        bmi_sepRNA$META_enr[13:16],  whr_sepRNA$META_enr[13:16],  whradj_sepRNA$META_enr[13:16],  hgb_sepRNA$META_enr[13:16], hba1c_sepRNA$META_enr[13:16] ),
                                enr_se = c(height_sepRNA$META_enr_SE[13:16], 
                                           bmi_sepRNA$META_enr_SE[13:16],  whr_sepRNA$META_enr_SE[13:16],  whradj_sepRNA$META_enr_SE[13:16],  hgb_sepRNA$META_enr_SE[13:16], hba1c_sepRNA$META_enr_SE[13:16] ),
                                pval = c(height_sepRNA$META_enr_p[13:16],
                                         bmi_sepRNA$META_enr_p[13:16],  whr_sepRNA$META_enr_p[13:16],  whradj_sepRNA$META_enr_p[13:16],  hgb_sepRNA$META_enr_p[13:16], hba1c_sepRNA$META_enr_p[13:16] ),
                                Q_P = c(  height_sepRNA$META_enr_Q_p[13:16],
                                         bmi_sepRNA$META_enr_Q_p[13:16], whr_sepRNA$META_enr_Q_p[13:16],  whradj_sepRNA$META_enr_Q_p[13:16],  hgb_sepRNA$META_enr_Q_p[13:16], hba1c_sepRNA$META_enr_Q_p[13:16] ),
                                ancestry = rep("META", 24), 
                                Trait = factor(c(rep("Height", 4), rep("BMI", 4), rep("WHR", 4), rep("WHRadjBMI", 4), rep("HGB", 4), rep("HbA1c", 4)), 
                                               levels=c("Height", "BMI", "WHR", "WHRadjBMI", "HGB", "HbA1c"))
                      )
meta_stack_sepRNA$annsmall <- c(rep(factor(c('rare\nlncRNA','common\nlncRNA', 'rare\nother ncRNA', 'common\nother ncRNA'), 
                                           levels=c('rare\nlncRNA','common\nlncRNA', 'rare\nother ncRNA', 'common\nother ncRNA')), 6))


all_sepRNA_pyramid <- ggplot(dplyr::filter(meta_stack_sepRNA, meta_stack_sepRNA$ancestry=="META"), 
                               aes(x = log10(enr), fill = Trait, y=annsmall)) +
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  labs(
    title = "Meta-analysis across genetic similarity groups: ncRNA",
    x="log10(fold enrichment)",
    y=NULL) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(limits=c(0, 3.5), n.breaks=6, oob = scales::squish) +
  scale_fill_manual(values = c("#68c0f2", "#0072B2", "#08bf8d","#096e52", "#f0d50a", "darkorange4")) +
  geom_errorbar(aes(xmin=(enr-1.96*enr_se), xmax=(enr+1.96*enr_se)),
                width=0, linewidth=0.3,
                position=position_dodge(.7, reverse=T)) + theme_linedraw(base_size=18) +
  geom_vline(xintercept = log10(1), col='black', linewidth=0.3) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(linetype=5, colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow=2, override.aes = list(size = 3)))
all_sepRNA_pyramid


ggarrange(plotlist=c(all_sepcod_pyramid, all_sepRNA_pyramid), 
          common.legend=T, legend="bottom",
          ncol=2,nrow=2, 
          widths=c(0.78,1), heights= c(1, 0.6))

ggsave('plots/Fig3AB_allsep.png', width=8, height=16)


############################################################################################################
# FIG3C, D: SIMULATIONS
############################################################################################################

### h2 simulation plot
sepcod_og <- read.csv("data/h2_results/sims/all_height_rintres_sepcod_foldenr.csv")
sepcod_og <- sepcod_og[1:12,1:7]
sepcod_og_full <- read.csv("data/h2_results/meta/afr_height_rintres_mph_newvep2_alpha_sepcod2.mq.vc.csv")
sepcod_og_full <- data.frame(annotation = c("0182.comm_syn_lo", "0183.comm_syn_hi"),
                             AFR_h2 = sepcod_og_full[15:16,7],
                             AFR_h2_SE = sepcod_og_full[15:16,8],
                             AFR_h2_p = c(0, 0),
                             AFR_enr_sum = sepcod_og_full[15:16,9],
                             AFR_enr_SE = sepcod_og_full[15:16,10],
                             AFR_enr_pval = c(0, 0)
)
sepcod_og_full$AFR_h2_p <- 2*pnorm(sepcod_og_full$AFR_h2/sepcod_og_full$AFR_h2_SE, lower.tail=F)
sepcod_og_full$AFR_enr_pval <- 2*pnorm(sepcod_og_full$AFR_enr_sum/sepcod_og_full$AFR_enr_SE, lower.tail=F)
sepcod_og_full <- rbind(sepcod_og[1:8,], sepcod_og_full)
sepcod_og_full <- rbind(sepcod_og_full, sepcod_og[9:12,], make.row.names = F)

sepcod_sim <- read.csv("data/h2_results/sims/all_sim2_height_rintres_seed.csv")
sepcod_sim <- sepcod_sim[1:14,]
sepcod_sim_nosyn <- read.csv("data/all_sim_nosyn_height_rintres_seed.csv")
sepcod_sim_nosyn <- sepcod_sim_nosyn[1:14,]
sepcod_sim_norna <- read.csv("data/all_sim2_norna_height_rintres_seed.csv")
sepcod_sim_norna <- sepcod_sim_norna[1:14,]
sepcod_sim_nomiss <- read.csv("data/all_sim_nomiss_seed.csv")
sepcod_sim_nomiss <- sepcod_sim_norna[1:14,]

annots_sepcod <- factor( c("011.rare_lof", "012.comm_lof", 
                           "013.rare_splice", "014.comm_splice",
                           "015.rare_miss", "016.comm_miss",
                           "017.rare_syn", "018.comm_syn",
                           "0182.comm_syn_lo", "0183.comm_syn_hi",
                           "019.rare_codoth", "020.comm_codoth",
                           "03.rare_RNA", "04.comm_RNA"))
annotation_sepcod <- c(rep(factor(annots_sepcod, levels=annots_sepcod), 4))
annots_sepcod_small <- c(rep(factor(c("rare lof", "common\n lof", 
                                      "rare\n splicing", "common\n splicing",
                                      "rare\n missense", "common\n missense",
                                      "rare\n synonymous", "common\n synonymous",
                                      "common\n synonymous\n (low LD)", "common\n synonymous\n (high LD)",
                                      "rare\n other", "common\n other",
                                      "rare ncRNA", "common\n ncRNA"), levels=c("rare lof", "common\n lof", 
                                                                                "rare\n splicing", "common\n splicing",
                                                                                "rare\n missense", "common\n missense",
                                                                                "rare\n synonymous", "common\n synonymous",
                                                                                "common\n synonymous\n (low LD)", "common\n synonymous\n (high LD)",
                                                                                "rare\n other", "common\n other",
                                                                                "rare ncRNA", "common\n ncRNA") ),5))

meta_stack_sim <- data.frame(annotsmall = annots_sepcod_small, 
                             h2 = c(sepcod_og_full$AFR_h2, sepcod_sim$meta_h2, sepcod_sim_nosyn$meta_h2, sepcod_sim_nomiss$meta_h2, sepcod_sim_norna$meta_h2),
                             h2_se = c(sepcod_og_full$AFR_h2_SE, sepcod_sim$meta_h2_se, sepcod_sim_nosyn$meta_h2_se, sepcod_sim_nomiss$meta_h2_se, sepcod_sim_norna$meta_h2_se),
                             h2_pval = c(sepcod_og_full$AFR_h2_p, sepcod_sim$meta_h2_p, sepcod_sim_nosyn$meta_h2_p, sepcod_sim_nomiss$meta_h2_p, sepcod_sim_norna$meta_h2_p),
                             enr = c(sepcod_og_full$AFR_enr_sum, sepcod_sim$meta_enr, sepcod_sim_nosyn$meta_enr, sepcod_sim_nomiss$meta_enr, sepcod_sim_norna$meta_enr), 
                             enr_se = c(sepcod_og_full$AFR_enr_SE, sepcod_sim$meta_enr_se, sepcod_sim_nosyn$meta_enr_se, sepcod_sim_nomiss$meta_enr_se, sepcod_sim_norna$meta_enr_se), 
                             enr_pval= c(sepcod_og_full$AFR_enr_pval, sepcod_sim$meta_enr_p, sepcod_sim_nosyn$meta_enr_p, sepcod_sim_nomiss$meta_enr_p, sepcod_sim_norna$meta_enr_p),
                             test=factor(c(rep("observed", 14), rep("sim:\nreplication", 14), rep("sim: no\nsynonymous", 14), rep("sim: no\nmissense", 14), rep("sim:\nno ncRNA", 14)), levels=c("observed", "sim:\nreplication", "sim: no\nsynonymous","sim: no\nmissense", "sim:\nno ncRNA")) 
)


####################################### plot scatter for simulations ###################################
sepcod_og <- read.csv("data/h2_results/sims/all_height_rintres_sepcod_gem_metafix.csv")
sepcod_sim <- read.csv("data/h2_results/sims/all_sim2_height_rintres_seed_foldenr_pluslog_gem.csv")

sim_reps <- data.frame(simall_h2 = sepcod_sim$meta_h2[c(1:8,11:24)],
                       empirical_h2 = sepcod_og$AFR_h2[1:22],
                       simall_se =  sepcod_sim$meta_h2[c(1:8,11:24)],
                       empirical_se = sepcod_og$AFR_h2_SE[1:22])
regression <- lm(sim_reps$simall_h2 ~ sim_reps$empirical_h2)

sim_regr_plot <- ggplot(sim_reps, aes(x=empirical_h2, y=simall_h2)) +
  geom_abline(intercept = regression$coefficients[1], slope = regression$coefficients[2], 
              color="grey50", linetype="dashed") +
  geom_errorbar(aes(y=simall_h2, xmin=empirical_h2-1.96*empirical_se, xmax=empirical_h2+1.96*empirical_se), 
                width=0, linewidth=0.2, col='grey50') +
  geom_errorbar(aes(x=empirical_h2, ymin=simall_h2-1.96*simall_se, ymax=simall_h2+1.96*simall_se), 
                width=0, linewidth=0.2, col='grey50') +
  geom_abline(slope=1, col='red') +
  geom_point(size=2, col='black') +
  labs(
    title = "Simulation heritability replication",
    x= bquote(paste("observed ", {h^2}[WGS])),
    y=bquote(paste("simulated ", {h^2}[WGS])) ) +
  scale_x_continuous(limits=c(-0.025,0.3), oob = scales::squish) +
  scale_y_continuous(limits=c(-0.025,0.3), oob = scales::squish) +
  theme_linedraw(base_size=18) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank())
sim_regr_plot

############################################## simulation fig ##########################################
meta_stack_sim <- dplyr::filter(meta_stack_sim, meta_stack_sim$annotsmall!= "rare\n other" & meta_stack_sim$annotsmall!= "common\n synonymous" & meta_stack_sim$annotsmall!= "common\n other")

sim_bar_plot <- ggplot(meta_stack_sim, aes(fill=test,  x=h2, y=reorder(annotsmall, desc(annotsmall)) ) ) +
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  scale_fill_manual(values = c("#88ccee",  "#1c3d50", "#fe6200", "brown4", "#ffc72d"), name="Analysis type") +
  guides(fill=guide_legend(nrow=2)) +
  geom_errorbar(aes(xmin=h2-1.96*h2_se, xmax=h2+1.96*h2_se),
                width=0, linewidth=0.3, position=position_dodge(.7, reverse=T)) +
  labs(
    title = bquote(paste("AFR height ", {h^2}[WGS] ," simulations")),
    y="",
    x=bquote({h^2}[WGS])) +
  scale_x_continuous(limits=c(-0.01,0.05), n.breaks=6, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y = element_blank()
  )
sim_bar_plot

ggsave('plots/sepcod_sims_enrichments_bar.png', width=8.5, height=11)

ggarrange(sim_regr_plot,sim_bar_plot, common.legend=F, ncol=1,nrow=2, heights=c(0.6,1), align="v")
ggsave('plots/Fig3CD_ggplot.png', width=8, height=15)


############################################### build fig 3 ###############################################

ggarrange(
  ggarrange(plotlist=c(all_sepcod_pyramid, all_sepRNA_pyramid), 
            common.legend=T, legend="bottom",
            ncol=2,nrow=2, 
            widths=c(0.73,1), heights= c(1, 0.6)),
  ggarrange(sim_regr_plot, sim_bar_plot, common.legend=F, ncol=1,nrow=2, heights=c(0.58,1), align="v")
)

ggsave('plots/Fig3.png', width=16, height=18)


############################################################################################################
# SUPPLEMENTARY 1: overlaps
############################################################################################################
library(corrplot)

anc_rc_overlaps <- read.csv("data/h2_results/ancestry_rc_overlaps.csv")
anc_rc_overlaps_X <- read.csv("data/h2_results/ancestry_rc_overlaps_X.csv")
anc_rc_overlaps[,2:7] <- anc_rc_overlaps[,2:7] + anc_rc_overlaps_X[,2:7]

rownames(anc_rc_overlaps) <- anc_rc_overlaps$anc
colnames(anc_rc_overlaps) <- c("", anc_rc_overlaps$anc)
anc_rc_overlaps <- anc_rc_overlaps[,2:7]


anc_rc_percentages <- anc_rc_overlaps
anc_rc_percentages[1,] <- anc_rc_percentages[1,]/anc_rc_percentages[1,1]
anc_rc_percentages[2,] <- anc_rc_percentages[2,]/anc_rc_percentages[2,2]
anc_rc_percentages[3,] <- anc_rc_percentages[3,]/anc_rc_percentages[3,3]
anc_rc_percentages[4,] <- anc_rc_percentages[4,]/anc_rc_percentages[4,4]
anc_rc_percentages[5,] <- anc_rc_percentages[5,]/anc_rc_percentages[5,5]
anc_rc_percentages[6,] <- anc_rc_percentages[6,]/anc_rc_percentages[6,6]

rownames(anc_rc_percentages) <- c("fraction of\nAFR rare", "fraction of\nAMR rare", "fraction of\nEUR rare",  "fraction of\nAFR common", "fraction of\nAMR common", "fraction of\nEUR common")

par(mfrow=c(1,2))
corrplot(as.matrix(anc_rc_overlaps), type="lower", method="circle", tl.col='black', bg="grey80", is.corr=F, col.lim=c(0,70000000), cl.length=8,tl.cex=0.7, cl.cex=0.6)
corrplot(as.matrix(anc_rc_percentages), method="shade", col.lim=c(0,1), tl.col='black', tl.cex=0.7, cl.cex=0.6)

############################################################################################################
# SUPPLEMENTARY 2: rare, common, total (per ancestry, per trait)
############################################################################################################
anc_cols <- c("#2596be","#7ad239","#a124bd", "#1c3d50")

ggplot(plot_df_rc_peranc, aes(fill=ancestry,y=h2, x=trait ) ) +
  geom_bar(position="dodge", width=.7, stat="identity") +
  scale_fill_manual(values = c("#88ccee", "#1c3d50", "#2596be", 
                               "#bde99c", "#3d691d", "#7ad239", 
                               "pink2", "#43086e", "#a124bd",
                               "grey", "grey10", "grey40"), name="MAF range") +
  scale_y_continuous(limits=c(-0.12,0.8), n.breaks=6, oob = scales::squish) +
  geom_errorbar(aes(ymin=h2-ci, ymax=h2+ci),
                width=0, linewidth=0.5, position=position_dodge(.7)) +
  labs(
    x="trait",
    y="h2WGS") +
  theme_linedraw(base_size=18) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.y = element_line(colour="grey30"),
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom",
        axis.text.y = element_text(hjust=0.5))

ggsave('plots/supp2.png', width=14, height=14)


############################################################################################################
# SUPPLEMENTARY 3: scatter rare and common (per ancestry, per trait)
############################################################################################################

idx_rare <- seq(1,14,2)
idx_comm <- seq(2,14,2)

annot_cols <- c("#fe6200", "#ffc72d", "#2596be", "#1c3d50", "#785ef0", "#6ebd33", "#3d691d")

ann_forplot <- factor(c("coding", "RNA exons", "proximal UTR", "proximal not UTR", 
                 "intronic", "distal regulatory", "un-annotated"), levels=c("coding", "RNA exons", "proximal UTR", "proximal not UTR", 
                                                                            "intronic", "distal regulatory", "un-annotated"))

supp3_master_df <- data.frame(
  trait = factor(c(rep("Height",14), rep("BMI", 14), rep("WHR", 14),rep("WHRadjBMI", 14), rep("HGB", 14), rep("T2D", 14), rep("HbA1c", 14)), 
                 levels = c("Height", "BMI", "WHR", "WHRadjBMI", "HGB", "T2D", "HbA1c") ), 
  annot = ann_forplot, 
  AFR_rare_h2 = c(height$AFR_h2[idx_rare], bmi$AFR_h2[idx_rare], whr$AFR_h2[idx_rare], whradj$AFR_h2[idx_rare], hgb$AFR_h2[idx_rare], t2d$AFR_h2[idx_rare], rep(NA,7) ),
  AFR_rare_h2_SE = c(height$AFR_h2_SE[idx_rare], bmi$AFR_h2_SE[idx_rare], whr$AFR_h2_SE[idx_rare], whradj$AFR_h2_SE[idx_rare], hgb$AFR_h2_SE[idx_rare], t2d$AFR_h2_SE[idx_rare], rep(NA,7) ),
  AFR_comm_h2 = c(height$AFR_h2[idx_comm], bmi$AFR_h2[idx_comm], whr$AFR_h2[idx_comm], whradj$AFR_h2[idx_comm], hgb$AFR_h2[idx_comm], t2d$AFR_h2[idx_comm], rep(NA,7) ),
  AFR_comm_h2_SE = c(height$AFR_h2_SE[idx_comm], bmi$AFR_h2_SE[idx_comm], whr$AFR_h2_SE[idx_comm], whradj$AFR_h2_SE[idx_comm], hgb$AFR_h2_SE[idx_comm], t2d$AFR_h2_SE[idx_comm], rep(NA,7) ),
  
  AMR_rare_h2 = c(height$AMR_h2_SE[idx_comm], bmi$AMR_h2_SE[idx_comm], whr$AMR_h2_SE[idx_comm], whradj$AMR_h2_SE[idx_comm], hgb$AMR_h2_SE[idx_comm], t2d$AMR_h2_SE[idx_comm], hba1c$AMR_h2_SE[idx_comm] ),
  AMR_rare_h2_SE = c(height$AMR_h2_SE[idx_rare], bmi$AMR_h2_SE[idx_rare], whr$AMR_h2_SE[idx_rare], whradj$AMR_h2_SE[idx_rare], hgb$AMR_h2_SE[idx_rare], t2d$AMR_h2_SE[idx_rare], hba1c$AMR_h2_SE[idx_rare] ),
  AMR_comm_h2 = c(height$AMR_h2[idx_comm], bmi$AMR_h2[idx_comm], whr$AMR_h2[idx_comm], whradj$AMR_h2[idx_comm], hgb$AMR_h2[idx_comm], t2d$AMR_h2[idx_comm], hba1c$AMR_h2[idx_comm] ),
  AMR_comm_h2_SE = c(height$AMR_h2_SE[idx_comm], bmi$AMR_h2_SE[idx_comm], whr$AMR_h2_SE[idx_comm], whradj$AMR_h2_SE[idx_comm], hgb$AMR_h2_SE[idx_comm], t2d$AMR_h2_SE[idx_comm], hba1c$AMR_h2_SE[idx_comm] ),
  
  EUR_rare_h2 = c(height$EURf_h2[idx_rare], bmi$EURf_h2[idx_rare], whr$EURf_h2[idx_rare], whradj$EURf_h2[idx_rare], hgb$EURf_h2[idx_rare], t2d$EURf_h2[idx_rare], hba1c$EURf_h2[idx_rare] ),
  EUR_rare_h2_SE = c(height$EURf_h2_SE[idx_rare], bmi$EURf_h2_SE[idx_rare], whr$EURf_h2_SE[idx_rare], whradj$EURf_h2_SE[idx_rare], hgb$EURf_h2_SE[idx_rare], t2d$EURf_h2_SE[idx_rare], hba1c$EURf_h2_SE[idx_rare] ),
  EUR_comm_h2 = c(height$EURf_h2[idx_comm], bmi$EURf_h2[idx_comm], whr$EURf_h2[idx_comm], whradj$EURf_h2[idx_comm], hgb$EURf_h2[idx_comm], t2d$EURf_h2[idx_comm], hba1c$EURf_h2[idx_comm] ),
  EUR_comm_h2_SE = c(height$EURf_h2_SE[idx_comm], bmi$EURf_h2_SE[idx_comm], whr$EURf_h2_SE[idx_comm], whradj$EURf_h2_SE[idx_comm], hgb$EURf_h2_SE[idx_comm], t2d$EURf_h2_SE[idx_comm], hba1c$EURf_h2_SE[idx_comm] )
)


plot_rc_peranc_func <- function(anc, trait) {
  df <- supp3_master_df[supp3_master_df$trait==trait,]
  ggplot(df, aes(col=annot, x=df[,paste0(anc, "_comm_h2")], y=df[,paste0(anc, "_rare_h2")] ) ) +
  geom_abline(slope=1, col='red', linetype=2, linewidth=0.5) +
    
   geom_errorbar(aes(y=.data[[paste0(anc, "_rare_h2")]], xmin=.data[[paste0(anc, "_comm_h2")]] - 1.96*.data[[paste0(anc, "_comm_h2_SE")]], xmax=.data[[paste0(anc, "_comm_h2")]] + 1.96*.data[[paste0(anc, "_comm_h2_SE")]]), 
                 width=0, linewidth=0.4, col='grey70') +
   geom_errorbar(aes(x=.data[[paste0(anc, "_comm_h2")]], ymin=.data[[paste0(anc, "_rare_h2")]] - 1.96*.data[[paste0(anc, "_rare_h2_SE")]], ymax=.data[[paste0(anc, "_rare_h2")]] + 1.96*.data[[paste0(anc, "_rare_h2_SE")]]), 
               width=0, linewidth=0.4, col='grey70') +
    
  geom_point(size=2, pch=16) +
  labs(y=bquote(paste("rare ", {h^2}[WGS])), x=bquote(paste("common ", {h^2}[WGS])), title= paste0(anc, " ", trait)) +
  scale_color_manual(values = rep(annot_cols,2), name="Functional annotation") +
  scale_y_continuous(labels=abs, limits=c(-0.05,0.25), n.breaks=3, oob = scales::squish) +
  scale_x_continuous(labels=abs, limits=c(0,0.25), n.breaks=3, oob = scales::squish) +
  guides(color = guide_legend(override.aes = list(size = 3))) + 
  theme_linedraw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10 ))
}

afr_height <- plot_rc_peranc_func("AFR", "Height")
amr_height <- plot_rc_peranc_func("AMR", "Height")
eur_height <- plot_rc_peranc_func("EUR", "Height")

afr_bmi <- plot_rc_peranc_func("AFR", "BMI")
amr_bmi <- plot_rc_peranc_func("AMR", "BMI")
eur_bmi <- plot_rc_peranc_func("EUR", "BMI")

afr_whr <- plot_rc_peranc_func("AFR", "WHR")
amr_whr <- plot_rc_peranc_func("AMR", "WHR")
eur_whr <- plot_rc_peranc_func("EUR", "WHR")

afr_whradj <- plot_rc_peranc_func("AFR", "WHRadjBMI")
amr_whradj <- plot_rc_peranc_func("AMR", "WHRadjBMI")
eur_whradj <- plot_rc_peranc_func("EUR", "WHRadjBMI")

afr_hgb <- plot_rc_peranc_func("AFR", "HGB")
amr_hgb <- plot_rc_peranc_func("AMR", "HGB")
eur_hgb <- plot_rc_peranc_func("EUR", "HGB")

afr_t2d <- plot_rc_peranc_func("AFR", "T2D")
amr_t2d <- plot_rc_peranc_func("AMR", "T2D")
eur_t2d <- plot_rc_peranc_func("EUR", "T2D")

afr_hba1c <- plot_rc_peranc_func("AFR", "HbA1c")
amr_hba1c <- plot_rc_peranc_func("AMR", "HbA1c")
eur_hba1c <- plot_rc_peranc_func("EUR", "HbA1c")


ggarrange(plotlist=c(afr_height, amr_height, eur_height,
                     afr_bmi, amr_bmi, eur_bmi,
                     afr_whr, amr_whr, eur_whr,
                     afr_whradj, amr_whradj, eur_whradj,
                     afr_hgb, amr_hgb, eur_hgb,
                     afr_t2d, amr_t2d, eur_t2d,
                     afr_hba1c, afr_hba1c, afr_hba1c), 
          common.legend=T, legend="bottom",
          ncol=3, nrow=7)

ggsave('plots/supp3_enr_h2_scatter.png', width=10, height=16)



############################################################################################################
# SUPPLEMENTARY 4: h2WGS versus enrichment scatter (all traits)
############################################################################################################
trait_cols <- c("#68c0f2", "#0072B2", "#08bf8d","#096e52", "#f0d50a", "darkorange", "darkorange4")
#trait_cols <- c("#FF9999", "brown1", "#990000","#330000", "#FFCC00", "darkorange", "darkorange4")

supp4_log10pvals_rare <- ggplot(rare_enr, aes(col=Trait, x=-log10(h2_pval), y=-log10(enr_pval)) ) +
  scale_color_manual(values = trait_cols, name="Trait") +
  geom_abline(slope=1, col='grey50', linewidth=0.5, linetype=2) +
  geom_point(size=1) +
  labs(title="", y=bquote(paste("-",log[10], "P (fold enrichment)")), x=bquote(paste("-", log[10], "P (", {h^2}[WGS], ")"))) +
  scale_x_continuous(limits=c(0,85), oob = scales::squish) +
  scale_y_continuous(limits=c(0, 85), oob = scales::squish) +
  theme_linedraw(base_size=16) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank())

supp4_log10pvals_comm <- ggplot(comm_enr, aes(col=Trait, x=-log10(h2_pval), y=-log10(enr_pval)) ) +
  scale_color_manual(values = trait_cols, name="Trait") +
  geom_abline(slope=1, col='grey50', linewidth=0.5, linetype=2) +
  geom_point(size=1) +
  labs(title="", y=bquote(paste("-",log[10], "P (fold enrichment)")), x=bquote(paste("-", log[10], "P (", {h^2}[WGS], ")"))) +
  scale_x_continuous(limits=c(0,85), oob = scales::squish) +
  scale_y_continuous(limits=c(0, 85), oob = scales::squish) +
  theme_linedraw(base_size=16) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank())

supp4_enr_h2_rare <- ggplot(rare_enr, aes(col=Trait, x=h2, y=enr)) +
  scale_color_manual(values = trait_cols, name="Trait") +
  geom_errorbar(aes(y=enr, xmin=h2-h2_se*1.96, xmax=h2+h2_se*1.96),
                width=0, linewidth=0.2, col='grey70') +
  geom_errorbar(aes(x=h2, ymin=enr-1.96*enr_se, ymax=enr+1.96*enr_se),
                width=0, linewidth=0.2, col='grey70') +
  geom_hline(yintercept=0, col='black', linewidth=0.2) +
  geom_vline(xintercept=0, col='black', linewidth=0.2) +
  geom_point(size=1) +
  labs(title="", y="fold enrichment", x=bquote({h^2}[WGS])) +
  scale_x_continuous(limits=c(-0.1,0.3), oob = scales::squish) +
  scale_y_continuous(limits=c(-25,65), oob = scales::squish) +
  theme_linedraw(base_size=16) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank())

supp4_enr_h2_comm <- ggplot(comm_enr, aes(col=Trait, x=h2, y=enr)) +
  scale_color_manual(values = trait_cols, name="Trait") +
  geom_errorbar(aes(y=enr, xmin=h2-h2_se*1.96, xmax=h2+h2_se*1.96),
                width=0, linewidth=0.2, col='grey70') +
  geom_errorbar(aes(x=h2, ymin=enr-1.96*enr_se, ymax=enr+1.96*enr_se),
                width=0, linewidth=0.2, col='grey70') +
  geom_hline(yintercept=0, col='black', linewidth=0.2) +
  geom_vline(xintercept=0, col='black', linewidth=0.2) +
  geom_point(size=1) +
  labs(title="", y="fold enrichment", x=bquote({h^2}[WGS])) +
  scale_x_continuous(limits=c(-0.1,0.3), oob = scales::squish) +
  scale_y_continuous(limits=c(-25,65), oob = scales::squish) +
  theme_linedraw(base_size=16) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank())

ggarrange(plotlist=c(supp4_log10pvals_rare, supp4_enr_h2_rare, supp4_log10pvals_comm, supp4_enr_h2_comm), 
          ncol=2, nrow=2, common.legend = T, legend="bottom")

ggsave("plots/supp4_pval_h2_enr_scatter.png", width=10, height=9)



############################################################################################################
# SUPPLEMENTARY 5: simulations for noncoding h2
############################################################################################################
annot <- c("rare\ncoding", "common\ncoding", "rare ncRNA", "common ncRNA", "rare\nproximal UTR", "common\nproximal UTR", "rare proximal\nnon-UTR", "common proximal\nnon-UTR", "rare intronic", "common intronic", "rare distal\nregulatory", "common distal\nregulatory", "rare non-\nannotated", "common non-\nannotated")

##### simulations for prox and reg
meta_stack_sim_prox_reg_intr <- read.csv("data/h2_results/sims/meta_stack_sim_prox_reg.csv")
meta_stack_sim_prox_reg_intr$test <- factor(c(
  rep("observed", 14), rep("sim:\nreplication", 14), 
  rep("sim: no\nproximal UTR", 14), rep("sim: no proximal\nnon-UTR", 14), 
  rep("sim: no\nintronic", 14), rep("sim: no\nregulatory", 14)),
  levels=c("observed", "sim:\nreplication", "sim: no\nproximal UTR", "sim: no proximal\nnon-UTR", "sim: no\nintronic", "sim: no\nregulatory"))
meta_stack_sim_prox_reg_intr$annsmall <- c(rep(factor(annot, levels=annot), 6))


ggplot(meta_stack_sim_prox_reg_intr, aes(fill=test,  x=h2, y=reorder(annsmall, desc(annsmall)) ) ) +
  geom_bar(position=position_dodge(reverse=T), width=.7, stat="identity") +
  scale_fill_manual(values = c("darkorange","#b22222","#2596be","#1c3d50","#785ef0", "#6ebd33" ), name="Analysis type") +
  guides(fill=guide_legend(nrow=2)) +
  geom_errorbar(aes(xmin=h2-1.96*h2_se, xmax=h2+1.96*h2_se),
                width=0, linewidth=0.3, position=position_dodge(.7, reverse=T)) +
  labs(
    title = "AFR height h2WGS simulations",
    y="",
    x="h2WGS") +
  scale_x_continuous(limits=c(-0.01,0.3), n.breaks=6, oob = scales::squish) +
  theme_linedraw(base_size=18) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        plot.background = element_blank(), panel.background = element_blank(),
        legend.background = element_blank(),
        axis.ticks.y = element_blank()
        #axis.text.x = element_text(angle = 45, hjust=1, vjust=0.8)
  )

ggsave('plots/supp5_annot_sims_h2_bar_newCIs.png', width=10, height=14)
