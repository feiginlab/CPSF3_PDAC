rm(list = ls())
library(tidyverse)
library(ggrepel)
library(scales)

## 397 is JTE607 vs DMSO
## 364 is sh1CPSF3 vs shNTC
## save(dp397,  expr397,
##      dp364,  expr364,
##      file = './Comp_Data.RData')
load(file = './Comp_Data.RData')

## OutPut0 <- read.delim(
##     '../output/DaPars/Results_All_Prediction_Results_GRCh38v22.txt')
## ## OutPut0 <- read.delim('../output/DaPars/Results_GRCh38v38_All_Prediction_Results.txt')

## head(OutPut0)
## dim(OutPut0)

################################################################################
################################################################################

## @knitr Visualization
## Group 1 is Control: Panc1 DMSO
## Group 2 is Experimental: Panc1 JTE607

ggOut <- dp397 %>%
    mutate(Alter = ifelse(PDUI_Group_diff > 0, 'Shortened', 'Lengthened')) %>%
    subset(!is.na(Alter)) %>%
    mutate(Alter = ifelse(abs(PDUI_Group_diff) > 0.1, Alter,  'ns')) %>%
    ggplot(aes(x=Group_A_Mean_PDUI, y=Group_B_Mean_PDUI, fill = Alter)) +
    geom_point(shape=21, alpha =0.6, size=1.7) +
    scale_fill_manual(values =c(Shortened = 'firebrick',
                      Lengthened = 'navy', ns = 'darkgrey'))+
    labs(x='Panc DMSO', y='Panc JTE607', ### Modify accordingly!
         title = 'PDUI Scattered\n0.1 cut-off')+
    theme_bw()+
    geom_abline(intercept = c(-0.1, 0.1), slope=1, size=1.3, linetype='dashed')
ggsave('./JTE607_vs_DMSO_PDUI_Scattered_01.pdf',
       width=8, height=7)

library(scales)
ggOut <- dp397 %>%
    mutate(PDUI_diff = -PDUI_Group_diff, ## Small adjustment
        Alter = ifelse(PDUI_Group_diff > 0, 'Shortened', 'Lengthened')) %>%
    subset(!is.na(Alter)) %>%
    mutate(Alter = ifelse(abs(PDUI_Group_diff) > 0.1, Alter,  'ns')) %>%
    ggplot(aes(x=PDUI_diff, y=-log10(adjusted.P_val), fill = Alter)) +
    geom_point(shape=21, alpha =0.6, size=1.7) +
    ## scale_y_log10()+
    scale_fill_manual(values =c(Shortened = 'firebrick',
                      Lengthened = 'navy', ns = 'darkgrey'))+
    labs(x='Change in PDUI (Panc1 JTE607 - Panc1_DMSO)', y='-log10(FDR)',
         title = 'PDUI Scattered\n0.1 cut-off')+
    theme_bw() +
    geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') +
    scale_y_continuous(trans = log2_trans(),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)))
ggsave('./JTE607_vs_DMSO_Volcano_01.pdf',
       width=8, height=7)

