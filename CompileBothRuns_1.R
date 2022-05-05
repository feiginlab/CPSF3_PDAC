## @knitr LibPrep

rm(list = ls())
library(tidyverse)

## Group1 is control and Group2 is experimental !!
## diff_PDUI should be PDUI_contr - PDUI_exp
## PDUI_Group has cleaner interpretation:
## Large pos values mean lengthening of experimental cond,
## negative are shortened in the experimental condition

## @knitr LoadDataJTE607

(fiIn <- list.files(path='./Results_All_Prediction_JTE607_vs_DMSO_GRCh38v22.txt',
                    full.names=TRUE))
dp397 <- read.table(fiIn, header=TRUE) %>%
    subset(!is.na(PDUI_Group_diff)) %>%
    subset(select = grep(pattern = '^[AB]_[0-9]+_.*',
                         x = colnames(.),
                         value=TRUE, invert = TRUE)) %>%
    mutate(PDUI_Group = - PDUI_Group_diff,  ## my own label!!
           APA_Label = ifelse(PDUI_Group > 0, 'Lengthend', 'Shortend'), 
           diff_cut0p1 = ifelse(abs(PDUI_Group) > 0.1,
                                '>abs(0.1_PDUI_diff)', '<abs(0.1_PDUI_diff)'),

           Symbol = mystrsplit(Gene, '\\|', 2),
           strand = mystrsplit(Gene, '\\|', 4),
           chr = mystrsplit(Loci, ':', 1),

           start = mystrsplit(Loci, ':', 2) %>%
           mystrsplit(delim = '-', sel = 1) %>% an,

           end = mystrsplit(Loci, ':', 2) %>%
           mystrsplit(delim = '-', sel = 2) %>% an) %>%
    select(chr, start, end, strand, Symbol,
           Group_A_Mean_PDUI, Group_B_Mean_PDUI,
           PDUI_Group_diff, P_val, adjusted.P_val, Pass_Filter,
           PDUI_Group, APA_Label, diff_cut0p1, PDUI_Group,
           Gene, fit_value, Predicted_Proximal_APA)

## Retrieve Expression
expr397 <- read.csv('./DEcounts_JTE607_vs_DMSO.csv',
                   check.names=FALSE, row.names=1) %>%
    subset(select = colnames(.) %>%
               grep(pattern='raw.*|norm.*', value=TRUE, invert=TRUE)) %>%
    mutate(DE = ifelse(log2FoldChange > log2(2) & padj < 0.05, 'Up',
                ifelse(log2FoldChange < -log2(2) & padj < 0.05, 'Down',
                       'NS')) %>%
               factor(levels = c('Up', 'Down', 'NS')) ) %>%
    dplyr::rename(Symbol = Gene) 

## subset(padj < 0.05 & abs(log2FoldChange) > log2(2))
head(expr397)
dim(expr397)

#########################################################################

## @knitr LoadData_siCPSF3
(fiIn <- list.files(
    path='Results_All_Prediction_JTE607_vs_DMSO_GRCh38v22.txt', full.names=TRUE))
dp364 <- read.table(fiIn, header=TRUE) %>%
    subset(!is.na(PDUI_Group_diff)) %>%
    subset(select = grep(pattern = '^[AB]_[0-9]+_.*',
                         x = colnames(.),
                         value=TRUE, invert = TRUE)) %>%
    mutate(PDUI_Group = - PDUI_Group_diff,  ## my Label: delta Exp - Control
           APA_Label = ifelse(PDUI_Group > 0, 'Lengthend', 'Shortend'),
           diff_cut0p1 = ifelse(abs(PDUI_Group) > 0.1,
                                '>abs(0.1_PDUI_diff)',
                                '<abs(0.1_PDUI_diff)'),

           Symbol = mystrsplit(Gene, '\\|', 2),
           strand = mystrsplit(Gene, '\\|', 4),
           chr = mystrsplit(Loci, ':', 1),

           start = mystrsplit(Loci, ':', 2) %>%
           mystrsplit(delim = '-', sel = 1) %>% an,

           end = mystrsplit(Loci, ':', 2) %>%
           mystrsplit(delim = '-', sel = 2) %>% an) %>%
    select(chr, start, end, strand, Symbol,
           Group_A_Mean_PDUI, Group_B_Mean_PDUI,
           PDUI_Group_diff, P_val, adjusted.P_val, Pass_Filter,
           PDUI_Group, APA_Label, diff_cut0p1, PDUI_Group,
           Gene, fit_value, Predicted_Proximal_APA)

dp364 %>% head
dp364 %>% dim
with(dp364, table(Pass_Filter))

expr364 <- read.csv(
    './DEcounts_sh1CPSF3_vs_shNTC.csv',
    check.names=FALSE, row.names=1) %>%
    subset(select = colnames(.) %>%
            grep(pattern='raw.*|norm.*', value=TRUE, invert=TRUE)) %>%
    mutate(DE = ifelse(log2FoldChange > log2(2) & padj < 0.05, 'Up',
                        ifelse(log2FoldChange < -log2(2) & padj < 0.05, 'Down',
                               'NS')) %>%
               factor(levels = c('Up', 'Down', 'NS')) ) %>%
    dplyr::rename(Symbol=Gene)  
dim(expr364)
head(expr364)

## @knitr CompileLoad
## subset(padj < 0.05 & abs(log2FoldChange) > log2(2))
dir.create('../output/APA_Expr_Int')

## APA filt: None. diff_cut0p1: 0.1 delta PDUI label
## Expr filt: None. Label: adjpval < 0.05 & |FC| > 2, DE label. 
save(dp397,  expr397,
     dp364,  expr364,
     file = '../output/APA_Expr_Int/Comp_Data.RData')

load(file = '../output/APA_Expr_Int/Comp_Data.RData')
