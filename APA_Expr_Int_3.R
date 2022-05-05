## This script integrates expression and APA data

## @knitr LibPrep

rm(list = ls())
library('tidyverse')
library('ggrepel')
library('scales')

res.out <- 'APA_Expr_Int'

## save(dp397, expr397,
##      dp364, expr364,
load(file = './Comp_Data.RData')

geneFeat <- read.table('./Type_GRCh38.txt') ## u3
colnames(geneFeat) <- c('ENS', 'Symbol', 'Feat')
head(geneFeat)

## @knitr VennDiagram
library('VennDiagram')
library('RColorBrewer')

OverlOut <- function(List, outFi){
    out <- List
    outList <- list(out1=setdiff(out[[1]], out[[2]]),
                    out2=setdiff(out[[2]], out[[1]]),
                    inter=intersect(out[[1]], out[[2]]))
    names(outList) <- c(names(out)[1], names(out)[2], 'intersect')
    outDF <- unrag(outList) %>% apply(., 2, unlist)
    
    write.csv(outDF,
              paste0(outFi, '.csv'),
              row.names=FALSE,
              quote=FALSE,
              na='')
}

myCol <- brewer.pal(3, "Dark2")

dp397 %<>% mutate(Symbol2 = ifelse(PDUI_Group > 0, paste0(Symbol, '+'), paste0(Symbol, '-')))
dp364 %<>% mutate(Symbol2 = ifelse(PDUI_Group > 0, paste0(Symbol, '+'), paste0(Symbol, '-')))

## Symbol was replaced by Symbol2!!!
curOv <- list(JTE607=subset(dp397, abs(PDUI_Group) > 0.1 & P_val < 0.05)$Symbol2 %>% unique,
              sh1CPSF3=subset(dp364, abs(PDUI_Group) > 0.1 & P_val < 0.05)$Symbol2 %>% unique)
OverlOut(curOv,
         outFi = paste0('./Venn_JTE607_vs_DMSO_AbsPDUI_0p1_pval0p05'))
venn.diagram(x = curOv,
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('./Venn_AbsPDUI_0p1_pval0p05.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty = 1,
             ## lty='blank',
             cex = 1.3,
             main='APA analysis comparison\n|PDUI| > 0.1 and pval < 0.05',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))

curOv <- list(JTE607=subset(dp397, PDUI_Group > 0.1 & P_val < 0.05)$Symbol %>% unique,
              sh1CPSF3=subset(dp364, PDUI_Group > 0.1 & P_val < 0.05)$Symbol %>% unique)
OverlOut(curOv,
         outFi = paste0('./Venn_PDUI_0p1_pval0p05_Length'))

venn.diagram(x = curOv,
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('./Venn_PDUI_0p1_pval0p05_Length.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty=1,
             cex = 1.3,
             main='APA analysis comparison\nPDUI > 0.1 & pval < 0.05',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))

curOv <- list(JTE607=subset(dp397, PDUI_Group < -0.1 & P_val < 0.05)$Symbol %>% unique,
              sh1CPSF3=subset(dp364, PDUI_Group < -0.1 & P_val < 0.05)$Symbol %>% unique)
OverlOut(curOv,
         outFi = paste0('./Venn_JTE607_vs_DMSO_pval0p05_Short'))

venn.diagram(x = curOv,
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('./Venn_PDUI_0p1_pval0p05_Short.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty=1,
             cex = 1.3,
             main='APA analysis comparison\nPDUI < -0.1 & pval < 0.05',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))

## Integration matching
#########################

tempCols <- c(OverExpr_Len = 'royalblue', UnderExpr_Len = 'cyan',
              OverExpr_Short = 'red3', UnderExpr_Short = 'indianred1',
              NS = 'ivory2')


## Expr Sign: 1.2 Fold & padj < 0.05; PDUI: abs(0.1)
##            ########

Mer397 <- merge(dp397, expr397)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                             PDUI_Group > 0.1 & P_val < 0.05,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 & P_val < 0.05,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group > 0.1 & P_val < 0.05,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 & P_val < 0.05,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, diff_cut0p1, DE)
dim(Mer397)
head(Mer397)

ggOut <- Mer397 %>%
    subset(abs(PDUI_Group) > 0.1) %>%
    ggplot(aes(x=PDUI_Group, y=log2FoldChange, fill = ExprPDUI)) +
    geom_hline(yintercept = 0, color = 'darkgrey') +
    geom_vline(xintercept = 0, color = 'darkgrey') +
    geom_point(shape=21, alpha = 0.7, size = 2) +
    scale_fill_manual(values = tempCols)+
    labs(x='Change in PDUI (JTE607 - DMSO)',
         y='Expression log2(JTE607 / DMSO)',
         ## title = 'Difference PDUI and Expression\n0.1 cut-off',
         title = '',
         fill = '') +
    theme_bw() +
    theme(legend.position = c(0.85, 0.15),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 3, shape = 21))) 
ggsave(paste0('./JTE607_vs_DMSO_dPDUI_0p1_pval0p05_DEgenes_1p2FC.pdf'),
       ggOut, 
       width=8, height=7)
write.csv(Mer397 %>% subset(ExprPDUI != 'NS') %>% arrange(ExprPDUI),
         paste0('./JTE607_vs_DMSO_dPDUI_0p1_pval0p05_DEgenes_1p2FC.csv'),
         quote = FALSE)
head(Mer397)
write.csv(with(Mer397 %>% arrange(ExprPDUI),
               table(ExprPDUI) %>% data.frame),
         paste0('./JTE607_vs_DMSO_dPDUI_0p1_pval0p05_DEgenes_1p2FC_table.csv'),
         quote = FALSE)
table(subset(Mer397, abs(PDUI_Group) > 0.1)$ExprPDUI)

### CPSF3 Experiment

## Expr Sign: 1.2 Fold & padj < 0.05; PDUI: abs(0.1)
Mer364 <- merge(dp364, expr364)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                             PDUI_Group > 0.1 & P_val < 0.05,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 & P_val < 0.05,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group > 0.1 & P_val < 0.05,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 & P_val < 0.05,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, diff_cut0p1, DE)

ggOut <- Mer364 %>%
    subset(abs(PDUI_Group) > 0.1) %>%
    ggplot(aes(x=PDUI_Group, y=log2FoldChange, fill = ExprPDUI)) +
    geom_hline(yintercept = 0, color = 'darkgrey') +
    geom_vline(xintercept = 0, color = 'darkgrey') +
    geom_point(shape=21, alpha = 0.8, size = 2) +
    scale_fill_manual(values = tempCols)+
    labs(x='Change in PDUI (sh1CPSF3 - shNTC)',
         y='Expression log2(sh1CPSF3 / shNTC)',
         title = '',
         fill = '') +
    theme_bw() +
    theme(legend.position = c(0.85, 0.15),
          legend.title = element_blank(),
          legend.background = element_blank()) +
    guides(fill = guide_legend(override.aes = list(size = 3, shape = 21))) 
ggsave(paste0('./sh1CPSF3_vs_shNTC_dPDUI_0p1_pval0p05_DEgenes_1p2FC.pdf'),
       ggOut, 
       width=8, height=7)
write.csv(Mer364 %>% subset(ExprPDUI != 'NS'),
     paste0('./sh1CPSF3_vs_shNTC_dPDUI_0p1_pval0p05_DEgenes_1p2FC.csv'),
         quote = FALSE)
write.csv(with(Mer364 %>% arrange(ExprPDUI),
               table(ExprPDUI) %>% data.frame),
          paste0('./sh1CPSF3_vs_shNTC_dPDUI_0p1_pval0p05_DEgenes_1p2FC_table.csv'),
         quote = FALSE)

table(subset(Mer364, abs(PDUI_Group) > 0.1)$ExprPDUI)
Mer364 %>% head

IntFin364 <- Mer364 %>% subset(ExprPDUI != 'NS') %>% .$Symbol
IntFin397 <- Mer397 %>% subset(ExprPDUI != 'NS') %>% .$Symbol

venn.diagram(x = list(IntFin364,
                      IntFin397),
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('/Venn_Integr_PDUI_0p1_pval0p05.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty=1,
             cex = 1.3,
             main='APA: abs(PDUI) > 0.1 & pval < 0.05 \n Expr: |FC| > 2 & padj < 0.05\nGene Symbols',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))

##### Without pval < 0.05 and |FC| > 1.2

Mer364 <- merge(dp364, expr364)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                             PDUI_Group > 0.1 ,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 ,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group > 0.1,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, -diff_cut0p1, -DE)


Mer397 <- merge(dp397, expr397)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                             PDUI_Group > 0.1,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group > 0.1 ,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.2) & padj < 0.05 &
                         PDUI_Group < -0.1 ,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, -diff_cut0p1, -DE)

IntFin364 <- Mer364 %>% subset(ExprPDUI != 'NS') %>% .$Symbol %>% unique
IntFin397 <- Mer397 %>% subset(ExprPDUI != 'NS') %>% .$Symbol %>% unique

curOv <- list(JTE607 = IntFin397,
              sh1CPSF3 = IntFin364)
OverlOut(curOv,
         outFi = paste0('./Venn_Integr_PDUI_0p1_1p2FC'))
venn.diagram(x = curOv,
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('../output/',
                    res.out, '/Venn_Integr_PDUI_0p1_1p2FC.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty=1,
             cex = 1.3,
             main='APA: abs(PDUI) > 0.1 \n Expr: |FC| > 1.2 & padj < 0.05\nGene Symbols',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))


##### Without pval < 0.05 and |FC| > 1.5

Mer364 <- merge(dp364, expr364)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.5) & padj < 0.05 &
                             PDUI_Group > 0.1 ,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.5) & padj < 0.05 &
                         PDUI_Group < -0.1 ,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.5) & padj < 0.05 &
                         PDUI_Group > 0.1,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.5) & padj < 0.05 &
                         PDUI_Group < -0.1,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, -diff_cut0p1, -DE)


Mer397 <- merge(dp397, expr397)  %>%
    ## subset(abs(log2FoldChange) > log2(1.5) & padj < 0.05) %>%
    mutate(ExprPDUI = ifelse(log2FoldChange > log2(1.5) & padj < 0.05 &
                             PDUI_Group > 0.1,
                         'OverExpr_Len',
                  ifelse(log2FoldChange > log2(1.5) & padj < 0.05 &
                         PDUI_Group < -0.1,
                         'OverExpr_Short',
                  ifelse(log2FoldChange < -log2(1.5) & padj < 0.05 &
                         PDUI_Group > 0.1 ,
                         'UnderExpr_Len',
                  ifelse(log2FoldChange < -log2(1.5) & padj < 0.05 &
                         PDUI_Group < -0.1 ,
                         'UnderExpr_Short',
                         'NS')
                  ) ) )
           ) %>%
    select(-Pass_Filter, -diff_cut0p1, -DE)

IntFin364 <- Mer364 %>% subset(ExprPDUI != 'NS') %>% .$Symbol %>% unique
IntFin397 <- Mer397 %>% subset(ExprPDUI != 'NS') %>% .$Symbol %>% unique

curOv <- list(JTE607 = IntFin397,
              sh1CPSF3 = IntFin364)
OverlOut(curOv,
         outFi = paste0('./Venn_Integr_PDUI_0p1_1p5FC'))
venn.diagram(x = curOv,
             category.names = c('JTE607\nvs\nDMSO', 'sh1CPSF3\nvs\nshNTC'),
             paste0('../output/',
                    res.out, '/Venn_Integr_PDUI_0p1_1p5FC.png'),
             alpha=0.3,
             fill=c('forestgreen', 'orange'),
             col = 'navy',
             lwd=1,
             lty=1,
             cex = 1.3,
             main='APA: abs(PDUI) > 0.1 \n Expr: |FC| > 1.5 & padj < 0.05\nGene Symbols',
             imagetype='png',
             cat.default.pos = 'outer',
             cat.pos=c(-30, 0),
             cat.cex = 1.3)
file.remove(list.files(paste0('../output/', res.out),
                       '.*log$', full.names=TRUE))
