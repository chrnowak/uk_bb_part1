############# ############# ############# ############# ############# ############# ############# 
############# Generate forest & scatter/funnel plots from summary-level results #################
############# ############# ############# ############# ############# ############# ############# 

dir <- "~/Desktop/UKBB_project/"
setwd(dir)

library(metafor)
library(gridExtra)
library(ggplot2)

results <- read.delim("MR_results_for_plot.txt")  ## MR results from MendelianRandomization  package - all as beta / se - hence for forestplot-OR transform = exp
MV_result <- read.delim("MV_MR_for_plots.txt")    ## MV MR results for T2D and FG

#### forest plot generics ####
names2 <- c("Mortality"," ","Circulatory disease"," ","IHD"," ","Cerebrovascular disease"," ","Ischemic stroke"," ", "Hemorrhagic stroke"," ")
names3 <- c("Mortality"," "," ", "Circulatory disease"," "," ","IHD"," "," ","Cerebrovascular disease"," "," ","Ischemic stroke"," "," ","Hemorrhagic stroke"," ","")

### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### T2D ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######

MV_b <- c(MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Mortality")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Circulatory Disease")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "IHD")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Cerebrovascular disease")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Ischaemic Stroke")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Haemorrhagic stroke")),][c(3,5,7)])
MV_se <- c(MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Mortality")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Circulatory Disease")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "IHD")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Cerebrovascular disease")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Ischaemic Stroke")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "T2D" & MV_result$Outcome == "Haemorrhagic stroke")),][c(4,6,8)])

x2 <- results$causal_b[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "t2d_grs_69"),
                         which(results$MR_model == "Penalized robust IVW" & results$GRS == "t2d_grs_53"))]
x2 <- x2[c(1,7,2,8,3,9,4,10,5,11,6,12)] ## re-order for plot

sei2 <- results$se[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "t2d_grs_69"),
                     which(results$MR_model == "Penalized robust IVW" & results$GRS == "t2d_grs_53"))]
sei2 <- sei2[c(1,7,2,8,3,9,4,10,5,11,6,12)]

## PDF 6.3 by 7.5
## pdf(w = 6.3, h = 7.5, "MR_plot_T2D.pdf")
forest.default(x = c(x2, unlist(MV_b)) , sei = c(sei2, unlist(MV_se)), psize = 1, refline = 1,
               col = c(rep(c("black","gray50"),6), rep("gray20",3), rep("gray40",3),
                       rep("gray20",3), rep("gray40",3), rep("gray20",3), rep("gray40",3)),
               xlim = c(-2.1,3), at = c(0,1,2), transf = exp,
               slab = c(names2, names3), 
               ilab = c(rep(c("broad GRS","strict GRS"),6), rep(c("full","obesity","BP"),6)),
               ilab.xpos = -1, ilab.pos = 4, cex = 0.5,
               ylim = c(-3,44),
               rows = c(39:38, 36:35, 33:32, 30:29, 27:26, 24:23,
                        19:17, 15:13, 11:9, 7:5, 3:1, -1:-3),
               xlab = "Causal OR (95% CI) per OR increase in genetic T2D risk")
text(-1.095,40.75,"Main MR analysis",cex = 0.55,font = 1,pos = 2)
abline(22,0,lty = 2,lwd = 0.5)
text(-1.145,20.5,"Multivariable MR",cex = 0.55,font = 1,pos = 2)
## dev.off()

### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### FG  ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######

x2 <- results$causal_b[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "fg_grs_33"),
                         which(results$MR_model == "Penalized robust IVW" & results$GRS == "fg_grs_25"))]
x2 <- x2[c(1,7,2,8,3,9,4,10,5,11,6,12)]

sei2 <- results$se[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "fg_grs_33"),
                     which(results$MR_model == "Penalized robust IVW" & results$GRS == "fg_grs_25"))]
sei2 <- sei2[c(1,7,2,8,3,9,4,10,5,11,6,12)]

MV_b <- c(MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Mortality")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Circulatory Disease")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "IHD")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Cerebrovascular disease")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Ischaemic Stroke")),][c(3,5,7)],
          MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Haemorrhagic stroke")),][c(3,5,7)])
MV_se <- c(MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Mortality")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Circulatory Disease")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "IHD")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Cerebrovascular disease")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Ischaemic Stroke")),][c(4,6,8)],
           MV_result[c(which(MV_result$trait == "FG" & MV_result$Outcome == "Haemorrhagic stroke")),][c(4,6,8)])

## PDF  9.5 by 7.5 ######
## pdf(w = 9.5, h = 7.5, "MR_plot_FG.pdf")
forest.default(x = c(x2, unlist(MV_b)) , sei = c(sei2, unlist(MV_se)), psize = 1, refline = 1,
               col = c(rep(c("black","gray50"),6), rep("gray20",3), rep("gray40",3),
                       rep("gray20",3), rep("gray40",3), rep("gray20",3), rep("gray40",3)),
               xlim = c(-3.1, 6), at = c(0:5), transf=exp,
               slab = c(names2, names3), 
               ilab = c(rep(c("broad GRS","strict GRS"),6), rep(c("full","obesity","BP"),6)),
               ilab.xpos = -1, ilab.pos = 4,cex = 0.5,
               ylim = c(-3,44),
               rows = c(39:38, 36:35, 33:32, 30:29, 27:26, 24:23,
                        19:17, 15:13, 11:9, 7:5, 3:1, -1:-3),
               xlab =  "Causal OR (95% CI) per mmol/l increase in fasting glucose", cex.lab = 0.55)
text(-1.935,40.75, "Main MR analysis", cex = 0.55, font = 1, pos = 2)
abline(22,0, lty = 2, lwd = 0.5)
text(-2,20.5,"Multivariable MR", cex = 0.55, font = 1, pos = 2)
## dev.off()

####### without MV_MR ######
#### 6.5 by 4.5 #####
forest.default(x = x2 , sei = sei2, psize = 1, refline = 1,
               col = rep(c("black","gray50"),6),
               xlim = c(-3.1, 5.5), at = c(0:4), transf=exp,
               slab = names2,  
               ilab = rep(c("broad GRS","strict GRS"),6),
               ilab.xpos = -1,cilab.pos = 4,ccex = 0.5,
               ylim = c(22,44),
               rows = c(39:38, 36:35, 33:32, 30:29, 27:26, 24:23),
               xlab =  "Causal OR (95% CI) per mmol/l increase in fasting glucose", cex.lab = 0.55)

### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### IR & IIS ## ###### ### ###### ### ###### ### ###### ### ###### ### ######
### ### ###### ### ###### ### ###### ### ### ### ### ### ###### ### ###### ### ###### ### ###### ### ###### ### ######

x2 <- results$causal_b[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "ir_grs"),
                         which(results$MR_model == "Penalized robust IVW" & results$GRS == "iis_grs"))]
x2 <- x2[c(1,7,2,8,3,9,4,10,5,11,6,12)]

sei2 <- results$se[c(which(results$MR_model == "Penalized robust IVW" & results$GRS == "ir_grs"),
                     which(results$MR_model == "Penalized robust IVW" & results$GRS == "iis_grs"))]
sei2 <- sei2[c(1,7,2,8,3,9,4,10,5,11,6,12)]

## PDF 6.00 x 4.50
## pdf(w = 6, h = 4.5, "MR_plot_IR_IIS.pdf")
forest.default(x = x2 , sei = sei2, psize = 1.2, refline = 1,
               col = rep(c("black","gray50"),6),
               xlim = c(-3.1, 5.5), at = c(0,1,2,3,4), transf=exp,
               slab = c(names2), 
               ilab = rep(c("IR","IIS"),6),
               ilab.xpos = -1, ilab.pos = 4, cex = 0.5,
               ylim = c(22,44),
               rows = c(39:38, 36:35, 33:32, 30:29, 27:26, 24:23),
               xlab = "Causal OR (95% CI) per SD-unit increase IR or IIS",cex.lab = 0.55)
# text(-1.113,40.75,"Main MR analysis",cex=0.55,font=1,pos=2)
# abline(22,0,lty=2,lwd=0.5)
# text(-1.158,20.5,"Multivariable MR",cex=0.55,font=1,pos=2)
## dev.off()

#### #### ######## #### ######## #### ######## #### ######## #### ######## #### ######## #### ####

