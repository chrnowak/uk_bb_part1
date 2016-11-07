#######################################################################################################################
##### Function to generate 95% CI-crosshair scatter plots and instrument-strength funnel plots based on ###############
##### summary-level results from MR analysis ##########################################################################

dir <- "~/Desktop/UKBB_project/"
setwd(dir)

library(ggplot2)
library(gridExtra)
 
## load results - res = MR results (for regression lines), res_2 = per SNP associations with risk and outcome

res <- read.delim("MR_results_for_plot.txt") 
res_2 <- read.delim("SUPP_TABLE_SNPS_Scatter.txt")

#### The MR_plot_function generates PDF-graphs in the working directory that feature the crosshair-scatter (left) and funnel plot (right) on a single page
#### the results are taken from the summary level output from the analysis script, so no potentially confidential data is used

MR_plot_function <- function (results, results_2, trait, risk , outcome ){  
  
  ## select risk trait
  temp <- results_2[which(results_2$Trait == trait),]
  
  ## select outcome
  if (outcome == "mortality") {
    beta <- temp$effect_mort            
    se <- temp$se_mort
    outcome <- "Mortality"
    } else if (outcome == "circulatory disease"){
      beta <- temp$effect_circ           
      se <- temp$se_circ
      outcome <- "Circulatory Disease"
      } else if (outcome == "ihd"){
        beta <- temp$effect_ihd           
        se <- temp$se_ihd
        outcome <- "IHD"
        } else if (outcome == "cerebrovascular disease"){
          beta <- temp$effect_cere           
          se <- temp$se_cere
          outcome <- "Cerebrovascular disease"
          } else if (outcome == "isch stroke"){
            beta <- temp$effect_isch           
            se <- temp$se_isch
            outcome <- "Ischaemic Stroke"
            } else if (outcome == "haem stroke"){
              beta <- temp$effect_haem           
              se <- temp$se_haem
              outcome <- "Haemorrhagic stroke"
            }
  
  ## select GRS - with/without obesity-SNPs excluded
  if (risk == "T2D_full"){
    grs_temp <- "t2d_grs_69"
    } else if (risk == "T2D_strict"){
      grs_temp <- "t2d_grs_53"
      }  else if (risk == "FG_full"){
        grs_temp <- "fg_grs_33"
        }  else if (risk == "FG_strict"){
          grs_temp <- "fg_grs_25"
          }  else if (risk == "IR"){
            grs_temp <- "ir_grs"
            }  else if (risk == "IIS"){
              grs_temp <- "iis_grs"
            }  
  
  ## define penalized-robust-IVW and penalized-MR-Egger values for regression lines
  IVW_beta <- results$causal_b[which(results$MR_model == "Penalized robust IVW" & results$GRS == grs_temp & results$Outcome == outcome )]
  Egger_intercept <- results$causal_b[which(results$MR_model == "Penalized MR-Egger (intercept)" & results$GRS == grs_temp  & results$Outcome == outcome)]
  Egger_beta <-results$causal_b[which(results$MR_model == "Penalized MR-Egger"  & results$GRS == grs_temp & results$Outcome == outcome )]
  
  ## format for plot
  df <- data.frame(x = temp$effect_risk, y = beta, 
                   ymin = beta - (1.96*se),ymax = beta + (1.96*se), 
                   xmin = temp$effect_risk - (1.96*(temp$se_risk)), xmax = temp$effect_risk + (1.96*(temp$se_risk)))
  
  ## plot 1: crosshair scatter (genetic associations with outcome v. genetic associations with risk)
  p1 <- ggplot(data = df, aes(x = x,y = y)) + geom_point() + 
    geom_errorbar(aes(ymin = ymin, ymax = ymax)) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    ggtitle(paste(risk, outcome, sep = "_")) + 
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
    labs(x = "Genetic association with risk factor", y = "Genetic association with outcome") +
    geom_abline(slope = Egger_beta, intercept = Egger_intercept, col="red") +
    geom_abline(slope = IVW_beta, intercept = 0, col="blue")  
  
  ## plot 2: funnel of instrument strength v. instrumental variable estimate
  seiv <- (se / temp$effect_risk)
  dfunnel <- data.frame(y = temp$effect_risk/se, x = beta/temp$effect_risk,
                        xmin = (beta/temp$effect_risk) - (1.96*seiv), xmax = (beta/temp$effect_risk) + (1.96*seiv))
  p2 <- ggplot(data = dfunnel,aes(x = x, y = abs(y))) + 
    geom_point(size = 2) + 
    geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = IVW_beta, col = "blue") +
    geom_vline(xintercept = Egger_beta, col = "red") +
    ggtitle(paste(risk, outcome,sep = "_")) + 
    labs(y = "IV strength",x = "IV estimate")  +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) 
  
  ## OUTPUT: both plots combined left/right on one page as PDF graph named as per risk/outcome
  pdf(w = 12, h = 5.5, paste(risk, outcome, "pdf", sep = "."))
  print(grid.arrange(p1, p2, ncol = 2))
  dev.off()
}

## e.g. to plot for the causal effect of T2D on IHD and display as differently-coloured regression lines the
## IVW/Egger estimates for either the full T2D-GRS or the T2D-GRS excluding obesity-SNP do:

MR_plot_function(results,results_2, "T2D", "T2D_full","ihd")

MR_plot_function(results,results_2, "T2D", "T2D_strict","ihd")

########################################################################################################################

