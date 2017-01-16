## Project Name: NationalMappingForestAttributes
## Authors: Giona Matasci (giona.matasci@gmail.com), Geordie Hoabart (ghobart@nrcan.gc.ca), Harold Zald (hsz16@humboldt.edu)       
## File Name: prog5_ecozone_stats_plots.R                           
## Objective: Produce boxplots based on the descriptive stats extracted with 'prog5_ecozone_stats.py'

#### TO DO -------------------------------------------------------------------

## STILL TO DO:
# Prior to actual run:
# - change stats.by.eco.dir to proper one

## SOLVED:
# -

#### READS/WRITES ------------------------------------------------------------

## READS:
# - 

## WRITES:
# - 

#### INIT --------------------------------------------------------------------

print('Prog4c: boxplots by ecozone') 

rm(list=ls())

param_file <- "D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/wkg/AllUTMzones_paramsGL.Rdata"
load(param_file)

param_file_prog2 = file.path(base_wkg_dir, 'AllUTMzones_params2.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog2)

param_file_prog3 = file.path(base_wkg_dir, 'AllUTMzones_params3.Rdata', fsep = .Platform$file.sep) 
load(param_file_prog3)

source("D:/Research/ANALYSES/NationalMappingForestAttributes/WKG_DIR_NationalMappingForestAttributes/code/Functions_NatMapping_R.R")


#### SCRIPT SPECIFIC PARAMETERS ---------------------------------------------

params4c <- list()

params4c$stats.vect <- c('mean', '10%', '25%', '50%', '75%', '90%')

# params4c$targ.names.lg <- list("elev_mean", "elev_stddev", "elev_cv", "elev_p95", "percentage_first_returns_above_2m", "percentage_first_returns_above_mean", "loreys_height", "basal_area", "gross_stem_volume", "total_biomass")
params4c$targ.names.lg <- list("gross_stem_volume", "total_biomass")

params4c$ylims <- list(c(0, 135), c(0, 270))

params4c$fill.color.vect <- c('chocolate1', 'olivedrab1')
params4c$mapped.ecozones <- list("Boreal Cordillera", "Boreal Plains", "Boreal Shield East", "Boreal Shield West", "Hudson Plains",
                                 "Taiga Cordillera", "Taiga Plains", "Taiga Shield East", "Taiga Shield West")
# params4c$mapped.ecozones <- list("Boreal Cordillera", "Boreal Plains") 

stats.by.eco.dir <- file.path(base_results_dir, "StatsByEcozone")
# stats.by.eco.dir <- file.path(base_results_dir, "StatsByEcozone_1stRun_Still_w_DEM_zeros_errors")


#### LOAD PACKAGES ----------------------------------------------------------

list.of.packages <- c("ggplot2",
                      "ggthemes",
                      "grid",
                      "gridExtra",
                      "GGally",
                      "psych",
                      "plyr",
                      "dplyr",   ## to be loaded before foreach to avoid "assertion failed" errors
                      "vegan",
                      "lubridate", 
                      "data.table",
                      "RColorBrewer",
                      "gplots"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]   ## named vector members whose name is "Package"
if(length(new.packages)) install.packages(new.packages)
for (pack in list.of.packages){
  library(pack, character.only=TRUE)
}

#### START ------------------------------------------------------------------

tic <- proc.time() ## start clocking global time

## subdirectory to save model assessment plots
StatsByEcozone.subdir <- file.path(base_figures_dir, "StatsByEcozone", fsep = .Platform$file.sep)
if (! file.exists(StatsByEcozone.subdir)){dir.create(StatsByEcozone.subdir, showWarnings = F, recursive = T)}

#### READ DATA --------------------------------------------------------------

idx.plt <- 1
plot.list <- list()
for (targ in params4c$targ.names.lg) {
  
  idx.targ <- which(params3$targ.names.lg==targ)

  boxpl.values.df <- data.frame(Ecozone= rep('', length(params4c$mapped.ecozones)), mean=rep(NA, length(params4c$mapped.ecozones)), c10=rep(NA, length(params4c$mapped.ecozones)), c25=rep(NA, length(params4c$mapped.ecozones)), 
                                c50=rep(NA, length(params4c$mapped.ecozones)), c75=rep(NA, length(params4c$mapped.ecozones)), c90=rep(NA, length(params4c$mapped.ecozones)), x=seq(1, length(params4c$mapped.ecozones)),
                                stringsAsFactors=FALSE)
  
  for (z in 1:length(params4c$mapped.ecozones)) {
    zone.name <- params4c$mapped.ecozones[[z]]
    desc.stat.table <- fread(file.path(stats.by.eco.dir, sprintf('%s_stats.csv', zone.name)), header=T)
    cmd <- sprintf('boxpl.values.df[z, seq(2, ncol(boxpl.values.df)-1)] <- desc.stat.table[V1 %%in%% params4c$stats.vect]$%s', targ)
    eval(parse(text=cmd))
    boxpl.values.df[z,1] <- zone.name
  }
  
  if (targ == params4c$targ.names.lg[length(params4c$targ.names.lg)]) {   ## if plot is the last one set specific values
    labels.x <- unlist(params4c$mapped.ecozones)
    margins.vect <- c(0,0,0,0)
  } else {
    labels.x <- rep('', length(params4c$mapped.ecozones))
    margins.vect <- c(0,0,-0.5,0)
  }
  
  plot.list[[idx.plt]] <- ggplot(boxpl.values.df, aes(x = as.factor(x))) + 
    geom_boxplot(aes(ymin=c10, lower=c25, middle=c50, upper=c75, ymax=c90), stat = "identity", fill=params4c$fill.color.vect[idx.plt]) +
    geom_point(aes(x=x, y=mean), color='black', size=1, shape=1, stroke=2) +
    labs(y = sprintf("%s [%s]", params3$targ.names.plots[idx.targ], params3$targ.units.plots[idx.targ]), x = '') +
    scale_x_discrete(labels=labels.x) +
    coord_cartesian(ylim = params4c$ylims[[idx.plt]]) +
    theme(panel.background = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA),
          panel.grid.major = element_line(colour="black", size=0.3, linetype="dashed"),
          panel.grid.minor = element_line(colour="black", size=0.3, linetype="dashed"),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin=unit(margins.vect, "cm"))
  
  idx.plt <- idx.plt+1
  
}

fig.name.str <- file.path(StatsByEcozone.subdir, sprintf("Boxplot_estimates_ECO.pdf"), sep='')
pdf(fig.name.str)   ## , width=6, height=4
  do.call("grid.arrange", c(plot.list, nrow=2, ncol=1))
dev.off()
  
# gr.idx <- 1
# for (pred.gr.names in params3$groups.to.plot) {
#   gr.name <- unlist(params3$groups.to.plot[gr.idx])
#   gr.idx <- gr.idx+1
#   
#     pred.names <- "Ch_attr"
#     X.data.to.plot <- as.data.frame(X.trn.val[, pred.names])  ## subset to keep only variable of interest
#     colnames(X.data.to.plot) <- pred.names
#     
#     Y.data.to.plot.Ch_attr <- Y.trn.val[, unlist(params3$targ.key.names.lg)]
#     colnames(Y.data.to.plot.Ch_attr) <- params3$targ.key.names.sh
#     
#     df.to.plot <- cbind(X.data.to.plot, Y.data.to.plot.Ch_attr)
#     
#     ## Canadian level (merged UTM zones) boxplots by Change_attr (categorical var)
#     plot1 <- ggplot(df.to.plot, aes(x=Ch_attr, y=p_1r_2m, fill=Ch_attr)) +         ## canopy cover above 2m
#       geom_boxplot(notch=F, outlier.shape = NA) + 
#       scale_x_discrete(labels=params3$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))
#     plot2 <- ggplot(df.to.plot, aes(x=Ch_attr, y=tot_biom, fill=Ch_attr)) +          ## biomass
#       geom_boxplot(notch=F, outlier.shape = NA) + 
#       scale_x_discrete(labels=params3$Ch_attr.classes[1+as.integer(levels(df.to.plot$Ch_attr))]) + 
#       theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#       coord_cartesian(ylim = c(0, 250000))
#     str <- file.path(CAN.subdir, sprintf("CAN_boxplots_keyTarg_VS_%s.pdf", gr.name), sep='')
#     pdf(str, width=10, height=5)
#     grid.arrange(plot1, plot2, nrow=1, ncol=2, top=sprintf("%s VS pct_1ret_ab_2m and tot_biom", gr.name))  ## arrange the two plots in same two subplots
#     dev.off()
# 
# }