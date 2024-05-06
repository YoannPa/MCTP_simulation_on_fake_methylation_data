#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Simulate_MCTP_on_fake_methylation_###########################################
Version = '0.0.1'
Date = '2024-05-04'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 4.3.3 (2024-02-29)',
'RStudio 2023.09.1 Build 494 © 2009-2023','!!Add dependencies here!!')
Description = 'script description here'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##IMPORTS
setwd("~/MCTP_simulation_on_fake_methylation_data/src/")
Imports = c("truncnorm", "nparcomp", "data.table", "parallel")
invisible(lapply(Imports, library, character.only = T))
source("random_generation_methylation_values.R")

# Simulate MCTP with Tukey contrast 10,000 times using different asymptotic
# approximation on random generation of bimodal
# methylation distributions
asy_method <- c("fisher", "mult.t", "log.odds", "normal")
ls_bimodal_asy <- lapply(
  X = asy_method, FUN = function(asy){
    ls_res_bimodal <- mclapply(X = seq(10000), mc.cores = 8, FUN = function(i){
      ls_grp <- list(
        data.table("groups" = 1, "betas" = rbimodal_meth(n = 1000)),
        data.table("groups" = 2, "betas" = rbimodal_meth(n = 1000)),
        data.table("groups" = 3, "betas" = rbimodal_meth(n = 1000))
      )
      dt_sim_bimodal <- rbindlist(ls_grp)

      res_mctp <- mctp(
        betas ~ groups, data = dt_sim_bimodal, asy.method = asy, type = "Tukey",
        alternative = "two.sided", plot.simci = FALSE,
        info = FALSE)$Analysis.Inf$p.Value
    })
    unlist(ls_res_bimodal)
  })
names(ls_bimodal_asy) <- asy_method

# Get percentages of p.values below 0.05 for each asymptotic approximation
res_bimodal <- vapply(X = names(ls_bimodal_asy), FUN = function(asy){
  paste0(round(length(ls_bimodal_asy[[asy]][ls_bimodal_asy[[asy]] <= 0.05])/
      length(ls_bimodal_asy[[asy]])*100, 2), "%")
}, FUN.VALUE = character(length = 1L))

# Same with trimodal methylation distributions
ls_trimodal_asy <- lapply(
  X = asy_method, FUN = function(asy){
    ls_res_trimodal <- mclapply(X = seq(10000), mc.cores = 10, FUN = function(i){
      ls_grp <- list(
        data.table("groups" = 1, "betas" = rtrimodal_meth(n = 1000)),
        data.table("groups" = 2, "betas" = rtrimodal_meth(n = 1000)),
        data.table("groups" = 3, "betas" = rtrimodal_meth(n = 1000))
      )
      dt_sim_trimodal <- rbindlist(ls_grp)

      res_mctp <- mctp(
        betas ~ groups, data = dt_sim_trimodal, asy.method = asy,
        type = "Tukey", alternative = "two.sided", plot.simci = FALSE,
        info = FALSE)$Analysis.Inf$p.Value
    })
    unlist(ls_res_trimodal)
  })
names(ls_trimodal_asy) <- asy_method

res_trimodal <- vapply(X = names(ls_trimodal_asy), FUN = function(asy){
  paste0(round(length(ls_trimodal_asy[[asy]][ls_trimodal_asy[[asy]] <= 0.05])/
                 length(ls_trimodal_asy[[asy]])*100, 2), "%")
}, FUN.VALUE = character(length = 1L))

# Put results together
dt_sim_res <- as.data.table(
  rbind(res_bimodal, res_trimodal), keep.rownames = "distributions")
dt_sim_res$distributions <- c("rbimodal_meth()", "rtrimodal_meth()")
write(
  x = knitr::kable(x = dt_sim_res),
  file = "MCTP_simulation_on_methylation_results.txt")
