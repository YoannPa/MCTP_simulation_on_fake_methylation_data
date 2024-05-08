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

##PARAMETERS
asy_method <- c("fisher", "mult.t", "log.odds", "normal")
sample_sizes <- c(10, 100, 1000)

##ANALYSIS
# Simulate MCTP with Tukey contrast 10,000 times using different asymptotic
# approximation on random generation of bimodal
# methylation distributions
ls_bimodal_asy <- lapply(X = asy_method, FUN = function(asy){
  cat(asy,"\n")
  ls_sample_sizes <- lapply(X = sample_sizes, FUN = function(sampl){
    cat("\tSample size=", sampl, "\n")
    if(Sys.info()[["sysname"]] == "Linux"){
      ls_res_bimodal <- mclapply(X = seq(1000), mc.cores = 15, FUN = function(i){
        ls_grp <- list(
          data.table("groups" = 1, "betas" = rbimodal_meth(n = sampl)),
          data.table("groups" = 2, "betas" = rbimodal_meth(n = sampl)),
          data.table("groups" = 3, "betas" = rbimodal_meth(n = sampl))
        )
        dt_sim_bimodal <- rbindlist(ls_grp)
        
        res_mctp <- mctp(
          betas ~ groups, data = dt_sim_bimodal, asy.method = asy,
          type = "Tukey", alternative = "two.sided", plot.simci = FALSE,
          info = FALSE)$Analysis.Inf$p.Value
      })
    } else {
      ls_res_bimodal <- lapply(X = seq(1000), FUN = function(i){
        ls_grp <- list(
          data.table("groups" = 1, "betas" = rbimodal_meth(n = sampl)),
          data.table("groups" = 2, "betas" = rbimodal_meth(n = sampl)),
          data.table("groups" = 3, "betas" = rbimodal_meth(n = sampl))
        )
        dt_sim_bimodal <- rbindlist(ls_grp)
        
        res_mctp <- mctp(
          betas ~ groups, data = dt_sim_bimodal, asy.method = asy,
          type = "Tukey", alternative = "two.sided", plot.simci = FALSE,
          info = FALSE)$Analysis.Inf$p.Value
      })
    }
    unlist(ls_res_bimodal)
  })
  names(ls_sample_sizes) <- as.character(sample_sizes)
  ls_sample_sizes
})
names(ls_bimodal_asy) <- asy_method

# Get percentages of p.values below 0.05 for each asymptotic approximation
# and sample sizes
res_bimodal <- lapply(X = names(ls_bimodal_asy), FUN = function(asy){
  ls_asy_res <- lapply(X = seq_along(ls_bimodal_asy[[asy]]), FUN = function(i){
    err_rate <- paste0(round(
      length(ls_bimodal_asy[[asy]][[i]][ls_bimodal_asy[[asy]][[i]] <= 0.05])/
        length(ls_bimodal_asy[[asy]][[i]])*100, 2), "%")
    data.table(
      "Asympt. approx." = asy,
      "Sample size" = names(ls_bimodal_asy[[asy]])[i],
      "Sim. err. rate" = err_rate)
  })
  rbindlist(l = ls_asy_res)
})
res_bimodal <- rbindlist(l = res_bimodal)

# Same with trimodal methylation distributions
ls_trimodal_asy <- lapply(X = asy_method, FUN = function(asy){
  cat(asy,"\n")
  ls_sample_sizes <- lapply(X = sample_sizes, FUN = function(sampl){
    cat("\tSample size=", sampl, "\n")
    if(Sys.info()[["sysname"]] == "Linux"){
      ls_res_trimodal <- mclapply(X = seq(1000), mc.cores = 15, FUN = function(i){
        ls_grp <- list(
          data.table("groups" = 1, "betas" = rtrimodal_meth(n = sampl)),
          data.table("groups" = 2, "betas" = rtrimodal_meth(n = sampl)),
          data.table("groups" = 3, "betas" = rtrimodal_meth(n = sampl))
        )
        dt_sim_trimodal <- rbindlist(ls_grp)
        
        res_mctp <- mctp(
          betas ~ groups, data = dt_sim_trimodal, asy.method = asy,
          type = "Tukey", alternative = "two.sided", plot.simci = FALSE,
          info = FALSE)$Analysis.Inf$p.Value
      })
    } else {
      ls_res_trimodal <- lapply(X = seq(1000), FUN = function(i){
        ls_grp <- list(
          data.table("groups" = 1, "betas" = rtrimodal_meth(n = sampl)),
          data.table("groups" = 2, "betas" = rtrimodal_meth(n = sampl)),
          data.table("groups" = 3, "betas" = rtrimodal_meth(n = sampl))
        )
        dt_sim_trimodal <- rbindlist(ls_grp)
        
        res_mctp <- mctp(
          betas ~ groups, data = dt_sim_trimodal, asy.method = asy,
          type = "Tukey", alternative = "two.sided", plot.simci = FALSE,
          info = FALSE)$Analysis.Inf$p.Value
      })
    }
    unlist(ls_res_trimodal)
  })
  names(ls_sample_sizes) <- as.character(sample_sizes)
  ls_sample_sizes
})
names(ls_trimodal_asy) <- asy_method

res_trimodal <- lapply(X = names(ls_trimodal_asy), FUN = function(asy){
  ls_asy_res <- lapply(X = seq_along(ls_trimodal_asy[[asy]]), FUN = function(i){
    err_rate <- paste0(round(
      length(ls_trimodal_asy[[asy]][[i]][ls_trimodal_asy[[asy]][[i]] <= 0.05])/
        length(ls_trimodal_asy[[asy]][[i]])*100, 2), "%")
    data.table(
      "Asympt. approx." = asy,
      "Sample size" = names(ls_trimodal_asy[[asy]])[i],
      "Sim. err. rate" = err_rate)
  })
  rbindlist(l = ls_asy_res)
})
res_trimodal <- rbindlist(l = res_trimodal)

# Put results together
ls_res <- list(
  "`rbimodal_meth()`" = res_bimodal, "`rtrimodal_meth()`" = res_trimodal)
dt_sim_res <- rbindlist(l = ls_res, idcol = "Distributions")
write(
  x = knitr::kable(x = dt_sim_res),
  file = "MCTP_simulation_on_methylation_results.txt")
