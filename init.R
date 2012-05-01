## Chose repository
local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.rakanu.com/"
  options(repos = r)
})

## Default packages
R_DEFAULT_PACKAGES='utils,grDevices,graphics,stats,ggplot2,reshape,plyr,sqldf'

## First and last things to do
.First <- function() cat("
   Welcome to R Ronert!

")
.Last <- function()  cat("
   Goodbye -.-!

")

## Check and install packages
library(utils)
list.of.packages <- c("ggplot2", "reshape", "reshape2", "lubridate", "plyr",
                      "stringr",  "digest", "mutatr", "brew",
                      "testthat", "sinartra", "helpr",  "RCurl",
                      "Rcpp", "XML", "MASS", "RTextTools", "tikzDevice",
                      "sqldf", "e1071", "foreign", "RMySQL", "RSQLite",
                      "xtable", "randomForest", "RColorBrewer",
                      "RODBC", "lattice", "Hmisc", "Sim.DiffProc",
                      "lme4", "data.table", "PerformanceAnalytics",
                      "boot", "MCMCglmm", "rpart", "survey",
                      "twitteR", "Matrix", "rgl", "zoo", "foreach",
                      "quantmod", "xts", "plm", "multicomp", "glmnet",
                      "tm", "coin", "caret", "plotrix", "reporttools",
                      "animation", "xlsReadWrite", "xlsx", "Cairo",
                      "iterators", "unknownR", "DPpackage",
                      "fitdistrplus", "IBrokers", "mlbench",
                      "heatmap.plus", "FNN", "vars", "memoise",
                      "MHadaptive", "RgoogleMaps", "klaR", "tree")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## tikzDevice for org-export
old <- getOption("defaultPackages")
options(defaultPackages = c(old, "tikzDevice"))

message("\n******************************\nSuccessfully loaded init.R\n******************************")
