## Chose repository
local({
  r <- getOption("repos")
  r["CRAN"] <- "http://cran.rakanu.com/"
  options(repos = r)
})

## Check and install packages
library(utils)
list.of.packages <- c("devtools", "ggplot2", "reshape", "reshape2", "lubridate", "plyr",
"stringr", "digest", "mutatr", "brew", "testthat", "sinartra",
"helpr", "RCurl", "Rcpp", "XML", "MASS", "RTextTools", "tikzDevice",
"sqldf", "e1071", "foreign", "RMySQL", "RSQLite", "xtable",
"randomForest", "RColorBrewer", "RODBC", "lattice", "Hmisc",
"Sim.DiffProc", "lme4", "data.table", "PerformanceAnalytics", "boot",
"bootstrap", "MCMCglmm", "rpart", "survey", "twitteR", "Matrix",
"rgl", "zoo", "foreach", "quantmod", "xts", "plm", "multicomp", "tm",
"coin", "caret", "plotrix", "reporttools", "animation",
"xlsReadWrite", "xlsx", "Cairo", "iterators", "unknownR", "DPpackage",
"fitdistrplus", "IBrokers", "mlbench", "heatmap.plus", "FNN", "vars",
"memoise", "MHadaptive", "RgoogleMaps", "klaR", "tree", "rggobi",
"classifly", "clusterfly", "meifly", "profr", "ProjectTemplate",
"makeR", "benchmark", "ElemStatLearn", "mclust", "flexmix", "magic",
"nnet", "maxLik", "TeachingDemos", "googleVis", "rainbow",
"orderbook", "segue", "partykit", "hints", "mvtnorm", "formatR",
"mboost", "ggmap", "MCMCpack", "ada", "Boruta", "bst", "caTools",
"class", "Cubist", "earth", "elasticnet", "ellipse", "evtree",
"fastICA", "foba", "gam", "GAMens", "gbm", "glmnet", "grid", "hda",
"HDclassif", "ipred", "kernlab", "klaR", "kohonen", "KRLS", "lars",
"leaps", "LogicForest", "logicFS", "LogicReg", "mboost", "mda",
"mgcv", "mlbench", "neuralnet", "nnet", "nodeHarvest", "obliqueRF",
"pamr", "partDSA", "party", "penalized", "penalizedLDA", "pls",
"pROC", "proxy", "qrnn", "quantregForest", "randomForest", "RANN",
"rda", "relaxo", "rFerns", "rocc", "rpart", "rrcov", "RRF", "rrlda",
"RSNNS", "RWeka", "sda", "SDDA", "sparseLDA", "spls", "stepPlr",
  "superpc", "vbmp", "arules", "arulesViz", "doMC", "doParallel", "gregmisc",
                      "rdatamrket", "datamart", "RUnit")

new.packages <- list.of.packages[!(list.of.packages %in%
                                   installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## Default Packages
# old <- getOption("defaultPackages")
# options(defaultPackages = c(old, list.of.packages))

# Update Packages
update.packages(ask=FALSE)

# Some useful aliases
cd <- setwd
pwd <- getwd
lss <- dir

# Override q() to not save by default.
# Same as saying q("no")
q <- function (save="no", ...) {
  quit(save=save, ...)
}

message("\n******************************\nSuccessfully loaded init.R\n******************************")
