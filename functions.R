####################
## Help Functions ##
####################

## Browse all Vignettes
man <- function() {
browseVignettes(package = NULL, lib.loc = NULL, all = TRUE)
}

## search cran for keywords by regex
cranSearch = function(q='')
{
  library(XML)
  u<-'http://cran.r-project.org/web/packages/'
  d<-readHTMLTable(u)
  d[[1]]$V1

  lib.names <- as.vector(d[[1]]$V1)
  lib.descriptions <- as.vector(d[[1]]$V2)
  lib.names[grep(q,lib.descriptions)]
}


## Listing objects in your global environment
## A simple ls() doesn’t really tell you enough useful information at
## a glance. Most often I just want to know what I named certain
## data.frames or functions. This handy little function, called as lsa() will do that for you:
lsa<-function()
  {
    obj_type<-function(x) { class(get(x)) }
    foo=data.frame(sapply(ls(envir=.GlobalEnv),obj_type))
    foo$object_name=rownames(foo)
    names(foo)[1]="class"
    names(foo)[2]="object"
    return(unrowname(foo))
  }

## Listing all functions in a certain package
## This can be called with lsp(). The pattern argument will allow you to quickly find the right function if you vaguely remember the name.
lsp <-function(package, all.names = FALSE, pattern)
{
  package <- deparse(substitute(package))
  ls(
     pos = paste("package", package, sep = ":"),
     all.names = all.names,
     pattern = pattern
     )
}

## Purpose: Package installation and loading
## from https://github.com/gimoya/theBioBucket-Archives/blob/master/R/Functions/instant_pkgs.R
instant_pkgs <- function(pkgs) {
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
        install.packages(pkgs_miss)
    }

        if (length(pkgs_miss) == 0) {
        message("\n ...Packages were already installed!\n")
        }

    # install packages not already loaded:
    pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
    if (length(pkgs_miss) > 0) {
        install.packages(pkgs_miss)
    }

    # load packages not already loaded:
    attached <- search()
    attached_pkgs <- attached[grepl("package", attached)]
    need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]

    if (length(need_to_attach) > 0) {
      for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
    }

        if (length(need_to_attach) == 0) {
        message("\n ...Packages were already loaded!\n")
        }
}

#######################
## Data Manipulation ##
#######################
## Transpose a numeric data frame with ID in first column
tdf <- function(d) {
  row.names(d) <- d[[1]]
  d[[1]] <- NULL
  d <- as.data.frame(t(d))
  d$id <- row.names(d)
  d <- cbind(d[ncol(d)], d[-ncol(d)])
  row.names(d) <- NULL
  d
}

## Convert selected columns of a data frame to factor variables
factorcols <- function(d, ...) lapply(d, function(x) factor(x, ...))

## Returns a logical vector TRUE for elements of X not in Y
"%nin%" <- function(x, y) !(x %in% y)

## Returns names(df) in single column, numbered matrix format.
n <- function(df) matrix(names(df))

## Single character shortcuts for summary() and head().
s <- base::summary
h <- utils::head

## ht==headtail, i.e., show the first and last 10 items of an object
ht <- function(d) rbind(head(d,10),tail(d,10))

## Show the first 5 rows and first 5 columns of a data frame or matrix
hh <- function(d) d[1:5,1:5]

## Takes a dataframe and a column name, and moves that column to the front of the DF.
moveColFront <- function(d=dataframe, colname="colname") {
  index <- match(colname, names(d))
  cbind(d[index],d[-index])
}

## Permutes a column in a data.frame, sets seed optionally
permute <- function (dataframe, columnToPermute="column", seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  colindex <- which(names(dataframe)==columnToPermute)
  permutedcol <- dataframe[ ,colindex][sample(1:nrow(dataframe))]
  dataframe[colindex] <- permutedcol
  return(dataframe)
}

## Summarize missing data in a data frame. Return a list (lpropmiss) or data frame (propmiss)
lpropmiss <- function(dataframe) lapply(dataframe,function(x) data.frame(nmiss=sum(is.na(x)), n=length(x), propmiss=sum(is.na(x))/length(x)))
propmiss <- function(dataframe) {
  m <- sapply(dataframe, function(x) {
    data.frame(
               nmiss=sum(is.na(x)),
               n=length(x),
               propmiss=sum(is.na(x))/length(x)
               )
  })
  d <- data.frame(t(m))
  d <- sapply(d, unlist)
  d <- as.data.frame(d)
  d$variable <- row.names(d)
  row.names(d) <- NULL
  d <- cbind(d[ncol(d)],d[-ncol(d)])
  return(d[order(d$propmiss), ])
}

############
## System ##
############

## Read data on clipboard.
read.cb <- function(...) {
  ismac <- Sys.info()[1]=="Darwin"
  if (!ismac) read.table(file="clipboard", ...)
  else read.table(pipe("pbpaste"), ...)
}

## Open current directory on mac
macopen <- function(...) system("open .")

## name("test.png") results in "C:/R/2010-04-20-test.png" if running this in C:/R on April 20 2010.
name <- function(filename="filename") paste(getwd(),"/",Sys.Date(),"-",filename,sep="")

###########
## Plots ##
###########

## Make a pretty QQ plot of p-values
qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
                                        #e = -log10( 1:length(o)/length(o) )
  e = -log10( ppoints(length(pvector) ))
  plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
  abline(0,1,col="red")
}

## Draw a histogram with normal overlay (From http://www.statmethods.net/graphs/density.html)
histnormal <- function(d, main=NULL, xlab=NULL, breaks="FD", ...) {
  if (any(is.na(d))) warning(paste(sum(is.na(d)), "missing values"));   d <- na.omit(d)
  h <- hist(d, plot=FALSE, breaks=breaks, ...)
  x <- seq(min(d), max(d), length=40)
  y <- dnorm(x, mean=mean(d), sd=sd(d))
  y <- y*diff(h$mids[1:2])*length(d)
  hist(d, col="gray50", main=main, xlab=xlab, ylim=c(0,max(y)), breaks=breaks,...)
  lines(x,y, col="blue", lwd=2)
  rug(x)
}

## Draw a histogram with density overlay
histdensity <- function(x, main=NULL, breaks="FD", ...) {
  if (any(is.na(x))) warning(paste(sum(is.na(x)), "missing values"));   x <- na.omit(x)
  hist(x, col="gray50", probability=TRUE,  breaks=breaks, main=main, ...)
  lines(density(x, na.rm = TRUE), col = "blue", lwd=2)
  rug(x)
}

## Plot scatterplot with trendline and confidence interval (From http://tinyurl.com/3bvrth7)
scatterci <- function(x, y, ...) {
  plot(x, y, ...)
  mylm <- lm(y~x)
  abline(mylm, col="blue")
  x=sort(x)
  prd<-predict(mylm,newdata=data.frame(x=x),interval = c("confidence"), level = 0.95)
  lines(x,prd[,2],col="blue",lty=3)
  lines(x,prd[,3],col="blue",lty=3)
}

## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

## Makes a better scatterplot matrix.
## Stolen from the PerformanceAnalytics package: http://cran.r-project.org/web/packages/PerformanceAnalytics/index.html
## Also see http://moderntoolmaking.blogspot.com/2011/08/graphically-analyzing-variable.html
## To color code points based on levels of a factor, use these args:
## pairs.perfan(d, bg=c("red","blue")[d$factor], pch=21)
betterpairs <- function (R, histogram = TRUE, ...)
{
  x=as.matrix(R) # in PerformanceAnalytics: x = checkData(R, method = "matrix")
  if (mode(x)!="numeric") stop("Must pass in only numeric values")
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use = use))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) cex <- 0.8/strwidth(txt)
    test <- cor.test(x, y)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * r)
    text(0.8, 0.8, Signif, cex = cex, col = 2)
  }
  f <- function(t) dnorm(t, mean = mean(x), sd = sd(x))
                                        # Useful function for histogram showing density overlay and rug
  hist.panel = function(x, ...) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
    rug(x)
  }
  if (histogram) pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, diag.panel = hist.panel, ...)
  else           pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, ...)
}

################
## Statistics ##
################

## Get the proportion variation explained. See this website for more details: http://goo.gl/jte8X
rsq <- function(predicted, actual) 1-sum((actual-predicted)^2)/sum((actual-mean(actual))^2)

## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function(X, dfr = nrow(X) - 2) {
  R <- cor(X)
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr / (1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R)==col(R)]<-NA
  R
}

## This function accepts a GLM object and does a LR chi-square test on the fit.
lrt <- function (modelobject) {
  lrtest.chi2 <- model$null.deviance - model$deviance # Difference in deviance between model with intercept only and full model.  This is the likelihood ratio test statistic (-2(log(L))).
  lrtest.df   <- model$df.null - model$df.residual # Difference in DF.  Make sure this equals the number of predictors in the model!
  fitpval     <- 1-pchisq(lrtest.chi2,lrtest.df)
  cat("Likelihood ratio test on model fit:\n\n")
  data.frame(lrtest.chi2=lrtest.chi2,lrtest.df=lrtest.df,fitpval=fitpval) #Output gives you the chisquare, df, and p-value.
}

## This function does the same thing as lrtest in the Design package, but doesn't do as much checking.
## Remember, the lrt has to test the same model (model fit on same observations)
## Also the drop1(fullmodel,test="Chisq") does something similar.
lrt2 <- function (full,reduced) {
  if (reduced$deviance<=full$deviance) stop ("Reduced model not worse than full.")
  if (reduced$df.residual<=full$df.residual) stop ("Reduced model doesn't have more degrees of freedom.")
  lrtest.chi2 <- reduced$deviance-full$deviance
  lrtest.df   <- reduced$df.residual - full$df.residual
  fitpval        <- 1-pchisq(lrtest.chi2,lrtest.df)
  cat("Likelihood ratio test on two models:\n\n")
  data.frame(lrtest.chi2=lrtest.chi2,lrtest.df=lrtest.df,fitpval=fitpval)
}

## This gets the overall anova p-value out of a linear model object
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

## Imputes the median value of a vector, matrix, or data frame.
## Stolen from na.roughfix function in the randomForest package.
na.roughfix <- function (object=NULL, ...) {
  if (class(object) == "data.frame") {
    isfac <- sapply(object, is.factor)
    isnum <- sapply(object, is.numeric)
    if (any(!(isfac | isnum)))
      stop("dfMedianImpute only works for numeric or factor")
    roughfix <- function(x) {
      if (any(is.na(x))) {
        if (is.factor(x)) {
          freq <- table(x)
          x[is.na(x)] <- names(freq)[which.max(freq)]
        }
        else {
          x[is.na(x)] <- median(x, na.rm = TRUE)
        }
      }
      x
    }
    object[] <- lapply(object, roughfix)
    return(object)
  }
  else if(is.atomic(object)) {
    d <- dim(object)
    if (length(d) > 2)
      stop("vectorMedianImpute can't handle objects with more than two dimensions")
    if (all(!is.na(object)))
      return(object)
    if (!is.numeric(object))
      stop("vectorMedianImpute can only deal with numeric data.")
    if (length(d) == 2) {
      hasNA <- which(apply(object, 2, function(x) any(is.na(x))))
      for (j in hasNA) object[is.na(object[, j]), j] <- median(object[,
                                                                      j], na.rm = TRUE)
    }
    else {
      object[is.na(object)] <- median(object, na.rm = TRUE)
    }
    return(object)
  }
  else stop("Object is not a data frame or atomic vector")
}

# Display/extract regression coefficients
coef.lme <- function(mod){
    res <- data.frame(
        "Beta.CI" = paste(round(summary(mod)$coefficients$fixed, 3), " (",round(summary(mod)$coefficients$fixed-1.96*sqrt(diag(mod$varFix)), 2), ",", round(summary(mod)$coefficients$fixed+1.96*sqrt(diag(mod$varFix)), 2),")", sep=""),
        "P.value" = round(2 * pt(-abs(summary(mod)$coefficients$fixed/sqrt(diag(mod$varFix))), summary(mod)$fixDF[[1]]), 3)
    )
    return(res)
}

#########
## End ##
#########
message("\n******************************\nSuccessfully loaded functions.R\n******************************\n")


#####################
## Display Fortune ##
#####################
library(fortunes);
message(fortune(sample(1:nrow(read.fortunes()),1)))
message("\n")
