setwd("/Users/lawrence/Desktop/stats790/Assignment1/")

X <- read.csv("ChettyEtAl_CZ_data.csv")
Y <- read.csv("ChettyEtAl_mobility.csv")

dim(X)  # Checking to see whether data has same number of rows; should be 741
dim(Y)

X[1,]   # looking for where X has an extra row
X[742,]
X <- X[1:741,]  # For some reason, X had a blank row at the bottom. Removed here

colnames(Y)  # Checking to see where the relevant mobility stat appears
             # turns out its in colum 7, "P.Child.in.Q5...Parent.in.Q1...80.85.Cohort"
Z <- Y[,7]  # This is our regressand, the variable we want to predict 

NAspots <- is.na(Z)
Z <- Z[!NAspots]
X <- X[!NAspots,]

#####################
#####################
##### PROBLEM A #####
#####################
#####################

colnames(X)
  # 4 - Population; 12 - household Income; 7 - Racial segregation
  # 14 - Top 1% income share; 21 - school expenditure per student
  # 37 - violent crime rate; 11 - % with short commute

for (i in c(4,12,7,14,21,37,11)) {
  myFit <- lm(Z ~ X[,i]); myFit
  png(filename = paste(colnames(X)[i],"_vs_mobility.png",sep=""), width = 400, height = 400)
  plot(X[,i],Z, xlab = colnames(X)[i], ylab = colnames(Y)[7])
  lines(range(X[,i], na.rm=TRUE), myFit$coefficients[1] + myFit$coefficients[2]*range(X[,i], na.rm=TRUE))
  dev.off()
}

myFit <- lm(Z ~ log(X[,4])); myFit
png(filename = paste("log.",colnames(X)[4],"_vs_mobility.png",sep=""), width = 400, height = 400)
plot(log(X[,4]),Z, xlab = paste("log.",colnames(X)[4],sep=""), ylab = colnames(Y)[7])
lines(range(log(X[,4]), na.rm=TRUE), myFit$coefficients[1] + myFit$coefficients[2]*range(log(X[,4]), na.rm=TRUE))
dev.off()

#####################
#####################
##### PROBLEM B #####
#####################
#####################

#####################################################################
# Regress against racial segregation, using Gaussian kernel smoothing
#####################################################################

myKSmooth <- function(x, y, kernel = c("box", "normal"), bandwidth = 0.5,
                      range.x = range(x), x.points)
{
  # performs ksmooth, and fills in NA predictions with 1-nearest neighbor projections
  # (which is a very good approximation of kernel smoothing in a region of sparse data,
  # which is where ksmooth can erroneously give NA)
  myFit <- ksmooth(x,y,kernel=kernel,bandwidth=bandwidth,range.x=range.x,x.points=x.points)
    # ... run ksmooth
  for ( j in c( 1:length(x.points) )[is.na(myFit$y)] )
      # Loop through all indices j where this fit gives an NA projection
    myFit$y[j] <- mean(y[which.min(abs(x-x.points[j]))])
      # "which.min" finds the index (or indices) i where x[i] gets closest to the x coord
      #   x.points[j] corresponding to the missing y value. 
      # "y[ ... ]" takes the input y-value corresponding to those x[i]. 
      # "mean( ... )" averages in case of an exact tie
  myFit
}

myBand <- sqrt(var(X[,7],na.rm=TRUE)); myBand
mySeq <- seq(from = min(X[,7],na.rm=TRUE),to = max(X[,7],na.rm=TRUE), length.out = 4000)
myFit <- lm(Z ~ X[,7])$coefficients

png(filename = paste("KSmooth_",colnames(X)[7],"_vs_mobility.png",sep=""), width = 400, height = 400)
plot(X[,7],Z, xlab = colnames(X)[7], ylab = colnames(Y)[7], col="grey")
lines(myKSmooth(X[,7], Z, "normal", bandwidth = myBand/3, x.points = mySeq), col = "green")
lines(myKSmooth(X[,7], Z, "normal", bandwidth = myBand, x.points = mySeq), col = "blue")
lines(myKSmooth(X[,7], Z, "normal", bandwidth = myBand*3, x.points = mySeq), col = "red")
lines(myKSmooth(X[,7], Z, "normal", bandwidth = 9*myBand, x.points = mySeq), col = "magenta")
lines(myKSmooth(X[,7], Z, "normal", bandwidth = 27*myBand, x.points = mySeq), col = "orange")
lines(range(X[,7], na.rm=TRUE), myFit[1] + myFit[2]*range(X[,7], na.rm=TRUE))
legend("topright",legend = c("linear", paste("h =",round(27*myBand,4)), paste("h =",round(9*myBand,4)), 
                             paste("h =",round(3*myBand,4)), paste("h =",round(myBand,4)), 
                             paste("h =",round(myBand/3,4))),
       lty = rep("solid",6), col = c("black", "orange", "magenta", "red", "blue", "green"))
dev.off()

###########################################################
# 10-fold cross-validation for kernel smoothing and lin reg
###########################################################

# Z has 729 entries
myAssign <- c(rep(1:10,72),1:9)
  # We represent each 10-fold split as a sequence of seventy
  #   three 1s, 2s, ..., 9s, and seventy two 10s
  # To generate a new split, we permute this sequence

numSplits <- 100
  # number of times to run cross-validation for each regression method
bWs <- myBand*c(27,9,3,1,1/3)
numBands <- length(bWs)
ksMSE <- rep(0,numBands)
  # MSE in ksmooth with h in bWs
lrMSE <- 0
  # MSE in linear regression

for (i in 1:numSplits) {
  myAssign <- sample(myAssign)
    # generate a new 10-fold split
    # to minimize the number of times we need to do this, we use the same split for each of our regressions
  for (m in 1:numBands) {
      # cycle through the chosen bandwidths for kernel smoothing
    for (j in 1:10) {
      ksMSE[m] = ksMSE[m] + 
        sum((myKSmooth(X[myAssign != j,7], Z, "normal", bandwidth = bWs[m], x.points = X[myAssign == j,7])$y
             - Z[myAssign == j])^2)
    }
  }
    # now do linear regression
  for (j in 1:10) {
    myFit <- lm( Z[myAssign!=j] ~ X[myAssign != j,7] )$coefficients
    lrMSE <- lrMSE + sum((myFit[1] + myFit[2]*X[myAssign==j,7] - Z[myAssign==j])^2)
  }
}
ksMSE <- ksMSE / (729*numSplits); ksMSE
lrMSE <- lrMSE / (729*numSplits); lrMSE


##################################################################
# Now let's try k nearest neighbors, again with racial segregation
##################################################################

library(FNN)
mySeq <- matrix(mySeq, byrow = TRUE)
myFit <- lm(Z ~ X[,7])$coefficients

png(filename = paste("KNN_",colnames(X)[7],"_vs_mobility.png",sep=""), width = 400, height = 400)
plot(X[,7],Z, xlab = colnames(X)[7], ylab = colnames(Y)[7], col="grey")
lines(mySeq, knn.reg(train = X[,7], test = mySeq, y = Z, k = 20)$pred,
      col = "green")
lines(mySeq, knn.reg(train = X[,7], test = mySeq, y = Z, k = 40)$pred,
      col = "magenta")
lines(mySeq, knn.reg(train = X[,7], test = mySeq, y = Z, k = 80)$pred,
      col = "blue")
lines(mySeq, knn.reg(train = X[,7], test = mySeq, y = Z, k = 160)$pred,
      col = "red")
lines(range(X[,7], na.rm=TRUE), myFit[1] + myFit[2]*range(X[,7], na.rm=TRUE))
legend("topright",legend = c("linear", "k = 160", "k = 80", "k = 40", "k = 20"),
       lty = rep("solid",6), col = c("black", "red", "blue", "magenta", "green"))
dev.off()

##########################################################
# 10-fold cross-validation for kNN and lin reg
##########################################################

numSplits <- 100
# number of times to run cross-validation for each regression method
k <- c(160,80,40,20)
numBands <- length(k)
knnMSE <- rep(0,numBands)
  # MSE in ksmooth with h in bWs
lrMSE <- 0
  # MSE in linear regression

for (i in 1:numSplits) {
  myAssign <- sample(myAssign)
  # generate a new 10-fold split
  # to minimize the number of times we need to do this, we use the same split for each of our regressions
  for (m in 1:numBands) {
    # cycle through the chosen bandwidths for kernel smoothing
    for (j in 1:10) {
      knnMSE[m] = knnMSE[m] + 
        sum( ( knn.reg(train = X[myAssign!=j,7], test = matrix(X[myAssign==j,7],byrow=TRUE), 
                       y = Z[myAssign!=j], k = k[m])$pred - Z[myAssign == j] )^2 )
    }
  }
  # now do linear regression
  for (j in 1:10) {
    myFit <- lm( Z[myAssign!=j] ~ X[myAssign != j,7] )$coefficients
    lrMSE <- lrMSE + sum((myFit[1] + myFit[2]*X[myAssign==j,7] - Z[myAssign==j])^2)
  }
}
knnMSE <- knnMSE / (729*numSplits); knnMSE
lrMSE <- lrMSE / (729*numSplits); lrMSE

#####################
#####################
##### PROBLEM C #####
#####################
#####################

##############
# ADA prob 3 #
##############

# We do not want to include CZ (identification number) in our regression, as this
# is unique for each census zone. Likewise, CZ Name is not useful.

# We could in principle include state as a factor, but with the large number of 
# distinct states (51, including DC) relative to the size of the data set and 
# several states with very few CZs, we runa serious risk of overfitting state 
# differences and washing out the signal in our other factors.

##################################################################################################
# Ignoring CZ, Name, and State, but otherwise keeping covariates and throwing away incomplete rows
##################################################################################################

dim(X)
  # 729 x 41
myFit <- lm(Z ~ as.matrix(X[,4:41]))
summary(myFit)
  # Residual std error: .02021

################################################################
# Alternatively, also throwing away highly incomplete covariates
################################################################

completeRows <- complete.cases(X)
sum(completeRows)
  # 418 rows are complete
naCount <- rep(0,41)
for(j in 1:41) naCount[j] <- sum(is.na(X[,j]))
naCount
  # there are four covariates in which lots (> 100) of data is missing, in columns 24:27
  # in all other covariates, at most 27 entries missing
sum(is.na(X[,24]) | is.na(X[,25]) | is.na(X[,26]) | is.na(X[,27]))
  # 258 out of 729 rows in which at least one of these is missing
colnames(X)[24:27]
  # HS dropout rate; num colleges per capita; college tuition; college grad rate
  # In the previous regression, only num colleges per capita was a significant regressor

newX <- as.matrix(X[,c(4:23,28:41)])
  # Exclude the first three columns
sum(complete.cases(newX))
  # counting complete rows now; up from 418 to 633
summary(lm(Z[completeRows] ~ newX[completeRows,]))
  # Before fitting all rows that are now complete, we try the regression on 
  #   the 418 rows that were complete previously, with the reduced set of covariates
  # residual std error: 0.02037
  # looks like a modest increase from when we included the education covariates
myFit <- lm(Z ~ newX)
summary(myFit)


##############
# ADA prob 4 #
##############

which(X[,2] == "Pittsburgh")
  # row 229
complete.cases(X[229,])
  # TRUE, so we could use either the regression including educational factors or
  # the one without
Z[229]
  # Find the mobility of Pittsburgh
PittPrediction <-  sum(myFit$coefficients * c(1,newX[229,])); PittPrediction
  # Find predicted mobility. Hey, not bad!
  # Could use myFit$fitted.values, but then would need to find which row in the 
  # reduced matrix, after removing incompletes, is Pittsburgh.
  # I prefer seeing the linear algebra

colnames(newX)
  # 30 - violent crime rate; 5 - income segregation; 11 - top 1% income share

PittPrediction + myFit$coefficients[31] * newX[229,30]
  # predicted mobility in Pittsburgh if violent crime doubled
PittPrediction - myFit$coefficients[31] * (newX[229,30] / 2)
  # predicted mobility in Pittsburgh if violent crime halved

((1 - PittPrediction) / myFit$coefficients[6]) + newX[229,5]
  # the expression inside parentheses computes how much income segregation
  # would need to change for mobility to reach 1; then we add how much
  # income segregation we already have to find the final value needed
newX[229,11] - (PittPrediction / myFit$coefficients[12])
  # the expression inside parentheses computes how much top 1% income share
  # would need to change for mobility to reach 0; then we add how much
  # income share they already have

##############
# ADA prob 5 #
##############

sum(Z*X[,4]) / sum(X[,4])
  # average mobility over census zones, weighted by population

# This problem is asking about changing college tuition, so we need to go 
#   back to the regression that includes that covariate
# College tuition is X[,26]

myFit <- lm(Z ~ as.matrix(X[,4:41]))
myFit$coefficients[24]
# Because myFit includes an intercept but excludes X[1:3,] terms, 
# myFit$coefficients[24] is the coeff associated with college tuition

PredChange <- -1 * myFit$coefficients[24] * X[,26]
  # predicted mobility change if tuition sent to 0
range(PredChange , na.rm = TRUE)
mean(PredChange , na.rm=TRUE)
median(PredChange, na.rm=TRUE)

sum(Z[!is.na(X[,26])]*X[!is.na(X[,26]),4]) / sum(X[!is.na(X[,26]),4])
  # pop-weighted average mobility over CZs where we have college tuition data
  # (slightly lower than overall pop-weighted average mobility)

sum(PredChange * X[,4], na.rm=TRUE) / sum(X[!is.na(X[,26]),4])
  # pop-weighted average predicted change in mobility over CZs where we have 
  # tuition data

# For the last part of the problem, we need the std error associated with 
# myFit$coefficients[24]. From summary(myFit), this is approx. 3.599e-07.

