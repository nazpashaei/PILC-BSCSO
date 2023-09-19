library(RSNNS)
library(glmnet)
library(ggpubr)
library(matlab)
library(monmlp)
library(rgl)
library(car)
library(MASS)
library(foreign)
library(caret)
library(kernlab)
library(preprocessCore)
library(farff)
library(statip)
library(irr)
library(naivebayes)
library(forcats)
library(caret)
library(randomForest)
library(e1071)
library(CMA)
library(tidyverse)
library(C50)
library(xgboost)
library(klaR)
library(CMA)
library(corpcor)
library(xgboost)
library(pROC)
library(mltools)
library(pls)
library(plsRglm)
library(MLmetrics)
library(Metrics)
library(yardstick)
library(mltools)
library(ROCit)
library(RSSL)
library(sparsediscrim)
library(pracma)
library(round)
library(ipred)

ptm <- proc.time()
set.seed(2)

#==================================================================================
data <- readARFF("~/Zscore-colon-DEGS-358.arff")
class <- data[, ncol(data)]

data <- data.frame(data, check.names = TRUE)
dim(data)
names(data)[1:4] 
iden <- names(data)
rownames(data)
iden[ncol(data)]

#----------------------------------------
class <-data[, ncol(data)]
#------------------------------------------------
levels(as.factor(class))

col = (ncol(data))


l = class
l = as.factor(l)

oout <- matrix(0, nrow = 1, ncol = 3)


#===================================================================================

fitness <- function(star) {
  index <- which(star == 1)
  if (length(index) == 0 | length(index) == 1) {
    return(0)
  }
  subData <- data[, index]
  n2 = ncol(subData)
  
  dimo <- dim(subData)
  subData$class <- data$class
  
  #----------------------------------
  #SETTING k-fold or LOOV
  type <- "k-fold"
  kf = 10
  
  score = list()
  for (zz in 1:3) {
    if (type == "k-fold") {
      a <- createFolds(subData$class, k = kf, list = FALSE)
      subData$id <- a
      t = kf
      
    } else{
      t <- nrow(subData)
      subData$id <- 1:t
    }
    #============================================================================
    
    for (i in 1:t) {
      train1 <- subData[subData$id != i,-ncol(subData)]# delete id
      test1 <-  subData[subData$id == i,-ncol(subData)]#
      test_lable <- test1$class
      test2 <- test1[, -ncol(test1)]
      
      trainyy <- train1$class
      testyy <- test1[, ncol(test1)]
      
      
      
      model = svm(
        class ~ .,
        data = train1,
        type = 'C-classification',
        kernel = "linear",
        scale = FALSE
      )
      pred <- predict(model, test2, type = "class")
      results <- kappa2(data.frame(test_lable, pred))$value;accuracy1 <-results*100
      score[[i]] = accuracy1

    }
    
    oout[zz] <- round(mean(unlist(score)), digits = 2)
  }
  return((mean(oout)))
  
}
#=======================================================

nfeature <- function(vec) {
  vec <- as.numeric(vec)
  a <- table(vec)
  out <- a[names(a) == 1]
  if (length(out) == 0) {
    out = 0
  }
  return(out)
}

#--------------
dim(data)
row = nrow(data)
#==============================================================================
RouletteWheelSelection <- function(p) {
  r = runif(1, 0)
  
  s = sum(p)
  
  p = p / s
  
  C = cumsum(p)
  
  j = which.max(r <= C)
  return(j)
}




#Best_PD,PDBest_P,PDConv
SearchAgents_no = 100

Max_iter = 50

LB = 0
UB = 1
nVar = col - 1
Dim = nVar



BestFit = matrix(0, nrow = 1, ncol = Dim)
fitness1 = matrix(0, nrow = 1, ncol = SearchAgents_no)
Best_Score = -Inf

Positions=matrix(data = sample(0:1, size = SearchAgents_no * Dim, replace = TRUE),nrow = SearchAgents_no,ncol = Dim)
Convergence_curve = matrix(0, nrow = 1, ncol = Max_iter)
Convergence_feature= matrix(0, nrow = 1, ncol = Max_iter)
p = matrix(1:nrow(Positions),
           nrow = 1,
           ncol = nrow(Positions))

t = 0

k = runif(1, 0)

#-----------------------------------------------
while (t < Max_iter) {
  for (i in 1:nrow(Positions)) {
    fitness1[1, i] = fitness(Positions[i, ])
    
    if (fitness1[1, i] > Best_Score) {
      Best_Score = fitness1[1, i]
      
      BestFit[1, ] = Positions[i, ]
      
    }
  }
  S = 1
  #floor(runif(1, min=1, max=10))   ### S is maximum Sensitivity range
  rg = S - ((S) * t / (Max_iter))
  ### guides R
  for (i in 1:nrow(Positions)) {
    r = runif(1, 0) * rg
    
    R = ((2 * rg) * runif(1, 0)) - rg
    ###   controls to transition phases
    for (j in 1:ncol(Positions)) {
      teta = RouletteWheelSelection(p)
      
      if (0 < R & R < 1) {
        ### R value is between -1 and 1
        Rand_position = abs(runif(1, 0) * BestFit[1, j] - Positions[i, j])
        
        Positions[i, j] = BestFit[1, j] - r * Rand_position * cos(teta)
        
      } else{
        cp = floor(SearchAgents_no * runif(1, 0) + 1)
        
        CandidatePosition = Positions[cp, ]
        
        Positions[i, j] = r * (CandidatePosition[j] - runif(1, 0) * Positions[i, j])
        
      }
      if (2 / (1 + exp(-2 * Positions[i, j])) - 1 > runif(1, 0)) {
        Positions[i, j] = 1
      } else{
        Positions[i, j] = 0
      }
      
    }
    #============================
    p = runif(1, 0)
    
    if (TRUE) {
      if (p < 0.5) {
        P = CrossoverX(BestFit[1, ], Positions[i, ])
        fp1 = fitness(P[1, ])
        fp2 = fitness(P[2, ])
        #
        if ((fp1 > fp2) &&
            (fp1 > Best_Score) ||
            ((fp1 == Best_Score) &&
             (nfeature(P[1, ])[[1]] < nfeature(BestFit)[[1]]))) {
          BestFit = P[1, ]
          Best_Score = fp1
          cat("crossover1", '\n')
          
          cat(Best_Score, '\n')
        } else if (fp2 > Best_Score) {
          BestFit = P[2, ]
          Best_Score = fp2
          cat("crossover2", '\n')
          
          cat(Best_Score, '\n')
        }
      }
    }
    #enhance part=============================
    if (TRUE) {
      if (p > 0.5) {
        k = 0.05
        GO = (LB + UB) / 2 + (LB + UB) / (2 * k) - Positions[i, ] / k
        sss = and(GO, BestFit[1,])
        solo = as.integer(sss)
        
        fsolo = fitness(solo)
        if ((fsolo > Best_Score) ||
            ((
              fsolo == Best_Score && nfeature(solo) < nfeature(BestFit[1,])
            ))) {
          cat("opposition-based", '\n')
          Best_Score = fsolo
          BestFit[1,] = solo
          
          Positions[i, ] = solo
          fitness1[i] = fsolo
        }
      }
    }
    #---------------------------------
  }
  t = t + 1
  
  Convergence_curve[1, t] = Best_Score
  Convergence_feature[1, t] = sum(BestFit[1, ])
  cat(
    "At iteration",
    (t),
    "the best solution fitness is:",
    (Best_Score),
    "Number of features:",
    sum(BestFit[1, ]),
    '\n'
  )
  
}
e = proc.time() - ptm
cputime = 0
cputime = e[1] + e[2]
print(BestFit[1, ])
q1 = which(BestFit[1, ] > 0)
print(q1)
print(cputime)
print(length(q1))
print(Convergence_curve)
print(Convergence_feature)