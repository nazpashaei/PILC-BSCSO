# PILC-BSCSO
An Efficient Binary Sand Cat Swarm Optimization for Feature Selection in High-dimensional Biomedical Data
a new version of Binary Sand Cat Swarm Optimization (called PILC-BSCSO), incorporating a pinhole-imaging-based learning strategy and crossover operator is presented for selecting the most informative features. First, the crossover operator is used to strengthen the search capability of BSCSO. Second, the pinhole-imaging learning strategy is utilized to effectively increase exploration capacity while avoiding premature convergence. The Support Vector Machine (SVM) classifier with a linear kernel is used to assess classification accuracy. 

##How to get started
Install R 4.1.3 and RStudio
Install R Packages
Download Microarray datasets in your main working directory.

Required dependencies (R packages)
library(MASS)
library(foreign)
library(caret)
library(kernlab)
library(farff)
library(irr)
library(naivebayes)
library(randomForest)
library(e1071)
library(klaR)

