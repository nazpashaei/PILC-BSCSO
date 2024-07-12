# PILC-BSCSO: An Efficient Binary Sand Cat Swarm Optimization for Feature Selection in High-dimensional Biomedical Data

PILC-BSCSO (Pinhole-Imaging Learning based Crossover Binary Sand Cat Swarm Optimization) is an enhanced version of Binary Sand Cat Swarm Optimization (BSCSO), designed for feature selection in high-dimensional biomedical data. It incorporates a pinhole-imaging-based learning strategy and a crossover operator to improve search capability and exploration capacity. The Support Vector Machine (SVM) classifier with a linear kernel is used to assess classification accuracy.

## How to Get Started

### Prerequisites

1. **Install R 4.1.3 and RStudio**
   - Download and install [R](https://www.r-project.org/) version 4.1.3 or later.
   - Download and install [RStudio](https://www.rstudio.com/) for a user-friendly interface.

2. **Install Required R Packages**
   - Open R or RStudio and install the following packages:

     ```R
     install.packages(c(
       "MASS", "foreign", "caret", "kernlab", "farff", "irr", "naivebayes",
       "randomForest", "e1071", "klaR"
     ))
     ```

3. **Download Microarray Datasets**
   - Download microarray datasets and place them in your main working directory.

### Usage

1. **Load Required Packages**

   ```R
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
##citation
Pashaei, E. An Efficient Binary Sand Cat Swarm Optimization for Feature Selection in High-Dimensional Biomedical Data. Bioengineering 2023, 10, 1123. https://doi.org/10.3390/bioengineering10101123.
