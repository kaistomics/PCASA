# PCASA
PCASA (Prioritization of combinatorial cancer-associated surface antigens) is a tool for predicting the best gene combination targets of surface antigens by classification of malignant and non-malignant cells using single-cell RNA-seq data.

## Requirements
> Required modules for RF (Random forest) in R.

* R (>= 4.0): https://cloud.r-project.org/bin/linux/ubuntu%bionic-cran40/
* ROCR (>= 1.0-11), caret (>= 6.0-88), e1071 (>= 1.7-8), ggplot2 (>= 3.3.5), gplots (>=3.1.1), randomForest (>= 4.6-14), rpart (>= 4.1-15)

> Required modules for CNN (Convolutional neural network) in Python.

* Python (>= 3.6.9): https://www.python.org/downloads/source/
* pip3 (>= 21.2.4), numpy (>= 1.19.5), matplotlib (>= 3.3.4), pandas (>= 1.1.5), scikit-learn (>= 0.24.1), 
* keras (>= 2.6.0), tensorflow (>= 2.6.0), tf-keras-vis (>= 0.8.0)

## Running codes
PCASA is a three-step program :
* Step 1. prioritization of important single-genes well-classifying tumor and normal cells
* Step 2. calculation of the expressing cell fraction (ECF) for each gene-combination
* Step 3. prioritization of gene-combinations well-classifying tumor and normal cells.
All of the steps can be run with the following command lines.

