# PCASA
#### PCASA (Prioritization of Combinatorial Cancer-Associated Surface Antigens) is a tool for predicting the best gene combination targets of surface antigens by classification of malignant and non-malignant cells using single-cell RNA-seq data.

## Requirements
> Required packages for RF (Random forest) in R.

* R (>= 3.6.0): https://cloud.r-project.org/bin/linux/ubuntu/
* ROCR (>= 1.0-11), caret (>= 6.0-88), e1071 (>= 1.7-8), ggplot2 (>= 3.3.5), gplots (>=3.1.1), randomForest (>= 4.6-14), rpart (>= 4.1-15)

> Required modules for CNN (Convolutional neural network) in Python.

* Python (>= 3.6.9): https://www.python.org/downloads/source/
* pip3 (>= 21.2.4), numpy (>= 1.19.5), matplotlib (>= 3.3.4), pandas (>= 1.1.5), scikit-learn (>= 0.24.1), keras (>= 2.6.0), tensorflow (>= 2.6.0), tf-keras-vis (>= 0.8.0)

## Running codes
### Overview
> PCASA is a three-step program.
* Step 1. Prioritization of the single genes well-classifying tumor and normal cells.
* Step 2. Calculation of the expressing cell fraction (ECF) for each gene combination.
* Step 3. Prioritization of the gene combinations well-classifying tumor and normal cells.

### Input files
> The first command for each step requires input data.
* For step 1, a tab-delimited file containing cell-by-gene sparse matrix composed of log-transformed counts with a column showing the binary class (Tumor, 1; Normal, 0).
* For steps 2 & 3, a two-column tab-delimited file containing cell-code and cell-type.
* For steps 2 & 3, a tab-delimited file containing gene-by-cell sparse matrix composed of log-transformed counts.

### Step1. RF - Cell classifier for single genes
```
cd code

Rscript step-1a__random_forest.R ../data/input-1__scrna_class.txt
# 'input-1__scrna_class.txt': Cell-by-Gene-Matrix-with-Class

python step-1b__random_forest.py
```
### Step2. Expression logic evaluator for ECF
```
python step-2__gate_coverage_calc.py \
../data/input-2a__scrna_annotation.txt \
../data/input-2b__scrna_gc-matrix.txt
# 'input-2a__scrna_annotation.txt' : Cell-Type-Annotation
# 'input-2b__scrna_gc-matrix.txt' : Gene-by-Cell-sparse-Matrix
```
### Step3. CNN - Cell classifier for gene combinations
```
python step-3a__cnn_gradcam.py \
../data/input-2a__scrna_annotation.txt \
../data/input-2b__scrna_gc-matrix.txt
# 'input-2a__scrna_annotation.txt' : Cell-Type-Annotation
# 'input-2b__scrna_gc-matrix.txt' : Gene-by-Cell-sparse-Matrix

python step-3b__cnn_gradcam.py
python step-3c__cnn_gradcam.py
```

## Pre-process
> The directory named 'preprocess' contains the codes for pre-processing.
* Initial filtering, normalization, and dimensional reduction by Seurat
* Reference-based cell-type annotation by SingleR
* Aneuploid cell prediction by CopyKat
* Volume-dependent sub-sampling by Geosketch
* Batch correction by BBKNN
* Geosketch and BBKNN codes are based on the module named 'scjp' that Dr.Park made manually.
