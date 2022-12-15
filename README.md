# COMP479 - GTEx Sample Classifier

## Introduction

While DNA is identical across all of the cells in an individual’s body and rarely changes, RNA levels are highly variable, changing with external factors over time and in different environments. Different human tissues are the result of variable RNA expression across an individual’s body; cells in different parts of the body need to perform different functions, and thus express genes in different amounts and will have different RNA levels. For this reason, we have decided to see if RNA expression is an effective predictor of the human organ that cells are from. 

For this analysis, we decided to use gene expression data produced by the Gentoype-Tissue Expression (GTEx) project. This is a public resource that has worked to collect cell samples from different parts of the human body and measure gene expression levels. These data were produced with a method called RNA-Seq. This method performs sequencing on all RNA transcripts collected from a cell sample, thus helping to identify exactly which genes were expressed in the sample. In order to make a classifier for the GTEx RNA-Seq dataset, we are using logistic regression. Logistic regression is a supervised machine learning method that, based on certain variables, tries to classify observations into distinct categorical classes. Often those outcome categorical classes (“labels”) are binary, but logistic regression can be modified to also work on multiclass classification. Thus, as the dataset we decided to work with has 30 distinct labels, we sought to use it to evaluate its applicability. 

This is a final project developed for the COMP479 Machine Learning course at Loyola University Chicago. 

## Requirements

In this GitHub repo, you can find the main script (`main_script.py`) and some test data. For the full data, you may download it directly from GTEx Portal in [here](https://gtexportal.org/home/datasets) - files `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz` and `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`. Besides those files, the general requirements to run the main script are:

* Python (3.6 or later)
* Pandas (1.3.4 or later)
* Matplotlib (3.4.3 or later)
* Numpy (1.20.3 or later)
* Scikit-learn (0.24.2 or later)
* Seaborn (0.11.2 or later)
* Umap-learn (0.5.3 or later)

## Implementation

The GTEx RNA-Seq dataset contains measured expression levels for different genes across samples from 30 human tissues. However, not all genes are differentially expressed across distinct tissues. Thus, our first step to reduce our dataset was to filter out genes with variance less or equal than 10 TPMs. This step guarantees that we would only be working with genes whose expression greatly vary among our samples, meaning that they could be optimal parameters for a classifier.

Next, we further sought to reduce the number of variables of our dataset. To do that, we computed average TPM levels for each gene across all samples, and selected the top 100. This step assures that the predictor variables in our classifier are highly expressed genes. This is important because, as with any data, background noise can be an issue. Consequently, by selecting the top 100 expressed genes, we are confident that the read counts (TPMs) are above background noise. 

Lastly, we randomly split our data into training and test datasets (85% and 15%, respectively). With the training dataset, we built and optimized a logistic regression classifier. Applying 5-fold cross validation, we implemented a grid search that evaluated performance of our classifier in order to identify the best values for two hyperparameters: i) using elastic-net, we wanted to identify the best mixing parameter for our model, and thus we tested distinct L1 ratio values ranging from 0.0 to 1.0, with a 0.1 step; ii) additionally, we tested different regularization strength values (C), ranging from 0.1 to 1.0, with a 0.1 step. We used the best hyperparameters in our classifier to evaluate its performance in the test dataset. 

## Results

For the results, please read the final project report included in this repo (`project_report.pdf`).

## Credits

Daniel Araujo, Olaf Garcia, Henry Wittich.
