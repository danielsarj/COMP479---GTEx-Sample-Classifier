import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix 
from sklearn.decomposition import PCA

# read annotation file and get only the column we need - sample and tissue name
anno_file = pd.read_csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')
anno_file = anno_file[['SAMPID', 'SMTS']]

# recode tissue ID to numeric
set_of_tissues = set(anno_file.loc[:,'SMTS'])
tissue_codes = pd.concat([pd.DataFrame({'SMTS':list(set_of_tissues)}), pd.DataFrame({'SMTS_ID':range(len(set_of_tissues))})], axis=1)
tissue_codes.to_csv('tissue_codes.txt', index = False)

# adds tissue code ID to annotation data frame
anno_file = pd.merge(tissue_codes, anno_file, how='inner', on='SMTS')

# read expression file and add a column for variance
exp_file = pd.read_csv('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', sep='\t', skiprows=2)
exp_file['variance'] = exp_file.iloc[:,2:exp_file.shape[1]].var(axis=1)

# filter by variance column
exp_file = exp_file.loc[exp_file['variance'] > 10.0]
exp_file = exp_file.drop(columns=['variance'])

# add average TPM column & filter
exp_file['avg_tpm'] = exp_file.iloc[:,2:exp_file.shape[1]].mean(axis=1)
exp_file = exp_file.nlargest(100,'avg_tpm')
exp_file = exp_file.drop(columns=['avg_tpm'])

exp_file = exp_file.drop(columns=['Name']).transpose()
exp_file.rename(columns=exp_file.iloc[0], inplace = True)
exp_file.drop(exp_file.index[0], inplace = True)

# turn row names into column so we can merge afterwards
exp_file = exp_file.rename_axis('SAMPID').reset_index()
 
# merging annotation and expression data frames by sample ID
merged_file = pd.merge(anno_file, exp_file, how='inner', on='SAMPID')

# dropping sample ID column
merged_file = merged_file.drop(columns=['SAMPID', 'SMTS'])

# shuffling the data frame and slicing it (70% for training, 15% for development, 15% for test)
merged_file = merged_file.sample(frac=1).reset_index(drop=True)
row_15 = 0
for n in range(0, merged_file.shape[0]+1):
    if round(n/merged_file.shape[0], 2) == 0.15:
        row_15 = n

training_data = merged_file.loc[row_15+1:merged_file.shape[0],]
test_data = merged_file.loc[0:row_15,]
training_data = training_data.astype('float64')
test_data = test_data.astype('float64')

# transforming datasets from pandas to numpy
training_data = training_data.to_numpy()
test_data = test_data.to_numpy()

# running GridSearchCV & and fitting data
params = {'l1_ratio': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
         'C': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]}

logistic_model = GridSearchCV(LogisticRegression(penalty='elasticnet', solver='saga', max_iter=1000), param_grid=params, scoring='accuracy', cv=5)
logistic_model.fit(training_data[:,1:], training_data[:,0].ravel())
print(logistic_model.best_params_)

# predicting on test data
predictions = logistic_model.predict(test_data[:,1:])
print(classification_report(test_data[:,0].ravel(), predictions))
predictions.to_csv('predicted_labels.txt', index = False)
test_data[:,0].to_csv('actual_labels.txt', index = False)

# running a PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(test_data[:,1:])
principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
principalDf.to_csv('testdata_firsttwoPCs.txt', index = False)
