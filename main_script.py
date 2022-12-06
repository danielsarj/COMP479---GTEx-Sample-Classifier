import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV

# read annotation file and get only the column we need - sample and tissue name
anno_file = pd.read_csv('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t')
anno_file = anno_file[['SAMPID', 'SMTS']]

# recode tissue ID to numeric
set_of_tissues = set(anno_file.loc[:,'SMTS'])
tissue_codes = pd.concat([pd.DataFrame({'SMTS':list(set_of_tissues)}), pd.DataFrame({'SMTS_ID':range(len(set_of_tissues))})], axis=1)
tissue_codes.to_csv('tissue_codes.txt', index = False)

# adds tissue code ID to annotation data frame
anno_file = pd.merge(tissue_codes, anno_file, how='inner', on='SMTS')

# read expression file, drop first column, and transpose it
exp_file = pd.read_csv('test_exp_file.txt.gz', sep='\t', skiprows=2)
exp_file = exp_file.drop(columns=['Name']).transpose()
exp_file.rename(columns=exp_file.iloc[0], inplace = True)
exp_file.drop(exp_file .index[0], inplace = True)

# remove genes with TPM <= 0.1
exp_file = exp_file[exp_file.columns[exp_file.mean(axis=0) > 0.1]]

# remove genes with TPM variance <= 0.1
exp_file = exp_file[exp_file.columns[exp_file.var(axis=0) > 0.1]]

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

# running GridSearchCV
params = {'l1_ratio': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
         'C': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]}
         
logistic_model = GridSearchCV(LogisticRegression(penalty='elasticnet', solver='saga', max_iter=1000), param_grid=params, scoring='accuracy', cv=5)

#logistic_model = LogisticRegression(penalty='elasticnet', l1_ratio=1/10, solver='saga', C=1, max_iter=10000)

logistic_model.fit(training_data[:,1:], training_data[:,0].ravel())
logistic_model.best_params_
print(logistic_model.best_params_)









