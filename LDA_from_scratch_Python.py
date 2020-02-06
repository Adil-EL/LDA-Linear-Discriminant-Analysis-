# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 11:10:27 2020

@author: Adil
"""

# The Linear Discriminant Analysis from scratch

# The required libraries

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#%matplotlib inline

# Import the data

columns = ["var","skewness","curtosis","entropy","class"]
df = pd.read_csv("http://archive.ics.uci.edu/ml/machine-learning-databases/00267/data_banknote_authentication.txt",index_col=False, names = columns)

# Checking the data 

f, ax = plt.subplots(1, 4, figsize=(10,3))
vis1 = sns.distplot(df["var"],bins=10, ax= ax[0])
vis2 = sns.distplot(df["skewness"],bins=10, ax=ax[1])
vis3 = sns.distplot(df["curtosis"],bins=10, ax= ax[2])
vis4 = sns.distplot(df["entropy"],bins=10, ax=ax[3])
f.savefig('subplot.png')

sns.pairplot(df, hue="class")

def My_LDA(X):
    
    # Data preprocessing
    n=X.shape[0]
    p=X.shape[1]
    
    Features=X[:,:p]
    classVector=X["class"]
    class1index=classVector=="class1"
    class2index=classVector=="class2"
    Class1=X[class1index,:]
    Class2=X[class2index,:]
    
    n1=Class1.shape[0]
    n2=Class2.shape[0]

    
    # Compute d-dimensional mean vectors for different classes
    m1=np.mean(Class1,axis=0)
    m2=np.mean(Class2,axis=0)
    m=np.mean(Features,axis=0)
    
    # Compute in-between class and with-in class scatter matrices
   
    B=np.dot((m1-m),(m1-m).T)*n1 +np.dot((m2-m),(m2-m).T)*n2
    B=B/n
    
    W1=np.cov(Class1)*(n1-1)/n1    
    W2=np.cov(Class2)*(n2-1)/n2   
  
    W=n1*W1+n2*W2
    W=W/n
    
    
    # Compute eigen values and correspanding eigen vectors for the scatter matrices
    
    lmd, Vect = np.linalg.eig(np.linalg.inv(W).dot(B))
    indices=np.argsort(lmd)[::-1]
    Vect=Vect[:,indices]
    
    
    # Transform the feature space
    X_lda = X.dot(Vect)
    
    return X_lda
   
    
    
    
    
    
    
    


    