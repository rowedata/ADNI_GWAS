# -*- coding: utf-8 -*-
"""
Created on Thu May 21 12:14:57 2020

@author: Benazir
"""
import numpy as np
import pandas as pd

import os
df = pd.read_csv('Documents/dissertation/data/ADNI/phenotypes.csv')


fam = pd.read_table('Documents/dissertation/data/ADNI/ADNI.fam',delimiter=' ',header=None)
fam.shape
fam.columns = ['index','Subject','bla', 'bla2','bla3','bla4']

new =fam.merge(df)
new.shape
skim =new.drop(columns=['index', 'bla','bla2','bla3', 'bla4'])
skim.to_csv('Documents/dissertation/data/ADNI/adni_subject_pheno.txt', index=False)

pure_pheno =skim.drop(columns='Subject')
pure_pheno.to_csv('Documents/dissertation/data/ADNI/adni_pheno.txt', index=False,header=None)

dementia = pure_pheno.replace({'CN':0,'Dementia':1,'MCI':np.nan })
mci =pure_pheno.replace({'CN':0,'Dementia':np.nan,'MCI':1 })
both = pure_pheno.replace({'CN':0,'Dementia':1,'MCI':1 })

dementia.to_csv('Documents/dissertation/data/ADNI/dementia.txt', index=False,header=None)
mci.to_csv('Documents/dissertation/data/ADNI/mci.txt', index=False,header=None)
both.to_csv('Documents/dissertation/data/ADNI/both.txt', index=False,header=None)
