#!/usr/bin/env python
# coding: utf-8

# In[62]:


import numpy as np
import pandas as pd
import os


# In[63]:


# pyseer output filtered by significance 
pyseer_p = 'pyseer_significant_genes_list_TOP_hits.csv'
upset_out = pyseer_p.replace('.csv', '_upset.csv')


# In[64]:


pyseer = pd.read_csv(pyseer_p, index_col=0)
pyseer['Genes'] = pyseer.index
pyseer


# In[65]:


def fix_traits(x):
    x = x.replace('pyseer_', '')
    x = x.replace('_gwas.txt', '')
    x = x.replace('_of_isolation_', ':')
    x = x.replace('_village_', '_village:')
    return x

pyseer['TRAITS'] = pyseer['TRAITS'].apply(fix_traits)
pyseer


# In[66]:


traits2remove = [
    'pyseer_Clade_of_isolation_Environmental_gwas.txt',
    'pyseer_host_genus_of_isolation_Capsicum_gwas.txt',
    'pyseer_host_spp_of_isolation_soil_banana_fallow_gwas.txt',
    'pyseer_host_spp_of_isolation_Zingiber_officinale_gwas.txt',
    'pyseer_host_spp_of_isolation_Surface_water_gwas.txt'
]

traits2remove = [fix_traits(x) for x in traits2remove]
#traits2remove
#traits2remove = []


# In[67]:


sets_for_upset = pyseer.groupby(['TRAITS'])['Genes'].apply(','.join).to_dict()
for t in sets_for_upset:
    sets_for_upset[t] = sets_for_upset[t].split(',')
len(sets_for_upset)


# In[68]:


seen = []
for t1 in sets_for_upset:
    seen.append(t1)
    for t2 in sets_for_upset:
        if t1!=t2 and not t2 in seen and not t1 in traits2remove and not t2 in traits2remove:
            set1 = set(sets_for_upset[t1])
            set2 = set(sets_for_upset[t2])
            intersec = set1.intersection(set2)
            if len(intersec) != 0 and (len(intersec) == len(set1) or len(intersec) == len(set2)):
                print(t1, t2)
                print(len(set1), len(set2), len(intersec))


# In[69]:



for t in traits2remove:
    del sets_for_upset[t]
len(sets_for_upset)


# In[70]:


dummies = pd.Series(sets_for_upset).str.join('|').str.get_dummies().T
pd.DataFrame.to_csv(dummies, upset_out)






