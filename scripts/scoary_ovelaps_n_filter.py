#!/usr/bin/env python
# coding: utf-8

# In[25]:


import numpy as np
import pandas as pd
import os


# In[26]:


sample = "rssc"
date_suffix2remove = '_30_01_2024_0947.results.csv'

scoary_results_path = 'scoary_analysis/rssc_scoary_hits/'
scoary_filtered_path = 'scoary_analysis/rssc_scoary_filtered/'

out = sample+'_olaps.xlsx'
upset_out = out.replace('_olaps.xlsx', '_upset.csv')


# In[27]:


# Nested traits! remove if identical (Identified manually)
traits2remove = [
    'host_spp:Nicotiana tabacum',
    'host_spp:Zingiber officinale',
    'host_spp:Cucumis sativus',
    'Clade:Environmental',
    'order:Environmental',
    'host_spp:soil_banana_fallow',
    'order:Solanales',
    'host_genus:Zingiber',
    'order:Alismatales',
    'host_genus:Anthurium',
    'host_genus:Water',
    'order:Solanales',
    'country:French Guiana'
]


# In[28]:


def file2index(f):
    return f.replace('_of_isolation_', ':').replace('_10_09_2024_1313.results.csv', '').replace('_or_nearest_village_', ':')


# # Filtering scoary results based on the sensititvity and specificity 
# 
# The main goal is to identify the genes important for an infection of a particular host or for survival in a particular region or environment. That's why the *sensitivity* is especially important paramenter: the genes which are required to infect tomato must be present in ALL tomato isolates. On the other hand specificity is not so relevant, because many rso strains are able to infect multiple hosts, so the genes important for tomato virulence can be present in musa isolate. For these reasons the following filtering thresholds are applied:
# 
# + sensitivity over 95%
# + specificity over 70%

# In[29]:


directory = os.fsencode(scoary_results_path)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".csv"):
        df = pd.read_csv(scoary_results_path + filename)
        before = len(df)
        df = df.loc[df['Sensitivity'] > 94]
        df = df.loc[df['Specificity'] > 69]
        if len(df)>0 and not file2index(filename) in traits2remove:
            print(filename.replace(date_suffix2remove, ''))
            print('Before', before)
            df.to_csv(scoary_filtered_path + filename)
            print('After', len(df))
            print()


# # Nested ovelaps 
# 
# Remove a trait from the table if overlap of the genes between 2 traits is equal to one of the traits genes set. The trait with the smaller number of genes gets removed. If the both traits have equal number of genes: remove the trait with the lower taxonomy rank

# In[30]:


sets_for_upset = {}
olaps = {}
directory = os.fsencode(scoary_filtered_path)

for file1 in os.listdir(directory):
    filename1 = os.fsdecode(file1)
    if filename1.endswith(".csv"):
        df1 = pd.read_csv(scoary_filtered_path + filename1)
        index1 = file2index(filename1)
        
        if len(df1) > 0 and not index1 in traits2remove:
            olaps[index1] = {}
            
            for file2 in os.listdir(directory):
                filename2 = os.fsdecode(file2)
                index2 = file2index(filename2)
                
                if filename2.endswith(".csv") and not index2 in traits2remove:
                    df2 = pd.read_csv(scoary_filtered_path + filename2)
        
                    if len(df2) > 0:
                        # compare pairs
                        genes1 = set(df1["Gene"].values)
                        genes2 = set(df2["Gene"].values)
                        intersec = genes1.intersection(genes2)
                        olaps[index1][index2] = len(intersec)

                        # check for nested overlap
                        if (intersec == genes1) and index1 != index2:
                            print('overlap', len(intersec))
                            print(index1, len(genes1))
                            print(index2, len(genes2))
                            print()
            sets_for_upset[index1] = set(df1["Gene"].values)
#sets_for_upset                


# In[31]:


olaps_df = pd.DataFrame.from_dict(olaps).sort_index().sort_index(axis=1)
olaps_df


# In[32]:


olaps_df.to_excel(out)


# In[33]:


sets_for_upset_updated = {}

for s in sets_for_upset:
    sets_for_upset_updated[s] = []
    for g in sets_for_upset[s]:
        if '--' in g:
            g = g.split('--')
            sets_for_upset_updated[s] = sets_for_upset_updated[s] + g
        else:
            sets_for_upset_updated[s] = sets_for_upset_updated[s] + [g]
#sets_for_upset_updated


# In[34]:


dummies = pd.Series(sets_for_upset_updated).str.join('|').str.get_dummies().T
pd.DataFrame.to_csv(dummies, upset_out)


# In[ ]:




