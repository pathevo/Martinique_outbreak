#!/usr/bin/env python
# coding: utf-8

# In[82]:


import pandas as pd
import numpy as np


# In[83]:


# num of RSSC strains outside of 4npd allowed to have the gene 
_CORE_FILTER_min = 10
_CORE_FILTER_max = 800

SAMPLE = '4npb'

# pangenome output 
panaroo_p = '/path/to/rssc_panaroo_gene_presence_absence.csv'

# strains lists in form of text files 
rssc = 'rssc.txt'
npb = SAMPLE+'.txt'
# output
out_core_info = SAMPLE+'_uniquely_present_genes.xlsx'
out_plot = SAMPLE+'_uniquely_present_genes_for_plot.csv'



# In[84]:


f_rssc = open(rssc, "r")
f_npb = open(npb, "r")

rssc = f_rssc.read().split("\n")[:-1]
npb = f_npb.read().split("\n")[:-1]

non_npb = list(set(rssc).difference(set(npb)))

print('RSSC strains number:',len(rssc), '4NPB strains:', len(npb), 'Non-4npb strains:', len(non_npb))
f_rssc.close()
f_npb.close()


# # load data
# 
# remove soft-core (genes present in 95%/825 of the sequences), remove cloud (present in <10 seq) 

# In[85]:


panaroo = pd.read_csv(panaroo_p, dtype=str)
#panaroo = panaroo[panaroo.columns[13:]]
panaroo = panaroo.astype({'No. isolates': 'int32'})
panaroo = panaroo[panaroo['No. isolates'] > _CORE_FILTER_min]
panaroo = panaroo[panaroo['No. isolates'] < _CORE_FILTER_max]
panaroo


# In[ ]:


# # compare presence and absence 
# 
# + present in all 4npb but absent in all rssc: all_npb - non_npb

# In[86]:


present_in_npb = panaroo.loc[panaroo.index][npb].dropna(thresh=int(len(npb)*0.95)).index

present_in_non_npb = panaroo.loc[panaroo.index][non_npb].dropna(thresh=_CORE_FILTER).index


# In[87]:


print(len(present_in_npb), len(present_in_non_npb))


# In[88]:


npb_only = list(set(present_in_npb).difference(set(present_in_non_npb)))


# In[89]:


panaroo.loc[npb_only].to_excel(out_core_info, index=False)






