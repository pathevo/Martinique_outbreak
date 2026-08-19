#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import re


# In[2]:


fna = snakemake.input['fna'] # genome fasta input 
bed_out = snakemake.output['bed'] # bed file with kmers positions 
kmers_p = snakemake.params['kmers'] # kmers list 


# In[3]:


fna_seqs = []
f = open(fna, 'r')
for l in f.readlines():
    if not l.startswith('>'):
        fna_seqs.append(l)

f.close()


# In[4]:


for l in fna_seqs:
    print(len(l))


# In[5]:


def traits_names(t):
    return t.replace('pyseer_', '').replace('_kmers.txt', '') + ' '


# In[6]:


kmers = pd.read_csv(kmers_p)
kmers = kmers[['variant', 'TRAITS']]
kmers['TRAITS'] = kmers['TRAITS'].apply(traits_names)
kmers = kmers.groupby('variant').sum()#.to_csv(kmers_p_out)


# In[7]:


f = open(bed_out, 'w')
contigs_counter = 0 
for l in fna_seqs:
    contigs_counter = contigs_counter + 1
    for k in kmers.index:
        if k in l:
            trait = kmers.loc[k]['TRAITS']
            occurences = [m.start() for m in re.finditer(k, l)]
            for o in occurences:
                end = str(o + len(k))
                f.write('contig_'+str(contigs_counter)+'\t'+str(o)+'\t'+str(end)+'\t'+trait+'\n')


f.close()





