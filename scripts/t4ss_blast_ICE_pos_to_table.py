#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import re 
from Bio import SeqIO
from Bio.Seq import Seq
import zipfile


# In[2]:


__SAVE_MOBS__ = False

insertion_site = 'TTTTTTAT'
insertion_site_RC = str(Seq(insertion_site).reverse_complement())

path2hits = 'rssc_ont_blast_ICE_regions_best_hits.csv' # ICE blast output 
complete_genomes = 'complete_MAR_NCBI.txt' # list of complete genomes names 
genomes_ann_p = 'annotations/' # directory containing annotations of all genomes 

mobscan_p = 'MOBscan/MOBscan_results/'
mobscan_out = mobscan_p.replace('MOBscan_results/', 'MOBscan_results_all.csv')
mobscan_tab = mobscan_p.replace('MOBscan_results/', 'MOBscan_results_TAB.csv')

out_tab = path2hits.replace("regions_best_hits.csv", "T4SS_regions.csv")


# In[3]:


# strains not shown because of the 2 t4ss
not_SHOWN_strains = ['RS10', 'UW386']

strains_MOB_and_T4SS = ['YC40-M', 'SL3175', 'T98']

to_reverse = ['RUN0820', 'RS24', '362200', 'HA4-1', 'PeaFJ1', 'Wj644', 'FJAT-91', 'RS10']

def reverese_hits_start(strain, slen, sstart, send):
    if not strain in to_reverse:
        return sstart
    else:
        return slen-send

def reverese_hits_end(strain, slen, sstart, send):
    if not strain in to_reverse:
        return send
    else:
        return slen-sstart


# In[4]:


motif_dict = {}

f = open(complete_genomes, 'r')
for l in f.readlines():
    strain = l.replace('\n', '')
    genome_p = genomes_ann_p + strain + '_annot/' + strain + '.fna'
    motif_dict[strain] = []
    
    for record in SeqIO.parse(genome_p, "fasta"):
        if "contig_1" in record.id:
            motif_occurences = [m.start() for m in re.finditer(insertion_site, str(record.seq))]
            motif_occurences = motif_occurences + [m.start() for m in re.finditer(insertion_site_RC, str(record.seq))]
            for i in motif_occurences:
                motif_dict[strain].append(i)
                #motif_dict['strain'].append(strain)
                #motif_dict['sseqid'].append('contig_1')
                #motif_dict['pos'].append(str(i))
                #motif_dict['hit_ID'].append(strain + 'motif_' + str(i))
                #motif_dict['region'].append('IME predicted by ICEfinder')
f.close()

for s in motif_dict:
    motif_dict[s].sort()
#motif_dict
#motif_df = pd.DataFrame.from_dict(motif_dict)
#motif_df


# In[5]:


mobscan_list = []
for file in os.listdir(mobscan_p):
    
    archive = zipfile.ZipFile(mobscan_p+file, 'r')
    res_file = archive.open('results_60.csv')

    d = pd.read_csv(res_file, sep='\t')
    if len(d) > 0:
        mobscan_list.append(d)
mobscan_df = pd.concat(mobscan_list, ignore_index=True).set_index('Query name')
mobscan_df.to_csv(mobscan_out)


# In[6]:


set(mobscan_df['Relaxase MOB family'])


# In[7]:


if __SAVE_MOBS__:
    mob_dict = {'gene_id':[], 'strain':[], 'sseqid':[], 'pos':[], 'hit_ID':[], 'region':[]}
    
    f = open(complete_genomes, 'r')
    for l1 in f.readlines():
        strain = l1.replace('\n', '')
        genome_p = genomes_ann_p + strain + '_annot/' + strain + '.gff3'
        
        # open the gff file 
        gen_f = open(genome_p, 'r')
        for l in gen_f.readlines():
            sl = l.split('\t')
            if l.startswith('contig_') and sl[2] == 'CDS':
    
                # get MOBs positions 
                gene_id = sl[-1].split(';')[0].replace('ID=', '')
                if gene_id in mobscan_df.index:
                    #for i in [sl[3], sl[4]]:
                    mob_dict['gene_id'].append(gene_id)
                    mob_dict['strain'].append(strain)
                    mob_dict['sseqid'].append(sl[0])
                    mob_dict['pos'].append(sl[3])
                    mob_dict['hit_ID'].append(strain + '_mob_start_' + sl[3])
                    if strain in not_SHOWN_strains + strains_MOB_and_T4SS:
                        mob_dict['region'].append('MOB relaxase with T4SS')
                    else:
                        mob_dict['region'].append('MOB relaxase')
    
        gen_f.close()     
    f.close()
    
    mob_df = pd.DataFrame.from_dict(mob_dict)
    mob_df = mob_df.set_index('gene_id')
    
    ##### Fix RS10 ###
    mob_df.at['RS10_05810', 'region'] = 'MOB relaxase'
    ##################
    
    mob_df    



# In[8]:


hits = pd.read_csv(path2hits)
hits = hits.loc[(hits['region'] == 'Soil associated genes') | (hits['region'] == 'Martinique associated genes')]
hits


# In[9]:


#abs(-10)


# In[10]:


t4ss_dict = {'molecule':[], 'start':[], 'end':[], 'gene_type':[], 'total_range':[], 'size': []}

f = open(complete_genomes, 'r')
for l1 in f.readlines():
    strain = l1.replace('\n', '')
    genome_p = genomes_ann_p + strain + '_annot/' + strain + '.gff3'

    #if not strain in not_SHOWN_strains:

    # get the genomes range including the T4SS blast hits 
    hits_for_strain = hits.loc[hits['strain'] == strain]
    if not len(hits_for_strain) == 0:
        t4ss_start = hits_for_strain['pos'].min()-1
        t4ss_end = hits_for_strain['pos'].max()+1

        motifs = motif_dict[strain]
        motif_start = max([i for i in motifs if (i < t4ss_start)])
        motif_end = min([i for i in motifs if (i> t4ss_end)])

        # search repeats position in 10kb range around t4ss
        if abs(t4ss_start-motif_start)>10000 or abs(t4ss_end-motif_end)>10000:
            start = t4ss_start
            end = t4ss_end
        else:
            start = motif_start
            end = motif_end
            offset_for_plotting = 300
            plot_i = 1
            for i in [start, end]:
                t4ss_dict['molecule'].append(strain)
                t4ss_dict['total_range'].append(end-start)
                t4ss_dict['start'].append(i - start - offset_for_plotting * plot_i)
                t4ss_dict['end'].append(i - start + len(insertion_site) - offset_for_plotting * plot_i)
                t4ss_dict['gene_type'].append('TTTTTTAT repeat')
                t4ss_dict['size'].append(1.3)
                plot_i = -1

        # set Mar genes positions
        hits_for_strain_mar = hits_for_strain.loc[hits_for_strain['region'] == 'Martinique associated genes']
        if not len(hits_for_strain_mar) == 0:
            mar_start = hits_for_strain_mar['pos'].min()-100
            mar_end = hits_for_strain_mar['pos'].max()+100
        else:
            mar_start = float("inf")
            mar_end = -float("inf")

        # open the gff file and search for the range 
        gen_f = open(genome_p, 'r')
        first = True
        for l in gen_f.readlines():
            if l.startswith('contig_1'):
                l = l.split('\t')

                # get T4SS blast genes
                #if int(l[4]) <= t4ss_end and int(l[3]) >= t4ss_start and l[2] == 'CDS':
                if int(l[4]) <= end and int(l[3]) >= start and l[2] == 'CDS':
                    #if strain == 'GMI1000':
                      #  print(motif_start, motif_end)
                    
                    t4ss_dict['molecule'].append(strain)
                    t4ss_dict['total_range'].append(end-start)
                    
                    if (l[6] == '+' and not strain in to_reverse) or (l[6] == '-' and strain in to_reverse):
                        t4ss_dict['start'].append(int(l[3])-start)
                        t4ss_dict['end'].append(int(l[4])-start)
                    else:
                        t4ss_dict['end'].append(int(l[3])-start)
                        t4ss_dict['start'].append(int(l[4])-start)

                    # get Gene name positions 
                    gene_id = l[-1].split(';')[0].replace('ID=', '')
                    if gene_id in mobscan_df.index:
                        print(mobscan_df.loc[gene_id]['Relaxase MOB family'])
                        t4ss_dict['gene_type'].append('MOBP relaxase')
                        if __SAVE_MOBS__:
                            mob_df.at[gene_id,'region'] = 'MOB relaxase with T4SS'
                            
                    elif 'Trb' in l[-1] or 'VirB' in l[-1]:
                        t4ss_dict['gene_type'].append('Trb proteins')
                    #elif 'ranscriptional regulator' in l[-1]:
                    #    t4ss_dict['gene_type'].append('Transcriptional regulator')
                    elif 'TraG' in l[-1]: 
                        t4ss_dict['gene_type'].append('TraG proteins')
                    elif 'integrase' in l[-1]:
                        t4ss_dict['gene_type'].append('Integrase') 
                    elif int(l[4]) <= mar_end and int(l[3]) >= mar_start:
                        t4ss_dict['gene_type'].append('Martinique associated genes')
                    else:
                        t4ss_dict['gene_type'].append('Other') 
                    
                    t4ss_dict['size'].append('NA')
        gen_f.close()
            
f.close()

t4ss_df = pd.DataFrame.from_dict(t4ss_dict)
t4ss_df


# In[11]:


t4ss_df["start_nonrev"] = t4ss_df["start"]
t4ss_df["end_nonrev"] = t4ss_df["end"]
t4ss_df["start"] = t4ss_df.apply(lambda x: reverese_hits_start(x.molecule, x.total_range, x.start_nonrev, x.end_nonrev), axis=1)
t4ss_df["end"] = t4ss_df.apply(lambda x: reverese_hits_end(x.molecule, x.total_range, x.start_nonrev, x.end_nonrev), axis=1)
t4ss_df = t4ss_df.drop(["start_nonrev","end_nonrev"], axis=1)


# In[12]:


if __SAVE_MOBS__:
    mob_df.to_csv(mobscan_tab, index=False)

#for s in not_SHOWN_strains:
t4ss_df = t4ss_df.loc[~(t4ss_df['molecule'] == 'UW386') | ~(t4ss_df['start'] <59200)]
t4ss_df = t4ss_df.loc[~(t4ss_df['molecule'] == 'RS10') | ~(t4ss_df['start'] <84000)]
print(t4ss_df)

t4ss_df.to_csv(out_tab, index=False)
t4ss_df



