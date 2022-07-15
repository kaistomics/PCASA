############################################################################
#### Notice: Check the information of author, organ, organ_spec, and batch
############################################################################
author = 'OH.Lee'
organ = 'Ovary'
organ_spec = 'Ovary'
############################################################################

f = open("./temp_result/scrna.gene.txt","r")
f.readline()
line_1st = 'No.\tGene_symbol'
line_ext = list()
for line in f:
    line_ext.append(line.strip())
f.close()
line_ext.insert(0, line_1st)
f_out = open("./temp_result/scrna.gene.mod.txt","w")
f_out.write('\n'.join(line_ext))
f_out.close()

####

import time
start = time.time()
print("#### ALL steps : 4 #### (Note. Pls modify the information above for annotation)")

####

cells_tmp = [line.strip().split('\t')[0] for line in open("./temp_result/test_copykat_prediction.txt")][1:]
cells_cpk = dict()
for cell in cells_tmp:
    cells_cpk[cell] = ''

####

def head_mod(f_name):
    f = open(f_name, "r")
    line_1st = f.readline().strip().split('\t')
    check = 0
    if 'Gene_symbol' not in line_1st:
        line_1st.insert(0,'Gene_symbol')
        check = 1
    if check == 1:
        temp = ['\t'.join(line_1st)]
        for line in f:
            temp.append(line.strip())
        f.close()
        f_out = open(f_name,"w")
        f_out.write('\n'.join(temp))
        f_out.close()

def head_cell_mod(f_name):    
    f = open(f_name, "r")
    line_1st = f.readline().strip().split('\t')
    check = 0
    if 'Gene_symbol' not in line_1st:
        line_1st.insert(0,'Gene_symbol')
        check = 1
    c_list = [0]
    cells = list()
    for i in range(len(line_1st)):
        if line_1st[i] in cells_cpk:
            c_list.append(i)
            cells.append(line_1st[i])
    line_ext = ['Gene_symbol\t'+'\t'.join(cells)]
    for line in f:
        comp = line.strip().split('\t')
        temp = list()
        for i in c_list:
            temp.append(comp[i])
        line_ext.append('\t'.join(temp))
    f.close()
    f_out = open(f_name,"w")
    f_out.write('\n'.join(line_ext))
    f_out.close()

f_name = "./temp_result/scrna.seurat-norm.log1p.tsv"
head_cell_mod(f_name)
f_name = "./temp_result/scrna.seurat-raw.counts.tsv"
head_cell_mod(f_name)
print("Step 1 : Done - Copykat cells only")
print("---- time :", time.time() - start)
start = time.time()

##

import pandas as pd
from pandas import DataFrame
import numpy as np
import scanpy as sc

##

f_name = "./temp_result/scrna.seurat-raw.counts.tsv"
f_df = pd.read_table(f_name, header=0)
f_df_tp = f_df.set_index('Gene_symbol').T
f_df_tp.to_csv('./temp_result/scanpy_input_cnt.tsv', sep='\t', na_rep='NaN')

f_name = "./temp_result/scanpy_input_cnt.tsv"
head_mod(f_name)
print("Step 2 : Done - Preparing scanpy input")
print("---- time :", time.time() - start)
start = time.time()

####
####

add = './temp_result/'

adata = sc.read_text(add+'scanpy_input_cnt.tsv', delimiter=None, first_column_names=None, dtype='float')

#anno = pd.read_csv(add+'scanpy_input_ann.tsv', sep='\t')
#adata.obs = anno

feat = pd.read_csv(add+'scrna.gene.mod.txt', sep='\t')
adata.var['gene_symbols'] = list(feat['Gene_symbol'])

##

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata.raw = adata

##

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

##

sc.tl.leiden(adata)
pd.DataFrame(data=adata.obs["leiden"], index=adata.obs_names).to_csv(add+'scrna.cluster.txt',sep='\t')

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
pd.DataFrame(
    {'C'+group + '_' + key[:-1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).to_csv(add+'scrna.cluster.markers.txt',sep='\t')

sc.tl.tsne(adata, perplexity=30, learning_rate=1000, random_state=0)
pd.DataFrame(data=adata.obsm["X_tsne"], index=adata.obs_names).to_csv(add+'scrna.tsne.txt',sep='\t')

##

sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')
sc.tl.umap(adata)

##

data_mod = adata.X
DF = pd.DataFrame(data=data_mod, index=adata.obs_names, columns=adata.raw.var_names)
DF.index.name = 'Gene_symbol'
DF.to_csv(add+'scanpy_input_cnt.tsv',sep='\t')

##

adata.write(add+'Processed_scRNA.h5ad')
print("Step 3 : Done - Scanpy reanalysis for Clustering, Marker-genes, PCA, UMAP and TSNE")
print("---- time :", time.time() - start)
start = time.time()

####
####

f = open("./temp_result/scrna.cluster.txt","r")
f.readline()
clus_d = dict()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    clus = comp[-1].strip()
    clus_d[cell] = clus
f.close()

####

f = open("./temp_result/scrna.celltype.txt","r")
f.readline()
ctype_d = dict()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    ctype = comp[-2].strip().split('.')[1]
    ctype = ctype.replace(" ","_")
    if ctype[-1] != 's':
        ctype = ctype+'s'
    ctype_d[cell] = ctype
f.close()

####

f = open("./temp_result/scrna.tsne.txt","r")
f.readline()
line_1st = 'CellCode\tCluster\ttSNE.1\ttSNE.2\tCellType\tAuthor\tOrgan\tOrgan_Specific\tBatch'
line_ext = list()
line_1st_2 = 'CellCode\tCellType\tAuthor\tOrgan\tOrgan_Specific\tBatch'
line_ext_2 = list()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    if cell not in cells_cpk:
        continue
    clus = clus_d.get(cell)
    tsne1 = comp[1].strip()
    tsne2 = comp[2].strip()
    ctype = ctype_d.get(cell)
    batch = cell.split('_')[0]
    line_ext.append(cell+'\t'+clus+'\t'+tsne1+'\t'+tsne2+'\t'+ctype+'\t'+author+'\t'+organ+'\t'+organ_spec+'\t'+batch)
    line_ext_2.append(cell+'\t'+ctype+'\t'+author+'\t'+organ+'\t'+organ_spec+'\t'+batch)
f.close()

line_ext.insert(0, line_1st)
f_out = open("./temp_result/scrna.all-annot.txt","w")
f_out.write('\n'.join(line_ext))
f_out.close()

line_ext_2.insert(0, line_1st_2)
f_out = open("./temp_result/scanpy_input_ann.tsv","w")
f_out.write('\n'.join(line_ext_2))
f_out.close()

print("Step 4 : Done - Rearrangement of annotation files")
print("---- time :", time.time() - start)

####

