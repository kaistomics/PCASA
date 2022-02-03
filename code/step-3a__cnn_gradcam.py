#########################################################
#### import ####

import datetime
import time
import sys
import random
import numpy as np
from itertools import combinations
from itertools import permutations

start = time.time()


#########################################################
#### Prepare input data ####

f = open("../results/result-1__rf_var-imp.txt","r")
f.readline()
rf_d = dict()
surf = list()
for line in f:
    comp = line.strip().split(' ')
    gene = comp[0].strip()[1:-1]
    rf_d[gene] = float(comp[3])
    surf.append(gene)
surf = surf[:100]
print(len(surf))
f.close()

input = sys.argv
input1 = input[1]
input2 = input[2]
print(input1,input2)

#f = open("input-2a__scrna_annotation.txt","r")
f = open(input1,"r")
f.readline()
cell_d = dict()
c_tum_d = dict()
c_nor_d = dict()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    cell_d[cell] = ''
    c_type = comp[1].strip()
    if 'Tumor' in c_type:
        c_tum_d[cell] = ''
    else:
        c_nor_d[cell] = ''
f.close()

#f = open("input-2b__scrna_gc-matrix.txt","r")
f = open(input2,"r")
line_1st = f.readline().strip().split('\t')

all_d = dict()

cells_t = sorted(c_tum_d.keys())
cells_n = sorted(c_nor_d.keys())

print('Tumor cells  : ',len(cells_t))
print('Normal cells : ',len(cells_n))

c_list = list()
for i in range(len(line_1st)):
    if line_1st[i] in cell_d:
        c_list.append(i)

genes = list()
genes2 = list()
all_vals = list()

for line in f:
    comp = line.strip().split('\t')
    gene = comp[0].strip()
    if gene not in surf:
        gene_mod = gene.replace('-','.')
        if gene_mod not in surf:
            continue
        else:
            genes2.append(gene_mod)
    else:
        genes2.append(gene)
    genes.append(gene)
    for i in c_list:
        if i != 0:
            cell = line_1st[i].strip()
            key = cell+gene
            val = comp[i].strip()
            all_d[key] = float(val)
            all_vals.append(float(val))
f.close()

for gene in surf:
    if gene not in genes2:
        print(gene)

print(len(surf))
print(len(genes))

val_max = np.max(all_vals)
#print(val_max)


#########################################################
#### Making training data ####

c2 = list(permutations(genes,2)) #****
c2 = random.sample(c2,len(c2))
c2_list = list()
for combi in c2:
    c2_list.append(combi[0]+':'+combi[1])

##

gene_row1 = list()
gene_row2 = list()
for gene_l in c2:
    gene_row1.append(gene_l[0])
    gene_row2.append(gene_l[1])

##

tum_all = list()

for cell in cells_t:

    row1 = list()
    row2 = list()

    rowA = list()
    rowB = list()

    for gene in gene_row1:
        key = cell+gene
        val = all_d.get(key)
        row1.append(val)
        if gene not in rf_d:
            gene = gene.replace('-','.')
        rowA.append(rf_d.get(gene))
    for gene in gene_row2:
        key = cell+gene
        val = all_d.get(key)
        row2.append(val)
        if gene not in rf_d:
            gene = gene.replace('-','.')
        rowB.append(rf_d.get(gene))

    row1_arr = np.array(row1)
    rowA_arr = np.array(rowA)
    row2_arr = np.array(row2)
    rowB_arr = np.array(rowB)

    tum_all.append([row1_arr,rowA_arr,row2_arr,rowB_arr])

tum_arr = np.array(tum_all)

##

nor_all = list()

for cell in cells_n:

    row1 = list()
    row2 = list()

    rowA = list()
    rowB = list()

    for gene in gene_row1:
        key = cell+gene
        val = all_d.get(key)
        row1.append(val)
        if gene not in rf_d:
            gene = gene.replace('-','.')
        rowA.append(rf_d.get(gene))
    for gene in gene_row2:
        key = cell+gene
        val = all_d.get(key)
        row2.append(val)
        if gene not in rf_d:
            gene = gene.replace('-','.')
        rowB.append(rf_d.get(gene))

    row1_arr = np.array(row1)
    rowA_arr = np.array(rowA)
    row2_arr = np.array(row2)
    rowB_arr = np.array(rowB)

    nor_all.append([row1_arr,rowA_arr,row2_arr,rowB_arr])

nor_arr = np.array(nor_all)

##

np.save("../results/training_tumor.npy",tum_arr)
np.save("../results/training_normal.npy",nor_arr)

f_out = open("../results/genes_shuffle.txt","w")
f_out.write('\n'.join(c2_list))
f_out.close()


####

end = time.time()
sec = (end - start)
result = datetime.timedelta(seconds=sec)
result_list = str(datetime.timedelta(seconds=sec)).split(".")
print(result_list[0])

