##########################################################
#### Pacakge ####

import datetime
import time
import sys
import os
import numpy as np

start = time.time()


##########################################################
#### Calculation of z-scores for each Grad-CAM result ####

add0 = '../results/combi_cnn/gcam/'
f_list0 = os.listdir(add0)

for j in range(10):
    check = str(j+1)

    if 'gcam_'+check not in f_list0:
        continue

    add = '../results/combi_cnn/gcam/gcam_'+check+'/'
    all_d = dict()
    f_list = os.listdir(add)
    out = 0
    for f_name in f_list:
        f = open(add+f_name,"r")
        genes = f.readline().strip().split('\t')
        comp = f.readline().strip().split('\t')[1:]
        if len(comp) < 3:
            out += 1
            continue
        f.close()
        ##
        val_all = list()
        for val in comp:
            val_all.append(float(val))
        val_max = np.max(val_all)
        val_min = np.min(val_all)
        gap = val_max-val_min
        if gap == 0:
            out += 1
            continue
        ##
        for i in range(len(genes)):
            gene = genes[i].strip()
            val = float(comp[i])
            gap = val_max-val_min
            val_mod = (val-val_min)/(val_max-val_min)
            if gene in all_d:
                all_d[gene] = all_d.get(gene)+val_mod
            else:
                all_d[gene] = val_mod
    total = len(f_list)
    print("total: ",total)
    print("out : ",out)
    total = total-out
    
    line_ext = list()
    for gene in genes:
        val = all_d.get(gene)/float(total)
        line_ext.append([gene+'\t'+str(val),gene,val])
    
    line_final = list()
    line_ext_sort = sorted(line_ext, key=lambda line_l: line_l[-1], reverse=True)
    for line_l in line_ext_sort:
        line_final.append(line_l[0])
    line_final.insert(0,'Gene_Combination\tWeight_for_CNN-prediction')
    
    print("Done",check)
    f_out_name = add0+"gcam_res_"+check+".txt"
    f_out = open(f_out_name,"w")
    f_out.write('\n'.join(line_final))
    f_out.close()


##########################################################
#### Calculation of mean values for each combination

add0 = '../results/combi_cnn/gcam/'
f_list0 = os.listdir(add0)

all_d = dict()
total = 0
for j in range(10):
    check = str(j+1)
    count = 2475

    f_in_name = "gcam_res_"+check+".txt"
    if f_in_name not in f_list0:
        continue

    total += 1
    genes = [line.strip().split('\t')[0] for line in open(add0+f_in_name)][1:count]

    for gene in genes:
        if gene in all_d:
            all_d[gene] = all_d.get(gene)+1
        else:
            all_d[gene] = 1

    all_genes = [line.strip().split('\t')[0] for line in open(add0+f_in_name)][1:]

final = list()
genes = sorted(all_d.keys())
genes_ext = list()
print("Total of GradCAM results : ",total)
for gene in genes:
    val = all_d.get(gene)
    if val == total:
        gene_mod = gene.replace(':','_')
        final.append(gene_mod)
        genes_ext.append(gene)

print('Final : ',len(genes_ext))
f_out = open(add0+"final_gcam_genes.txt","w")
f_out.write('\n'.join(final))
f_out.close()


##########################################################
#### Arrangement of results ####

all_d = dict()
rank_d = dict()
for gene in all_genes:
    all_d[gene] = []
    rank_d[gene] = []

for j in range(total):
    check = str(j+1)
    f_in_name = add0+"gcam_res_"+check+".txt"
    f = open(f_in_name,"r")
    f.readline()
    cnt = 0
    for line in f:
        cnt += 1
        comp = line.strip().split('\t')
        gene = comp[0]
        all_d[gene] = all_d.get(gene)+[float(comp[1])]
        rank_d[gene] = rank_d.get(gene)+[float(cnt)]
    f.close()

line_ext = list()
for gene in all_genes:
    vals = all_d.get(gene)
    val_mean = np.mean(vals)
    ranks = rank_d.get(gene)
    rank_mean = np.mean(ranks)
    gene_mod = gene.replace(':','_')
    line_ext.append([gene_mod+'\t'+str(val_mean)+'\t'+str(rank_mean),gene,rank_mean,val_mean])

line_final = list()
line_ext_sort = sorted(line_ext, key=lambda line_l: line_l[-1], reverse=True)
for line_l in line_ext_sort:
    line_final.append(line_l[0])

line_final2 = list()
temp_d = dict()
for line in line_final:
    comp = line.strip().split('\t')
    combi = '_'.join(sorted(comp[0].split('_')))
    weight = comp[1]
    if combi not in temp_d:
        line_final2.append(combi+'\t'+weight)
        temp_d[combi] = ''
line_final2.insert(0, 'Gene_Symbol\tWeight_Mean')

f_out = open("../results/result-3__gradcam_weight.txt","w")
f_out.write('\n'.join(line_final2))
f_out.close()


####

end = time.time()
sec = (end - start)
result = datetime.timedelta(seconds=sec)
result_list = str(datetime.timedelta(seconds=sec)).split(".")
print(result_list[0])

