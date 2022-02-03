########################################################
#### Packages ####

import datetime
import time
import sys
from itertools import combinations

start = time.time()


########################################################
#### Prepare RF-result and 2 inputs ####

input = sys.argv
input1 = input[1]
input2 = input[2]
print(input1,input2)

surf_goi = [line.strip().split(' ')[0][1:-1] for line in open("../results/result-1__rf_var-imp.txt")][1:101]

#f = open("input-2a__scrna_annotation.txt","r")
f = open(input1,"r")
f.readline()
d_all_tum = dict()
d_all_nor = dict()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    anno = comp[1].strip()
    if 'Tumor' in anno:
        d_all_tum[cell] = ''
    if 'Normal' in anno:
        d_all_nor[cell] = ''
f.close()

#f = open("input-2b__scrna_gc-matrix.txt","r")
f = open(input2,"r")
line_1st = f.readline().strip().split('\t')
c_tum = list()
c_nor = list()
for i in range(len(line_1st)):
    if i != 0:
        if line_1st[i] in d_all_tum:
            c_tum.append(i)
        if line_1st[i] in d_all_nor:
            c_nor.append(i)

total_tum = float(len(c_tum))
total_nor = float(len(c_nor))

d_all = dict()

for line in f:
    comp = line.strip().split('\t')
    gene = comp[0].strip()
    if gene in surf_goi:
        temp_tum = list()
        temp_nor = list()
        cnt_tum = 0.0
        cnt_nor = 0.0
        for i in c_tum:
            val = float(comp[i])
            temp_tum.append(val)
            if val != 0.0:
                cnt_tum += 1.0
        for i in c_nor:
            val = float(comp[i])
            temp_nor.append(val)
            if val != 0.0:
                cnt_nor += 1.0
        d_all[gene] = [temp_tum,temp_nor]
f.close()

print("Job1 : Done\n")


########################################################
#### Calculation of Coverages for AND- and OR-gates ####

surf_goi.sort()
c2 = list(combinations(surf_goi,2))
print("all combi count : ",len(c2),"\n")

line_ext_1st = 'Gene_combination\tCoverage_tumor\tCoverage_normal'
line_ext_and = list()
line_ext_or  = list()
all_count = 0
for pair in c2:
    all_count += 1
    if all_count % 2000 == 0:
        print(all_count)
    gene1 = pair[0].strip()
    gene2 = pair[1].strip()
    p_name = gene1+'_'+gene2
    ##
    if gene1 not in d_all or gene2 not in d_all:
        continue
    temp_tum1 = d_all.get(gene1)[0]
    temp_nor1 = d_all.get(gene1)[1]
    temp_tum2 = d_all.get(gene2)[0]
    temp_nor2 = d_all.get(gene2)[1]
    ##
    cnt_tum_and = 0.0
    cnt_nor_and = 0.0
    cnt_tum_or = 0.0
    cnt_nor_or = 0.0
    for i in range(int(total_tum)):
        if temp_tum1[i] != 0.0 and temp_tum2[i] != 0.0:
            cnt_tum_and += 1.0
            cnt_tum_or  += 1.0
        if temp_tum1[i] != 0.0 or temp_tum2[i] != 0.0:
            cnt_tum_or  += 1.0
    for i in range(int(total_nor)):
        if temp_nor1[i] != 0.0 and temp_nor2[i] != 0.0:
            cnt_nor_and += 1.0
        if temp_nor1[i] != 0.0 or temp_nor2[i] != 0.0:
            cnt_nor_or += 1.0
    prop_tum_and = str(cnt_tum_and/total_tum)
    prop_tum_or  = str(cnt_tum_or/total_tum)
    prop_nor_and = str(cnt_nor_and/total_nor)
    prop_nor_or  = str(cnt_nor_or/total_nor)
    ##
    line_ext_and.append([p_name+'\t'+prop_tum_and+'\t'+prop_nor_and, p_name, float(prop_nor_and), float(prop_tum_and)])
    line_ext_or.append([p_name+'\t'+prop_tum_or+'\t'+prop_nor_or, p_name, float(prop_nor_or), float(prop_tum_or)])

print("Job2 : Done\n")

line_final = list()
line_ext_sort = sorted(line_ext_and, key=lambda line_l: line_l[-1], reverse=True)
for line_l in line_ext_sort:
    line_final.append(line_l[0])

line_final.insert(0, line_ext_1st)
f_out = open("../results/result-2__coverage_gate-and.txt","w")
f_out.write('\n'.join(line_final))
f_out.close()

line_final = list()
line_ext_sort = sorted(line_ext_or, key=lambda line_l: line_l[-1], reverse=True)
for line_l in line_ext_sort:
    line_final.append(line_l[0])

line_final.insert(0, line_ext_1st)
f_out = open("../results/result-2__coverage_gate-or.txt","w")
f_out.write('\n'.join(line_final))
f_out.close()
        

########################################################
#### Calculation of Coverages for NOT-gates ####

from itertools import permutations
c2 = list(permutations(surf_goi,2))
print("all combi count : ",len(c2),"\n")

line_ext_1st = 'Gene_combination\tCoverage_tumor\tCoverage_normal'
line_ext_not  = list()
all_count = 0
for pair in c2:
    all_count += 1
    if all_count % 2000 == 0:
        print(all_count)
    gene1 = pair[0].strip()
    gene2 = pair[1].strip()
    p_name = gene1+'_'+gene2
    ##
    if gene1 not in d_all or gene2 not in d_all:
        continue
    temp_tum1 = d_all.get(gene1)[0]
    temp_nor1 = d_all.get(gene1)[1]
    temp_tum2 = d_all.get(gene2)[0]
    temp_nor2 = d_all.get(gene2)[1]
    ##
    cnt_tum_not = 0.0
    cnt_nor_not = 0.0
    for i in range(int(total_tum)):
        if temp_tum1[i] != 0.0 and temp_tum2[i] == 0.0:
            cnt_tum_not  += 1.0
    for i in range(int(total_nor)):
        if temp_nor1[i] != 0.0 and temp_nor2[i] == 0.0:
            cnt_nor_not += 1.0
    prop_tum_not = str(cnt_tum_not/total_tum)
    prop_nor_not = str(cnt_nor_not/total_nor)
    ##
    line_ext_not.append([p_name+'\t'+prop_tum_not+'\t'+prop_nor_not, p_name, float(prop_nor_not), float(prop_tum_not)])

print("Job3 : Done\n")

line_final = list()
line_ext_sort = sorted(line_ext_not, key=lambda line_l: line_l[-1], reverse=True)
for line_l in line_ext_sort:
    line_final.append(line_l[0])

line_final.insert(0, line_ext_1st)
f_out = open("../results/result-2__coverage_gate-not.txt","w")
f_out.write('\n'.join(line_final))
f_out.close()


####

end = time.time()
sec = (end - start)
result = datetime.timedelta(seconds=sec)
result_list = str(datetime.timedelta(seconds=sec)).split(".")
print(result_list[0])

