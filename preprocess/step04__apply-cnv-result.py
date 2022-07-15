####

# CNV result

epi_cnv_d = dict()

f = open("./temp_result/test_copykat_prediction.txt","r")
f.readline()
for line in f:
    comp = line.strip().split('\t')
    cell = comp[0].strip()
    clas = comp[1].strip()
    if clas == 'aneuploid':
        epi_cnv_d[cell] = 'tumor'
f.close()

####

f = open("./temp_result/scrna.all-annot.txt","r")

line_1st = f.readline().strip()
line_ext = list()
line_ext2 = list()

temp_rem = dict()
temp_tum = dict()

for line in f:
    comp = line.strip().split('\t')
    cell = comp[0]
    c_type = comp[4]
    ##
    if c_type == 'Epithelial_cells':
        if cell in epi_cnv_d:
            temp_tum[cell] = ''
            line_ext.append('\t'.join(comp[:4]+['Tumor_cells']+comp[5:]))
            line_ext2.append('\t'.join(comp[:4]+['Tumor']+comp[5:]))
        else:
            line_ext.append(line.strip())
            line_ext2.append('\t'.join(comp[:4]+['Normal']+comp[5:]))
    else:
        line_ext.append(line.strip())
        line_ext2.append('\t'.join(comp[:4]+['Normal']+comp[5:]))
f.close()

line_ext.insert(0, line_1st)
f_out = open("./temp_result/scrna.all-annot.t-lab.txt","w")
f_out.write('\n'.join(line_ext))
f_out.close()

line_ext2.insert(0, line_1st)
f_out = open("./temp_result/scrna.all-annot.t-lab.bi.txt","w")
f_out.write('\n'.join(line_ext2))
f_out.close()

##

f = open("./temp_result/scanpy_input_ann.tsv","r")

line_1st = f.readline().strip()
line_ext1 = list()
line_ext2 = list()
line_ext3 = list()

for line in f:
    comp = line.strip().split('\t')
    cell = comp[0]
    if cell in temp_rem: ####!####
        continue
    elif cell in temp_tum:
        line_ext1.append(cell+'\t'+'Tumor_cells'+'\t'+'\t'.join(comp[2:]))
        line_ext3.append(cell+'\t'+'Tumor_cells'+'\t'+'\t'.join(comp[2:]))
    else:
        line_ext2.append(line.strip())
        line_ext3.append(line.strip())
f.close()
        
line_ext1.insert(0, line_1st)
f_out = open("./temp_result/scanpy_input_ann_tum.tsv","w")
f_out.write('\n'.join(line_ext1))
f_out.close()

line_ext2.insert(0, line_1st)
f_out = open("./temp_result/scanpy_input_ann_nor.tsv","w")
f_out.write('\n'.join(line_ext2))
f_out.close()

####

f = open("./temp_result/scanpy_input_cnt.tsv","r")

line_1st = f.readline().strip()
line_ext1 = list()
line_ext2 = list()
line_ext3 = list()

for line in f:
    comp = line.strip().split('\t')
    cell = comp[0]
    if cell in temp_rem:
        continue
    elif cell in temp_tum:
        line_ext1.append(line.strip())
        line_ext3.append(line.strip())
    else:
        line_ext2.append(line.strip())
        line_ext3.append(line.strip())
f.close()

line_ext1.insert(0, line_1st)
f_out = open("./temp_result/scanpy_input_cnt_tum.tsv","w")
f_out.write('\n'.join(line_ext1))
f_out.close()

line_ext2.insert(0, line_1st)
f_out = open("./temp_result/scanpy_input_cnt_nor.tsv","w")
f_out.write('\n'.join(line_ext2))
f_out.close()

####

