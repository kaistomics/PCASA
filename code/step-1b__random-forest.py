f_name = '../results/result-1__rf_log.txt'
f = open(f_name,"r")

line_ext = list()
for line in f:
    if 'Kappa : ' in line:
        kappa = 'KAPPA\t'+line.strip().split(' : ')[1]
    if 'Accuracy : ' in line and 'Balanced' not in line:
        accu = 'ACCURACY\t'+line.strip().split(' : ')[1]
    if 'Sensitivity : ' in line:
        sens = 'SENSITIVITY\t'+line.strip().split(' : ')[1]
    if 'Specificity : ' in line:
        spec = 'SPECIFICITY\t'+line.strip().split(' : ')[1]
    if 'error rate:' in line:
        err = 'ERROR\t'+str(float(line.strip().split(': ')[1][:-1])/100.0)
    if '[1]' in line and '[[1]]' not in line and 'Done' not in line:
        auc = 'AUC\t'+line.strip().replace('[1] ','')
f.close()

line_ext.append(kappa)
line_ext.append(accu)
line_ext.append(spec)
line_ext.append(sens)
line_ext.append(err)
line_ext.append(auc)

f_out = open("../results/result-1__rf_performance.txt","w")
f_out.write('\n'.join(line_ext))
f_out.close()

