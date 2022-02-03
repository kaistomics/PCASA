#################################################################################
#### Package
#################################################################################

import datetime
import time
import sys

####

from tensorflow import keras
from keras.utils import np_utils
from keras.models import Sequential
from keras.layers import Dense, Conv2D, MaxPooling2D, Dropout, Flatten, pooling
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split

####

import tensorflow as tf
from tensorflow.keras import backend as K
from tf_keras_vis.utils import normalize
from tf_keras_vis.gradcam import GradcamPlusPlus

####

from sklearn.metrics import roc_curve, auc
from numpy import interp
from itertools import cycle
import os

####

#np_load_old = np.load
#np.load = lambda *a,**k: np_load_old(*a, allow_pickle=True, **k)

####

f_list = os.listdir('../results/')
dir_name = 'combi_cnn'
if dir_name not in f_list:
    os.system("mkdir ../results/"+dir_name)
    os.system("mkdir ../results/"+dir_name+'/gcam/')
else:
    os.system("rm -r ../results/"+dir_name)
    os.system("mkdir ../results/"+dir_name)
    os.system("mkdir ../results/"+dir_name+'/gcam/')

start = time.time()


#################################################################################
#### CNN - Classification of tumor & normal cells using gene-combinations
#################################################################################

X_nor = np.load("../results/training_normal.npy")
X_abn = np.load("../results/training_tumor.npy")

Y_nor = np.ndarray(shape=(X_nor.shape[0], 1))
for i in range(X_nor.shape[0]):
    Y_nor[i] = 0

Y_abn = np.ndarray(shape=(X_abn.shape[0], 1))
for i in range(Y_abn.shape[0]):
    Y_abn[i] = 1

X_nor = X_nor.reshape(X_nor.shape[0], 4, 9900, 1)
X_abn = X_abn.reshape(X_abn.shape[0], 4, 9900, 1)
X_nor = X_nor.astype('float')
X_abn = X_abn.astype('float')

print(X_nor.shape)
print(X_abn.shape)

X_all = np.concatenate((X_nor, X_abn), axis=0)
Y_all = np.concatenate((Y_nor, Y_abn), axis=0).astype('int')

print (X_all.shape)
print (Y_all.shape)

tmp = np.zeros((Y_all.shape[0], 2))
for i in range(Y_all.shape[0]):
   loci = int(Y_all[i])
   tmp[i][loci] = 1

Y_all = tmp

####

def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

def f1_m(y_true, y_pred):
    precision = precision_m(y_true, y_pred)
    recall = recall_m(y_true, y_pred)
    return 2*((precision*recall)/(precision+recall+K.epsilon()))

####

line_final = list()

for j in range(10):

    count = str(j+1)

    X_tmp, X_test, Y_tmp, Y_test = train_test_split(X_all, Y_all, test_size=0.2)
    X_train, X_validation, Y_train, Y_validation = train_test_split(X_tmp, Y_tmp, test_size=0.25)
    
    ####
    
    model = Sequential()
    model.add(Conv2D(32, (2,1), strides=(2,1), activation='relu', input_shape=(4,9900,1)))
    model.add(Conv2D(filters=32, kernel_size=(2,1), strides=(1,1), activation='relu'))
    model.add(pooling.MaxPooling2D(pool_size=(1,11)))
    model.add(Flatten())
    model.add(Dense(450, activation='relu'))
    model.add(Dropout(0.25))
    model.add(Dense(45, activation='relu'))
    model.add(Dense(2, activation='sigmoid', name="visualized_layer"))
    
    #### Compile the model ####
    model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy',f1_m, precision_m, recall_m])
    
    #### Training & Test - Fit the model ####
    history = model.fit(X_train, Y_train, validation_data=(X_validation, Y_validation), epochs=10, batch_size=500, verbose=1)
    
    #### Evaluate the model ####
    loss, accuracy, f1_score, precision, recall = model.evaluate(X_test, Y_test, verbose=1)

    #### Save stat and model ####
    stat = 'loss\t'+str(loss)+'\n'+'accuracy\t'+str(accuracy)+'\n'+'f1_score\t'+\
            str(f1_score)+'\n'+'precision\t'+str(precision)+'\n'+'recall\t'+str(recall)
    f_out = open("../results/combi_cnn/combi-cnn_"+count+"_stat.txt","w")
    f_out.write(stat)
    f_out.close()
    model.save('../results/combi_cnn/combi-cnn_'+count+'.h5')
    
   
    #################################################################################
    #### Save performance result.
    #################################################################################

    print('CNN Round : ',count)
    line_final.append('\n# Round No.'+str(count))

    Y_pred_keras = model.predict(X_test)
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(2):
        fpr[i], tpr[i], _ = roc_curve(Y_test[:, i], Y_pred_keras[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    #### Micro average ####
    fpr["micro"], tpr["micro"], _ = roc_curve(Y_test.ravel(), Y_pred_keras.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    
    #### Macro average ####
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(2)]))
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(2):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    mean_tpr /= 2
    
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    AUC_val = roc_auc["micro"]
    cnt_nor = len(X_nor)
    cnt_tum = len(X_abn)
    print(AUC_val,cnt_nor,cnt_tum)

    line_ext = [line.strip() for line in open('../results/combi_cnn/combi-cnn_'+count+'_stat.txt')]
    line_ext.append('AUC\t'+str(AUC_val))
    line_ext.append('Tumor/Normal\t'+str(cnt_tum)+'/'+str(cnt_nor))
    f_out = open('../results/combi_cnn/combi-cnn_'+count+'_stat.txt',"w")
    f_out.write('\n'.join(line_ext))
    f_out.close()
    line_final = line_final+line_ext

    
    #################################################################################
    #### GradCAM++ visualization
    #################################################################################
    
    f_list = os.listdir('../results/combi_cnn/gcam/')
    dir_name = 'gcam_'+count
    if dir_name not in f_list:
        os.system("mkdir ../results/combi_cnn/gcam/"+dir_name)
    else:
        os.system("rm -r ../results/combi_cnn/gcam/"+dir_name)
        os.system("mkdir ../results/combi_cnn/gcam/"+dir_name)

    def loss(output):
        return output[0][1]

    def model_modifier(m):
        m.layers[-1].activation = tf.keras.activations.linear
        return m

    col_names = [line.strip() for line in open("../results/genes_shuffle.txt")]
    for i in range(X_abn.shape[0]):
        if i%1000 == 0:
            print(i)
        input_image = X_abn[i]
        gradcam = GradcamPlusPlus(model, model_modifier, clone=False)
        result = normalize(gradcam(loss, input_image, penultimate_layer=-1))
        result = pd.DataFrame(result[0])
        result.columns = col_names
        result.to_csv("../results/combi_cnn/gcam/"+dir_name+"/gcam_temp_"+str(i+1)+".txt", sep="\t")    


line_final2 = list()
for line in line_final:
    line_mod = line.upper()
    line_final2.append(line_mod)

f_out = open('../results/result-3__cnn_performance.txt',"w")
f_out.write('\n'.join(line_final2[1:]))
f_out.close()


####

end = time.time()
sec = (end - start)
result = datetime.timedelta(seconds=sec)
result_list = str(datetime.timedelta(seconds=sec)).split(".")
print(result_list[0])

