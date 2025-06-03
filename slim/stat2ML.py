#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 17:32:26 2022

@author: xinjunzhang
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 00:32:59 2021

@author: xinjunzhang
"""


import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
from inspect import signature
import seaborn as sns
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.feature_selection import RFE
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
import os,argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from itertools import cycle
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from itertools import cycle
import pickle
from sklearn.svm import SVC
from sklearn.decomposition import PCA
import gzip
import random



def prec_recall_mod (cm, num_classes): #this is to get precision recall for multiple flanking classes
    prec_list = []
    reca_list = []    
    for i in range(0,num_classes):
        prec_i = round(cm[i,i]/sum(cm[:,i]),4)
        reca_i = round(cm[i,i]/sum(cm[i,:]),4)
        prec_list.append(prec_i)
        reca_list.append(reca_i)
    prec_list = [0 if np.isnan(x) else x for x in prec_list]
    reca_list = [0 if np.isnan(x) else x for x in reca_list]    
    print(prec_list)
    print(reca_list)    
    return prec_list,reca_list

def plot_prcurve_average_mod(testY,yscore,savepath,num_classes):
    y_score = yscore
    ydata = testY #binarized
    precision = dict()
    recall = dict()
    average_precision = dict()
    thresh = dict()
    for i in range(num_classes):
        precision[i], recall[i], thresh[i] = precision_recall_curve(ydata[:, i],y_score[:, i],pos_label=1)
        average_precision[i] = average_precision_score(ydata[:, i], y_score[:, i])
    precision["micro"], recall["micro"], thresh["micro"] = precision_recall_curve(ydata.ravel(),y_score.ravel())
    average_precision["micro"] = average_precision_score(ydata, y_score,average="micro")
    # setup plot details
    #colors = cycle(['blue', 'red', 'yellow'])
    colors = ["b","r"] + sns.color_palette("Greens_r",num_classes)[2:num_classes]
    #sns.palplot(colors)
    #colors = sns.palplot('blue', 'red',sns.color_palette("Oranges_r",7)[2:7])
    #choose several shades of yellows
    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        #plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))
    #lines.append(l)
    #labels.append('f1 score curves')
    #l, = plt.plot(recall["micro"], precision["micro"], color='black', lw=2)
    #lines.append(l)
    #labels.append('Average Precision-Recall (area = {0:0.2f})'''.format(average_precision["micro"]))
    for i, color in zip(range(num_classes), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-Recall for class {0} (area = {1:0.2f})'
                  ''.format(i, average_precision[i]))
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall curve of all classes')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=10))
    fig.savefig(savepath+'.png')
    plt.show()



def plot_roc_mod(testY,yscore,savepath,num_classes):
    y_score = yscore
    ydata = testY #binarized
    fpr = dict()
    tpr = dict()
    thresh = dict()
    for i in range(num_classes):
        fpr[i], tpr[i], thresh[i] = roc_curve(ydata[:, i],y_score[:, i],pos_label=1)
    #fpr["micro"], tpr["micro"], thresh["micro"] = roc_curve(ydata.ravel(),y_score.ravel())
# setup plot details
    colors = ["b","r"] + sns.color_palette("Greens_r",num_classes)[2:num_classes]
    plt.figure(figsize=(7, 8))
    lines = []
    labels = []
    for i, color in zip(range(num_classes), colors):
        plt.plot([0, 1], [0, 1], 'k--')
        l, = plt.plot(fpr[i], tpr[i], color=color, lw=2)
        lines.append(l)
        labels.append('ROC for class {0}'
                  ''.format(i))
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC curve of all classes')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=10))    
    fig.savefig(savepath+'.png')    
    plt.show()



def write_importance (outfile,feature_list,score_list):
    with open(outfile, "w") as feature_file:
        feature_name = feature_list
        score = score_list
        feature_file.writelines(i+"\t" for i in feature_name)
        feature_file.write("\n")
        feature_file.writelines(str(i)+"\t" for i in score)
        feature_file.write("\n")

def plot_featurescore (scores_file,savefile):
    #scores_file = DIR_save+whichML+"/cluster/"+"feature_importance_recessive+partial.txt.gzrecessive+partial_.txt"
    #scores_file = "feature-feature_importance_rec+part_set6-all.txt"
    with open(scores_file,"rt") as file:
        names = file.readline()
        scores = file.readline()        
    names = names.split("\t")
    names = names[0:len(names)-1]
    scores = scores.split("\t")
    scores = scores[0:len(scores)-1]  
    scores   = [float(i) for i in scores]    
    scores, names = zip(*sorted(zip(scores, names)))    
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25,left=0.25)
    plt.barh(range(len(names)),scores, align="center")
    plt.yticks(range(len(names)),names,fontsize=8)
    plt.xlabel("Feature Importance",fontsize=16)
    plt.ylabel("Features",fontsize=16)
    plt.title("Feature Importance Ranking")
    plt.show()
    fig.savefig(savefile,orientation='landscape') 


def ML_cm(dataframe,whichML,savename,to_scale,use_thresh):
    ETC = ExtraTreesClassifier(n_estimators=100, random_state=0,max_features="sqrt",min_samples_split=10,min_samples_leaf=10) #saga, sag, lbfgs,liblinear
    RF = RandomForestClassifier(n_estimators=100, random_state=0,max_depth=2) #saga, sag, lbfgs,liblinear
    L0LR = LogisticRegression(solver = 'saga', multi_class='multinomial',penalty='l1',tol=0.01,C=1e10)
    L1LR = LogisticRegression(solver = 'saga', multi_class='multinomial',penalty='l1')
    L2LR = LogisticRegression(solver = 'lbfgs', multi_class='multinomial',penalty='l2')
    MLP = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(512, 256, 128), random_state=1)    
    RBF = SVC(kernel='rbf',class_weight='balanced',probability=True) #kernal, penalized    
    if (whichML=="ETC"):
        ML = ETC
    elif (whichML=="RF"):
        ML = RF
    elif (whichML=="L0LR"):
        ML = L0LR
    elif (whichML=="L1LR"):
        ML = L1LR       
    elif (whichML=="L2LR"):
        ML = L2LR        
    elif (whichML=='MLP'):
        ML = MLP   
    elif (whichML=='RBF'):
        ML = RBF
        

    dataframe = dataframe.replace(float("inf"), 100) 
    dataframe = dataframe.fillna(-100)
    array = dataframe.values

    X = array[:,:]
    
    if to_scale==True:
        scaler = StandardScaler()
        scaler.fit(X[:,5:dataframe.shape[1]])
        X[:,5:dataframe.shape[1]] = scaler.transform(X[:,5:dataframe.shape[1]])
    
    
    Y = array[:,4]
    #X=X.astype(float)
    Y=Y.astype(int)
    train_X, test_X, train_Y, test_Y = train_test_split( X, Y, test_size=1/9.0, random_state=0)
    #train = train_X[:,5:X.shape[1]].astype(float)    
    #np.argwhere(np.isnan(train))
    
    ML.fit(train_X[:,5:X.shape[1]].astype(float), train_Y)
    #filename = DIR_save+whichML+'/feature-'+savename+'_finalized_model.sav'
    filename = DIR_save+whichML+"/"+savename+'_finalized_model.sav'
    pickle.dump(ML, open(filename, 'wb'))    

    y_score = ML.predict_proba(test_X[:,5:X.shape[1]].astype(float))
    pred = ML.predict(test_X[:,5:X.shape[1]].astype(float))   
    ydata = label_binarize(test_Y, classes=[*range(3)])[:,0:2]  
    cm = metrics.confusion_matrix(test_Y, pred)
    print(cm)
    np.savetxt(DIR_save+whichML+"/cm_"+savename+".txt", cm.astype(int), fmt='%i',delimiter="\t")
    prec,reca = prec_recall_mod(cm,2)
    prec_reca_mat = np.row_stack((prec,reca))
    np.savetxt(DIR_save+whichML+"/pr_"+savename+".txt", prec_reca_mat.astype(int), fmt='%i',delimiter="\t")  
        
    plot_prcurve_average_mod(ydata,y_score,DIR_save+whichML+"/pr-curve_"+savename,2)      
    plot_roc_mod(ydata,y_score,DIR_save+whichML+"/roc-curve_"+savename+"_01",2)
 
# =============================================================================
    if use_thresh==True:
        for prob_thresh in [0.5,0.55,0.6]:
            y_score = ML.predict_proba(test_X[:,5:X.shape[1]].astype(float))
            pred = ML.predict(test_X[:,5:X.shape[1]].astype(float))   
            ydata = label_binarize(test_Y, classes=[*range(3)])[:,0:2]  
    
            print(prob_thresh)
            y_score_mod = y_score
            pred_mod = pred
            for i in range(0,len(y_score_mod)):
                if y_score_mod[i][1]>=prob_thresh:
                    y_score_mod[i][0]=0
                    y_score_mod[i][1]=1
                    pred_mod[i]=1
                else:
                    y_score_mod[i][0]=1
                    y_score_mod[i][1]=0
                    pred_mod[i]=0
                
            savename_mod = savename+str(prob_thresh)
            cm = metrics.confusion_matrix(test_Y, pred_mod)
            print(cm)
            np.savetxt(DIR_save+whichML+"/segs_bythresh/cm_"+savename_mod+".txt", cm.astype(int), fmt='%i',delimiter="\t")
            prec,reca = prec_recall_mod(cm,2)
            prec_reca_mat = np.row_stack((prec,reca))
            np.savetxt(DIR_save+whichML+"/segs_bythresh/pr_"+savename_mod+".txt", prec_reca_mat.astype(int), fmt='%i',delimiter="\t")  
                
            plot_prcurve_average_mod(ydata,y_score_mod,DIR_save+whichML+"/segs_bythresh/pr-curve_"+savename_mod,2)      
            plot_roc_mod(ydata,y_score_mod,DIR_save+whichML+"/segs_bythresh/roc-curve_"+savename_mod+"_01",2)
    
# =============================================================================
    if whichML == "ETC":
        print(ML.feature_importances_)
        score = ML.feature_importances_.tolist()
        feature_names = list(dataframe.columns.values)[5:X.shape[1]]
        features = feature_names#[1:]
        write_importance (DIR_save+whichML+"/feature_importance_"+savename+".txt",features,score)        
        plot_featurescore(DIR_save+whichML+"/feature_importance_"+savename+".txt",DIR_save+whichML+"/feature_importance_"+savename+".png")

def ML_test(test,testname,whichML,savename,to_scale,use_thresh):


    test = test.replace(float("inf"), 100) 
    test = test.fillna(-100)
    test = test.sample(frac=1)
    
    
    array_test = test.values
    test_X = array_test[:,:] 
    
    if to_scale==True:
        scaler = StandardScaler()
        test_X[:,5:test.shape[1]] = scaler.transform(test_X[:,5:test.shape[1]])
    
    test_Y = array_test[:,4].astype(int)
 
    
    
    filename = DIR_save+whichML+"/segs/"+savename+'_finalized_model.sav'
    ML = pickle.load(open(filename, "rb"))
    

    y_score = ML.predict_proba(test_X[:,5:test_X.shape[1]].astype(float))
    pred = ML.predict(test_X[:,5:test_X.shape[1]].astype(float))   
    ydata = label_binarize(test_Y, classes=[*range(3)])[:,0:2]  
    cm = metrics.confusion_matrix(test_Y, pred)
    print(cm)
    np.savetxt(DIR_save+whichML+"/segs/cm_"+savename+"-"+testname+".txt", cm.astype(int), fmt='%i',delimiter="\t")
    prec,reca = prec_recall_mod(cm,2)
    prec_reca_mat = np.row_stack((prec,reca))
    np.savetxt(DIR_save+whichML+"/segs/pr_"+savename+"-"+testname+".txt", prec_reca_mat.astype(int), fmt='%i',delimiter="\t")  
        
    plot_prcurve_average_mod(ydata,y_score,DIR_save+whichML+"/segs/pr-curve_"+savename+"-"+testname,2)      
    plot_roc_mod(ydata,y_score,DIR_save+whichML+"/segs/roc-curve_"+savename+"_01"+"-"+testname,2)
 
# =============================================================================
    if use_thresh==True:
        for prob_thresh in [0.5,0.55,0.6]:
            y_score = ML.predict_proba(test_X[:,5:test_X.shape[1]].astype(float))
            pred = ML.predict(test_X[:,5:test_X.shape[1]].astype(float))   
            ydata = label_binarize(test_Y, classes=[*range(3)])[:,0:2]  
    
            print(prob_thresh)
            y_score_mod = y_score
            pred_mod = pred
            for i in range(0,len(y_score_mod)):
                if y_score_mod[i][1]>=prob_thresh:
                    y_score_mod[i][0]=0
                    y_score_mod[i][1]=1
                    pred_mod[i]=1
                else:
                    y_score_mod[i][0]=1
                    y_score_mod[i][1]=0
                    pred_mod[i]=0
                
            savename_mod = savename+str(prob_thresh)
            cm = metrics.confusion_matrix(test_Y, pred_mod)
            print(cm)
            np.savetxt(DIR_save+whichML+"/segs_bythresh/cm_"+savename_mod+".txt", cm.astype(int), fmt='%i',delimiter="\t")
            prec,reca = prec_recall_mod(cm,2)
            prec_reca_mat = np.row_stack((prec,reca))
            np.savetxt(DIR_save+whichML+"/segs_bythresh/pr_"+savename_mod+".txt", prec_reca_mat.astype(int), fmt='%i',delimiter="\t")  
                
            plot_prcurve_average_mod(ydata,y_score_mod,DIR_save+whichML+"/segs_bythresh/pr-curve_"+savename_mod,2)      
            plot_roc_mod(ydata,y_score_mod,DIR_save+whichML+"/segs_bythresh/roc-curve_"+savename_mod+"_01",2)
    

def pre_process_csv (dataframe):
    dataframe = dataframe.fillna(-100)
    dataframe = dataframe.replace(float("inf"), 100) 
    dataframe = dataframe.sample(frac=1)
    return dataframe

if __name__ == "__main__":

    DIR_master = "/Users/xinjunzhang/Desktop/Dominance_project/simulation/stat_1kseg/"
    DIR_save = "/Users/xinjunzhang/Desktop/Dominance_project/simulation/stat_1kseg/ML/"
    
    os.chdir(DIR_master)
        
    #first, process training data to divide into different groups
    #actual data starts from the 4th column index
    training_data_path = '220725combined_stats_1MB_1ksegbyexon.csv'
    training = pd.read_csv(training_data_path, sep=",")
    
    training = pre_process_csv (training)
    #names = list(training.columns.values)
    
    train_all = training
    train_exon = training.loc[training["exon_window"]!=0]
    
    train_single = training.loc[training["segment"]=="seg1"]
    train_minus1 = training.loc[training["segment"]!="seg1"]
    train_highexon = training.loc[training["exon_density"]>=700]
    train_lowexon = training.loc[training["exon_density"]<700]
    train_half = training.iloc[0:int(train_all.shape[0]/2)]
    train_otherhalf = training.iloc[int(train_all.shape[0]/2):train_all.shape[0]]
    
    train_all.to_csv("ML/train_all_stat.csv", index=False)
    train_exon.to_csv("ML/train_exon_stat.csv", index=False)
    train_single.to_csv("ML/train_single_stat.csv", index=False)
    train_minus1.to_csv("ML/train_minus1_stat.csv", index=False)
    train_highexon.to_csv("ML/train_highexon_stat.csv", index=False)
    train_lowexon.to_csv("ML/train_lowexon_stat.csv", index=False)
    train_half.to_csv("ML/train_half_stat.csv", index=False)
    train_otherhalf.to_csv("ML/train_otherhalf_stat.csv", index=False)
    
    ##########################################################
    #test_all = pd.read_csv('220725combined_testing_stats.csv', sep=",")
    #test_varym = pd.read_csv('210818combined_testing-varym_stats.csv', sep=",")
    #test_varymu = pd.read_csv('210818combined_testing-varymu_stats.csv', sep=",")
    #test_varyt = pd.read_csv('210818combined_testing-varyt_stats.csv', sep=",")
    
    
    ######
    #get different types of training
    all_files = os.listdir("ML/")
    all_files = [x for x in all_files if "train_" in x]
    
    for n in range(0,len(all_files)):
        file = all_files[n]
        train = pd.read_csv(DIR_save+file, sep=",")
        savename=file.split("_")[1]
        to_scale=False
        use_thresh=False
        print(savename)
        for whichML in ["ETC"]: #["ETC","RF","L0LR","L1LR","L2LR"]     
            print(whichML)
            ML_cm(train,whichML,savename,to_scale,use_thresh=False)
        
    
    #####turned out ETC single is doing somehting, so break them down by segment
    train_all = training
    all_segs = list(set(train_all["segment"]))
    
    exons = []
    recrates = []
    introg_anc = []
    
    for seg in all_segs:
        train_single = training.loc[training["segment"]==seg]
        #train_single = train_exon.loc[training["segment"]==seg]
        print(seg)
        recrates.append(train_single["mean_recrate"].values[0])
        exons.append(train_single["exon_density"].values[0])
        introg_anc.append(train_single["mean_introg_anc"].values[0])
        
        #train_single.to_csv("ML/train_"+seg+"_stat.csv", index=False)
        
        use_thresh=False
        to_scale=False
        #to_scale=True
        #ML_cm(train_single,"ETC",seg,to_scale,use_thresh)
        #ML_cm(train_single,"RF",seg,to_scale,use_thresh)
        #ML_cm(train_single,"L0LR",seg,to_scale,use_thresh)

    info = pd.DataFrame({"segment":all_segs,"exon_density":exons,"recombination_rate":recrates,"mean_introgressed_ancestry":introg_anc})
    info.to_csv("segment_information.csv", index=False)
    
    testing = pd.read_csv(DIR_save+"train-all40_highexon800_stat_1MB.csv", sep=",")
    ML_test(testing,"highexon_200-on-800_1MB","ETC","all40-2batches_highexon2001MB-scaled",True,False)
    
    #####find the best performing seg=seg1, try the model on test data
    
    #test_varym = test_varym.loc[test_varym["segment"]=="seg1"]
    #test_varymu = test_varymu.loc[test_varymu["segment"]=="seg1"]
    #test_varyt = test_varyt.loc[test_varyt["segment"]=="seg1"]
    
    #ML_test(test_varym,"varym","ETC","seg1",False,False)
    #ML_test(test_varymu,"varymu","ETC","seg1",False,False)
    #ML_test(test_varyt,"varyt","ETC","seg1",False,False)
