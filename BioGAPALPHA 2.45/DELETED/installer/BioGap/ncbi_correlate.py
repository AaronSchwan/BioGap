from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re

import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


import ncbi_data_pull

def best_fit_slope_and_intercept(xs,ys):
    m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
         ((mean(xs)*mean(xs)) - mean(xs*xs)))
    b = mean(ys) - m*mean(xs)
    return m, b

def squared_error(ys_orig,ys_line):
    return sum((ys_line - ys_orig) * (ys_line - ys_orig))

def coefficient_of_determination(ys_orig,ys_line):
    y_mean_line = [mean(ys_orig) for y in ys_orig]
    squared_error_regr = squared_error(ys_orig, ys_line)
    squared_error_y_mean = squared_error(ys_orig, y_mean_line)
    return 1 - (squared_error_regr/squared_error_y_mean)


   
def correlate(prob,probC,dataS,outliers):
    if(dataS == 'Small Cell Lung Cancer'):
        with open('Data Sources\_NCBI_\smallcell.csv','r') as file:
            urls = file.read().split('
')
            file.close()

    if(dataS == 'Lung Adenocarcinoma'):
        with open('Data Sources\_NCBI_\lungadeno.csv','r') as file:
            urls = file.read().split('
')
            file.close()

    if(dataS == 'Lung Squamous Carcinoma'):
        with open('Data Sources\_NCBI_\lungsquamous.csv','r') as file:
            urls = file.read().split('
')
            file.close()

    
    
    if(len(urls) > 0 and prob != 'Search Prob Number'):      

        dataT = ncbi_data_pull.search_ncbi(urls,prob)
        
    rawDataT = dataT;

    if(dataS == 'Small Cell Lung Cancer'):
        with open('Data Sources\_NCBI_\smallcell.csv','r') as file:
            urls = file.read().split('
')
            file.close()

    if(dataS == 'Lung Adenocarcinoma'):
        with open('Data Sources\_NCBI_\lungadeno.csv','r') as file:
            urls = file.read().split('
')
            file.close()

    if(dataS == 'Lung Squamous Carcinoma'):
        with open('Data Sources\_NCBI_\lungsquamous.csv','r') as file:
            urls = file.read().split('
')
            file.close()


    if(len(urls) > 0 and probC != 'Correlate to Prob'):      

        dataTC = ncbi_data_pull.search_ncbi(urls,probC)
        
    rawDataTC = dataTC;

    rawData = [];
    rawDataC = [];
    
    if(outliers == "Include Outliers"):
        rawData = rawDataT;
        rawDataC = rawDataTC;


    if(outliers == "Exclude Outliers"):
        
        
        
        q75, q25 = np.percentile(rawDataT,[75,25])
        iqrnum = q75-q25

        lowout = q25-1.5*iqrnum;
        highout = q75+1.5*iqrnum;

        
        
        breaker = 0;
        x=0;

        xx = len(rawDataT)-1;
        while(xx>0):
            
            if(highout>=rawDataT[x]>=lowout ):
                rawData.insert(breaker,rawDataT[x])
                rawDataC.insert(breaker,rawDataTC[x])

                
                
                breaker+=1


            xx-=1;
            x+=1;

       
    ys = np.array(rawDataC)
    xs = np.array(rawData)

   
             
    m, b = best_fit_slope_and_intercept(xs,ys)
    regression_line = [(m*x)+b for x in xs]

    r_squared = coefficient_of_determination(ys,regression_line)
    
    
    

    idsforline = [];
    rawMax=max(rawData)+100
      
    rawMin=min(rawData)
       
    insert=0;
        
    while(rawMax>rawMin):
        idsforline.insert(insert,insert)                                                                                                 
            
        insert+=1
        rawMax-=1
        
    
    regline = [];
    
    

    
    rawMax=max(rawData)
      
    rawMin=min(rawData)
    
    regline.insert(0,m*rawMin+b)
    regline.insert(1,m*rawMax+b)
        



    MR = mean(rawData)
    MC = mean(rawDataC)
    ttop=MR - MC;
    n1 = np.std(rawData)/len(rawData)
    n2 = np.std(rawDataC)/len(rawDataC)
    n3 = n1+n2                        
    tbottom=np.sqrt(n3)
    ttest = ttop/tbottom

    return r_squared,ttest,rawDataC,rawData,regline ;


    




        
