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

from lifelines.statistics import logrank_test

    
global gene
global rawData
global dataS
with open("LISTOFGENES.csv","r") as file:
    GENE = file.read().split('\n');
    file.close
    
g = len(GENE)-2;

print(g)
while(g>0):
    dataS = "Lung Adenocarcinoma (TCGA, Provisional)"
    outliers = 2;
    gene = GENE[g]
    print(gene)
    try:
        

        

      
        
        #main url
        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=luad_tcga&cancer_study_id=luad_tcga&genetic_profile_ids=luad_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=luad_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
        lengths = 230
        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

       



        
        with open(os.path.join("Data Sources\Data "+dataS+"\Kaplan Meiers\Paitent Data\Survival In Days.csv"),"r") as survival:

            SurvivalT = survival.read().split('\n')
            survival.close()
        
        global Survival
        Survival = [];
        del Survival[:];
        
        with open(os.path.join("Data Sources\Data "+dataS+"\Kaplan Meiers\Paitent Data\Censord Data.csv"),"r") as censord:

            CensordT = censord.read().split('\n')
            censord.close()
        
        global Censord
        Censord = [];
        del Censord[:];



        
        urlmid = gene
        url =urlstart+urlmid+urlend
        
        #main numbers
        resp = urllib.request.urlopen(url)
        respData = resp.read()
       
        if(str(respData).find("# Warning:  Unknown gene:")==2):

           print("no",gene)
            
        
        if(str(respData).find("# Warning:  Unknown gene:")!=2):
                
            rawDataT = str(respData).replace(gene,"1")    
            newData = re.findall(r'[-+]?\d*\.\d+|\d+',str(rawDataT))
            


            try:
                #Varibles
                z = 0;
                x = 5;
                y = 0;
                rawDataT = [];
                del rawDataT[:];
                #insert data to rawData
                while(y<lengths):
                
                    rawDataT.insert(y,float(newData[x]))

                    x+=4
                    y+=1

                
                if(len(rawDataT) != lengths):
                    print("Error: Something Went Wrong With The Datas List It Isnt the Right Length")
                    sys.exit()
                


              
                
                rawData = [];
                del rawData[:];
                
                

                if(outliers == 2):

                    medi = median(rawDataT)


                    qu1 = [];
                    qu3 = [];

                    xx=len(rawDataT)
                    x=0;
                    bre1 = 0;
                    bre3 =0;
                    while(xx>0):
                        if(rawDataT[x]<medi):
                            qu1.insert(bre1,rawDataT[x])
                            bre1+=1

                        if(rawDataT[x]>medi):
                            qu3.insert(bre3,rawDataT[x])
                            bre3+=1
                        
                        x+=1
                        xx-=1



                    medq1 = median(qu1)

                    medq3 = median(qu3)

                    iqr = medq3-medq1;
                    xx=len(rawDataT)

                    lowout = medq1-1.5*iqr;
                    highout = medq3+1.5*iqr;
                    breaker = 0;
                    x=0;
                    
                    while(xx>0):
                        
                        if(highout>=rawDataT[x]>=lowout ):
                            rawData.insert(breaker,rawDataT[x])
                            Survival.insert(breaker,int(SurvivalT[x]))
                            Censord.insert(breaker,int(CensordT[x]))
                            breaker+=1


                        
                        
                        x+=1
                        xx-=1
                    
                    

            except:
                pass;
    except:
        pass;
    if(2 == 2):
        step = int(max(rawData)/100);
        num = 1;
        pvalues = [];
        hs = [];
        ls=[];
        z = 0;
        l1 = 0;
        h1 = step*num;
        

        while(num<100):
            
            l2 = max(rawData);
            l1 = 0;
            h1 = step*num;
            
            high = 0;
            kaplandata1 =[];
                

                
            xxx = int(max(Survival));
            breaker1=0;
                
                
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h1 >= rawData[x] >= l1):
                    kaplandata1.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                        
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            if(len(kaplandata1) >= 15):

                #second line
                while(l2>h1):
                    
                    lh = high*step;
                    
                    mr = max(rawData);
                    
                    l2 = mr-lh;
                    h2 = max(rawData);
                    kaplandata2 =[];
                    

                        
                    xxx = int(max(Survival));
                    breaker1=0;
                        
                        
                    x = 0;
                    xx = len(rawData)
                    while(xx>0):
                        if(h2 >= rawData[x] >= l2):
                            kaplandata2.insert(breaker1,int(Survival[x]))
                            breaker1+=1;
                                
                        x+=1;
                        xx-=1;
                    x=0;
                    xxx-=1;

                    if(len(kaplandata1) >= 15 and len(kaplandata2) >= 15):
                        lowkap = sorted(kaplandata1);
                        highkap = sorted(kaplandata2);

                        summary, p_value, results = logrank_test(lowkap, highkap, alpha=.95)
                        pvalues.insert(z,p_value);
                        hs.insert(z,l2);
                        ls.insert(z,h1);
                        z+=1;
                   
                    high+=1;
                
                
   
            num+=1;
            



        kapmin = min(pvalues)
        
        
        if(2 == 2):
            i = pvalues.index(kapmin)
            
            l1 =0;
            h1 = ls[i]

            h2= max(rawData);
            l2 = hs[i];

            with open("bestkaplan.csv","a") as savefile:
                savefile.write("\n");
                savefile.write(gene);
                savefile.write(',');
                savefile.write(str(h1));
                savefile.write(',');
                savefile.write(str(l2));
                savefile.write(',');
                savefile.write(str(kapmin));
                savefile.close()
    g-=1;
