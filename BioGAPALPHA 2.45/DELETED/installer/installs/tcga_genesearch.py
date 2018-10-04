from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re

def geneSearch(entry,outliers,dataS,STAGES,AGES,SEXES,RACES):
    
    
    pdata = {};
    
    rawData = [];
    AliveDead = [];
    Censored = [];
    StageData = [];
    Survival = [];
    SexData = [];
    AgeData = [];
    RaceData = [];

    rawDataT = [];
    AliveDeadT = [];
    CensoredT = [];
    StageDataT = [];
    SurvivalT = [];
    SexDataT = [];
    AgeDataT = [];
    RaceDataT = [];

    rawDataT_ = [];
    AliveDeadT_ = [];
    CensoredT_ = [];
    StageDataT_ = [];
    SurvivalT_ = [];
    SexDataT_ = [];
    AgeDataT_ = [];
    RaceDataT_ = [];


    if(dataS == "Data Source:"):
        
        return;
            

    if(dataS == "Lung Squamous Cell Carcinoma (TCGA, Provisional)"):
        #main url
        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=lusc_tcga&cancer_study_id=lusc_tcga&genetic_profile_ids=lusc_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=lusc_tcga_all&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
        lengths = 177
        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

        with open("Data Sources\Data Lung Squamous Cell Carcinoma (TCGA, Provisional)\PDATALUSC.csv",'r') as file:
            data = file.read().split('\n');
            file.close();
        
    if(dataS == "Lung Adenocarcinoma (TCGA, Provisional)"):
        #main url
        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=luad_tcga&cancer_study_id=luad_tcga&genetic_profile_ids=luad_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=luad_tcga_all&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
        lengths = 230
        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

        with open("Data Sources\Data Lung Adenocarcinoma (TCGA, Provisional)\PDATALUAD.csv",'r') as file:
            data = file.read().split('\n');
            file.close();

        
        
    xx = len(data)-2 
    x = 0;
    while(xx>=0):
        stuff = data[x].split(',')
        
        stuff1 = [int(stuff[1]),int(stuff[2]),int(stuff[3]),float(stuff[4]),int(stuff[5]),stuff[6],int(stuff[7])]
        pdata[stuff[0]] = stuff1

            
            
        x+=1;
        xx-=1;
       



        
 
    
    urlmid = entry
    url =urlstart+urlmid+urlend
    
    #main numbers
    resp = urllib.request.urlopen(url)
    respData = str(str(str(resp.read()).split('\\t')).split('\\n')).split(',')
    xx = len(respData)-1
    x = 6;
    breaker = 0;

    #print(pdata)
    while(xx>0):
        try:
            
            tcgaval = str(str(respData[x]).replace("'",'').replace('\\','').replace('"','').replace(' ',''))
            
            data = pdata[tcgaval];
            if(str(str(respData[x+1]).replace("'",'').replace('\\','').replace('"','').replace(' ','')) != 'NaN'):
                data.insert(8,float(str(respData[x+1]).replace("'",'').replace('\\','').replace('"','').replace(' ','')))
                
                StageDataT.insert(breaker,data[0])
                rawDataT.insert(breaker,data[7])
                CensoredT.insert(breaker,data[6])
                AliveDeadT.insert(breaker,data[1])
                SexDataT.insert(breaker,data[4])
                SurvivalT.insert(breaker,data[2])
                RaceDataT.insert(breaker,data[5])
                AgeDataT.insert(breaker,data[3])
                
                breaker+=1;


                                                               
        except:
            pass;
        
        x+=1;
        xx-=1;

    if(outliers == "Include Outliers"):
        rawDataT_ = rawDataT;
        AliveDeadT_ = AliveDeadT;
        CensoredT_ = CensoredT;
        StageDataT_ = StageDataT;
        SurvivalT_ = SurvivalT;
        SexDataT_ = SexDataT;
        AgeDataT_ = AgeDataT;
        RaceDataT_ = RaceDataT;

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
                rawDataT_.insert(breaker,rawDataT[x])
                StageDataT_.insert(breaker,float(StageDataT[x]))
                SurvivalT_.insert(breaker,int(SurvivalT[x]))
                CensoredT_.insert(breaker,int(CensoredT[x]))
                AliveDeadT_.insert(breaker,int(AliveDeadT[x]))
                SexDataT_.insert(breaker,int(SexDataT[x]))
                AgeDataT_.insert(breaker,float(AgeDataT[x]))
                RaceDataT_.insert(breaker,RaceDataT[x])
                
                
                breaker+=1


            xx-=1;
            x+=1;
    xx = len(rawDataT_)-1;
    x = 0;
    lowage = AGES[0];
    highage = AGES[1];
    braker = 0;
    
    while(xx>0):
        if(StageDataT_[x] in STAGES and RaceDataT_[x] in RACES and SexDataT_[x] in SEXES and lowage<=AgeDataT_[x]<=highage):
            rawData.insert(breaker,rawDataT_[x]);
            AliveDead.insert(breaker,AliveDeadT_[x]);
            Censored.insert(breaker,CensoredT_[x]);
            StageData.insert(breaker,StageDataT_[x]);
            Survival.insert(breaker,SurvivalT_[x]);
            SexData.insert(breaker,SexDataT_[x]);
            AgeData.insert(breaker,AgeDataT_[x]);
            RaceData.insert(breaker,RaceDataT_[x]);

            breaker+=1;
            
        xx-=1;
        x+=1;

   
    return StageData,rawData,RaceData,SexData,Survival,Censored,AliveDead,AgeData;

    


        
