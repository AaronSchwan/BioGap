
from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re


def searchcellline(cellline):
    
    url = "https://www.ncbi.nlm.nih.gov/geoprofiles/?term=CONTROL+"+cellline


    resp = urllib.request.urlopen(url)
    respData = resp.read()

    listall = respData.split()


    xx = len(listall)
    x = 0;
    ids = [];
    breaker = 0;


    while(xx>0):

        if("ID=G" in str(listall[x])):
            idval = (str(listall[x]).replace("b'href=",'').replace("/geo/tools/profileGraph.cgi?",'').replace('"><img','').replace("'",'').replace('"',''))
            ids.insert(breaker,idval)
            breaker+=1;
            
        x+=1;
        xx-=1;

    z = len(ids)-1
    ##print(ids)
    idcontrols = [];
    idcontrolsdata = [];
    breaker = 0;
    while(z>1):
       
        url = "https://www.ncbi.nlm.nih.gov/geo/tools/profileGraph.cgi?"+str(ids[z])

        resp = urllib.request.urlopen(url)
        respData = resp.read()

        listall = str(respData).split()

        controlurls = [];
        breker = 0;

        xx = len(listall)
        x = 0;

        while(xx>0):

            if('href="/geo/query/acc.cgi?acc=' in str(listall[x])):
                cont = (str(listall[x]).replace('href=','https://www.ncbi.nlm.nih.gov').replace('"',''))
               
                controlurls.insert(breaker,cont)

                
                breaker+=1;
                
            x+=1;
            xx-=1;
        
       
        ss = int(len(controlurls)/2)-1



        while(ss>0):
            url = controlurls[ss]
           
            resp = urllib.request.urlopen(url)
            respData = resp.read()


            listall = str(respData).split('justify">')


        
            
            title_and_sample_type  = str(listall[1]).split('</td>\\n</tr>\\n<tr valign="top"><td nowrap>Sample type')
            
            
            title = title_and_sample_type[0]
            sampletype = str(title_and_sample_type[1]).replace('</td>\\n</tr>\\n<tr valign="top"><td nowrap>&nbsp;</td>\\n<td></td>\\n</tr>\\n<tr valign="top"><td nowrap>Source name</td>\\n<td style="text-align: ','').replace('</td>\\n<td>','')

            
            source_name_organism  = str(listall[2]).split('<br></td>\\n</tr>\\n<tr valign="top"><td nowrap>Organism</td>\\n<td>')
           
            
            sourcename = source_name_organism[0]
           



            xx = len(listall)-1

            while(xx>0):
                if('>Description</td>\\n<td style="text-align:' in str(listall[xx])):
                    description_ = str(listall[xx+1]).split('</td>')
                    description = str(description_[0]).replace('<br>','\n')
                xx-=1;
            


            listall = [str(title).replace(',',' '),str(sampletype).replace(',',' '),str(sourcename).replace(',',' '),str(description).replace(',',' '),controlurls[ss]]
            
          
            #str(listall).count('control') > 0 and 
            if(cellline in str(listall)):
                
                idcontrols.insert(breaker,listall)
               

                breaker+=1;
                

            ss-=1;


        z-=1;

    idcontrols.sort()

    i = len(idcontrols) - 1
    while i > 0:  
        if idcontrols[i] == idcontrols[i - 1]:
            idcontrols.pop(i)
        i -= 1

    URLTABLE = [];
    FINAL = [];
    breaker =0;
    breaker2 = 0;

    x = len(idcontrols)-1;
    while(x>=0):
        idurl = idcontrols[x]
        url = idurl[4]
        resp = urllib.request.urlopen(url)
        respData = resp.read()
        listall = str(respData).split()

        xx = len(listall)-1;

        GSMNUM = url.replace("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",'')

        while(xx>0):
            
            if('&amp;id=' in str(listall[xx]) and GSMNUM in str(listall[xx]) and "&amp;db=" in str(listall[xx])):
                table = str(listall[xx]).split("&amp;")
                

                urltable = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?view=data&"+str(table[1])+'&'+str(table[2])+"&"+str(table[3]).replace("\\',",'')
                URLTABLE.insert(breaker,urltable)
                idurl.insert(5,urltable)
                FINAL.insert(breaker,idurl);
                breaker+=1;
            if("Contact" in str(listall[xx])):
                first = str(listall[xx+1]).replace("name</td>\\n<td>",'')
                last = str(listall[xx+2]).replace('</td>\\n</tr>\\n<tr','')
                
                name = first+' '+last;
                idurl.insert(6,name)
                
                FINAL.insert(breaker,idurl);
                breaker+=1;
                
                
            xx-=1;


        x-=1;

    FINAL.sort()

    i = len(FINAL) - 1
    while i > 0:  
        if FINAL[i] == FINAL[i - 1]:
            FINAL.pop(i)
        i -= 1

    x = len(FINAL);

    
    if(x == 0):
        print('no data found for ',cellline)
        
    if(x>=1):
        with open('Data Sources\Celllines\cell_lines.csv','a') as file:
            
            file.write(str(cellline+'1'))
            file.write(',')
            file.write(str(FINAL[0]).replace('[','').replace(']',''))
            file.write('\n');
            file.close();

    if(x>=2):
        with open('Data Sources\Celllines\cell_lines.csv','a') as file:
            
            file.write(str(cellline+'2'))
            file.write(',')
            file.write(str(FINAL[1]).replace('[','').replace(']',''))
            file.write('\n');
            file.close();

    if(x>=3):
        with open('Data Sources\Celllines\cell_lines.csv','a') as file:

            file.write(str(cellline+'3'))
            file.write(',')
            file.write(str(FINAL[2]).replace('[','').replace(']',''))
            file.write('\n');
            file.close();

    if(x>=4):
        with open('Data Sources\Celllines\cell_lines.csv','a') as file:

            file.write(str(cellline+'4'))
            file.write(',')
            file.write(str(FINAL[3]).replace('[','').replace(']',''))
            file.write('\n');
            file.close();
                
    return FINAL



