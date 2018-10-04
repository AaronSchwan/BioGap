
import csv
import os
import urllib.request
import urllib.parse
import re




def search_ncbi(urls,prob):
    data = [];
    breaker = 0;
    
    urls = list(filter(None, urls))

    
    z = len(urls)-1;
    while(z>=0):
        
        resp = urllib.request.urlopen(urls[z])
        respData = str(str(str(resp.read()).split('\\t')).split('\\n')).split(',')

        x = 21;
        found = 0;
        
        while(x<len(respData) and found == 0):

            if(prob in str(respData[x]).replace('"','').replace("'",'')):
                data.insert(breaker,float(str(respData[x+1]).replace("'",'').replace('\\','').replace('"','').replace(' ','')))
                
                breaker+=1;
                found = 1;

           
            x+=1;

        z-=1;
        
    return data;    


