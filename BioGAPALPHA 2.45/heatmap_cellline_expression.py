import matplotlib.pyplot as plt
import numpy as np
import ncbi_correlate
import time
import ncbi_data_pull

def heatmap(prob):
    global FINAL

    pdata = []

    FINAL = []
    del FINAL[:]

    

    with open('Data Sources\Celllines\cell_lines.csv','r') as file:
        data = file.read().split('\n');
        file.close();
        


    data = list(filter(None, data))

    xx = len(data)-1; 
    x = 0;
    while(xx>=0):
        stuff = data[x].split(',')

        
        stuff1 = [str(stuff[1]),str(stuff[2]),str(stuff[3]),str(stuff[4]),str(stuff[5]),str(stuff[6]),str(stuff[7])]
        pdata.insert(x,stuff1)



            
            
        x+=1;
        xx-=1;
    
    xx = len(pdata)
    x = 0;

    row = [];
    labels = [];

    breaker = 0;
    
    while(xx>0):

        try:
            data = pdata[x]
            title = str(data[0]).replace("'","");
            urls = [str(data[5]).replace("'","")];
            
        
            
            dataT = ncbi_data_pull.search_ncbi(urls,prob)
            

            if(dataT[0] != ''):
                row.insert(breaker,dataT[0])
                labels.insert(breaker,title)


                breaker+=1;
                
        except:
            pass;
        
        xx-=1;
        x+=1;

    rows = [row];

    plt.matshow(rows)



    x_pos = np.arange(len(labels))
    plt.xticks(x_pos,labels, rotation='vertical')


    plt.show()


