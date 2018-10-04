import matplotlib.pyplot as plt
import numpy as np
import ncbi_correlate
import time


def heatmap(labels,dataS,outliers):


    xx = len(labels)
    x = 0;

    rows = [[] for i in range(1, xx+1)]
    while(xx>0):

        entry = labels[x];

        ss = len(labels);
        s = 0;
        
        while(ss>0):

            Rsquared,ttest,rawDataC,rawData,regline = ncbi_correlate.correlate(entry,labels[s],dataS,outliers);
            rows[x].append(Rsquared)



            ss-=1;
            s+=1;
        
        xx-=1;
        x+=1;

        

    plt.matshow(rows)



    x_pos = np.arange(len(labels))
    plt.xticks(x_pos,labels)

    y_pos = np.arange(len(labels))
    plt.yticks(y_pos,labels)

    plt.show()


