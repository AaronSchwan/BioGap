
import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Kaplan_Meier_Graphing_expected
import Kaplan_Meier_Graphing


def kaplanmeiers(entry,rawData,Survival,AliveDead):


    
    
    P.figure()

    step = 0;
    stepval = int(max(rawData)/100);

    pvals = [];
    L2 = [];
    H1 = [];
    
    pbreaker = 0;

    minpat = len(rawData)/8;

    while(step<100):
        

       
        l1 = 0;
        h1 = step*stepval;

        

        kaplandata1 = [];
        kaplancen1 = [];
        

        
        xxx = int(max(Survival));
        breaker1=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            
            if(h1 >= rawData[x] >= l1):
                
                kaplandata1.insert(breaker1,Survival[x])
                kaplancen1.insert(breaker1,AliveDead[x])
                breaker1+=1;
                
                
            x+=1;
            
            xx-=1;

        if(len(kaplandata1) > minpat):
            highstep = 0;

            h2 = max(rawData)
            l2 = h2-stepval*highstep;
            
            while(h1 < l2):
                
                h2 = max(rawData)
                l2 = h2-stepval*highstep;

                kaplandata2 = [];
                kaplancen2 = [];
                xxx = int(max(Survival));
                breaker2=0;
                x = 0;
                xx = len(rawData)
                while(xx>0):
                    if(h2 >= rawData[x] >= l2):
                        kaplandata2.insert(breaker2,Survival[x])
                        kaplancen2.insert(breaker2,AliveDead[x])
                        breaker2+=1;
                        
                    x+=1;
                    xx-=1;
                if(len(kaplandata1) > minpat and len(kaplandata2) > minpat):

                    


                    vals = Kaplan_Meier_Graphing_expected.kaplan(kaplandata1,kaplancen1,kaplandata2,kaplancen2)
                    

                    pvals.insert(pbreaker,vals)
                    H1.insert(pbreaker,h1);
                    L2.insert(pbreaker,l2);
                    pbreaker+=1;

                highstep+=1;


        

        step+=1;
        

    minp = min(pvals)

    minpval = pvals.index(minp)

    h1 = H1[minpval]
    l1 = 0;
    h2 = max(rawData);
    l2 = L2[minpval]

    kaplandata = [];
    kaplancen = [];
    

    
    xxx = int(max(Survival));
    breaker1=0;
    
    
    x = 0;
    xx = len(rawData)
    while(xx>0):
        if(h1 >= rawData[x] >= l1):
            
            kaplandata.insert(breaker1,Survival[x])
            kaplancen.insert(breaker1,AliveDead[x])
            breaker1+=1;
            
        x+=1;
        
        xx-=1;


    x=0;
    xxx-=1;

    
    

    kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
    idsdays1 = kapper[1]
    kap1per = kapper[0]
    censoreds1 = kapper[2]
    censoredids1 = kapper[3]


    kaplandata = [];
    kaplancen = [];
    

    
    xxx = int(max(Survival));
    breaker1=0;
    
    
    x = 0;
    xx = len(rawData)
    while(xx>0):
        if(h2 >= rawData[x] >= l2):
            
            kaplandata.insert(breaker1,Survival[x])
            kaplancen.insert(breaker1,AliveDead[x])
            breaker1+=1;
            
        x+=1;
        
        xx-=1;


    x=0;
    xxx-=1;

    
    

    kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
    idsdays2 = kapper[1]
    kap2per = kapper[0]
    censoreds2 = kapper[2]
    censoredids2 = kapper[3]

   



    kaplab1 = entry+' Expression from ' + str(l1) +' to ' +str(h1)
    P.plot(idsdays1,kap1per,linestyle="-",color="green",label= kaplab1)
    P.scatter(censoredids1,censoreds1,marker ='x',color='black')
    kaplab2 = entry+' Expression from ' + str(l2) +' to ' +str(h2)   
    P.plot(idsdays2,kap2per,linestyle="-",color="blue",label= kaplab2)
    P.scatter(censoredids2,censoreds2,marker ='x',color='black')


    P.grid(True)

    P.legend()

    
    P.title('Survival with ' + entry +"  p-value:"+str(minp),fontsize=24);
    P.xlabel("Survival Duration (days)",fontsize=20);
    P.ylabel("Percent Survival",fontsize=20);
   
        
        
            
    P.show()




