
import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import Kaplan_Meier_Graphing


def kaplanmeiers(entry,dataS,rawData,AliveDead,Survival):
    with open("Program Setup\KaplanMeier1-10\KaplanMeierLabels.csv",'r') as file:
        labels = file.read().split('
')
        file.close()

    kaplantitle = labels[0]
    kaplanxaxis = labels[1]
    kaplanyaxis = labels[2]

    kaplabel1 = labels[3]
    kaplabel2 = labels[4]
    kaplabel3 = labels[5]
    kaplabel4 = labels[6]
    kaplabel5 = labels[7]
    kaplabel6 = labels[8]
    kaplabel7 = labels[9]
    kaplabel8 = labels[10]
    kaplabel9 = labels[11]
    kaplabel10 = labels[12]

    with open("Program Setup\KaplanMeier1-10\KaplanMeierLinesColors.csv",'r') as file:
        Colors = file.read().split('
')
        file.close()

    kaplan1color = Colors[0]
    kaplan2color = Colors[1]
    kaplan3color = Colors[2]
    kaplan4color = Colors[3]
    kaplan5color = Colors[4]
    kaplan6color = Colors[5]
    kaplan7color = Colors[6]
    kaplan8color = Colors[7]
    kaplan9color = Colors[8]
    kaplan10color = Colors[9]
    

    with open("Program Setup\KaplanMeier1-10\KaplanMeierLines.csv",'r') as file:
        lines = file.read().split('
')
        file.close()
    
    with open("Program Setup\KaplanMeier1-10\KaplanMeierLinesLow.csv",'r') as file:
        LOWST = file.read().split('
')
        file.close()

    with open("Program Setup\KaplanMeier1-10\KaplanMeierLinesHigh.csv",'r') as file:
        HIGHST = file.read().split('
')
        file.close()



    LOWS =[];
    HIGHS = [];

    x =0
    xx = len(LOWST)-1;
    while(xx>0):
       
        if(LOWST[x] != ''):
            LOWS.insert(x,float(LOWST[x]))
        if(LOWST[x] == ''):
            LOWS.insert(x,0)
        x+=1;
        xx-=1;

    x =0
    xx = len(HIGHST)-1;
    while(xx>0):
        if(HIGHST[x] != ''):
            HIGHS.insert(x,float(HIGHST[x]))
        if(HIGHST[x] == ''):
            HIGHS.insert(x,0)
        x+=1;
        xx-=1;

    

    
    
    P.figure()
    if(lines[0] == '1'):
        

       
        l1 = float(LOWS[0])
        h1 = float(HIGHS[0])

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
        idsdays = kapper[1]
        kap1per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

        if(kaplabel1 == "DEFAULT"):
            kaplab1 = entry+' Expression from ' + str(l1) +' to ' +str(h1)

        if(kaplabel1 != "DEFAULT"):
            kaplab1 = kaplabel1;
            
        P.plot(idsdays,kap1per,linestyle="-",color=kaplan1color,label= kaplab1)
        P.scatter(censoredids,censoreds,marker ='x',color='black')

       
    
    if(lines[1] == '1'):
        
        l2 = float(LOWS[1])
        h2 = float(HIGHS[1])
        kaplandata = [];
        kaplancen = [];
        xxx = int(max(Survival));
        breaker2=0;
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h2 >= rawData[x] >= l2):
                kaplandata.insert(breaker2,Survival[x])
                kaplancen.insert(breaker2,AliveDead[x])
                breaker2+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays2 = kapper[1]
        kap2per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

        
        if(kaplabel2 == "DEFAULT"):
            kaplab2 = entry+' Expression from ' + str(l2) +' to ' +str(h2)

        if(kaplabel2 != "DEFAULT"):
            kaplab2 = kaplabel2;
            
        P.plot(idsdays2,kap2per,linestyle="-",color=kaplan2color,label= kaplab2)
        P.scatter(censoredids,censoreds,marker ='x',color='black')


        
    if(lines[2] == '1'):
        
        l3 = float(LOWS[2])
        h3 = float(HIGHS[2])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker3=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h3 >= rawData[x] >= l3):
                kaplandata.insert(breaker3,int(Survival[x]))
                kaplancen.insert(breaker3,int(AliveDead[x]))
                breaker3+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays3 = kapper[1]
        kap3per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

            
        if(kaplabel3 == "DEFAULT"):
            kaplab3 = entry+' Expression from ' + str(l3) +' to ' +str(h3)

        if(kaplabel3 != "DEFAULT"):
            kaplab3 = kaplabel3;
            
        P.plot(idsdays3,kap3per,linestyle="-",color=kaplan3color,label= kaplab3)
        P.scatter(censoredids,censoreds,marker ='+',color='black')

    if(lines[3] == '1'):
        
        l4 = float(LOWS[3])
        h4 = float(HIGHS[3])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker4=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h4 >= rawData[x] >= l4):
                kaplandata.insert(breaker4,int(Survival[x]))
                kaplancen.insert(breaker4,int(AliveDead[x]))
                breaker4+=1;
                
            x+=1;
            xx-=1;

        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays4 = kapper[1]
        kap4per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

        if(kaplabel4 == "DEFAULT"):
            kaplab4 = entry+' Expression from ' + str(l4) +' to ' +str(h4)

        if(kaplabel4 != "DEFAULT"):
            kaplab4 = kaplabel4;
            
        P.plot(idsdays4,kap4per,linestyle="-",color=kaplan4color,label= kaplab4)
        P.scatter(censoredids,censoreds,marker ='+',color='black')
        
    if(lines[4] == '1'):
        
        l5 = float(LOWS[4])
        h5 = float(HIGHS[4])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker5=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h5 >= rawData[x] >= l5):
                kaplandata.insert(breaker5,int(Survival[x]))
                kaplancen.insert(breaker5,int(AliveDead[x]))
                breaker5+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays5 = kapper[1]
        kap5per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]


            
        if(kaplabel5 == "DEFAULT"):
            kaplab5 = entry+' Expression from ' + str(l5) +' to ' +str(h5)

        if(kaplabel5 != "DEFAULT"):
            kaplab5 = kaplabel5;
            
        P.plot(idsdays5,kap5per,linestyle="-",color=kaplan5color,label= kaplab5)
        P.scatter(censoredids,censoreds,marker ='+',color='black')


    if(lines[5] == '1'):
        
        l6 = float(LOWS[5])
        h6 = float(HIGHS[5])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker6=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h6 >= rawData[x] >= l6):
                kaplandata.insert(breaker6,int(Survival[x]))
                kaplancen.insert(breaker6,int(AliveDead[x]))
                breaker6+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays6 = kapper[1]
        kap6per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

            
        if(kaplabel6 == "DEFAULT"):
            kaplab6 = entry+' Expression from ' + str(l6) +' to ' +str(h6)

        if(kaplabel6 != "DEFAULT"):
            kaplab6 = kaplabel6;
            
        P.plot(idsdays6,kap6per,linestyle="-",color=kaplan6color,label= kaplab6)
        P.scatter(censoredids,censoreds,marker ='+',color='black')

    if(lines[6] == '1'):
        
        l7 = float(LOWS[6])
        h7 = float(HIGHS[6])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker7=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h7 >= rawData[x] >= l7):
                kaplandata.insert(breaker7,int(Survival[x]))
                kaplancen.insert(breaker7,int(AliveDead[x]))
                breaker7+=1;
                
            x+=1;
            xx-=1;

        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays7 = kapper[1]
        kap7per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]

            
        if(kaplabel7 == "DEFAULT"):
            kaplab7 = entry+' Expression from ' + str(l7) +' to ' +str(h7)

        if(kaplabel7 != "DEFAULT"):
            kaplab7 = kaplabel7;
            
        P.plot(idsdays7,kap7per,linestyle="-",color=kaplan7color,label= kaplab7)
        P.scatter(censoredids,censoreds,marker ='+',color='black')

    if(lines[7] == '1'):
        
        l8 = float(LOWS[7])
        h8 = float(HIGHS[7])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker8=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h8 >= rawData[x] >= l8):
                kaplandata.insert(breaker8,int(Survival[x]))
                kaplancen.insert(breaker8,int(AliveDead[x]))
                breaker8+=1;
                
            x+=1;
            xx-=1;

        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays8 = kapper[1]
        kap8per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]


        if(kaplabel8 == "DEFAULT"):
            kaplab8 = entry+' Expression from ' + str(l8) +' to ' +str(h8)

        if(kaplabel8 != "DEFAULT"):
            kaplab8 = kaplabel8;
            
        P.plot(idsdays8,kap8per,linestyle="-",color=kaplan8color,label= kaplab8)
        P.scatter(censoredids,censoreds,marker ='+',color='black')

    if(lines[8] == '1'):
        
        l9 = float(LOWS[8])
        h9 = float(HIGHS[8])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker9=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h9 >= rawData[x] >= l9):
                kaplandata.insert(breaker9,int(Survival[x]))
                kaplancen.insert(breaker9,int(AliveDead[x]))
                breaker9+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays9 = kapper[1]
        kap9per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]
        
        if(kaplabel9 == "DEFAULT"):
            kaplab9 = entry+' Expression from ' + str(l9) +' to ' +str(h9)

        if(kaplabel9 != "DEFAULT"):
            kaplab9 = kaplabel9;
            
        P.plot(idsdays9,kap9per,linestyle="-",color=kaplan9color,label= kaplab9)
        P.scatter(censoredids,censoreds,marker ='+',color='black')

    if(lines[9] == '1'):
        
        l10 = float(LOWS[9])
        h10 = float(HIGHS[9])

        kaplandata = [];
        kaplancen = [];
        

        
        xxx = int(max(Survival));
        breaker10=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h10 >= rawData[x] >= l10):
                kaplandata.insert(breaker10,int(Survival[x]))
                kaplancen.insert(breaker10,int(AliveDead[x]))
                breaker10+=1;
                
            x+=1;
            xx-=1;


        kapper = Kaplan_Meier_Graphing.kaplan(kaplandata,kaplancen)
        idsdays10 = kapper[1]
        kap10per = kapper[0]
        censoreds = kapper[2]
        censoredids = kapper[3]
            
        if(kaplabel10 == "DEFAULT"):
            kaplab10 = entry+' Expression from ' + str(l10) +' to ' +str(h10)

        if(kaplabel10 != "DEFAULT"):
            kaplab10 = kaplabel10;
            
        P.plot(idsdays10,kap10per,linestyle="-",color=kaplan10color,label= kaplab10)
        P.scatter(censoredids,censoreds,marker ='+',color='black')


    P.grid(True)

    P.legend()

    if(kaplantitle == "DEFAULT"):
        P.title('Survival with ' + entry,fontsize=24);
    if(kaplantitle != "DEFAULT"):
        P.title(kaplantitle,fontsize=24);
    if(kaplanxaxis == "DEFAULT"):
        P.xlabel("Survival Duration (days)",fontsize=20);
    if(kaplanxaxis != "DEFAULT"):
        P.xlabel(kaplanxaxis,fontsize=20);
    if(kaplanyaxis == "DEFAULT"):
        P.ylabel("Percent Survival",fontsize=20);
    if(kaplanyaxis != "DEFAULT"):
        P.ylabel(kaplanyaxis,fontsize=20);
        
        
            
    P.show()
