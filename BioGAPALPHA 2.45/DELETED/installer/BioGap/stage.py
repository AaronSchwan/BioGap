from statistics import mean
from statistics import median
import numpy as np
import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

def stage_data(entry,dataS,rawData,StageData,StageColors,medians,StageTitles,lengths):
    #stages
    stage1color = StageColors[0];
    stage2color = StageColors[1];
    stage3color = StageColors[2];
    stage4color = StageColors[3];

    stagetitle = StageTitles[0];
    stagexaxis = StageTitles[1];
    stageyaxis = StageTitles[2];
    
    


    StageDataLength = len(StageData)-1
    dinsert =0;

    StageOne = [];
    StageTwo = [];
    StageThree = [];
    StageFour = [];


    onebreaker = 0;
    twobreaker = 0;
    threebreaker = 0;
    fourbreaker = 0;
        
    while(StageDataLength>0):
        if(StageData[dinsert] == 1.0):
            StageOne.insert(onebreaker,rawData[dinsert])
            onebreaker+=1;
        if(StageData[dinsert] == 1.5):
            StageOne.insert(onebreaker,rawData[dinsert])
            onebreaker+=1;
        if(StageData[dinsert] == 2.0):
            StageTwo.insert(twobreaker,rawData[dinsert])
            twobreaker+=1;
        if(StageData[dinsert] == 2.5):
            StageTwo.insert(twobreaker,rawData[dinsert])
            twobreaker+=1;
        if(StageData[dinsert] == 3.0):
            StageThree.insert(threebreaker,rawData[dinsert])
            threebreaker+=1;
        if(StageData[dinsert] == 3.5):
            StageThree.insert(threebreaker,rawData[dinsert])
            threebreaker+=1;
        if(StageData[dinsert] == 4.0):
            StageFour.insert(fourbreaker,rawData[dinsert])
            fourbreaker+=1;

        StageDataLength-=1;
        dinsert+=1;
    P.figure()
    length = lengths*.5
    if(len(StageOne) != 0):
        idsone = [x for x in range(len(StageOne))]
        P.scatter(idsone,sorted(StageOne),color=stage1color,label="Stage One")
        if(medians == 1):
            midlinex = [0,length]
            med1 = [median(StageOne),median(StageOne)];
            P.plot(midlinex,med1,color=stage1color,linestyle="--")

    if(len(StageTwo) != 0):
        idstwo = [x for x in range(len(StageTwo))]
        P.scatter(idstwo,sorted(StageTwo),color=stage2color,label="Stage Two")
        if(medians == 1):
            midlinex = [0,length]
            med2 = [median(StageTwo),median(StageTwo)];
            P.plot(midlinex,med2,color=stage2color,linestyle="--")

    if(len(StageThree) != 0):
        idsthree = [x for x in range(len(StageThree))]
        P.scatter(idsthree,sorted(StageThree),color=stage3color,label="Stage Three")
        if(medians == 1):
            midlinex = [0,length]
            med3 = [median(StageThree),median(StageThree)];
            P.plot(midlinex,med3,color=stage3color,linestyle="--")

    if(len(StageFour) != 0):
        idsfour = [x for x in range(len(StageFour))]
        P.scatter(idsfour,sorted(StageFour),color=stage4color,label="Stage Four")
        if(medians == 1):
            midlinex = [0,length]
            med4 = [median(StageFour),median(StageFour)];
            P.plot(midlinex,med4,color=stage4color,linestyle="--")
        
    

    
    
    if(stagetitle == "DEFAULT"):
        titl =  " Stages of "+ dataS +" With the Expression of " +entry +"
  "
        P.title(titl)
        
    if(stagetitle != "DEFAULT"):
        P.title(stagetitle)
        

    if(stagexaxis == "DEFAULT"):
         P.xlabel('Number of Patients Per Stage')
    if(stagexaxis != "DEFAULT"):
        P.xlabel(stagexaxis)

    if(stageyaxis == "DEFAULT"):
        P.ylabel('Expression of '+ entry+' in RNA Seq V2')
    if(stageyaxis != "DEFAULT"):
        P.ylabel(stageyaxis)
  
    
    P.grid(True)
    P.legend()
    
    P.show()                
