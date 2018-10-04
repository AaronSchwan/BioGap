from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re


from numpy import arange, sin, pi
import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler

import Kaplan_Meier_Graphing_expected
import bestkaplansalltwoline
import bestkaplanstwoline
import stageFR
import tcga_genesearch


normalbackground = "#cce6ff"
sexbackground = "#d9f2d9"
agebackground = "#ffffcc"
stagebackground = "#ffeecc"

entry = "TNFSF9"
organ = "Lung"
outliers = "Include Outliers"

##if(organ == "Lung"):
##    #normal
##    STAGES=[1,2,3,4]
##    AGES=[0,10000]
##    SEXES=[1,0]
##    RACES=['white','not reported','black or african american','asian']
##    StageData,rawData,RaceData,SexData,Survival,Censored,AliveDead,AgeData = tcga_genesearch.geneSearch(entry,outliers,"Lung Squamous Cell Carcinoma (TCGA, Provisional)",STAGES,AGES,SEXES,RACES)
##    
##    medians = 1;
##
##    with open("Program Setup\Colors\Stage.csv",'r') as file:
##        StageColors = file.read().split('\n')
##        file.close();
##
##    with open("Program Setup\Titles\Stage.csv",'r') as file:
##        StageTitles = file.read().split('\n')
##        file.close();
##        
##    lengths = int(len(rawData)*2);
##    
##    P = stageFR.stage_data(entry,"Lung Squamous Cell Carcinoma (TCGA, Provisional)",rawData,StageData,StageColors,medians,StageTitles,lengths)
##    canvas = FigureCanvasTkAgg(P, master=root)
##    canvas.show()
    

class Window(Frame):
    

    def __init__(self, master=None):

        Frame.__init__(self, master)
        self.master = master
        
        
        self.init_window()
        
        

    def init_window(self):



        self.master.title("Biomarker Graphical Anaylisis Program")

        maincanvas=Canvas(bg="#e6ffff",width=200,height=200,scrollregion=(0,0,1300,1300))
        
        


        

        vbar=Scrollbar(root,orient=VERTICAL)
        vbar.pack(side=RIGHT,fill=Y)
        vbar.config(command=maincanvas.yview)

        maincanvas.config(width=200,height=200)
        LABEL=Canvas(bg="WHITE",width=10,height=20,scrollregion=(0,0,1300,1300))
        TITLES = Label(text="Aaron Schwan alpha 1.10")
        TITLES.configure(bg="WHITE")
        TITLES_window = LABEL.create_window(5,5,anchor=NW, window=TITLES)
        LABEL.pack(side=BOTTOM,fill=BOTH)
        

        maincanvas.config(yscrollcommand=vbar.set)
        maincanvas.pack(side=BOTTOM,expand=True,fill=BOTH)

        if(organ == "Lung"):
            #normal
            STAGES=[1]
            AGES=[0,10000]
            SEXES=[1,0]
            RACES=['white','not reported','black or african american','asian']
            
            StageData,rawData,RaceData,SexData,Survival,Censored,AliveDead,AgeData = tcga_genesearch.geneSearch(entry,outliers,"Lung Squamous Cell Carcinoma (TCGA, Provisional)",STAGES,AGES,SEXES,RACES)
            
            medians = 1;
            
            with open("Program Setup\Colors\Stage.csv",'r') as file:
                StageColors = file.read().split('\n')
                file.close();
            '''
            with open("Program SetupSetup\Titles\Stage.csv",'r') as file:
                StageTitles = file.read().split('\n')
                file.close();
            '''
            lengths = int(len(rawData)*2);

            P = stageFR.stage_data(entry,"Lung Squamous Cell Carcinoma (TCGA, Provisional)",rawData,StageData,StageColors,medians,lengths)
            



            

                        
            
            

        def Normal():
            #canvas.delete('all')

##            P = stageFR.stage_data(entry,"Lung Squamous Cell Carcinoma (TCGA, Provisional)",rawData,StageData,StageColors,medians,lengths)
##            canvas = FigureCanvasTkAgg(P, master=root)
##            canvas.show()
##            
            f = P
            a = f.add_subplot(111)
            t = arange(0.0,3.0,0.01)
            s = sin(2*pi*t)
            a.plot(P)


            dataPlot = FigureCanvasTkAgg(f, master=root)
            dataPlot.show()
##             dataPlot.get_tk_widget()




            maincanvas.create_window(5,5,anchor=NW ,tag = 'all',window=dataPlot.get_tk_widget())


            maincanvas.config(bg=normalbackground)
            NormalB.configure(bg = "#33B5E5")
            SexB.configure(bg = "#edefef")
            AgeB.configure(bg = "#edefef")
            StageB.configure(bg = "#edefef")
            
        def Sex():
            canvas.delete('all')


            canvas.config(bg=sexbackground)
            NormalB.configure(bg = "#edefef")
            SexB.configure(bg = "#33B5E5")
            AgeB.configure(bg = "#edefef")
            StageB.configure(bg = "#edefef")

        def Age():
            canvas.delete('all')

            canvas.config(bg=agebackground)
            NormalB.configure(bg = "#edefef")
            SexB.configure(bg = "#edefef")
            AgeB.configure(bg = "#33B5E5")
            StageB.configure(bg = "#edefef")

        def Stage():
            canvas.delete('all')

            canvas.config(bg=stagebackground)
            NormalB.configure(bg = "#edefef")
            SexB.configure(bg = "#edefef")
            AgeB.configure(bg = "#edefef")
            StageB.configure(bg = "#33B5E5")
            
        
        #Choosing Database
        choosebuttons=Canvas(bg="GREY",width=200,height=50)

        NormalB = Button(text = "Normal",command = Normal ,anchor = W)
        NormalB.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        NormalB_window = choosebuttons.create_window(10, 10, anchor=NW, window=NormalB)

        SexB = Button(text = "Sex Seperation",command = Sex, anchor = W)
        SexB.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        SexB_window = choosebuttons.create_window(100, 10, anchor=NW, window=SexB)

        AgeB = Button(text = "Age Seperation",command = Age, anchor = W)
        AgeB.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        AgeB_window = choosebuttons.create_window(190, 10, anchor=NW, window=AgeB)

        StageB = Button(text = "Stage Seperation",command = Stage, anchor = W)
        StageB.configure(width = 13, activebackground = "#33B5E5", relief = FLAT)
        StageB_window = choosebuttons.create_window(280, 10, anchor=NW, window=StageB)

        choosebuttons.pack(side=TOP,fill=BOTH)





        
        



root = Tk()

root.geometry("1200x900")
app = Window(root)   
root.mainloop()
