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



class Window(Frame):


    def __init__(self, master=None):

        Frame.__init__(self, master)
        self.master = master
        self.frame = Scrollbar(self)
        self.frame.grid(row=1,column=100)
        
        self.init_window()
        
        

    def init_window(self):

        self.master.title("Biomarker Graphical Anaylisis Program")
        self.pack(fill = BOTH, expand=1)
        
        
     

        

        #defaults
        global stagetitle
        global stagexaxis
        global stageyaxis
        global stage1color
        global stage2color
        global stage3color
        global stage4color
        global stagemark1
        global stagemark2
        global stagemark3
        global stagemark4
        
        
        stagemark1 ="o";
        stagemark2 ="o";
        stagemark3 ="o";
        stagemark4 ="o";
     
        stagetitle = "DEFAULT";
        stagexaxis = "DEFAULT";
        stageyaxis = "DEFAULT";
        
        stage1color = "lightgreen";
        stage2color = "Green";
        stage3color = "lightblue";
        stage4color = "Blue";

        global kaplantitle
        global kaplanxaxis
        global kaplanyaxis
        global kaplan1color
        global kaplan2color
        global kaplan3color
        global kaplan4color
        global kaplan5color
        global kaplan6color
        global kaplan7color
        global kaplan8color
        global kaplan9color
        global kaplan10color

        
        
       
     
        kaplantitle = "DEFAULT";
        kaplanxaxis = "DEFAULT";
        kaplanyaxis = "DEFAULT";
        
        kaplan1color = "lightgreen";
        kaplan2color = "Green";
        kaplan3color = "lightblue";
        kaplan4color = "Blue";
        kaplan5color = "Red";
        kaplan6color = "Purple";
        kaplan7color = "Yellow";
        kaplan8color = "Orange";
        kaplan9color = "Magenta";
        kaplan10color = "Black";
        


        global kaplabel1
        global kaplabel2
        global kaplabel3
        global kaplabel4
        global kaplabel5
        global kaplabel6
        global kaplabel7
        global kaplabel8
        global kaplabel9
        global kaplabel10
        

        kaplabel1 = "DEFAULT";
        kaplabel2 = "DEFAULT";
        kaplabel3 = "DEFAULT";
        kaplabel4 = "DEFAULT";
        kaplabel5 = "DEFAULT";
        kaplabel6 = "DEFAULT";
        kaplabel7 = "DEFAULT";
        kaplabel8 = "DEFAULT";
        kaplabel9 = "DEFAULT";
        kaplabel10 = "DEFAULT";
        









        #menu bar
        
        menu = Menu(self.master)
        self.master.config(menu=menu)

        
        file = Menu(menu)
        file.add_command(label="Help", command=self.client_help)
        file.add_command(label="Exit", command=self.client_exit)
        menu.add_cascade(label="File", menu=file)

        graphs = Menu(menu)
        graphs.add_command(label="Stage Graph", command=self.stage_graph_change)
        graphs.add_command(label="Kaplan Meier", command=self.kaplan_graph_change)
        
        
        menu.add_cascade(label="Graph Options", menu=graphs)
      
        
        
       

        #data source choose
        global var
        global out
        
        var = StringVar(self)
        var.set("Data Source:") # initial value
        
        option = OptionMenu(self, var, "Lung Squamous Cell Carcinoma (TCGA, Provisional)", "Lung Adenocarcinoma (TCGA, Provisional)")
        option.grid(row=1,column=1)
        
        out = StringVar(self)
        out.set("Include Outliers") # initial value

        optionout = OptionMenu(self, out, "Include Outliers", "Exclude Outliers")
        optionout.grid(row=1,column=2)
        
        #gene Searching        
        Label(self,text="Search Gene:",font=("Helvetica", 10)).grid(row=1,column=3)
        self.geneSearch = Entry(self)
        self.geneSearch.grid(row=1,column=4)
        geneSearchB = Button(self,text = "Search Gene",font=("Helvetica", 10),command = self.gene_search).grid(row=1,column=5)
        
        self.tell = Entry(self,width=40)
        self.tell.grid(row=1,column=6)
        self.tell.config(state=DISABLED)
        # spacer
        Label(self,text="",font=("Helvetica", 10)).grid(row=135,column=1)
        Label(self,text="Graphs And Data:",font= "Helvetica 10 underline").grid(row=136,column=1)
        Label(self,text="",font=("Helvetica", 10)).grid(row=137,column=1)

        #raw data display in entry
        rawdatadisplyret = Button(self,text = "Raw Data Download",command = self.raw_data_ret,font=("Helvetica", 10)).grid(row=138,column=2)
        rawdatadisply = Button(self,text = "Patient Data",command = self.raw_data,font=("Helvetica", 10)).grid(row=138,column=3)
        stagedatadisplay = Button(self,text = "Stage Data",command = self.stage_data,font=("Helvetica", 10)).grid(row=138,column=4)


        #Kaplan Meiers
        Label(self,text="",font="Helvetica 10 underline").grid(row=21,column=1)
        Label(self,text="Kaplan Meier Graph",font="Helvetica 10 underline").grid(row=21,column=1)
        Label(self,text="Line One",font="Helvetica 10").grid(row=23,column=1)
        Label(self,text="Line Two",font="Helvetica 10").grid(row=24,column=1)
        Label(self,text="Line Three",font="Helvetica 10").grid(row=25,column=1)
        Label(self,text="Line Four",font="Helvetica 10").grid(row=26,column=1)
        Label(self,text="Line Five",font="Helvetica 10").grid(row=27,column=1)
        Label(self,text="Line Six",font="Helvetica 10").grid(row=28,column=1)
        Label(self,text="Line Seven",font="Helvetica 10").grid(row=29,column=1)
        Label(self,text="Line Eight",font="Helvetica 10").grid(row=30,column=1)
        Label(self,text="Line Nine",font="Helvetica 10").grid(row=31,column=1)
        Label(self,text="Line Ten",font="Helvetica 10").grid(row=32,column=1)
        #Label(self,text="P Value",font="Helvetica 10").grid(row=33,column=1)
        
        Label(self,text="",font="Helvetica 10 underline").grid(row=22,column=1)
                        
        kaplan = Button(self,text = "Kaplan Meier",command = self.kaplan_data,font=("Helvetica", 10)).grid(row=23,column=4)

        chiget = Button(self,text = "p-value For",command = self.chi_squared,font=("Helvetica", 10)).grid(row=33,column=1)
        bestall = Button(self,text = "Include All Patients",command = self.mid_kap,font=("Helvetica", 10)).grid(row=24,column=4)
        bestgap = Button(self,text = "Best p-value",command = self.best_kap,font=("Helvetica", 10)).grid(row=25,column=4)
        
        global kap
        kap = StringVar(self)
        kap.set("One Line") # initial value

        optionkap = OptionMenu(self, kap, "One Line", "Two Lines", "Three Lines", "Four Lines", "Five Lines", "Six Lines", "Seven Lines", "Eight Lines", "Nine Lines", "Ten Lines")
        optionkap.grid(row=22,column=1)

        global chi1
        chi1 = StringVar(self)
        chi1.set("Line One And") # initial value

        optionchi1 = OptionMenu(self, chi1, "Line One And", "Line Two And", "Line Three And", "Line Four And", "Line Five And", "Line Six And", "Line Seven And", "Line Eight And", "Line Nine And", "Line Ten And")
        optionchi1.grid(row=33,column=2)

        global chi2
        chi2 = StringVar(self)
        chi2.set("Line One") # initial value

        optionchi2 = OptionMenu(self, chi2, "Line One", "Line Two", "Line Three", "Line Four", "Line Five", "Line Six", "Line Seven", "Line Eight", "Line Nine", "Line Ten")
        optionchi2.grid(row=33,column=3)

        self.chisq = Entry(self,width=40)
        self.chisq.grid(row=33,column=4)
        self.chisq.config(state=DISABLED)

        global kaplowone
        global kaphighone
        
        Label(self,text="Low Range").grid(row=22,column=2)
        kaplowone = Entry(self)
        kaplowone.grid(row=23,column=2)
        Label(self,text="High Range").grid(row=22,column=3)
        kaphighone = Entry(self)
        kaphighone.grid(row=23,column=3)

        global kaplowtwo
        global kaphightwo
        
        
        kaplowtwo = Entry(self)
        kaplowtwo.grid(row=24,column=2)
        
        kaphightwo = Entry(self)
        kaphightwo.grid(row=24,column=3)

        global kaplowthree
        global kaphighthree
        
        
        kaplowthree = Entry(self)
        kaplowthree.grid(row=25,column=2)
        
        kaphighthree = Entry(self)
        kaphighthree.grid(row=25,column=3)

        global kaplowfour
        global kaphighfour
        
        
        kaplowfour = Entry(self)
        kaplowfour.grid(row=26,column=2)
        
        kaphighfour = Entry(self)
        kaphighfour.grid(row=26,column=3)

        global kaplowfive
        global kaphighfive
        
        
        kaplowfive = Entry(self)
        kaplowfive.grid(row=27,column=2)
        
        kaphighfive = Entry(self)
        kaphighfive.grid(row=27,column=3)

        global kaplowsix
        global kaphighsix
        
        
        kaplowsix = Entry(self)
        kaplowsix.grid(row=28,column=2)
        
        kaphighsix = Entry(self)
        kaphighsix.grid(row=28,column=3)

        global kaplowseven
        global kaphighseven
        
        
        kaplowseven = Entry(self)
        kaplowseven.grid(row=29,column=2)
        
        kaphighseven = Entry(self)
        kaphighseven.grid(row=29,column=3)

        global kaploweight
        global kaphigheight
        
        
        kaploweight= Entry(self)
        kaploweight.grid(row=30,column=2)
        
        kaphigheight = Entry(self)
        kaphigheight.grid(row=30,column=3)

        global kaplownine
        global kaphighnine
        
        
        kaplownine = Entry(self)
        kaplownine.grid(row=31,column=2)
        
        kaphighnine = Entry(self)
        kaphighnine.grid(row=31,column=3)

        global kaplowten
        global kaphighten
        
        
        kaplowten = Entry(self)
        kaplowten.grid(row=32,column=2)
        
        kaphighten = Entry(self)
        kaphighten.grid(row=32,column=3)




        # spacer
        Label(self,text="",font=("Helvetica", 10)).grid(row=3,column=1)
        Label(self,text="Normal Statistical Measures:",font= "Helvetica 10 underline").grid(row=4,column=1)
        Label(self,text="",font=("Helvetica", 10)).grid(row=5,column=1)
        
        #average
        Label(self,text="Average:",font=("Helvetica", 10)).grid(row=6,column=2)
        self.average = Entry(self)
        self.average.grid(row=6,column=3)
        self.average.config(state=DISABLED)

        #standerd dev
        Label(self,text="Standard Deviation:",font=("Helvetica", 10)).grid(row=7,column=2)
        self.standerddev = Entry(self)
        self.standerddev.grid(row=7,column=3)
        self.standerddev.config(state=DISABLED)

        #IQR
        Label(self,text="Interquartile Range:",font=("Helvetica", 10)).grid(row=8,column=2)
        self.iqr = Entry(self)
        self.iqr.grid(row=8,column=3)
        self.iqr.config(state=DISABLED)

        #range
        Label(self,text="Range:",font=("Helvetica", 10)).grid(row=9,column=2)
        self.range = Entry(self)
        self.range.grid(row=9,column=3)
        self.range.config(state=DISABLED)

        #Min
        Label(self,text="Minimum:",font=("Helvetica", 10)).grid(row=10,column=2)
        self.minimum = Entry(self)
        self.minimum.grid(row=10,column=3)
        self.minimum.config(state=DISABLED)

        #Q1
        Label(self,text="First Quartile Median:",font=("Helvetica", 10)).grid(row=11,column=2)
        self.q1 = Entry(self)
        self.q1.grid(row=11,column=3)
        self.q1.config(state=DISABLED)
        #Med
        Label(self,text="Median:",font=("Helvetica", 10)).grid(row=12,column=2)
        self.median = Entry(self)
        self.median.grid(row=12,column=3)
        self.median.config(state=DISABLED)
        
        #Q3
        Label(self,text="Third Quartile Median:",font=("Helvetica", 10)).grid(row=13,column=2)
        self.q3 = Entry(self)
        self.q3.grid(row=13,column=3)
        self.q3.config(state=DISABLED)
        
        #Max
        Label(self,text="Maximum:",font=("Helvetica", 10)).grid(row=14,column=2)
        self.maximum = Entry(self)
        self.maximum.grid(row=14,column=3)
        self.maximum.config(state=DISABLED)

        # spacer
        Label(self,text="",font=("Helvetica", 10)).grid(row=15,column=1)
        Label(self,text="Correlation Data:",font= "Helvetica 10 underline").grid(row=16,column=1)
        Label(self,text="",font=("Helvetica", 10)).grid(row=17,column=1)

        #correlation
        Label(self,text="Correlate to Gene:",font=("Helvetica", 10)).grid(row=18,column=2)
        self.correlateGene = Entry(self)
        self.correlateGene.grid(row=18,column=3)
        geneSearchBCorrelate = Button(self,text = "Correlate",command = self.r_squared,font=("Helvetica", 10)).grid(row=18,column=4)
        geneCorrelate = Button(self,text = "Correlate To All",command = self.r_squared_all,font=("Helvetica", 10)).grid(row=19,column=4)

        Label(self,text="Linear Regression(r^2):",font=("Helvetica", 10)).grid(row=19,column=2)
        self.r_squaredE = Entry(self)
        self.r_squaredE.grid(row=19,column=3)
        self.r_squaredE.config(state=DISABLED)

        Label(self,text="t-test:",font=("Helvetica", 10)).grid(row=20,column=2)
        self.tTest = Entry(self)
        self.tTest.grid(row=20,column=3)
        self.tTest.config(state=DISABLED)


        
    def client_exit(self):
        exit()
    def client_help(self):
        os.system("Help.pdf");

    def stage_graph_change(self):
        
        t = Toplevel(self)
        t.wm_title("Stage Graph")
        t.wm_geometry("400x250")
        Label(t, text="Chart Title:").grid(row=1,column=1)

        global stagecharttitle
        stagecharttitle = Entry(t)
        stagecharttitle.delete(0,END)
        stagecharttitle.insert(0,stagetitle)
        stagecharttitle.grid(row=1,column=2)

        global stagechartxaxis
        Label(t, text="X Axis:").grid(row=2,column=1)
        stagechartxaxis = Entry(t)
        stagechartxaxis.delete(0,END)
        stagechartxaxis.insert(0,stagexaxis)
        stagechartxaxis.grid(row=2,column=2)

        global stagechartyaxis
        Label(t, text="Y Axis:").grid(row=3,column=1)
        stagechartyaxis = Entry(t)
        stagechartyaxis.delete(0,END)
        stagechartyaxis.insert(0,stageyaxis)
        stagechartyaxis.grid(row=3,column=2)

        
        Label(t, text="Stage One Color:").grid(row=4,column=1)
        global color1
        color1 = StringVar(self)
        color1.set(stage1color) # initial value

        optioncolor1 = OptionMenu(t, color1, "Blue", "Light Blue", "Green", "Light Green","Red","Orange", "Purple","Yellow")
        optioncolor1.grid(row=4,column=2)

        Label(t, text="Stage Two Color:").grid(row=5,column=1)
        global color2
        color2 = StringVar(self)
        color2.set(stage2color) # initial value

        optioncolor2 = OptionMenu(t, color2, "Blue", "Light Blue", "Green", "Light Green","Red","Orange", "Purple","Yellow")
        optioncolor2.grid(row=5,column=2)

        Label(t, text="Stage Three Color:").grid(row=6,column=1)
        global color3
        color3 = StringVar(self)
        color3.set(stage3color) # initial value

        optioncolor3 = OptionMenu(t, color3, "Blue", "Light Blue", "Green", "Light Green","Red","Orange", "Purple","Yellow")
        optioncolor3.grid(row=6,column=2)

        Label(t, text="Stage Four Color:").grid(row=7,column=1)
        global color4
        color4 = StringVar(self)
        color4.set(stage4color) # initial value

        optioncolor4 = OptionMenu(t, color4, "Blue", "Light Blue", "Green", "Light Green","Red","Orange", "Purple","Yellow")
        optioncolor4.grid(row=7,column=2)
        
        
        
        stagetitleB = Button(t,text="Apply",command=self.stage_graph_aesthetics).grid(row=20,column=1)
        
    def stage_graph_aesthetics(self):
        
        global stagetitle
        global stagexaxis
        global stageyaxis
        
        stagetitle = stagecharttitle.get()
        stagexaxis = stagechartxaxis.get()
        stageyaxis = stagechartyaxis.get()
        
        Scolor1 = color1.get()
        Scolor2 = color2.get()
        Scolor3 = color3.get()
        Scolor4 = color4.get()

        
        global stage1color
        
        
        
        
        if(Scolor1 == "Blue"):
            stage1color = "blue";
        if(Scolor1 == "Light Blue"):
            stage1color = "lightblue"; 
        if(Scolor1 == "Green"):
            stage1color = "green";
        if(Scolor1 == "Light Green"):
            stage1color = "lightgreen";
        if(Scolor1 == "Red"):
            stage1color = "red";
        if(Scolor1 == "Light Red"):
            stage1color = "lightred";
        if(Scolor1 == "Black"):
            stage1color = "black";
        if(Scolor1 == "Yellow"):
            stage1color = "yellow";
        if(Scolor1 == "Orange"):
            stage1color = "orange";
        if(Scolor1 == "Purple"):
            stage1color = "Purple";
        

        global stage2color
        
        if(Scolor2 == "Blue"):
            stage2color = "blue";
        if(Scolor2 == "Light Blue"):
            stage2color = "lightblue"; 
        if(Scolor2 == "Green"):
            stage2color = "green";
        if(Scolor2 == "Light Green"):
            stage2color = "lightgreen";
        if(Scolor2 == "Red"):
            stage2color = "red";
        if(Scolor2 == "Light Red"):
            stage2color = "lightred";
        if(Scolor2 == "Black"):
            stage2color = "black";
        if(Scolor2 == "Yellow"):
            stage2color = "yellow";
        if(Scolor2 == "Light Yellow"):
            stage2color = "light Yellow";
        if(Scolor2 == "Purple"):
            stage2color = "Purple";


        global stage3color
        
        if(Scolor3 == "Blue"):
            stage3color = "blue";
        if(Scolor3 == "Light Blue"):
            stage3color = "lightblue"; 
        if(Scolor3 == "Green"):
            stage3color = "green";
        if(Scolor3 == "Light Green"):
            stage3color = "lightgreen";
        if(Scolor3 == "Red"):
            stage3color = "red";
        if(Scolor3 == "Light Red"):
            stage3color = "lightred";
        if(Scolor3 == "Black"):
            stage3color = "black";
        if(Scolor3 == "Yellow"):
            stage3color = "yellow";
        if(Scolor3 == "Light Yellow"):
            stage3color = "light Yellow";
        if(Scolor3 == "Purple"):
            stage3color = "Purple";

        global stage4color
        
        if(Scolor4 == "Blue"):
            stage4color = "blue";
        if(Scolor4 == "Light Blue"):
            stage4color = "lightblue"; 
        if(Scolor4 == "Green"):
            stage4color = "green";
        if(Scolor4 == "Light Green"):
            stage4color = "lightgreen";
        if(Scolor4 == "Red"):
            stage4color = "red";
        if(Scolor4 == "Light Red"):
            stage4color = "lightred";
        if(Scolor4 == "Black"):
            stage4color = "black";
        if(Scolor4 == "Yellow"):
            stage4color = "yellow";
        if(Scolor4 == "Light Yellow"):
            stage4color = "light Yellow";
        if(Scolor4 == "Purple"):
            stage4color = "Purple";



    def kaplan_graph_change(self):
        
        t = Toplevel(self)
        t.wm_title("Kaplan Meier Graph")
        t.wm_geometry("460x500")
        Label(t, text="Chart Title:").grid(row=1,column=1)

        global kaplancharttitle
        kaplancharttitle = Entry(t)
        kaplancharttitle.delete(0,END)
        kaplancharttitle.insert(0,kaplantitle)
        kaplancharttitle.grid(row=1,column=2)

        global kaplanchartxaxis
        Label(t, text="X Axis:").grid(row=2,column=1)
        kaplanchartxaxis = Entry(t)
        kaplanchartxaxis.delete(0,END)
        kaplanchartxaxis.insert(0,kaplanxaxis)
        kaplanchartxaxis.grid(row=2,column=2)

        global kaplanchartyaxis
        Label(t, text="Y Axis:").grid(row=3,column=1)
        kaplanchartyaxis = Entry(t)
        kaplanchartyaxis.delete(0,END)
        kaplanchartyaxis.insert(0,kaplanyaxis)
        kaplanchartyaxis.grid(row=3,column=2)

        Label(t, text="Color").grid(row=4,column=2)
        Label(t, text="Label").grid(row=4,column=3)
        
        Label(t, text="Line One:").grid(row=5,column=1)
        global kaplancolor1
        kaplancolor1 = StringVar(self)
        kaplancolor1.set(kaplan1color) # initial value

        optionkaplancolor1 = OptionMenu(t, kaplancolor1, "Blue", "Light Blue","Cyan", "Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor1.grid(row=5,column=2)

        Label(t, text="Line Two:").grid(row=6,column=1)
        global kaplancolor2
        kaplancolor2 = StringVar(self)
        kaplancolor2.set(kaplan2color) # initial value

        optionkaplancolor2 = OptionMenu(t, kaplancolor2, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor2.grid(row=6,column=2)

        Label(t, text="Line Three:").grid(row=7,column=1)
        global kaplancolor3
        kaplancolor3 = StringVar(self)
        kaplancolor3.set(kaplan3color) # initial value

        optionkaplancolor3 = OptionMenu(t, kaplancolor3, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor3.grid(row=7,column=2)

        Label(t, text="Line Four:").grid(row=8,column=1)
        global kaplancolor4
        kaplancolor4 = StringVar(self)
        kaplancolor4.set(kaplan4color) # initial value

        optionkaplancolor4 = OptionMenu(t, kaplancolor4, "Blue", "Light Blue","Cyan", "Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor4.grid(row=8,column=2)

        Label(t, text="Line Five:").grid(row=9,column=1)
        global kaplancolor5
        kaplancolor5 = StringVar(self)
        kaplancolor5.set(kaplan5color) # initial value

        optionkaplancolor5 = OptionMenu(t, kaplancolor5, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor5.grid(row=9,column=2)

        Label(t, text="Line Six:").grid(row=10,column=1)
        global kaplancolor6
        kaplancolor6 = StringVar(self)
        kaplancolor6.set(kaplan6color) # initial value

        optionkaplancolor6 = OptionMenu(t, kaplancolor6, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor6.grid(row=10,column=2)

        Label(t, text="Line Seven:").grid(row=11,column=1)
        global kaplancolor7
        kaplancolor7 = StringVar(self)
        kaplancolor7.set(kaplan7color) # initial value

        optionkaplancolor7 = OptionMenu(t, kaplancolor7, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor7.grid(row=11,column=2)

        Label(t, text="Line Eight:").grid(row=12,column=1)
        global kaplancolor8
        kaplancolor8 = StringVar(self)
        kaplancolor8.set(kaplan8color) # initial value

        optionkaplancolor8 = OptionMenu(t, kaplancolor8, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor8.grid(row=12,column=2)
        Label(t, text="Line Nine:").grid(row=13,column=1)
        global kaplancolor9
        kaplancolor9 = StringVar(self)
        kaplancolor9.set(kaplan9color) # initial value

        optionkaplancolor9 = OptionMenu(t, kaplancolor9, "Blue", "Light Blue","Cyan", "Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor9.grid(row=13,column=2)

        Label(t, text="Line Ten:").grid(row=14,column=1)
        global kaplancolor10
        kaplancolor10 = StringVar(self)
        kaplancolor10.set(kaplan10color) # initial value

        optionkaplancolor10 = OptionMenu(t, kaplancolor10, "Blue", "Light Blue","Cyan" ,"Green", "Light Green","Red","Orange", "Purple","Magenta","Yellow","Black")
        optionkaplancolor10.grid(row=14,column=2)

        
        

        
        global kaplabelE1
        
        kaplabelE1 = Entry(t)
        kaplabelE1.delete(0,END)
        kaplabelE1.insert(0,kaplabel1)
        kaplabelE1.grid(row=5,column=3)

        global kaplabelE2
        
        kaplabelE2 = Entry(t)
        kaplabelE2.delete(0,END)
        kaplabelE2.insert(0,kaplabel2)
        kaplabelE2.grid(row=6,column=3)

        global kaplabelE3
        
        kaplabelE3 = Entry(t)
        kaplabelE3.delete(0,END)
        kaplabelE3.insert(0,kaplabel3)
        kaplabelE3.grid(row=7,column=3)

        global kaplabelE4
        
        kaplabelE4 = Entry(t)
        kaplabelE4.delete(0,END)
        kaplabelE4.insert(0,kaplabel4)
        kaplabelE4.grid(row=8,column=3)

        global kaplabelE5
        
        kaplabelE5 = Entry(t)
        kaplabelE5.delete(0,END)
        kaplabelE5.insert(0,kaplabel5)
        kaplabelE5.grid(row=9,column=3)

        global kaplabelE6
        
        kaplabelE6 = Entry(t)
        kaplabelE6.delete(0,END)
        kaplabelE6.insert(0,kaplabel6)
        kaplabelE6.grid(row=10,column=3)

        global kaplabelE7
        
        kaplabelE7 = Entry(t)
        kaplabelE7.delete(0,END)
        kaplabelE7.insert(0,kaplabel7)
        kaplabelE7.grid(row=11,column=3)

        global kaplabelE8
        
        kaplabelE8 = Entry(t)
        kaplabelE8.delete(0,END)
        kaplabelE8.insert(0,kaplabel8)
        kaplabelE8.grid(row=12,column=3)

        global kaplabelE9
        
        kaplabelE9 = Entry(t)
        kaplabelE9.delete(0,END)
        kaplabelE9.insert(0,kaplabel9)
        kaplabelE9.grid(row=13,column=3)

        global kaplabelE10
        
        kaplabelE10 = Entry(t)
        kaplabelE10.delete(0,END)
        kaplabelE10.insert(0,kaplabel10)
        kaplabelE10.grid(row=14,column=3)

        with open("minpatients.csv",'r') as file:
            PATS = file.read().split(',');
            file.close()

        minpat = int(PATS[0]);


        global kapminpat
        Label(t, text="Kaplan Meier Min Patients For Auto:").grid(row=15,column=1)
        kapminpat = Entry(t)
        kapminpat.delete(0,END)
        kapminpat.insert(0,minpat)
        kapminpat.grid(row=15,column=3)
        
        stagetitleB = Button(t,text="Apply",command=self.kaplan_graph_aesthetics).grid(row=16,column=1)
        
    def kaplan_graph_aesthetics(self):
        
        global kaplantitle
        global kaplanxaxis
        global kaplanyaxis
        
        kaplantitle = kaplancharttitle.get()
        kaplanxaxis = kaplanchartxaxis.get()
        kaplanyaxis = kaplanchartyaxis.get()

        global kaplabel1
        global kaplabel2
        global kaplabel3
        global kaplabel4
        global kaplabel5
        global kaplabel6
        global kaplabel7
        global kaplabel8
        global kaplabel9
        global kaplabel10
        
        
        kaplabel1 = kaplabelE1.get()
        kaplabel2 = kaplabelE1.get()
        kaplabel3 = kaplabelE1.get()
        kaplabel4 = kaplabelE1.get()
        kaplabel5 = kaplabelE1.get()
        kaplabel6 = kaplabelE1.get()
        kaplabel7 = kaplabelE1.get()
        kaplabel8 = kaplabelE1.get()
        kaplabel9 = kaplabelE1.get()
        kaplabel10 = kaplabelE1.get()
        
        
        Kcolor1 = kaplancolor1.get()
        Kcolor2 = kaplancolor2.get()
        Kcolor3 = kaplancolor3.get()
        Kcolor4 = kaplancolor4.get()
        Kcolor5 = kaplancolor5.get()
        Kcolor6 = kaplancolor6.get()
        Kcolor7 = kaplancolor7.get()
        Kcolor8 = kaplancolor8.get()
        Kcolor9 = kaplancolor9.get()
        Kcolor10 = kaplancolor10.get()

        global kaplan1color
        global kaplan2color
        global kaplan3color
        global kaplan4color
        global kaplan5color
        global kaplan6color
        global kaplan7color
        global kaplan8color
        global kaplan9color
        global kaplan10color
        
        kaplan1color = Kcolor1;
        kaplan2color = Kcolor2;
        kaplan3color = Kcolor3;
        kaplan4color = Kcolor4;
        kaplan5color = Kcolor5;
        kaplan6color = Kcolor6;
        kaplan7color = Kcolor7;
        kaplan8color = Kcolor8;
        kaplan9color = Kcolor9;
        kaplan10color = Kcolor10;

        
        
        
        
        
        if(Kcolor1 == "Light Blue"):
            kaplan1color = "lightblue"; 
        
        if(Kcolor1 == "Light Green"):
            kaplan1color = "lightgreen";

        if(Kcolor2 == "Light Blue"):
            kaplan2color = "lightblue"; 
        
        if(Kcolor2 == "Light Green"):
            kaplan2color = "lightgreen";


        if(Kcolor3 == "Light Blue"):
            kaplan3color = "lightblue"; 
        
        if(Kcolor3 == "Light Green"):
            kaplan3color = "lightgreen";
        

        if(Kcolor4 == "Light Blue"):
            kaplan4color = "lightblue"; 
        
        if(Kcolor4 == "Light Green"):
            kaplan4color = "lightgreen";
       
       

        if(Kcolor5 == "Light Blue"):
            kaplan5color = "lightblue"; 
       
        if(Kcolor5 == "Light Green"):
            kaplan5color = "lightgreen";
        
        if(Kcolor6 == "Light Blue"):
            kaplan6color = "lightblue"; 
        
        if(Kcolor6 == "Light Green"):
            kaplan6color = "lightgreen";
        

        if(Kcolor7 == "Light Blue"):
            kaplan7color = "lightblue"; 
       
        if(Kcolor7 == "Light Green"):
            kaplan7color = "lightgreen";
       

        if(Kcolor8 == "Light Blue"):
            kaplan8color = "lightblue"; 
       
        if(Kcolor8 == "Light Green"):
            kaplan8color = "lightgreen";
        
        if(Kcolor9 == "Light Blue"):
            kaplan9color = "lightblue"; 
        
        if(Kcolor9 == "Light Green"):
            kaplan9color = "lightgreen";
        
        if(Kcolor10 == "Light Blue"):
            kaplan10color = "lightblue"; 
        
        if(Kcolor10 == "Light Green"):
            kaplan10color = "lightgreen";
        with open("minpatients.csv",'w') as file:
            file.write(kapminpat.get())
            file.close()

        
        
            
    def own_data(self):
        print("own")



        
    def normal_data(self):

        #raw data stdev and mean
        ave = mean(rawData)
            
        self.average.config(state=NORMAL)
        self.average.delete(0,END)
        self.average.insert(0,ave)
        self.average.config(state=DISABLED)

        stan = np.std([rawData])
            
            
        self.standerddev.config(state=NORMAL)
        self.standerddev.delete(0,END)
        self.standerddev.insert(0,stan)
        self.standerddev.config(state=DISABLED)

        medi = median(rawData)
        
        self.median.config(state=NORMAL)
        self.median.delete(0,END)
        self.median.insert(0,medi)
        self.median.config(state=DISABLED)

        ran = max(rawData)-min(rawData)
        
        self.range.config(state=NORMAL)
        self.range.delete(0,END)
        self.range.insert(0,ran)
        self.range.config(state=DISABLED)

        global quart1
        global quart3
        
        quart1 = [];
        quart3 = [];

        xx=len(rawData)
        x=0;
        bre1 = 0;
        bre3 =0;
        while(xx>0):
            if(rawData[x]<medi):
                quart1.insert(bre1,rawData[x])
                bre1+=1

            if(rawData[x]>medi):
                quart3.insert(bre3,rawData[x])
                bre3+=1
            
            x+=1
            xx-=1
        
        

        medq1 = median(quart1)
        
        self.q1.config(state=NORMAL)
        self.q1.delete(0,END)
        self.q1.insert(0,medq1)
        self.q1.config(state=DISABLED)

        medq3 = median(quart3)
        
        self.q3.config(state=NORMAL)
        self.q3.delete(0,END)
        self.q3.insert(0,medq3)
        self.q3.config(state=DISABLED)


        IQR =  (medq3 - medq1);

        self.iqr.config(state=NORMAL)
        self.iqr.delete(0,END)
        self.iqr.insert(0,IQR)
        self.iqr.config(state=DISABLED)



        self.minimum.config(state=NORMAL)
        self.minimum.delete(0,END)
        self.minimum.insert(0,min(rawData))
        self.minimum.config(state=DISABLED)

        self.maximum.config(state=NORMAL)
        self.maximum.delete(0,END)
        self.maximum.insert(0,max(rawData))
        self.maximum.config(state=DISABLED)


        

        
    def gene_search(self):
        self.range.config(state=NORMAL)
        self.standerddev.config(state=NORMAL)
        self.median.config(state=NORMAL)
        self.average.config(state=NORMAL)
        self.tell.config(state=NORMAL)
        self.q1.config(state=NORMAL)
        self.q3.config(state=NORMAL)
        self.iqr.config(state=NORMAL)
        self.minimum.config(state=NORMAL)
        self.maximum.config(state=NORMAL)
        
        self.standerddev.delete(0,END)
        self.average.delete(0,END)
        self.median.delete(0,END)
        self.range.delete(0,END)
        self.tell.delete(0,END)
        self.q3.delete(0,END)
        self.q1.delete(0,END)
        self.iqr.delete(0,END)
        self.minimum.delete(0,END)
        self.maximum.delete(0,END)

        self.range.config(state=DISABLED)
        self.standerddev.config(state=DISABLED)
        self.median.config(state=DISABLED)
        self.average.config(state=DISABLED)
        self.tell.config(state=DISABLED)
        self.q1.config(state=DISABLED)
        self.q3.config(state=DISABLED)
        self.iqr.config(state=DISABLED)
        self.minimum.config(state=DISABLED)
        self.maximum.config(state=DISABLED)
        
        global entry
        entry = self.geneSearch.get()
           
        global rawData
        global dataS
        dataS = var.get()
        try:
            if(dataS == "Data Source:"):
                self.tell.config(state=NORMAL)
                self.tell.insert(0,"You must pick a data source")
                self.tell.config(state=DISABLED)

            if(dataS == "Own Data"):
                urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=lusc_tcga&cancer_study_id=lusc_tcga&genetic_profile_ids=lusc_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=lusc_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                lengths = 177
                urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

                return;
                    

            if(dataS == "Lung Squamous Cell Carcinoma (TCGA, Provisional)"):
                #main url
                urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=lusc_tcga&cancer_study_id=lusc_tcga&genetic_profile_ids=lusc_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=lusc_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                lengths = 177
                urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"
                
            if(dataS == "Lung Adenocarcinoma (TCGA, Provisional)"):
                #main url
                urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=luad_tcga&cancer_study_id=luad_tcga&genetic_profile_ids=luad_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=luad_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                lengths = 230
                urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

                           
                

            

            with open(os.path.join("Data Sources\Data "+dataS+"\Patient Data\Stage Data.csv"),"r") as stageData:
                
                StageDataT = stageData.read().split('\n')
                stageData.close()
            
            global StageData
            StageData = [];
            del StageData[:];

            
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

            with open(os.path.join("Data Sources\Data "+dataS+"\Kaplan Meiers\Paitent Data\Alive and Dead.csv"),"r") as alive:
                
                ADT = alive.read().split('\n')
                alive.close()
            global AliveDead
            AliveDead = [];
            del AliveDead[:];



            
            urlmid = self.geneSearch.get()
            url =urlstart+urlmid+urlend
            
            #main numbers
            resp = urllib.request.urlopen(url)
            respData = resp.read()
            
            if(str(respData).find("# Warning:  Unknown gene:")==2):

                self.tell.config(state=NORMAL)
                self.tell.insert(0,"# Warning:  Unknown gene:")
                self.tell.config(state=DISABLED)
            
            if(str(respData).find("# Warning:  Unknown gene:")!=2):
              
                rawDataT = str(respData).replace(entry,"1")    
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
                    


                  
                    outliers = out.get()
                    rawData = [];
                    del rawData[:];
                    
                    
                    if(outliers == "Include Outliers"):
                        rawData = rawDataT;
                        
                        xx = lengths;
                        x = 0;
                        while(xx>0):
                            StageData.insert(x,float(StageDataT[x]));
                            Survival.insert(x,int(SurvivalT[x]));
                            Censord.insert(x,int(CensordT[x]));
                            AliveDead.insert(x,int(ADT[x]));

                            xx-=1;
                            x+=1;
                        
                        self.normal_data()

                    if(outliers == "Exclude Outliers"):
                        
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
                                StageData.insert(breaker,float(StageDataT[x]))
                                Survival.insert(breaker,int(SurvivalT[x]))
                                Censord.insert(breaker,int(CensordT[x]))
                                AliveDead.insert(breaker,int(ADT[x]))
                                breaker+=1
                                

                            
                            
                            x+=1
                            xx-=1
                        
                        self.normal_data()
                            

                        
                        
                    
                except:
                    pass;
        except:
            pass;
        

        
            
    def stage_data(self):
        #stages
        

        

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
        idsone = [x for x in range(len(StageOne))]
        idstwo = [x for x in range(len(StageTwo))]
        
        idsthree = [x for x in range(len(StageThree))]
        idsfour = [x for x in range(len(StageFour))]

        P.scatter(idsone,sorted(StageOne),color=stage1color,label="Stage One")
        P.scatter(idstwo,sorted(StageTwo),color=stage2color,label="Stage Two")
        P.scatter(idsthree,sorted(StageThree),color=stage3color,label="Stage Three")
        P.scatter(idsfour,sorted(StageFour),color=stage4color,label="Stage Four")
        
        
        if(stagetitle == "DEFAULT"):
            titl =  " Stages of "+ dataS +" With the Expression of " + self.geneSearch.get() +"\n  "
            P.title(titl)
            
        if(stagetitle != "DEFAULT"):
            P.title(stagetitle)
            

        if(stagexaxis == "DEFAULT"):
             P.xlabel('Number of Patients Per Stage')
        if(stagexaxis != "DEFAULT"):
            P.xlabel(stagexaxis)

        if(stageyaxis == "DEFAULT"):
            P.ylabel('Expression of '+ self.geneSearch.get()+' in RNA Seq V2')
        if(stageyaxis != "DEFAULT"):
            P.ylabel(stageyaxis)
        midlinex = [0,len(StageOne)]
        med1 = [median(StageOne),median(StageOne)];
        med2 = [median(StageTwo),median(StageTwo)];
        med3 = [median(StageThree),median(StageThree)];
        med4 = [median(StageFour),median(StageFour)];
        
        P.plot(midlinex,med1,color=stage1color,linestyle="--")
        P.plot(midlinex,med2,color=stage2color,linestyle="--")
        P.plot(midlinex,med3,color=stage3color,linestyle="--")
        P.plot(midlinex,med4,color=stage4color,linestyle="--")
        
        P.grid(True)
        P.legend()
        
        P.show()                
                            



    def raw_data(self):
        #raw data graph
        dataS = var.get()
        P.figure()
        rawtitle = dataS+" Data for "+ self.geneSearch.get() +" Expression\n  "
        expression = rawData
        ids = [x for x in range(len(rawData))]
        P.bar(ids,sorted(rawData))
        P.ylabel('Expression of '+ self.geneSearch.get()+'in RNA Seq V2')
        P.xlabel('Patient Identification Number')
        P.title(rawtitle)
        P.grid(True)
        P.show()

    #retrive raw data
    def raw_data_ret(self):
        raw_ = str(rawData).replace("[",'')
        raw__=raw_.replace(']','')

        sur_ = str(Survival).replace("[",'')
        sur__=sur_.replace(']','')

        cen_ = str(Censord).replace("[",'')
        cen__=cen_.replace(']','')

        sta_ = str(StageData).replace("[",'')
        sta__=sta_.replace(']','')
        
        with open(os.path.join("Data For "+self.geneSearch.get()+".csv"),"w") as saveData:
            saveData.write('Expression')
            saveData.write(',')
            saveData.write(raw__)
            saveData.write('\n')
            saveData.write('Survival')
            saveData.write(',')
            saveData.write(sur__)
            saveData.write('\n')
            saveData.write('Censored')
            saveData.write(',')
            saveData.write(cen__)
            saveData.write('\n')
            saveData.write('Stage')
            saveData.write(',')
            saveData.write(sta__)
            saveData.close()
                      
                      

    #correlate to one gene
    def r_squared(self):

                          
                          
        def best_fit_slope_and_intercept(xs,ys):
            m = (((mean(xs)*mean(ys)) - mean(xs*ys)) /
                 ((mean(xs)*mean(xs)) - mean(xs*xs)))
            b = mean(ys) - m*mean(xs)
            return m, b

        def squared_error(ys_orig,ys_line):
            return sum((ys_line - ys_orig) * (ys_line - ys_orig))

        def coefficient_of_determination(ys_orig,ys_line):
            y_mean_line = [mean(ys_orig) for y in ys_orig]
            squared_error_regr = squared_error(ys_orig, ys_line)
            squared_error_y_mean = squared_error(ys_orig, y_mean_line)
            return 1 - (squared_error_regr/squared_error_y_mean)
        try:
                

            entryCorr = self.correlateGene.get()
            
            if(self.correlateGene.get() != ''):
            
                dataS = var.get()
                try:
                    if(dataS == "Data Source:"):
                        self.tell.config(state=NORMAL)
                        self.tell.insert(0,"You must pick a data source")
                        self.tell.config(state=DISABLED)
                            

                    if(dataS == "Lung Squamous Cell Carcinoma (TCGA, Provisional)"):
                        #main url
                        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=lusc_tcga&cancer_study_id=lusc_tcga&genetic_profile_ids=lusc_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=lusc_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                        lengths = 177
                        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"
                        
                    if(dataS == "Lung Adenocarcinoma (TCGA, Provisional)"):
                        #main url
                        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=luad_tcga&cancer_study_id=luad_tcga&genetic_profile_ids=luad_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=luad_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                        lengths = 230
                        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

                        

                    

                   


                    
                    urlmid = self.geneSearch.get()
                    url =urlstart+urlmid+urlend

                    #main numbers
                    resp = urllib.request.urlopen(url)
                    respData = resp.read()
                    if(str(respData).find("# Warning:  Unknown gene:")==2):

                        self.tell.config(state=NORMAL)
                        self.tell.insert(0,"# Warning:  Unknown gene:")
                        self.tell.config(state=DISABLED)
                    
                    if(str(respData).find("# Warning:  Unknown gene:")!=2):
                      
                        rawDataT = str(respData).replace(self.geneSearch.get(),"1")    
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
                            
                            


                          
                        except:
                            pass;
                except:
                    pass;
                

                
                try:
                    if(dataS == "Data Source:"):
                        self.tell.config(state=NORMAL)
                        self.tell.insert(0,"You must pick a data source")
                        self.tell.config(state=DISABLED)
                            

                    if(dataS == "Lung Squamous Cell Carcinoma (TCGA, Provisional)"):
                        #main url
                        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=lusc_tcga&cancer_study_id=lusc_tcga&genetic_profile_ids=lusc_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=lusc_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                        lengths = 177
                        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"
                        
                    if(dataS == "Lung Adenocarcinoma (TCGA, Provisional)"):
                        #main url
                        urlstart = "http://www.cbioportal.org/index.do?cancer_study_list=luad_tcga&cancer_study_id=luad_tcga&genetic_profile_ids=luad_tcga_rna_seq_v2_mrna&data_priority=0&case_set_id=luad_tcga_3way_complete&case_ids=&patient_case_select=sample&gene_set_choice=user-defined-list&gene_list="
                        lengths = 230
                        urlend = "&clinical_param_selection=null&tab_index=tab_download&transpose_matrix=on&Action=Submit"

                        

                    

                    



                    
                    urlmid = self.correlateGene.get()
                    url =urlstart+urlmid+urlend

                    #main numbers
                    respC = urllib.request.urlopen(url)
                    respDataC = respC.read()


                    if(str(respDataC).find("# Warning:  Unknown gene:")==2):

                        self.tell.config(state=NORMAL)
                        self.tell.insert(0,"# Warning:  Unknown gene:")
                        self.tell.config(state=DISABLED)
                    
                    if(str(respData).find("# Warning:  Unknown gene:")!=2):
                      
                        rawDataTC = str(respDataC).replace(self.correlateGene.get(),"1")    
                        newDataC = re.findall(r'[-+]?\d*\.\d+|\d+',str(rawDataTC))
                        
                        

                        try:
                            #Varibles
                            z = 0;
                            x = 5;
                            y = 0;
                            rawDataTC = [];
                            del rawDataTC[:];
                            #insert data to rawData
                            while(y<lengths):
                                    
                                rawDataTC.insert(y,float(newDataC[x]))

                                x+=4
                                y+=1
                                
                                    
                            if(len(rawDataTC) != lengths):
                                print("Error: Something Went Wrong With The Datas List It Isnt the Right Length")
                                sys.exit()
                            


                          
                            outliers = out.get()
                            rawDataC = [];
                            del rawDataC[:];
                            if(outliers == "Include Outliers"):

                                
                                rawDataC = rawDataTC;
                                
                                
                               

                            if(outliers == "Exclude Outliers"):
                                
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
                                        rawDataC.insert(breaker,rawDataTC[x])
                                        
                                        breaker+=1
                                        

                                    
                                    
                                    x+=1
                                    xx-=1
                                
                                
                                    

                                
                                
                            
                        except:
                            pass;
                except:
                    pass;
                

                
       
               
                
                ys = np.array(rawDataC)
                xs = np.array(rawData)

               
                            
                m, b = best_fit_slope_and_intercept(xs,ys)
                regression_line = [(m*x)+b for x in xs]

                r_squared = coefficient_of_determination(ys,regression_line)
                
                self.r_squaredE.config(state=NORMAL)
                self.r_squaredE.delete(0,END)
                self.r_squaredE.insert(0,str(r_squared))
                self.r_squaredE.config(state=DISABLED)
                P.figure()
                rawtitle = dataS+" Correlation of "+ self.geneSearch.get() +" and "+self.correlateGene.get()

                idsforline = [];
                rawMax=max(rawData)+100
                  
                rawMin=min(rawData)
                   
                insert=0;
                    
                while(rawMax>rawMin):
                    idsforline.insert(insert,insert)                                                                                                 
                        
                    insert+=1
                    rawMax-=1
                    

                regline = [];
                xxx=0;
                xx = len(idsforline)
                while(xx>0):
                    po = m*idsforline[xxx]
                    regline.insert(xxx,po)
                    
                    xxx+=1
                    xx-=1



                MR = mean(rawData)
                MC = mean(rawDataC)
                ttop=MR - MC;
                n1 =np.std(rawData)/len(rawData)
                n2 = np.std(rawDataC)/len(rawDataC)
                n3 = n1+n2                        
                tbottom=np.sqrt(n3)
                ttest = ttop/tbottom

                self.tTest.config(state=NORMAL)
                self.tTest.delete(0,END)
                self.tTest.insert(0,ttest)
                self.tTest.config(state=DISABLED)

                
                
                P.scatter(rawData,rawDataC)
                P.plot(idsforline,regline,linestyle="--",color="grey")
                P.ylabel('Expression of '+ self.correlateGene.get()+' in RNA Seq V2')
                P.xlabel('Expression of '+ self.geneSearch.get()+' in RNA Seq V2')
                P.title(rawtitle)
                P.grid(True)
                P.show()

        except:
            pass;
        if(self.correlateGene.get() == ''):
            self.tell.config(state=NORMAL)
            self.tell.delete(0,END)
            self.tell.insert(0,"You Must Pick A Gene To Correlate")
            self.tell.config(state=DISABLED)
    def r_squared_all(self):
        with open("entry.csv",'w') as file:
            file.write(entry);
            file.write('\n');
            file.write(dataS);
            file.write('\n');
            file.write(out.get());
            file.close();
            
        os.system("correlateall.py 1")
                    
            





        
    def kaplan_data(self):
        amtoflines = kap.get()
        if(amtoflines == "One Line"):
            kplots = 1;
            
            if(float(kaphighone.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan One is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighone.get())<float(kaplowone.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan One is Too Low")
                self.tell.config(state=DISABLED)
               
            
            
        if(amtoflines == "Two Lines"):
            kplots = 2;
            if(float(kaphightwo.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Two is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphightwo.get())<float(kaplowtwo.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Two is Too Low")
                self.tell.config(state=DISABLED)
                
        if(amtoflines == "Three Lines"):
            kplots = 3;

            if(float(kaphighthree.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Three is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighthree.get())<float(kaplowthree.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Three is Too Low")
                self.tell.config(state=DISABLED)
                
        if(amtoflines == "Four Lines"):
            kplots = 4;

            if(float(kaphighfour.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Four is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighfour.get())<float(kaplowfour.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Four is Too Low")
                self.tell.config(state=DISABLED)
                
        if(amtoflines == "Five Lines"):
            kplots = 5;

            if(float(kaphighfive.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Five is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighfive.get())<float(kaplowfive.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Five is Too Low")
                self.tell.config(state=DISABLED)

                
        if(amtoflines == "Six Lines"):
            kplots = 6;

            if(float(kaphighsix.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Six is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighsix.get())<float(kaplowsix.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Six is Too Low")
                self.tell.config(state=DISABLED)

                
        if(amtoflines == "Seven Lines"):
            kplots = 7;

            if(float(kaphighseven.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Seven is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighseven.get())<float(kaplowseven.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Seven is Too Low")
                self.tell.config(state=DISABLED)

                
        if(amtoflines == "Eight Lines"):
            kplots = 8;

            if(float(kaphigheight.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Eight is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphigheight.get())<float(kaploweight.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Eight is Too Low")
                self.tell.config(state=DISABLED)

                
        if(amtoflines == "Nine Lines"):
            kplots = 9;

            if(float(kaphighnine.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Nine is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighnine.get())<float(kaplownine.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Nine is Too Low")
                self.tell.config(state=DISABLED)

                
        if(amtoflines == "Ten Lines"):
            kplots = 10;

            if(float(kaphighten.get())<min(rawData)):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Ten is Too Low")
                self.tell.config(state=DISABLED)
            if(float(kaphighten.get())<float(kaplowten.get())):
                self.tell.config(state=NORMAL)
                self.tell.delete(0,END)
                self.tell.insert(0,"Your High Value For Kaplan Ten is Too Low")
                self.tell.config(state=DISABLED)
       
        global kaplandata1
        global kaplandata2
        global kaplandata3
        global kaplandata4
        global kaplandata5
        global kaplandata6
        global kaplandata7
        global kaplandata8
        global kaplandata9
        global kaplandata10
      
        P.figure()
        
        l1 = float(kaplowone.get())
        h1 = float(kaphighone.get())

        kaplandata1 = [];
        kaplancen1 = [];
        

        
        xxx = int(max(Survival));
        breaker1=0;
        
        
        x = 0;
        xx = len(rawData)
        while(xx>0):
            if(h1 >= rawData[x] >= l1):
                kaplandata1.insert(breaker1,int(Survival[x]))
                kaplancen1.insert(breaker1,int(AliveDead[x]))
                breaker1+=1;
                
            x+=1;
            xx-=1;


        x=0;
        xxx-=1;

        xx= max(kaplandata1)
        tot  = len(kaplandata1)
        division = len(kaplandata1)

        kap1per = [];

        
        kaplandata1, kaplancen1 = (list(t) for t in zip(*sorted(zip(kaplandata1, kaplancen1))))

        idsdays = [];

        print(kaplandata1)
        print(kaplancen1)

        

        breaker = 0;

        deaths = 0;
        censored = 0;

        censoreds = [];
        censoredids = [];

        breakercen = 0;
        
        while(xx>0):
            thisroundcensored = 0;
            
            howmany = 0;
            
            if(kaplandata1.count(x) >= 1):
                

                howmany = kaplandata1.count(x)-1;
                

                while(howmany >= 0):

                    
                    

                    if(kaplancen1[breaker] == 0):
                        deaths+=1;
                        
                    

                    if(kaplancen1[breaker] == 1):
                        thisroundcensored+=1;
                        censored+=1;
                        

                    howmany-=1;
                    
                    breaker+=1;
                    
            

            updiv = division-deaths-censored;
            
            downdiv = tot-censored;

            percent = updiv/downdiv;
            

            kap1per.insert(x,percent);
            idsdays.insert(x,x)

            if(thisroundcensored > 0):
                censoreds.insert(breakercen,kap1per[x]);
                censoredids.insert(breakercen,idsdays[x]);
                breakercen+=1;

            
                
                

            x+=1;   
            xx-=1;

            
        if(kaplabel1 == "DEFAULT"):
            kaplab1 = entry+' Expression from ' + str(l1) +' to ' +str(h1)

        if(kaplabel1 != "DEFAULT"):
            kaplab1 = kaplabel1;
            
        P.plot(idsdays,kap1per,linestyle="-",color=kaplan1color,label= kaplab1)
        P.scatter(censoredids,censoreds,marker ='+',color='black')


        
        
        if(kplots >= 2):
            l2 = float(kaplowtwo.get())
            h2 = float(kaphightwo.get())

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
            

                

            longestsurvival=max(kaplandata2)+1;

            x=1;
            breaker1=1;
            days = 0;
            riskset1 = len(kaplandata2)


            kap2 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata2) if x == days]
                riskset1-=len(indices)
                kap2.insert(days,riskset1)
                    
                days+=1;
                longestsurvival-=1;

            

            kap2per = [];
            x = 0;
            xx = len(kap2)
            tot = len(kaplandata2)

            while(xx>0):
                kap2per.insert(x,kap2[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays2=[];
            daysMax = int(len(kap2per))
            
            x=0;
            while(daysMax>0):
                idsdays2.insert(x,x)
                x+=1;
                daysMax-=1;
            if(kaplabel2 == "DEFAULT"):
                kaplab2 = entry+' Expression from ' + str(l2) +' to ' +str(h2)

            if(kaplabel2 != "DEFAULT"):
                kaplab12 = kaplabel2;
                
            P.plot(idsdays2,kap2per,linestyle="-",color=kaplan2color,label= kaplab2)


        if(kplots >= 3):
            l3 = float(kaplowthree.get())
            h3 = float(kaphighthree.get())

            kaplandata3 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h3 >= rawData[x] >= l3):
                    kaplandata3.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata3)+1;

            x=1;
            breaker3=1;
            days = 0;
            riskset3 = len(kaplandata3)


            kap3 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata3) if x == days]
                riskset3-=len(indices)
                kap3.insert(days,riskset3)
                    
                days+=1;
                longestsurvival-=1;

            

            kap3per = [];
            x = 0;
            xx = len(kap3)
            tot = len(kaplandata3)

            while(xx>0):
                kap3per.insert(x,kap3[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays3=[];
            daysMax3 = int(len(kap3per))
            
            x=0;
            while(daysMax3>0):
                idsdays3.insert(x,x)
                x+=1;
                daysMax3-=1;

            if(kaplabel3 == "DEFAULT"):
                kaplab3 = entry+' Expression from ' + str(l3) +' to ' +str(h3)

            if(kaplabel3 != "DEFAULT"):
                kaplab13 = kaplabel3;
           
            P.plot(idsdays3,kap3per,linestyle="-",color=kaplan3color,label=kaplab3)

                  

        if(kplots >= 4):
            l4 = float(kaplowfour.get())
            h4 = float(kaphighfour.get())

            kaplandata4 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h4 >= rawData[x] >= l4):
                    kaplandata4.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata4)+1;

            x=1;
            breaker4=1;
            days = 0;
            riskset4 = len(kaplandata4)


            kap4 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata4) if x == days]
                riskset4-=len(indices)
                kap4.insert(days,riskset4)
                    
                days+=1;
                longestsurvival-=1;

            

            kap4per = [];
            x = 0;
            xx = len(kap4)
            tot = len(kaplandata4)

            while(xx>0):
                kap4per.insert(x,kap4[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays4=[];
            daysMax4 = int(len(kap4per))
            
            x=0;
            while(daysMax4>0):
                idsdays4.insert(x,x)
                x+=1;
                daysMax4-=1;

            if(kaplabel4 == "DEFAULT"):
                kaplab4 = entry+' Expression from ' + str(l4) +' to ' +str(h4)

            if(kaplabel4 != "DEFAULT"):
                kaplab14 = kaplabel4;
           
            P.plot(idsdays4,kap4per,linestyle="-",color=kaplan4color,label=kaplab4)


        if(kplots >= 5):
            l5 = float(kaplowfive.get())
            h5 = float(kaphighfive.get())

            kaplandata5 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h5 >= rawData[x] >= l5):
                    kaplandata5.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata5)+1;

            x=1;
            breaker5=1;
            days = 0;
            riskset5 = len(kaplandata5)


            kap5 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata5) if x == days]
                riskset5-=len(indices)
                kap5.insert(days,riskset5)
                    
                days+=1;
                longestsurvival-=1;

            

            kap5per = [];
            x = 0;
            xx = len(kap5)
            tot = len(kaplandata5)

            while(xx>0):
                kap5per.insert(x,kap5[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays5=[];
            daysMax5 = int(len(kap5per))
            
            x=0;
            while(daysMax5>0):
                idsdays5.insert(x,x)
                x+=1;
                daysMax5-=1;

            if(kaplabel5 == "DEFAULT"):
                kaplab5 = entry+' Expression from ' + str(l5) +' to ' +str(h5)

            if(kaplabel5 != "DEFAULT"):
                kaplab15 = kaplabel5;
           
            P.plot(idsdays5,kap5per,linestyle="-",color=kaplan5color,label=kaplab5)


            
        if(kplots >= 6):
            l6 = float(kaplowsix.get())
            h6 = float(kaphighsix.get())

            kaplandata6 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h6 >= rawData[x] >= l6):
                    kaplandata6.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata6)+1;

            x=1;
            breaker6=1;
            days = 0;
            riskset6= len(kaplandata6)


            kap6 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata6) if x == days]
                riskset6-=len(indices)
                kap6.insert(days,riskset6)
                    
                days+=1;
                longestsurvival-=1;

            

            kap6per = [];
            x = 0;
            xx = len(kap6)
            tot = len(kaplandata6)

            while(xx>0):
                kap6per.insert(x,kap6[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays6=[];
            daysMax6 = int(len(kap6per))
            
            x=0;
            while(daysMax6>0):
                idsdays6.insert(x,x)
                x+=1;
                daysMax6-=1;

            if(kaplabel6 == "DEFAULT"):
                kaplab6 = entry+' Expression from ' + str(l6) +' to ' +str(h6)

            if(kaplabel6 != "DEFAULT"):
                kaplab16 = kaplabel6;
           
            P.plot(idsdays6,kap6per,linestyle="-",color=kaplan6color,label=kaplab6)

            
        if(kplots >= 7):
            l7 = float(kaplowseven.get())
            h7 = float(kaphighseven.get())

            kaplandata7 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h7 >= rawData[x] >= l7):
                    kaplandata7.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata7)+1;

            x=1;
            breaker7=1;
            days = 0;
            riskset7 = len(kaplandata7)


            kap7 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata7) if x == days]
                riskset7-=len(indices)
                kap7.insert(days,riskset7)
                    
                days+=1;
                longestsurvival-=1;

            

            kap7per = [];
            x = 0;
            xx = len(kap7)
            tot = len(kaplandata7)

            while(xx>0):
                kap7per.insert(x,kap7[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays7=[];
            daysMax7 = int(len(kap7per))
            
            x=0;
            while(daysMax7>0):
                idsdays7.insert(x,x)
                x+=1;
                daysMax7-=1;
            if(kaplabel7 == "DEFAULT"):
                kaplab7 = entry+' Expression from ' + str(l7) +' to ' +str(h7)

            if(kaplabel7 != "DEFAULT"):
                kaplab17 = kaplabel7;
           
            P.plot(idsdays7,kap7per,linestyle="-",color=kaplan7color,label=kaplab7)

        if(kplots >= 8):
            l8 = float(kaploweight.get())
            h8 = float(kaphigheight.get())

            kaplandata8 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h8 >= rawData[x] >= l8):
                    kaplandata8.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata8)+1;

            x=1;
            breaker8=1;
            days = 0;
            riskset8 = len(kaplandata8)


            kap8 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata8) if x == days]
                riskset8-=len(indices)
                kap8.insert(days,riskset8)
                    
                days+=1;
                longestsurvival-=1;

            

            kap8per = [];
            x = 0;
            xx = len(kap8)
            tot = len(kaplandata8)

            while(xx>0):
                kap8per.insert(x,kap8[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays8=[];
            daysMax8 = int(len(kap8per))
            
            x=0;
            while(daysMax8>0):
                idsdays8.insert(x,x)
                x+=1;
                daysMax8-=1;

            if(kaplabel8 == "DEFAULT"):
                kaplab8 = entry+' Expression from ' + str(l8) +' to ' +str(h8)

            if(kaplabel8 != "DEFAULT"):
                kaplab18 = kaplabel8;
           
            P.plot(idsdays8,kap8per,linestyle="-",color=kaplan8color,label=kaplab8)


        if(kplots >= 9):
            l9 = float(kaplownine.get())
            h9 = float(kaphighnine.get())

            kaplandata9 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h9 >= rawData[x] >= l9):
                    kaplandata9.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata9)+1;

            x=1;
            breaker9=1;
            days = 0;
            riskset9 = len(kaplandata9)


            kap9 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata9) if x == days]
                riskset9-=len(indices)
                kap9.insert(days,riskset9)
                    
                days+=1;
                longestsurvival-=1;

            

            kap9per = [];
            x = 0;
            xx = len(kap9)
            tot = len(kaplandata9)

            while(xx>0):
                kap9per.insert(x,kap9[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays9=[];
            daysMax9 = int(len(kap9per))
            
            x=0;
            while(daysMax9>0):
                idsdays9.insert(x,x)
                x+=1;
                daysMax9-=1;
            if(kaplabel9 == "DEFAULT"):
                kaplab9 = entry+' Expression from ' + str(l9) +' to ' +str(h9)

            if(kaplabel9 != "DEFAULT"):
                kaplab19 = kaplabel9;
           
            P.plot(idsdays9,kap9per,linestyle="-",color=kaplan9color,label=kaplab9)

        if(kplots >= 10):
            l10 = float(kaplowten.get())
            h10 = float(kaphighten.get())

            kaplandata10 =[];
            

            
            xxx = int(max(Survival));
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h10 >= rawData[x] >= l10):
                    kaplandata10.insert(breaker1,int(Survival[x]))
                    breaker1+=1;
                    
                x+=1;
                xx-=1;
            x=0;
            xxx-=1;
            

                

            longestsurvival=max(kaplandata10)+1;

            x=1;
            breaker10=1;
            days = 0;
            riskset10 = len(kaplandata10)


            kap10 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata10) if x == days]
                riskset10-=len(indices)
                kap10.insert(days,riskset10)
                    
                days+=1;
                longestsurvival-=1;

            

            kap10per = [];
            x = 0;
            xx = len(kap10)
            tot = len(kaplandata10)

            while(xx>0):
                kap10per.insert(x,kap10[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays10=[];
            daysMax10 = int(len(kap10per))
            
            x=0;
            while(daysMax10>0):
                idsdays10.insert(x,x)
                x+=1;
                daysMax10-=1;
            if(kaplabel10 == "DEFAULT"):
                kaplab10 = entry+' Expression from ' + str(l10) +' to ' +str(h10)

            if(kaplabel10 != "DEFAULT"):
                kaplab110 = kaplabel10;
           
            P.plot(idsdays10,kap10per,linestyle="-",color=kaplan10color,label=kaplab10)


                    
        P.grid(True)

        P.legend()
        if(kaplantitle == "DEFAULT"):
            P.title('Survival with ' + entry);
        if(kaplantitle != "DEFAULT"):
            P.title(kaplantitle);
        if(kaplanxaxis == "DEFAULT"):
           P.xlabel("Survival Duration (days)");
        if(kaplanxaxis != "DEFAULT"):
            P.xlabel(kaplanxaxis);
        if(kaplanyaxis == "DEFAULT"):
            P.ylabel("Percent Survival");
        if(kaplanyaxis != "DEFAULT"):
            P.ylabel(kaplanyaxis);
        
        
            
        P.show()

    def chi_squared(self):
        if(chi1.get() == "Line One And"):
            lowkap = kaplandata1;
        if(chi1.get() == "Line Two And"):
            lowkap = kaplandata2;
        if(chi1.get() == "Line Three And"):
            lowkap = kaplandata3;
        if(chi1.get() == "Line Four And"):
            lowkap = kaplandata4;
        if(chi1.get() == "Line Five And"):
            lowkap = kaplandata5;
        if(chi1.get() == "Line Six And"):
            lowkap = kaplandata6;
        if(chi1.get() == "Line Seven And"):
            lowkap = kaplandata7;
        if(chi1.get() == "Line Eight And"):
            lowkap = kaplandata8;
        if(chi1.get() == "Line Nine And"):
            lowkap = kaplandata9;
        if(chi1.get() == "Line Ten And"):
            lowkap = kaplandata10;
       


        if(chi2.get() == "Line One"):
            highkap = kaplandata1;
        if(chi2.get() == "Line Two"):
            highkap = kaplandata2;
        if(chi2.get() == "Line Three"):
            lowkap = kaplandata3;
        if(chi2.get() == "Line Four"):
            highkap = kaplandata4;
        if(chi2.get() == "Line Five"):
            highkap = kaplandata5;
        if(chi2.get() == "Line Six"):
            highkap = kaplandata6;
        if(chi2.get() == "Line Seven"):
            highkap = kaplandata7;
        if(chi2.get() == "Line Eight"):
            highkap = kaplandata8;
        if(chi2.get() == "Line Nine"):
            highkap = kaplandata9;
        if(chi2.get() == "Line Ten"):
            highkap = kaplandata10;
                 
       

        
        summary, p_value, results = logrank_test(lowkap, highkap, alpha=.95)

        
        self.chisq.config(state=NORMAL)
        self.chisq.delete(0,END)
        self.chisq.insert(0,p_value)
        self.chisq.config(state=DISABLED)

    def mid_kap(self):
        with open("minpatients.csv",'r') as file:
            PATS = file.read().split(',');
            file.close()

        minpat = int(PATS[0]);
        
        step = max(rawData)/100;
        num = 1;
        pvalues = [];
        hs = [];
        z = 0;

        while(num<100):
            l1 = 0;
            h1 = step*num;
            l2 = h1+1;
            h2 = max(rawData);

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

            #second line
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

            if(len(kaplandata1) >= minpat and len(kaplandata2) >= minpat):
                lowkap = sorted(kaplandata1);
                highkap = sorted(kaplandata2);

                summary, p_value, results = logrank_test(lowkap, highkap, alpha=.95)
                pvalues.insert(z,p_value);
                hs.insert(z,h1);
                z+=1;
   
            num+=1;



        kapmin = min(pvalues)
        
        if(2 == 2):
            i = pvalues.index(kapmin)
            
            l1 =0;
            h1 = hs[i]

            h2= max(rawData);
            l2 = h1+1;
            
            

            P.figure()
            
            

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
            

                
            
            longestsurvival=max(kaplandata1)+1;

            x=1;
            breaker1=1;
            days = 0;
            riskset1 = len(kaplandata1)


            kap1 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata1) if x == days]
                riskset1-=len(indices)
                kap1.insert(days,riskset1)
                    
                days+=1;
                longestsurvival-=1;

            

            kap1per = [];
            x = 0;
            xx = len(kap1)
            tot = len(kaplandata1)

            while(xx>0):
                kap1per.insert(x,kap1[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays=[];
            daysMax = int(len(kap1per))
            
            x=0;
            while(daysMax>0):
                idsdays.insert(x,x)
                x+=1;
                daysMax-=1;
            
            kaplab1 = entry+' Expression from ' + str(l1) +' to ' +str(h1)

            
                
            P.plot(idsdays,kap1per,linestyle="-",color=kaplan1color,label= kaplab1)

            #second
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
            

                
            
            longestsurvival=max(kaplandata2)+1;

            x=1;
            breaker1=1;
            days = 0;
            riskset1 = len(kaplandata2)


            kap2 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata2) if x == days]
                riskset1-=len(indices)
                kap2.insert(days,riskset1)
                    
                days+=1;
                longestsurvival-=1;

            

            kap2per = [];
            x = 0;
            xx = len(kap2)
            tot = len(kaplandata2)

            while(xx>0):
                kap2per.insert(x,kap2[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays2=[];
            daysMax = int(len(kap2per))
            
            x=0;
            while(daysMax>0):
                idsdays2.insert(x,x)
                x+=1;
                daysMax-=1;
            
            kaplab2 = entry+' Expression from ' + str(l2) +' to ' +str(h2)

            
                
            P.plot(idsdays2,kap2per,linestyle="-",color=kaplan2color,label= kaplab2)

            
            P.grid(True)

            P.legend()

            P.title('Survival with ' + entry + '     ' + "P Value:" + str(kapmin));
            P.xlabel("Survival Duration (days)");
            P.ylabel("Percent Survival");

            P.show()


    def best_kap(self):
        with open("minpatients.csv",'r') as file:
            PATS = file.read().split(',');
            file.close()

        minpat = int(PATS[0]);

        
        step = max(rawData)/100;
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
            if(len(kaplandata1) >= minpat):

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

                    if(len(kaplandata1) >= minpat and len(kaplandata2) >= minpat):
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
            
            

            P.figure()
            
            

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
            

                
            
            longestsurvival=max(kaplandata1)+1;

            x=1;
            breaker1=1;
            days = 0;
            riskset1 = len(kaplandata1)


            kap1 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata1) if x == days]
                riskset1-=len(indices)
                kap1.insert(days,riskset1)
                    
                days+=1;
                longestsurvival-=1;

            

            kap1per = [];
            x = 0;
            xx = len(kap1)
            tot = len(kaplandata1)

            while(xx>0):
                kap1per.insert(x,kap1[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays=[];
            daysMax = int(len(kap1per))
            
            x=0;
            while(daysMax>0):
                idsdays.insert(x,x)
                x+=1;
                daysMax-=1;
            
            kaplab1 = entry+' Expression from ' + str(l1) +' to ' +str(h1)

            
                
            P.plot(idsdays,kap1per,linestyle="-",color=kaplan1color,label= kaplab1)

            #second
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
            

                
            
            longestsurvival=max(kaplandata2)+1;

            x=1;
            breaker1=1;
            days = 0;
            riskset1 = len(kaplandata2)


            kap2 = [];

            while(longestsurvival>0):
                indices = [i for i, x in enumerate(kaplandata2) if x == days]
                riskset1-=len(indices)
                kap2.insert(days,riskset1)
                    
                days+=1;
                longestsurvival-=1;

            

            kap2per = [];
            x = 0;
            xx = len(kap2)
            tot = len(kaplandata2)

            while(xx>0):
                kap2per.insert(x,kap2[x]/tot)

                
                x+=1;
                xx-=1;
            idsdays2=[];
            daysMax = int(len(kap2per))
            
            x=0;
            while(daysMax>0):
                idsdays2.insert(x,x)
                x+=1;
                daysMax-=1;
            
            kaplab2 = entry+' Expression from ' + str(l2) +' to ' +str(h2)

            
                
            P.plot(idsdays2,kap2per,linestyle="-",color=kaplan2color,label= kaplab2)

            
            P.grid(True)

            P.legend()

            P.title('Survival with ' + entry +'     ' +"P Value:" + str(kapmin));
            P.xlabel("Survival Duration (days)");
            P.ylabel("Percent Survival");

            P.show()

    



        

root = Tk()

root.geometry("1200x900")
app = Window(root)
root.mainloop()
