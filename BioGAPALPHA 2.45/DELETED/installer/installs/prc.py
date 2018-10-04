from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re
from PIL import ImageTk,Image

import matplotlib.mlab as mlab
import pylab as P
import matplotlib.colors as colors
import matplotlib
from matplotlib import pyplot as plt
matplotlib.use("tkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure


import correlatetcga
import Kaplan_Meier_Graphing_1_10_lines
import Kaplan_Meier_Graphing_expected
import kaplan_change
import bestkaplansalltwoline
import bestkaplanstwoline
import stage
import tcga_genesearch
import heatmapforexpression

import ncbi_cellline_search
import heatmap_cellline_expression

import ncbi_data_pull

import ncbi_correlate
import heatmapforexpressionncbi


tcgabackground = "#cce6ff"
ncbitbackground = "#d9f2d9"
ncbicbackground = "#ffeecc"

GENE_TCGA = "Search Gene";
GENEC_TCGA = "Correlate to Gene";

global STAGES, SEXES, AGES ;


class Window(Frame):
    

    def __init__(self, master=None):

        Frame.__init__(self, master)
        self.master = master
        
        
        self.init_window()
        
        

    def init_window(self):


            
        menu = Menu(self.master)
        self.master.config(menu=menu)

        def client_exit():
            exit()
        def client_help():
            os.system("Help.pdf");
            
        file = Menu(menu)
        file.add_command(label="Help", command=client_help)
        file.add_command(label="Exit", command=client_exit)
        menu.add_cascade(label="File", menu=file)


##        def kaplanchange():
##            kaplan_change.kaplan_graph_change()
##
##        graphs = Menu(menu)
##        graphs.add_command(label="Stage Graph")
##        graphs.add_command(label="Kaplan Meier",command = kaplanchange)
##        
##        
##        menu.add_cascade(label="Graph Options", menu=graphs)

        self.master.title("Biomarker Graphical Anaylisis Program")

        canvas=Canvas(bg="#e6ffff",width=200,height=200,scrollregion=(0,0,1300,1300))
       

        load = Image.open("Program Setup\Pictures\Title.png")
        render = ImageTk.PhotoImage(load)
        img = Label(image = render)
        img.image = render
        backgroundLabel_window = canvas.create_window(600,450,anchor=CENTER,tag = 'all',window=img)
        


        

        vbar=Scrollbar(root,orient=VERTICAL)
        vbar.pack(side=RIGHT,fill=Y)
        vbar.config(command=canvas.yview)

        canvas.config(width=200,height=200)
        LABEL=Canvas(bg="WHITE",width=10,height=20,scrollregion=(0,0,1300,1300))
        TITLES = Label(text="Aaron Schwan alpha 2.45")
        TITLES.configure(bg="WHITE")
        TITLES_window = LABEL.create_window(5,5,anchor=NW, window=TITLES)
        LABEL.pack(side=BOTTOM,fill=BOTH)
        

        canvas.config(yscrollcommand=vbar.set)
        canvas.pack(side=BOTTOM,expand=True,fill=BOTH)

        

            
        def normal_data():
            Average.config(state = NORMAL)
            STDEV.config(state = NORMAL)
            IQR.config(state = NORMAL)
            Range.config(state = NORMAL)
            Minimum.config(state = NORMAL)
            FQM.config(state = NORMAL)
            Median.config(state = NORMAL)
            TQM.config(state = NORMAL)
            Maximum.config(state = NORMAL)

            Average.delete(0,END)
            STDEV.delete(0,END)
            IQR.delete(0,END)
            Range.delete(0,END)
            Minimum.delete(0,END)
            FQM.delete(0,END)
            Median.delete(0,END)
            TQM.delete(0,END)
            Maximum.delete(0,END)

            Average.insert(0,mean(rawData))
            STDEV.insert(0,np.std([rawData]))
            rang = max(rawData)-min(rawData)
            Range.insert(0,rang)
            q75, q25 = np.percentile(rawData,[75,25])
            iqrnum = q75-q25
            IQR.insert(0,iqrnum)
            Minimum.insert(0,min(rawData))
            FQM.insert(0,q25)
            Median.insert(0,median(rawData))
            TQM.insert(0,q75)
            Maximum.insert(0,max(rawData))


            Average.config(state = 'readonly')
            STDEV.config(state = 'readonly')
            IQR.config(state = 'readonly')
            Range.config(state = 'readonly')
            Minimum.config(state = 'readonly')
            FQM.config(state = 'readonly')
            Median.config(state = 'readonly')
            TQM.config(state = 'readonly')
            Maximum.config(state = 'readonly')

        def cellline_search_ncbi():
            global FINAL
            
            pdata = {}
            
            FINAL = []
            del FINAL[:]
            
            cellline = CellLineE.get()

            with open('Data Sources\Celllines\cell_lines.csv','r') as file:
                data = file.read().split('\n');
                file.close();
                


            data = list(filter(None, data))

            xx = len(data)-1; 
            x = 0;
            while(xx>=0):
                stuff = data[x].split(',')

                
                stuff1 = [str(stuff[1]),str(stuff[2]),str(stuff[3]),str(stuff[4]),str(stuff[5]),str(stuff[6]),str(stuff[7])]
                pdata[stuff[0]] = stuff1



                    
                    
                x+=1;
                xx-=1;
            
            try:
                data = pdata[str(cellline+'1')];
                FINAL.insert(0,data)
                
            except:
                pass;

            try:
                
                data = pdata[str(cellline+'2')];
                FINAL.insert(1,data)
            except:
                pass;

            try:
                data = pdata[str(cellline+'3')];
                FINAL.insert(2,data)
            except:
                pass;

            try:
                data = pdata[str(cellline+'4')];
                FINAL.insert(3,data)
            except:
                pass;

            
            
            if(len(FINAL) == 0):
                FINAL = ncbi_cellline_search.searchcellline(cellline)

            x = len(FINAL)

            global varc1
            global varc2
            global varc3
            global varc4
            global c1
            global c2
            global c3
            global c4

            varc1 = IntVar()
            varc2= IntVar()
            varc3 = IntVar()
            varc4 = IntVar()
          
            canvas.delete('check')

                
           
            

            
            if(x>=1):
              
                    
                info__ = FINAL[0]
                info = "Title:", info__[0],"Sample Type:", info__[1],"Source Name:", info__[2],"Description:", info__[3]
                info = str(info).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'')
                
                c1 = Checkbutton(variable=varc1)
                c1.configure(relief = FLAT,text = info, bg=ncbicbackground, activebackground = ncbicbackground)
                c1_window = canvas.create_window(100, 50, anchor=NW,tag = 'check', window=c1)

            if(x>=2):
                info__ = FINAL[1]
                info = "Title:", info__[0],"Sample Type:", info__[1],"Source Name:", info__[2],"Description:", info__[3]
                info = str(info).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'')

                c2 = Checkbutton(variable=varc2)
                c2.configure(relief = FLAT,text = info, bg=ncbicbackground, activebackground = ncbicbackground)
                c2_window = canvas.create_window(100, 100, anchor=NW,tag = 'check', window=c2)

            if(x>=3):
                info__ = FINAL[2]
                info = "Title:", info__[0],"Sample Type:", info__[1],"Source Name:", info__[2],"Description:", info__[3]
                info = str(info).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'')

                c3 = Checkbutton(variable=varc3)
                c3.configure(relief = FLAT,text = info, bg=ncbicbackground, activebackground = ncbicbackground)
                c3_window = canvas.create_window(100, 150, anchor=NW,tag = 'check', window=c3)

            if(x>=4):
                info__ = FINAL[3]
                info = "Title:", info__[0],"Sample Type:", info__[1],"Source Name:", info__[2],"Description:", info__[3]
                info = str(info).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'')

                c4 = Checkbutton(variable=varc4)
                c4.configure(relief = FLAT,text = info, bg=ncbicbackground, activebackground = ncbicbackground)
                c4_window = canvas.create_window(100, 200, anchor=NW,tag = 'check', window=c4)

            
            
        



        def search_prob():
            global urls
            
            urls = [];
            
            breaker = 0;
            
            if(varc1.get() == 1):
                data = str(FINAL[0]).split(',')
                urls.insert(breaker,str(data[5]).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'').replace(" ",''))

                breaker+=1;
                
            if(varc2.get() == 1):
                data = str(FINAL[1]).split(',')
                urls.insert(breaker,str(data[5]).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'').replace(" ",''))

                breaker+=1;
                
            if(varc3.get() == 1):
                data = str(FINAL[2]).split(',')
                urls.insert(breaker,str(data[5]).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'').replace(" ",''))

                breaker+=1;
                
            if(varc4.get() == 1):
                data = str(FINAL[3]).split(',')
                urls.insert(breaker,str(data[5]).replace('\n',';').replace('{',',').replace('}','').replace('"','').replace("'",'').replace('(','').replace(")",'').replace(" ",''))

                breaker+=1;

            prob = SearchProbE.get();

            dataT = ncbi_data_pull.search_ncbi(urls,prob)

                
            if(out.get() == "Include Outliers"):
                data = dataT
            if(out.get() == "Exclude Outliers"):

                
                q75, q25 = np.percentile(dataT,[75,25])

                iqrnum = q75-q25;

                lowout = q25-1.5*iqrnum;
                highout = q75+1.5*iqrnum;

                

                
                
                breaker = 0;
                x=0;

                xx = len(dataT)-1;
                while(xx>0):
                    
                    if(highout>=dataT[x]>=lowout ):
                        data.insert(breaker,dataT[x])

                        
                        
                        breaker+=1


                    xx-=1;
                    x+=1;

            global rawData
            rawData = data;

            normal_data();
        def cellline_data_graph():
          
            P.figure()
            rawtitle = " Data for "+ CellLineE.get() +" Expression\n  "
            expression = rawData
            ids = [x for x in range(len(rawData))]
            P.bar(ids,sorted(rawData))
            P.ylabel('Expression of '+ CellLineE.get()+'in RNA Seq V2')
            P.xlabel('Identification Number')
            P.title(rawtitle)
            P.grid(True)
            P.show()

        def ncbic_heatmap():
            probH = heatmapSearch.get();
            heatmap_cellline_expression.heatmap(probH);
            
        def ncbiT_search():
            data = [];
            
            if(var.get() == 'Small Cell Lung Cancer'):
                with open('Data Sources\_NCBI_\smallcell.csv','r') as file:
                    urls = file.read().split('\n')
                    file.close()

            if(var.get() == 'Lung Adenocarcinoma'):
                with open('Data Sources\_NCBI_\lungadeno.csv','r') as file:
                    urls = file.read().split('\n')
                    file.close()

            if(var.get() == 'Lung Squamous Carcinoma'):
                with open('Data Sources\_NCBI_\lungsquamous.csv','r') as file:
                    urls = file.read().split('\n')
                    file.close()

            prob =  geneSearch.get()
            
            if(len(urls) > 0 and prob != 'Search Prob Number'):      

                dataT = ncbi_data_pull.search_ncbi(urls,prob)
                
            if(out.get() == "Include Outliers"):
                data = dataT
            if(out.get() == "Exclude Outliers"):

                
                q75, q25 = np.percentile(dataT,[75,25])

                iqrnum = q75-q25;

                lowout = q25-1.5*iqrnum;
                highout = q75+1.5*iqrnum;

                

                
                
                breaker = 0;
                x=0;

                xx = len(dataT)-1;
                while(xx>0):
                    
                    if(highout>=dataT[x]>=lowout ):
                        data.insert(breaker,dataT[x])

                        
                        
                        breaker+=1


                    xx-=1;
                    x+=1;

            global rawData
            rawData = data;

            normal_data();
        def gene_correlate_ncbiT():
            outliers = out.get();
            dataS = var.get();
            prob = geneSearch.get();
            probC = geneCSearch.get();

            r_squared,ttest,rawDataC,rawData,regline = ncbi_correlate.correlate(prob,probC,dataS,outliers);

            RSQ.config(state = NORMAL)
            TT.config(state = NORMAL)

            RSQ.delete(0,END)
            TT.delete(0,END)

            RSQ.insert(0, r_squared)
            TT.insert(0, ttest)

            RSQ.config(state = 'readonly')
            
            TT.config(state = 'readonly')

            
            idsforline = [];
            idsforline.insert(0,min(rawData))
            idsforline.insert(1,max(rawData))
            
                           
            

              
            P.figure()
            
            P.scatter(rawData,rawDataC)

            if(YNREG.get() == 1):
                
                P.plot(idsforline,regline,linestyle="--",color="grey")
                
            P.ylabel('Expression of '+ geneCSearch.get() +' in RNA Seq V2')
            P.xlabel('Expression of '+ geneSearch.get() +' in RNA Seq V2')
                    
            P.title("Correlation of " + geneSearch.get() +" to " + geneCSearch.get())
            P.grid(True)
            P.show()

        def ncbi_heatmap():
            dataS = var.get();
            outliers = out.get();
            
            label_ = heatmapSearch.get();
            labels = label_.split(',');
            heatmapforexpressionncbi.heatmap(labels,dataS,outliers)
            
        def gene_search_special_tcga():

            Average.config(state = NORMAL)
            STDEV.config(state = NORMAL)
            IQR.config(state = NORMAL)
            Range.config(state = NORMAL)
            Minimum.config(state = NORMAL)
            FQM.config(state = NORMAL)
            Median.config(state = NORMAL)
            TQM.config(state = NORMAL)
            Maximum.config(state = NORMAL)
            RSQ.config(state = NORMAL)
            TT.config(state = NORMAL)

            Average.delete(0,END)
            STDEV.delete(0,END)
            IQR.delete(0,END)
            Range.delete(0,END)
            Minimum.delete(0,END)
            FQM.delete(0,END)
            Median.delete(0,END)
            TQM.delete(0,END)
            Maximum.delete(0,END)
            RSQ.delete(0,END)
            TT.delete(0,END)

            Average.config(state = 'readonly')
            STDEV.config(state = 'readonly')
            IQR.config(state = 'readonly')
            Range.config(state = 'readonly')
            Minimum.config(state = 'readonly')
            FQM.config(state = 'readonly')
            Median.config(state = 'readonly')
            TQM.config(state = 'readonly')
            Maximum.config(state = 'readonly')
            RSQ.config(state = 'readonly')
            TT.config(state = 'readonly')

            
            
            global entry 
            entry = geneSearch.get()
            
            global rawData
            global dataS
            dataS = var.get()
            
            global outliers

            outliers = out.get()

            global rawData
            global AliveDead
            global Censored
            global StageData
            global Survival

            global STAGES,AGES,SEXES,RACES

            STAGES = [];
            SEXES = [];
            RACES = [];
            AGES = [];
            sx = 0;
            sexx = 0;
            rx = 0;
        

            if(varS1.get() == 1):
                STAGES.insert(sx,1)
                sx+=1;
            if(varS2.get() == 1):
                STAGES.insert(sx,2)
                sx+=1;
            if(varS3.get() == 1):
                STAGES.insert(sx,3)
                sx+=1;
            if(varS4.get() == 1):
                STAGES.insert(sx,4)
                sx+=1;

            if(varR1.get() == 1):
                RACES.insert(rx,'white')
                rx+=1;
            if(varR2.get() == 1):
                RACES.insert(rx,'black or african american')
                rx+=1;
            if(varR3.get() == 1):
                RACES.insert(rx,'asian')
                rx+=1;
            if(varR4.get() == 1):
                RACES.insert(rx,'not reported')
                rx+=1;

            if(varMC.get() == 1):
                SEXES.insert(sexx,1)
                sexx+=1;
            if(varFC.get() == 1):
                SEXES.insert(sexx,0)
                sexx+=1;

            if(AEL.get() != '' and AEH.get() != ''):
                AGES.insert(0,int(AEL.get()))
                AGES.insert(1,int(AEH.get()))

            if(AEL.get() == '' or AEH.get() == ''):
                AGES=[1,1000000000000000];
                
                
                

            
            
            
            StageData,rawData,RaceData,SexData,Survival,Censored,AliveDead,AgeData = tcga_genesearch.geneSearch(entry,outliers,dataS,STAGES,AGES,SEXES,RACES)
            
            

        
            normal_data();
            
        def gene_search_tcga():
            
            Average.config(state = NORMAL)
            STDEV.config(state = NORMAL)
            IQR.config(state = NORMAL)
            Range.config(state = NORMAL)
            Minimum.config(state = NORMAL)
            FQM.config(state = NORMAL)
            Median.config(state = NORMAL)
            TQM.config(state = NORMAL)
            Maximum.config(state = NORMAL)
            RSQ.config(state = NORMAL)
            TT.config(state = NORMAL)

            Average.delete(0,END)
            STDEV.delete(0,END)
            IQR.delete(0,END)
            Range.delete(0,END)
            Minimum.delete(0,END)
            FQM.delete(0,END)
            Median.delete(0,END)

            TQM.delete(0,END)
            Maximum.delete(0,END)
            RSQ.delete(0,END)
            TT.delete(0,END)

            Average.config(state = 'readonly')
            STDEV.config(state = 'readonly')
            IQR.config(state = 'readonly')
            Range.config(state = 'readonly')
            Minimum.config(state = 'readonly')
            FQM.config(state = 'readonly')
            Median.config(state = 'readonly')
            TQM.config(state = 'readonly')
            Maximum.config(state = 'readonly')
            RSQ.config(state = 'readonly')
            TT.config(state = 'readonly')

            
            
            global entry 
            entry = geneSearch.get()
            
            global rawData
            global dataS
            dataS = var.get()
            
            global outliers

            outliers = out.get()

            global rawData
            global AliveDead
            global Censored
            global StageData
            global Survival

            global STAGES,AGES,SEXES,RACES
            

            STAGES=[1,2,3,4]
            AGES=[0,10000]
            SEXES=[1,0]
            RACES=['white','not reported','black or african american','asian']

  

            
            StageData,rawData,RaceData,SexData,Survival,Censored,AliveDead,AgeData = tcga_genesearch.geneSearch(entry,outliers,dataS,STAGES,AGES,SEXES,RACES)
            
            

        
            normal_data();



        def gene_correlate_tcga():
            corre = geneCSearch.get()
            
            Rsquared,ttest,rawDataC,rawData,regline = correlatetcga.correlate(corre,entry,dataS,outliers,STAGES,AGES,SEXES,RACES)
            


            

            RSQ.config(state = NORMAL)
            TT.config(state = NORMAL)

            RSQ.delete(0,END)
            TT.delete(0,END)

            RSQ.insert(0, Rsquared)
            TT.insert(0, ttest)

            RSQ.config(state = 'readonly')
            
            TT.config(state = 'readonly')

            
            idsforline = [];
            idsforline.insert(0,min(rawData))
            idsforline.insert(1,max(rawData))
            
                           
            

              
            P.figure()
            
            P.scatter(rawData,rawDataC)

            if(YNREG.get() == 1):
                
                P.plot(idsforline,regline,linestyle="--",color="grey")
                
            P.ylabel('Expression of '+ geneCSearch.get() +' in RNA Seq V2')
            P.xlabel('Expression of '+ geneSearch.get() +' in RNA Seq V2')
                    
            P.title("Correlation of " + geneSearch.get() +" to " + geneCSearch.get())
            P.grid(True)
            P.show()

        def tcga_heatmap():
            label_ = heatmapSearch.get();
            labels = label_.split(',');
            heatmapforexpression.heatmap(labels,dataS,outliers,STAGES,AGES,SEXES,RACES)
            
        def tcga_kaplan():
            Low1 = L1.get()
            Low2 = L2.get()
            Low3 = L3.get()
            Low4 = L4.get()
            Low5 = L5.get()
            Low6 = L6.get()
            Low7 = L7.get()
            Low8 = L8.get()
            Low9 = L9.get()
            Low10 = L10.get()

            High1 = H1.get()
            High2 = H2.get()
            High3 = H3.get()
            High4 = H4.get()
            High5 = H5.get()
            High6 = H6.get()
            High7 = H7.get()
            High8 = H8.get()
            High9 = H9.get()
            High10 = H10.get()

            with open('Program Setup\KaplanMeier1-10\KaplanMeierLinesLow.csv','w') as file:
                file.write(Low1)
                file.write('\n')
                file.write(Low2)
                file.write('\n')
                file.write(Low3)
                file.write('\n')
                file.write(Low4)
                file.write('\n')
                file.write(Low5)
                file.write('\n')
                file.write(Low6)
                file.write('\n')
                file.write(Low7)
                file.write('\n')
                file.write(Low8)
                file.write('\n')
                file.write(Low9)
                file.write('\n')
                file.write(Low10)
                file.write('\n')
                file.close()

            with open('Program Setup\KaplanMeier1-10\KaplanMeierLinesHigh.csv','w') as file:
                file.write(High1)
                file.write('\n')
                file.write(High2)
                file.write('\n')
                file.write(High3)
                file.write('\n')
                file.write(High4)
                file.write('\n')
                file.write(High5)
                file.write('\n')
                file.write(High6)
                file.write('\n')
                file.write(High7)
                file.write('\n')
                file.write(High8)
                file.write('\n')
                file.write(High9)
                file.write('\n')
                file.write(High10)
                file.write('\n')
                file.close()

            with open('Program Setup\KaplanMeier1-10\KaplanMeierLines.csv','w') as file:
                file.write(str(varg1.get()))
                file.write('\n')
                file.write(str(varg2.get()))
                file.write('\n')
                file.write(str(varg3.get()))
                file.write('\n')
                
                file.write(str(varg4.get()))
                file.write('\n')
                file.write(str(varg5.get()))
                file.write('\n')
                file.write(str(varg6.get()))
                file.write('\n')
                file.write(str(varg7.get()))
                file.write('\n')
                file.write(str(varg8.get()))
                file.write('\n')
                file.write(str(varg9.get()))
                file.write('\n')
                file.write(str(varg10.get()))
                file.write('\n')
                file.close()
            

            
            Kaplan_Meier_Graphing_1_10_lines.kaplanmeiers(entry,dataS,rawData,AliveDead,Survival)


            
        def bestkaps():
            if(all_pat.get() == 1):

                if(kaplineamount.get() == "2 Line"):
                    bestkaplansalltwoline.kaplanmeiers(entry,rawData,Survival,AliveDead)

            if(all_pat.get() == 0):

                if(kaplineamount.get() == "2 Line"):
                    bestkaplanstwoline.kaplanmeiers(entry,rawData,Survival,AliveDead)
                
            else:
                pass;


        def pval_kaps():
            dispatcherLow = {'1':L1,'2':L2,'3':L3,'4':L4,'5':L5,'6':L6,'7':L7,'8':L8,'9':L9,'10':L10}
            dispatcherHigh = {'1':H1,'2':H2,'3':H3,'4':H4,'5':H5,'6':H6,'7':H7,'8':H8,'9':H9,'10':H10}
            D1 = line1.get()
            l1 = float(dispatcherLow[D1].get())
            h1 = float(dispatcherHigh[D1].get())

            

            kaplandata1 = [];
            kaplancen1 = [];
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
            
            D2 = line2.get()
            l2 = float(dispatcherLow[D2].get())
            h2 = float(dispatcherHigh[D2].get())

            

            kaplandata2 = [];
            kaplancen2 = [];
            breaker1=0;
            
            
            x = 0;
            xx = len(rawData)
            while(xx>0):
                if(h2 >= rawData[x] >= l2):
                    
                    kaplandata2.insert(breaker1,Survival[x])
                    kaplancen2.insert(breaker1,AliveDead[x])
                    breaker1+=1;
                    
                x+=1;
                
                xx-=1;

            pvalforlines = Kaplan_Meier_Graphing_expected.kaplan(kaplandata1,kaplancen1,kaplandata2,kaplancen2)

            pvalE.config(state = 'normal')
            pvalE.delete(0,END)
            pvalE.insert(0,pvalforlines)
            pvalE.config(state='readonly')
            

    
        
        

        def stage_data_graph():
            medians = varMC1.get()

            with open("Program Setup\Colors\Stage.csv",'r') as file:
                StageColors = file.read().split('\n')
                file.close();

            with open("Program Setup\Titles\Stage.csv",'r') as file:
                StageTitles = file.read().split('\n')
                file.close();
                
            lengths = int(len(rawData)*2);
            
            stage.stage_data(entry,dataS,rawData,StageData,StageColors,medians,StageTitles,lengths)

        def patient_data_graph():
          
            P.figure()
            rawtitle = var.get()+" Data for "+ geneSearch.get() +" Expression\n  "
            expression = rawData
            ids = [x for x in range(len(rawData))]
            P.bar(ids,sorted(rawData))
            P.ylabel('Expression of '+ geneSearch.get()+'in RNA Seq V2')
            P.xlabel('Patient Identification Number')
            P.title(rawtitle)
            P.grid(True)
            P.show()

        def full_report_tcga():
            organ = ORGAN.get()
            outlier = IncExFR.get()

            if(organ != "Organ:" and entry != "Search Gene"):
                os.system("tcgafullreport.py 1");

        def gene_search_tcga_less_options():
            Average_window = canvas.create_window(250, 75, anchor=NW,tag = 'all', window=Average)
            STDEV_window = canvas.create_window(250, 95, anchor=NW,tag = 'all', window=STDEV)
            IQR_window = canvas.create_window(250, 115, anchor=NW,tag = 'all', window=IQR)
            Range_window = canvas.create_window(250, 135, anchor=NW,tag = 'all', window=Range)
            Minimum_window = canvas.create_window(250, 155, anchor=NW,tag = 'all', window=Minimum)
            FQM_window = canvas.create_window(250, 175, anchor=NW,tag = 'all', window=FQM)
            Median_window = canvas.create_window(250, 195, anchor=NW,tag = 'all', window=Median)
            TQM_window = canvas.create_window(250, 215, anchor=NW,tag = 'all', window=TQM)
            Maximum_window = canvas.create_window(250, 235, anchor=NW,tag = 'all', window=Maximum)
            stan_window = canvas.create_window(30,55, anchor=NW ,tag = 'all',window=standard)
            ave_window = canvas.create_window(147,75, anchor=NW ,tag = 'all',window=ave)
            stdev_window = canvas.create_window(120,95, anchor=NW ,tag = 'all',window=stdev)
            iqr_window = canvas.create_window(120,115, anchor=NW ,tag = 'all',window=iqr)
            ran_window = canvas.create_window(150,135, anchor=NW ,tag = 'all',window=ran)
            minimum_window = canvas.create_window(140,155, anchor=NW ,tag = 'all',window=minimum)
            fqm_window = canvas.create_window(115,175, anchor=NW ,tag = 'all',window=fqm)
            med_window = canvas.create_window(145,195, anchor=NW ,tag = 'all',window=med)
            tqm_window = canvas.create_window(115,215, anchor=NW ,tag = 'all',window=tqm)
            maximum_window = canvas.create_window(137,235, anchor=NW ,tag = 'all',window=maximum)
            corre_window = canvas.create_window(30,275, anchor=NW ,tag = 'all',window=corre)
            reglineYN_window = canvas.create_window(250, 360, anchor=NW,tag = 'all', window=reglineYN)
            rsq_window = canvas.create_window(175,320, anchor=NW ,tag = 'all',window=rsq)
            tt_window = canvas.create_window(185,340, anchor=NW ,tag = 'all',window=tt)
            geneCSearch_window = canvas.create_window(100, 295, anchor=NW, window=geneCSearch)
            RSQ_window = canvas.create_window(250, 320, anchor=NW,tag = 'all', window=RSQ)
            TT_window = canvas.create_window(250, 340, anchor=NW,tag = 'all', window=TT)
            buttonC_window = canvas.create_window(250, 290, anchor=NW,tag = 'all', window=ButtonC)
            kapti_window = canvas.create_window(30,390, anchor=NW ,tag = 'all',window=kapti)
            kaptiL_window = canvas.create_window(295,410, anchor=NW ,tag = 'all',window=kaptiL)
            kaptiH_window = canvas.create_window(445,410, anchor=NW ,tag = 'all',window=kaptiH)
            kaptiL1_window = canvas.create_window(170,430, anchor=NW ,tag = 'all',window=kaptiL1)
            kaptiL2_window = canvas.create_window(170,450, anchor=NW ,tag = 'all',window=kaptiL2)
            kaptiL3_window = canvas.create_window(170,470, anchor=NW ,tag = 'all',window=kaptiL3)
            kaptiL4_window = canvas.create_window(170,490, anchor=NW ,tag = 'all',window=kaptiL4)
            kaptiL5_window = canvas.create_window(170,510, anchor=NW ,tag = 'all',window=kaptiL5)
            kaptiL6_window = canvas.create_window(170,530, anchor=NW ,tag = 'all',window=kaptiL6)
            kaptiL7_window = canvas.create_window(170,550, anchor=NW ,tag = 'all',window=kaptiL7)
            kaptiL8_window = canvas.create_window(170,570, anchor=NW ,tag = 'all',window=kaptiL8)
            kaptiL9_window = canvas.create_window(170,590, anchor=NW ,tag = 'all',window=kaptiL9)
            kaptiL10_window = canvas.create_window(170,610, anchor=NW ,tag = 'all',window=kaptiL10)
            L1_window = canvas.create_window(250, 430, anchor=NW,tag = 'all', window=L1)
            L2_window = canvas.create_window(250, 450, anchor=NW,tag = 'all', window=L2)
            L3_window = canvas.create_window(250, 470, anchor=NW,tag = 'all', window=L3)
            L4_window = canvas.create_window(250, 490, anchor=NW,tag = 'all', window=L4)
            L5_window = canvas.create_window(250, 510, anchor=NW,tag = 'all', window=L5)
            L6_window = canvas.create_window(250, 530, anchor=NW,tag = 'all', window=L6)
            L7_window = canvas.create_window(250, 550, anchor=NW,tag = 'all', window=L7)
            L8_window = canvas.create_window(250, 570, anchor=NW,tag = 'all', window=L8)
            L9_window = canvas.create_window(250, 590, anchor=NW,tag = 'all', window=L9)
            L10_window = canvas.create_window(250, 610, anchor=NW,tag = 'all', window=L10)
            H1_window = canvas.create_window(405, 430, anchor=NW,tag = 'all', window=H1)
            H2_window = canvas.create_window(405, 450, anchor=NW,tag = 'all', window=H2)
            H3_window = canvas.create_window(405, 470, anchor=NW,tag = 'all', window=H3)
            H4_window = canvas.create_window(405, 490, anchor=NW,tag = 'all', window=H4)
            H5_window = canvas.create_window(405, 510, anchor=NW,tag = 'all', window=H5)
            H6_window = canvas.create_window(405, 530, anchor=NW,tag = 'all', window=H6)
            H7_window = canvas.create_window(405, 550, anchor=NW,tag = 'all', window=H7)
            H8_window = canvas.create_window(405, 570, anchor=NW,tag = 'all', window=H8)
            H9_window = canvas.create_window(405, 590, anchor=NW,tag = 'all', window=H9)
            H10_window = canvas.create_window(405, 610, anchor=NW,tag = 'all', window=H10)
            g1_window = canvas.create_window(560, 427, anchor=NW,tag = 'all', window=g1)
            g2_window = canvas.create_window(560, 447, anchor=NW,tag = 'all', window=g2)
            g3_window = canvas.create_window(560, 467, anchor=NW,tag = 'all', window=g3)
            g4_window = canvas.create_window(560, 487, anchor=NW,tag = 'all', window=g4)
            g5_window = canvas.create_window(560, 507, anchor=NW,tag = 'all', window=g5)
            g6_window = canvas.create_window(560, 527, anchor=NW,tag = 'all', window=g6)
            g7_window = canvas.create_window(560, 547, anchor=NW,tag = 'all', window=g7)
            g8_window = canvas.create_window(560, 567, anchor=NW,tag = 'all', window=g8)
            g9_window = canvas.create_window(560, 587, anchor=NW,tag = 'all', window=g9)
            g10_window = canvas.create_window(560, 607, anchor=NW,tag = 'all', window=g10)
            KAPGRAPH_window = canvas.create_window(250, 630, anchor=NW,tag = 'all', window=KAPGRAPH)
            kaplineamountoption_window = canvas.create_window(170, 750, anchor=NW,tag = 'all',  window=kaplineamountoption)
            kapL_window = canvas.create_window(30,720, anchor=NW ,tag = 'all',window=kapL)
            allpat_window = canvas.create_window(435, 752, anchor=NW,tag = 'all', window=allpat)
            KAPGRAPHALL_window = canvas.create_window(250, 752, anchor=NW,tag = 'all', window=KAPGRAPHALL)
            option1_window = canvas.create_window(250, 660, anchor=NW,tag = 'all',  window=option1)
            option2_window = canvas.create_window(330, 660, anchor=NW,tag = 'all',  window=option2)
            pvalB_window = canvas.create_window(420, 662, anchor=NW,tag = 'all', window=pvalB)
            pvalE_window = canvas.create_window(480, 664, anchor=NW,tag = 'all', window=pvalE)
            stageL_window = canvas.create_window(30,790, anchor=NW ,tag = 'all',window=stageL)
            StageB_window = canvas.create_window(170, 830, anchor=NW,tag = 'all', window=StageB)
            medianC_window = canvas.create_window(400, 830, anchor=NW,tag = 'all', window=medianC)
            pDataL_window = canvas.create_window(30,890, anchor=NW ,tag = 'all',window=pDataL)
            pDataB_window = canvas.create_window(170, 930, anchor=NW,tag = 'all', window=pDataB)
##            FRL_window = canvas.create_window(30,970, anchor=NW ,tag = 'all',window=FRL)
##            FRB_window = canvas.create_window(270, 1003, anchor=NW,tag = 'all', window=FRB)
##            option1FR_window = canvas.create_window(170, 1000, anchor=NW,tag = 'all',  window=option1FR)
##            option2FR_window = canvas.create_window(378, 1000, anchor=NW,tag = 'all',  window=option2FR)
            ButtonHeat_window = canvas.create_window(500, 315, anchor=NW,tag = 'all', window=ButtonHeat)
            heatmapSearch_window = canvas.create_window(500, 295, anchor=NW, window=heatmapSearch)

            canvas.delete(Less_Options)
            canvas.delete(sexL_window)
            canvas.delete(FC_window)
            canvas.delete(MC_window)
            canvas.delete(stageLL_window)
            canvas.delete(S1_window)
            canvas.delete(S2_window)
            canvas.delete(S3_window)
            canvas.delete(S4_window)
            canvas.delete(raceLL_window)
            canvas.delete(ageLL_window)
            canvas.delete(ageLLTO_window)
            canvas.delete(R1_window)
            canvas.delete(R2_window)
            canvas.delete(R3_window)
            canvas.delete(R4_window)
            canvas.delete(AEH_window)
            canvas.delete(AEL_window)

            global button1_window
            
            canvas.delete(button1_window)

            
            
            Button1 = Button(text = "Search", command = gene_search_tcga, anchor = W)
            Button1.configure(width = 10, activebackground = "#33B5E5",relief = FLAT)
            button1_window = canvas.create_window(700, 10, anchor=NW,tag = 'all', window=Button1)

            
            More_Options = Button(text = " More Options ",fg = 'blue', command = gene_search_tcga_more_options, anchor = W)
            More_Options.configure(width = 18, activebackground = "#33B5E5",bg=tcgabackground, relief = FLAT)
            More_Options_window = canvas.create_window(800, 10, anchor=NW,tag = 'all', window=More_Options)


            

            
            

        def gene_search_tcga_more_options():
            
            Average_window = canvas.create_window(250, 175, anchor=NW,tag = 'all', window=Average)
            STDEV_window = canvas.create_window(250, 195, anchor=NW,tag = 'all', window=STDEV)
            IQR_window = canvas.create_window(250, 215, anchor=NW,tag = 'all', window=IQR)
            Range_window = canvas.create_window(250, 235, anchor=NW,tag = 'all', window=Range)
            Minimum_window = canvas.create_window(250, 255, anchor=NW,tag = 'all', window=Minimum)
            FQM_window = canvas.create_window(250, 275, anchor=NW,tag = 'all', window=FQM)
            Median_window = canvas.create_window(250, 295, anchor=NW,tag = 'all', window=Median)
            TQM_window = canvas.create_window(250, 315, anchor=NW,tag = 'all', window=TQM)
            Maximum_window = canvas.create_window(250, 335, anchor=NW,tag = 'all', window=Maximum)
            stan_window = canvas.create_window(30,155, anchor=NW ,tag = 'all',window=standard)
            ave_window = canvas.create_window(147,175, anchor=NW ,tag = 'all',window=ave)
            stdev_window = canvas.create_window(120,195, anchor=NW ,tag = 'all',window=stdev)
            iqr_window = canvas.create_window(120,215, anchor=NW ,tag = 'all',window=iqr)
            ran_window = canvas.create_window(150,235, anchor=NW ,tag = 'all',window=ran)
            minimum_window = canvas.create_window(140,255, anchor=NW ,tag = 'all',window=minimum)
            fqm_window = canvas.create_window(115,275, anchor=NW ,tag = 'all',window=fqm)
            med_window = canvas.create_window(145,295, anchor=NW ,tag = 'all',window=med)
            tqm_window = canvas.create_window(115,315, anchor=NW ,tag = 'all',window=tqm)
            maximum_window = canvas.create_window(137,335, anchor=NW ,tag = 'all',window=maximum)
            corre_window = canvas.create_window(30,375, anchor=NW ,tag = 'all',window=corre)
            reglineYN_window = canvas.create_window(250, 460, anchor=NW,tag = 'all', window=reglineYN)
            rsq_window = canvas.create_window(175,420, anchor=NW ,tag = 'all',window=rsq)
            tt_window = canvas.create_window(185,440, anchor=NW ,tag = 'all',window=tt)
            geneCSearch_window = canvas.create_window(100, 395, anchor=NW, window=geneCSearch)
            RSQ_window = canvas.create_window(250, 420, anchor=NW,tag = 'all', window=RSQ)
            TT_window = canvas.create_window(250,440, anchor=NW,tag = 'all', window=TT)
            buttonC_window = canvas.create_window(250, 390, anchor=NW,tag = 'all', window=ButtonC)
            kapti_window = canvas.create_window(30,490, anchor=NW ,tag = 'all',window=kapti)
            kaptiL_window = canvas.create_window(295,510, anchor=NW ,tag = 'all',window=kaptiL)
            kaptiH_window = canvas.create_window(445,510, anchor=NW ,tag = 'all',window=kaptiH)
            kaptiL1_window = canvas.create_window(170,530, anchor=NW ,tag = 'all',window=kaptiL1)
            kaptiL2_window = canvas.create_window(170,550, anchor=NW ,tag = 'all',window=kaptiL2)
            kaptiL3_window = canvas.create_window(170,570, anchor=NW ,tag = 'all',window=kaptiL3)
            kaptiL4_window = canvas.create_window(170,590, anchor=NW ,tag = 'all',window=kaptiL4)
            kaptiL5_window = canvas.create_window(170,610, anchor=NW ,tag = 'all',window=kaptiL5)
            kaptiL6_window = canvas.create_window(170,630, anchor=NW ,tag = 'all',window=kaptiL6)
            kaptiL7_window = canvas.create_window(170,650, anchor=NW ,tag = 'all',window=kaptiL7)
            kaptiL8_window = canvas.create_window(170,670, anchor=NW ,tag = 'all',window=kaptiL8)
            kaptiL9_window = canvas.create_window(170,690, anchor=NW ,tag = 'all',window=kaptiL9)
            kaptiL10_window = canvas.create_window(170,710, anchor=NW ,tag = 'all',window=kaptiL10)
            L1_window = canvas.create_window(250, 530, anchor=NW,tag = 'all', window=L1)
            L2_window = canvas.create_window(250, 550, anchor=NW,tag = 'all', window=L2)
            L3_window = canvas.create_window(250, 570, anchor=NW,tag = 'all', window=L3)
            L4_window = canvas.create_window(250, 590, anchor=NW,tag = 'all', window=L4)
            L5_window = canvas.create_window(250, 610, anchor=NW,tag = 'all', window=L5)
            L6_window = canvas.create_window(250, 630, anchor=NW,tag = 'all', window=L6)
            L7_window = canvas.create_window(250, 650, anchor=NW,tag = 'all', window=L7)
            L8_window = canvas.create_window(250, 670, anchor=NW,tag = 'all', window=L8)
            L9_window = canvas.create_window(250, 690, anchor=NW,tag = 'all', window=L9)
            L10_window = canvas.create_window(250, 710, anchor=NW,tag = 'all', window=L10)
            H1_window = canvas.create_window(405, 530, anchor=NW,tag = 'all', window=H1)
            H2_window = canvas.create_window(405, 550, anchor=NW,tag = 'all', window=H2)
            H3_window = canvas.create_window(405, 570, anchor=NW,tag = 'all', window=H3)
            H4_window = canvas.create_window(405, 590, anchor=NW,tag = 'all', window=H4)
            H5_window = canvas.create_window(405, 610, anchor=NW,tag = 'all', window=H5)
            H6_window = canvas.create_window(405, 630, anchor=NW,tag = 'all', window=H6)
            H7_window = canvas.create_window(405, 650, anchor=NW,tag = 'all', window=H7)
            H8_window = canvas.create_window(405, 670, anchor=NW,tag = 'all', window=H8)
            H9_window = canvas.create_window(405, 690, anchor=NW,tag = 'all', window=H9)
            H10_window = canvas.create_window(405, 710, anchor=NW,tag = 'all', window=H10)
            g1_window = canvas.create_window(560, 527, anchor=NW,tag = 'all', window=g1)
            g2_window = canvas.create_window(560, 547, anchor=NW,tag = 'all', window=g2)
            g3_window = canvas.create_window(560, 567, anchor=NW,tag = 'all', window=g3)
            g4_window = canvas.create_window(560, 587, anchor=NW,tag = 'all', window=g4)
            g5_window = canvas.create_window(560, 607, anchor=NW,tag = 'all', window=g5)
            g6_window = canvas.create_window(560, 627, anchor=NW,tag = 'all', window=g6)
            g7_window = canvas.create_window(560, 647, anchor=NW,tag = 'all', window=g7)
            g8_window = canvas.create_window(560, 667, anchor=NW,tag = 'all', window=g8)
            g9_window = canvas.create_window(560, 687, anchor=NW,tag = 'all', window=g9)
            g10_window = canvas.create_window(560, 707, anchor=NW,tag = 'all', window=g10)
            KAPGRAPH_window = canvas.create_window(250, 730, anchor=NW,tag = 'all', window=KAPGRAPH)
            kaplineamountoption_window = canvas.create_window(170, 850, anchor=NW,tag = 'all',  window=kaplineamountoption)
            kapL_window = canvas.create_window(30,820, anchor=NW ,tag = 'all',window=kapL)
            allpat_window = canvas.create_window(435, 852, anchor=NW,tag = 'all', window=allpat)
            KAPGRAPHALL_window = canvas.create_window(250, 852, anchor=NW,tag = 'all', window=KAPGRAPHALL)
            option1_window = canvas.create_window(250, 760, anchor=NW,tag = 'all',  window=option1)
            option2_window = canvas.create_window(330, 760, anchor=NW,tag = 'all',  window=option2)
            pvalB_window = canvas.create_window(420, 762, anchor=NW,tag = 'all', window=pvalB)
            pvalE_window = canvas.create_window(480, 764, anchor=NW,tag = 'all', window=pvalE)
            stageL_window = canvas.create_window(30,890, anchor=NW ,tag = 'all',window=stageL)
            StageB_window = canvas.create_window(170, 930, anchor=NW,tag = 'all', window=StageB)
            medianC_window = canvas.create_window(400, 930, anchor=NW,tag = 'all', window=medianC)
            pDataL_window = canvas.create_window(30,990, anchor=NW ,tag = 'all',window=pDataL)
            pDataB_window = canvas.create_window(170, 1030, anchor=NW,tag = 'all', window=pDataB)
##            FRL_window = canvas.create_window(30,1070, anchor=NW ,tag = 'all',window=FRL)
##            FRB_window = canvas.create_window(270, 1103, anchor=NW,tag = 'all', window=FRB)
##            option1FR_window = canvas.create_window(170, 1100, anchor=NW,tag = 'all',  window=option1FR)
##            option2FR_window = canvas.create_window(378, 1100, anchor=NW,tag = 'all',  window=option2FR)
            ButtonHeat_window = canvas.create_window(500, 415, anchor=NW,tag = 'all', window=ButtonHeat)
            heatmapSearch_window = canvas.create_window(500, 395, anchor=NW, window=heatmapSearch)

            global button1_window

            canvas.delete(More_Options)
            canvas.delete(button1_window)

            
            Button1 = Button(text = "Search", command = gene_search_special_tcga, anchor = W)
            Button1.configure(width = 10, activebackground = "#33B5E5",relief = FLAT)
            button1_window = canvas.create_window(700, 10, anchor=NW,tag = 'all', window=Button1)
            

            global Less_Options
            
            Less_Options = Button(text = " Less Options ",fg = 'blue', command = gene_search_tcga_less_options, anchor = W)
            Less_Options.configure(width = 18, activebackground = "#33B5E5",bg=tcgabackground, relief = FLAT)
            Less_Options_window = canvas.create_window(800, 10, anchor=NW,tag = 'all', window=Less_Options)

            global varMC, varFC, MC, FC, sexL, stageLL,S1,S2,S3,S4,varS1,varS2,varS3,varS4,R1,R2,R3,R4,varR1,varR2,varR3,varR4,raceLL,ageLL,AEL,AEH,ageLLTOvarMC, MC_window, FC_window, sexL_window, stageLL_window,S1_window,S2_window,S3_window,S4_window,R1_window,R2_window,R3_window,R4_window,raceLL_window,ageLL_window,AEL_window,AEH_window,ageLLTO_window;
            
            varMC = IntVar()
            varFC = IntVar()
            varS1 = IntVar()
            varS2 = IntVar()
            varS3 = IntVar()
            varS4 = IntVar()
            varR1 = IntVar()
            varR2 = IntVar()
            varR3 = IntVar()
            varR4 = IntVar()

          

            MC = Checkbutton(variable=varMC)
            MC.configure(relief = FLAT,text='Male', bg=tcgabackground, activebackground = tcgabackground)
            MC_window = canvas.create_window(50, 50, anchor=NW,tag = 'all', window=MC)

            FC = Checkbutton(variable=varFC)
            FC.configure(relief = FLAT,text='Female', bg=tcgabackground, activebackground = tcgabackground)
            FC_window = canvas.create_window(130, 50, anchor=NW,tag = 'all', window=FC)

            sexL = Label(text = "Sex:",bg = tcgabackground)
            sexL_window = canvas.create_window(10,50, anchor=NW ,tag = 'all',window=sexL)

            stageLL = Label(text = "Stages:",bg = tcgabackground)
            stageLL_window = canvas.create_window(10,70, anchor=NW ,tag = 'all',window=stageLL)

            S1 = Checkbutton(variable=varS1)
            S1.configure(relief = FLAT,text='Stage 1', bg=tcgabackground, activebackground = tcgabackground)
            S1_window = canvas.create_window(50, 70, anchor=NW,tag = 'all', window=S1)

            S2 = Checkbutton(variable=varS2)
            S2.configure(relief = FLAT,text='Stage 2', bg=tcgabackground, activebackground = tcgabackground)
            S2_window = canvas.create_window(130, 70, anchor=NW,tag = 'all', window=S2)

            S3 = Checkbutton(variable=varS3)
            S3.configure(relief = FLAT,text='Stage 3', bg=tcgabackground, activebackground = tcgabackground)
            S3_window = canvas.create_window(210, 70, anchor=NW,tag = 'all', window=S3)

            S4 = Checkbutton(variable=varS4)
            S4.configure(relief = FLAT,text='Stage 4', bg=tcgabackground, activebackground = tcgabackground)
            S4_window = canvas.create_window(290, 70, anchor=NW,tag = 'all', window=S4)


            raceLL = Label(text = "Races:",bg = tcgabackground)
            raceLL_window = canvas.create_window(10,90, anchor=NW ,tag = 'all',window=raceLL)

            R1 = Checkbutton(variable=varR1)
            R1.configure(relief = FLAT,text='White', bg=tcgabackground, activebackground = tcgabackground)
            R1_window = canvas.create_window(50, 90, anchor=NW,tag = 'all', window=R1)

            R2 = Checkbutton(variable=varR2)
            R2.configure(relief = FLAT,text='Black or African American', bg=tcgabackground, activebackground = tcgabackground)
            R2_window = canvas.create_window(130, 90, anchor=NW,tag = 'all', window=R2)

            R3 = Checkbutton(variable=varR3)
            R3.configure(relief = FLAT,text='Asian', bg=tcgabackground, activebackground = tcgabackground)
            R3_window = canvas.create_window(310, 90, anchor=NW,tag = 'all', window=R3)

            R4 = Checkbutton(variable=varR4)
            R4.configure(relief = FLAT,text='Not Reported', bg=tcgabackground, activebackground = tcgabackground)
            R4_window = canvas.create_window(370, 90, anchor=NW,tag = 'all', window=R4)

            ageLL = Label(text = "Age Range (years):",bg = tcgabackground)
            ageLL_window = canvas.create_window(10,120, anchor=NW ,tag = 'all',window=ageLL)

            ageLLTO = Label(text = "TO",bg = tcgabackground)
            ageLLTO_window = canvas.create_window(155,120, anchor=NW ,tag = 'all',window=ageLLTO)

            AEL = Entry()
            AEL.configure(width = 3, relief = FLAT)
            AEL_window = canvas.create_window(130, 120, anchor=NW,tag = 'all',  window=AEL)

            AEH = Entry()
            AEH .configure(width = 3, relief = FLAT)
            AEH_window = canvas.create_window(180, 120, anchor=NW,tag = 'all',  window=AEH)


            

            

            

                        
            
            

        def TCGATUMOR():
            canvas.delete('all')
            canvas.delete('check')

            global button1_window
            
            Button1 = Button(text = "Search", command = gene_search_tcga, anchor = W)
            Button1.configure(width = 10, activebackground = "#33B5E5",relief = FLAT)
            button1_window = canvas.create_window(700, 10, anchor=NW,tag = 'all', window=Button1)

            global More_Options

            More_Options = Button(text = " More Options ",fg = 'blue', command = gene_search_tcga_more_options, anchor = W)
            More_Options.configure(width = 18, activebackground = "#33B5E5",bg=tcgabackground, relief = FLAT)
            More_Options_window = canvas.create_window(800, 10, anchor=NW,tag = 'all', window=More_Options)
            

            global geneSearch
            geneSearch = Entry()
            geneSearch.insert(0,GENE_TCGA)
            geneSearch .configure(width = 25, relief = FLAT)
            geneSearch_window = canvas.create_window(525, 15, anchor=NW,tag = 'all',  window=geneSearch)

            global var
            var = StringVar(root)
            var.set("Data Source:") # initial value
            
            option = OptionMenu(root,var, "Lung Squamous Cell Carcinoma (TCGA, Provisional)", "Lung Adenocarcinoma (TCGA, Provisional)")
            option.configure(width = 50, activebackground = "#33B5E5", relief = FLAT)
            option_window = canvas.create_window(10, 10, anchor=NW,tag = 'all',  window=option)

            global out
            
            out = StringVar(root)
            out.set("Include Outliers") # initial value

            optionout = OptionMenu(root, out, "Include Outliers", "Exclude Outliers")
            optionout.configure(width = 15, activebackground = "#33B5E5", relief = FLAT)
            optionout_window = canvas.create_window(370, 10, anchor=NW,tag = 'all', window=optionout)

            #Normal Data
            
            global Average
            Average = Entry()
            Average.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Average_window = canvas.create_window(250, 75, anchor=NW,tag = 'all', window=Average)
            
            global STDEV
            STDEV = Entry()
            STDEV.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            STDEV_window = canvas.create_window(250, 95, anchor=NW,tag = 'all', window=STDEV)
            
            global IQR
            IQR = Entry()
            IQR.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            IQR_window = canvas.create_window(250, 115, anchor=NW,tag = 'all', window=IQR)
            
            global Range
            Range = Entry()
            Range.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Range_window = canvas.create_window(250, 135, anchor=NW,tag = 'all', window=Range)
            
            global Minimum
            Minimum = Entry()
            Minimum .configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Minimum_window = canvas.create_window(250, 155, anchor=NW,tag = 'all', window=Minimum)
            
            global FQM
            FQM = Entry()
            FQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            FQM_window = canvas.create_window(250, 175, anchor=NW,tag = 'all', window=FQM)
            
            global Median
            Median = Entry()
            Median.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Median_window = canvas.create_window(250, 195, anchor=NW,tag = 'all', window=Median)
            
            global TQM
            TQM = Entry()
            TQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            TQM_window = canvas.create_window(250, 215, anchor=NW,tag = 'all', window=TQM)

            global Maximum
            Maximum = Entry()
            Maximum.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Maximum_window = canvas.create_window(250, 235, anchor=NW,tag = 'all', window=Maximum)
            #Normal Labels

            global standard
            
            standard = Label(text = "Standard Statistics:",bg = tcgabackground)
            stan_window = canvas.create_window(30,55, anchor=NW ,tag = 'all',window=standard)

            global ave

            ave = Label(text = "Average:",bg = tcgabackground)
            ave_window = canvas.create_window(147,75, anchor=NW ,tag = 'all',window=ave)

            global stdev

            stdev = Label(text = "Standard Deviation:",bg = tcgabackground)
            stdev_window = canvas.create_window(120,95, anchor=NW ,tag = 'all',window=stdev)

            global iqr

            iqr = Label(text = "Interquartile Range:",bg = tcgabackground)
            iqr_window = canvas.create_window(120,115, anchor=NW ,tag = 'all',window=iqr)

            global ran

            ran = Label(text = "Range:",bg = tcgabackground)
            ran_window = canvas.create_window(150,135, anchor=NW ,tag = 'all',window=ran)

            global minimum

            minimum = Label(text = "Minimum:",bg = tcgabackground)
            minimum_window = canvas.create_window(140,155, anchor=NW ,tag = 'all',window=minimum)

            global fqm

            fqm = Label(text = "First Quartile Median:",bg = tcgabackground)
            fqm_window = canvas.create_window(115,175, anchor=NW ,tag = 'all',window=fqm)

            global med

            med = Label(text = "Median:",bg = tcgabackground)
            med_window = canvas.create_window(145,195, anchor=NW ,tag = 'all',window=med)

            global tqm

            tqm = Label(text = "Third Quartile Median:",bg = tcgabackground)
            tqm_window = canvas.create_window(115,215, anchor=NW ,tag = 'all',window=tqm)

            global maximum

            maximum = Label(text = "Maximum:",bg = tcgabackground)
            maximum_window = canvas.create_window(137,235, anchor=NW ,tag = 'all',window=maximum)

            
            
            #correlation
            global corre
            
            corre = Label(text = "Correlations:",bg = tcgabackground)
            corre_window = canvas.create_window(30,275, anchor=NW ,tag = 'all',window=corre)

            global YNREG
            global reglineYN
            
            YNREG = IntVar()

            reglineYN = Checkbutton(variable=YNREG)
            reglineYN.configure(relief = FLAT,text="Graph Regression Line", bg=tcgabackground, activebackground = tcgabackground)
            reglineYN_window = canvas.create_window(250, 360, anchor=NW,tag = 'all', window=reglineYN)

            global rsq

            rsq = Label(text = "R Squared:",bg = tcgabackground)
            rsq_window = canvas.create_window(175,320, anchor=NW ,tag = 'all',window=rsq)

            global tt

            tt = Label(text = "t-test:",bg = tcgabackground)
            tt_window = canvas.create_window(185,340, anchor=NW ,tag = 'all',window=tt)
            
            global geneCSearch
            geneCSearch = Entry()
            geneCSearch.insert(0,GENEC_TCGA)
            geneCSearch .configure(width = 24)
            geneCSearch_window = canvas.create_window(100, 295, anchor=NW, window=geneCSearch)

            global RSQ
            RSQ = Entry()
            RSQ.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            RSQ_window = canvas.create_window(250, 320, anchor=NW,tag = 'all', window=RSQ)

            global TT
            TT = Entry()
            TT.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            TT_window = canvas.create_window(250, 340, anchor=NW,tag = 'all', window=TT)

            global ButtonC

            ButtonC = Button(text = "Correlate to Gene", command = gene_correlate_tcga, anchor = W)
            ButtonC.configure(width = 21, activebackground = "#33B5E5", relief = FLAT)
            buttonC_window = canvas.create_window(250, 290, anchor=NW,tag = 'all', window=ButtonC)
            #heatmaps

            global heatmapSearch
            heatmapSearch = Entry()
            heatmapSearch.insert(0,'Heat Map Genes')
            heatmapSearch .configure(width = 40)
            heatmapSearch_window = canvas.create_window(500, 295, anchor=NW, window=heatmapSearch)

            global ButtonHeat

            ButtonHeat = Button(text = "Start Heat Map", command = tcga_heatmap, anchor = W)
            ButtonHeat.configure(width = 21, activebackground = "#33B5E5", relief = FLAT)
            ButtonHeat_window = canvas.create_window(500, 315, anchor=NW,tag = 'all', window=ButtonHeat)

            #Kaplan Meiers
            global kapti
            global kaptiL
            global kaptiH
            global kaptiL1
            global kaptiL2
            global kaptiL3
            global kaptiL4
            global kaptiL5
            global kaptiL6
            global kaptiL7
            global kaptiL8
            global kaptiL9
            global kaptiL10
            
            
            kapti = Label(text = "Kaplan Meiers:",bg = tcgabackground)
            kapti_window = canvas.create_window(30,390, anchor=NW ,tag = 'all',window=kapti)

            kaptiL = Label(text = "Low Value",bg = tcgabackground)
            kaptiL_window = canvas.create_window(295,410, anchor=NW ,tag = 'all',window=kaptiL)

            kaptiH = Label(text = "High Value",bg = tcgabackground)
            kaptiH_window = canvas.create_window(445,410, anchor=NW ,tag = 'all',window=kaptiH)

            kaptiL1 = Label(text = "First Line:",bg = tcgabackground)
            kaptiL1_window = canvas.create_window(170,430, anchor=NW ,tag = 'all',window=kaptiL1)

            kaptiL2 = Label(text = "Second Line:",bg = tcgabackground)
            kaptiL2_window = canvas.create_window(170,450, anchor=NW ,tag = 'all',window=kaptiL2)

            kaptiL3 = Label(text = "Third Line:",bg = tcgabackground)
            kaptiL3_window = canvas.create_window(170,470, anchor=NW ,tag = 'all',window=kaptiL3)

            kaptiL4 = Label(text = "Fourth Line:",bg = tcgabackground)
            kaptiL4_window = canvas.create_window(170,490, anchor=NW ,tag = 'all',window=kaptiL4)

            kaptiL5 = Label(text = "Fifth Line:",bg = tcgabackground)
            kaptiL5_window = canvas.create_window(170,510, anchor=NW ,tag = 'all',window=kaptiL5)

            kaptiL6 = Label(text = "Sixth Line:",bg = tcgabackground)
            kaptiL6_window = canvas.create_window(170,530, anchor=NW ,tag = 'all',window=kaptiL6)

            kaptiL7 = Label(text = "Seventh Line:",bg = tcgabackground)
            kaptiL7_window = canvas.create_window(170,550, anchor=NW ,tag = 'all',window=kaptiL7)

            kaptiL8 = Label(text = "Eighth Line:",bg = tcgabackground)
            kaptiL8_window = canvas.create_window(170,570, anchor=NW ,tag = 'all',window=kaptiL8)

            kaptiL9 = Label(text = "Ninth Line:",bg = tcgabackground)
            kaptiL9_window = canvas.create_window(170,590, anchor=NW ,tag = 'all',window=kaptiL9)

            kaptiL10 = Label(text = "Tenth Line:",bg = tcgabackground)
            kaptiL10_window = canvas.create_window(170,610, anchor=NW ,tag = 'all',window=kaptiL10)
            

            global L1
            L1 = Entry()
            L1.configure(width = 25, relief = FLAT)
            L1_window = canvas.create_window(250, 430, anchor=NW,tag = 'all', window=L1)

            global L2
            L2 = Entry()
            L2.configure(width = 25, relief = FLAT)
            L2_window = canvas.create_window(250, 450, anchor=NW,tag = 'all', window=L2)

            global L3
            L3 = Entry()
            L3.configure(width = 25, relief = FLAT)
            L3_window = canvas.create_window(250, 470, anchor=NW,tag = 'all', window=L3)

            global L4
            L4 = Entry()
            L4.configure(width = 25, relief = FLAT)
            L4_window = canvas.create_window(250, 490, anchor=NW,tag = 'all', window=L4)

            global L5
            L5 = Entry()
            L5.configure(width = 25, relief = FLAT)
            L5_window = canvas.create_window(250, 510, anchor=NW,tag = 'all', window=L5)

            global L6
            L6 = Entry()
            L6.configure(width = 25, relief = FLAT)
            L6_window = canvas.create_window(250, 530, anchor=NW,tag = 'all', window=L6)

            global L7
            L7 = Entry()
            L7.configure(width = 25, relief = FLAT)
            L7_window = canvas.create_window(250, 550, anchor=NW,tag = 'all', window=L7)

            global L8
            L8 = Entry()
            L8.configure(width = 25, relief = FLAT)
            L8_window = canvas.create_window(250, 570, anchor=NW,tag = 'all', window=L8)

            global L9
            L9 = Entry()
            L9.configure(width = 25, relief = FLAT)
            L9_window = canvas.create_window(250, 590, anchor=NW,tag = 'all', window=L9)

            global L10
            L10 = Entry()
            L10.configure(width = 25, relief = FLAT)
            L10_window = canvas.create_window(250, 610, anchor=NW,tag = 'all', window=L10)

         
            global H1
            H1 = Entry()
            H1.configure(width = 25, relief = FLAT)
            H1_window = canvas.create_window(405, 430, anchor=NW,tag = 'all', window=H1)

            global H2
            H2 = Entry()
            H2.configure(width = 25, relief = FLAT)
            H2_window = canvas.create_window(405, 450, anchor=NW,tag = 'all', window=H2)

            global H3
            H3 = Entry()
            H3.configure(width = 25, relief = FLAT)
            H3_window = canvas.create_window(405, 470, anchor=NW,tag = 'all', window=H3)

            global H4
            H4 = Entry()
            H4.configure(width = 25, relief = FLAT)
            H4_window = canvas.create_window(405, 490, anchor=NW,tag = 'all', window=H4)

            global H5
            H5 = Entry()
            H5.configure(width = 25, relief = FLAT)
            H5_window = canvas.create_window(405, 510, anchor=NW,tag = 'all', window=H5)

            global H6
            H6 = Entry()
            H6.configure(width = 25, relief = FLAT)
            H6_window = canvas.create_window(405, 530, anchor=NW,tag = 'all', window=H6)

            global H7
            H7 = Entry()
            H7.configure(width = 25, relief = FLAT)
            H7_window = canvas.create_window(405, 550, anchor=NW,tag = 'all', window=H7)

            global H8
            H8 = Entry()
            H8.configure(width = 25, relief = FLAT)
            H8_window = canvas.create_window(405, 570, anchor=NW,tag = 'all', window=H8)

            global H9
            H9 = Entry()
            H9.configure(width = 25, relief = FLAT)
            H9_window = canvas.create_window(405, 590, anchor=NW,tag = 'all', window=H9)

            global H10
            H10 = Entry()
            H10.configure(width = 25, relief = FLAT)
            H10_window = canvas.create_window(405, 610, anchor=NW,tag = 'all', window=H10)

            global varg1
            global varg2
            global varg3
            global varg4
            global varg5
            global varg6
            global varg7
            global varg8
            global varg9
            global varg10

            global g1
            global g2
            global g3
            global g4
            global g5
            global g6
            global g7
            global g8
            global g9
            global g10
            

            
            varg1 = IntVar()
            varg2= IntVar()
            varg3 = IntVar()
            varg4 = IntVar()
            varg5 = IntVar()
            varg6 = IntVar()
            varg7 = IntVar()
            varg8 = IntVar()
            varg9 = IntVar()
            varg10 = IntVar()
            

            g1 = Checkbutton(variable=varg1)
            g1.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g1_window = canvas.create_window(560, 427, anchor=NW,tag = 'all', window=g1)

            g2 = Checkbutton(variable=varg2)
            g2.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g2_window = canvas.create_window(560, 447, anchor=NW,tag = 'all', window=g2)

            g3 = Checkbutton(variable=varg3)
            g3.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g3_window = canvas.create_window(560, 467, anchor=NW,tag = 'all', window=g3)

            g4 = Checkbutton(variable=varg4)
            g4.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g4_window = canvas.create_window(560, 487, anchor=NW,tag = 'all', window=g4)

            g5 = Checkbutton(variable=varg5)
            g5.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g5_window = canvas.create_window(560, 507, anchor=NW,tag = 'all', window=g5)

            g6 = Checkbutton(variable=varg6)
            g6.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g6_window = canvas.create_window(560, 527, anchor=NW,tag = 'all', window=g6)

            g7 = Checkbutton(variable=varg7)
            g7.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g7_window = canvas.create_window(560, 547, anchor=NW,tag = 'all', window=g7)

            g8 = Checkbutton(variable=varg8)
            g8.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g8_window = canvas.create_window(560, 567, anchor=NW,tag = 'all', window=g8)

            g9 = Checkbutton(variable=varg9)
            g9.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g9_window = canvas.create_window(560, 587, anchor=NW,tag = 'all', window=g9)

            g10 = Checkbutton(variable=varg10)
            g10.configure(relief = FLAT, bg=tcgabackground, activebackground = tcgabackground)
            g10_window = canvas.create_window(560, 607, anchor=NW,tag = 'all', window=g10)

            global KAPGRAPH

            KAPGRAPH = Button(text = "                                Graph Kaplan Meier's", anchor = W, command = tcga_kaplan)
            KAPGRAPH.configure(width = 43, relief = FLAT)
            KAPGRAPH_window = canvas.create_window(250, 630, anchor=NW,tag = 'all', window=KAPGRAPH)

            #best kaps
            global kaplineamount
            kaplineamount = StringVar(root)
            kaplineamount.set("2 Line") # initial value

            global kaplineamountoption
            
            kaplineamountoption = OptionMenu(root,kaplineamount,'2 Line')
            kaplineamountoption.configure(width = 5,height = 1, activebackground = "#33B5E5", relief = FLAT)
            kaplineamountoption_window = canvas.create_window(170, 750, anchor=NW,tag = 'all',  window=kaplineamountoption)

            global kapL

            kapL = Label(text = "Auto Kaplan Meier Graphs:",bg = tcgabackground)
            kapL_window = canvas.create_window(30,720, anchor=NW ,tag = 'all',window=kapL)
            
            global all_pat
            all_pat = IntVar()

            global allpat

            allpat = Checkbutton(variable=all_pat)
            allpat.configure(relief = FLAT,text="Include All Patients", bg=tcgabackground, activebackground = tcgabackground)
            allpat_window = canvas.create_window(435, 752, anchor=NW,tag = 'all', window=allpat)

            global KAPGRAPHALL

            KAPGRAPHALL = Button(text = "Best Kaplan Meiers Graph", anchor = W, command = bestkaps)
            KAPGRAPHALL.configure(width = 25, relief = FLAT)
            KAPGRAPHALL_window = canvas.create_window(250, 752, anchor=NW,tag = 'all', window=KAPGRAPHALL)

            #pvals
            global line1
            line1 = StringVar(root)
            line1.set("Line:") # initial value

            global option1
            
            option1 = OptionMenu(root,line1, "1",'2','3','4','5','6','7','8','9','10')
            option1.configure(width = 5,height = 1, activebackground = "#33B5E5", relief = FLAT)
            option1_window = canvas.create_window(250, 660, anchor=NW,tag = 'all',  window=option1)

            global line2
            line2 = StringVar(root)
            line2.set("Line:") # initial value

            global option2
            
            option2 = OptionMenu(root,line2,"1",'2','3','4','5','6','7','8','9','10')
            option2.configure(width = 5,height = 1, activebackground = "#33B5E5", relief = FLAT)
            option2_window = canvas.create_window(330, 660, anchor=NW,tag = 'all',  window=option2)

            global pvalB

            pvalB = Button(text = "p-value", anchor = W, command = pval_kaps)
            pvalB.configure(width = 6,height = 1, relief = FLAT)
            pvalB_window = canvas.create_window(420, 662, anchor=NW,tag = 'all', window=pvalB)

            global pvalE
            pvalE = Entry()
            pvalE.configure(width = 25, relief = FLAT,state="readonly",readonlybackground="White")
            pvalE_window = canvas.create_window(480, 664, anchor=NW,tag = 'all', window=pvalE)

            #stages
            global stageL
            global StageB
            
            stageL = Label(text = "Stage Graphs:",bg = tcgabackground)
            stageL_window = canvas.create_window(30,790, anchor=NW ,tag = 'all',window=stageL)

            StageB = Button(text = "Graph Stage Data", anchor = W, command = stage_data_graph)
            StageB.configure(width = 30, relief = FLAT)
            StageB_window = canvas.create_window(170, 830, anchor=NW,tag = 'all', window=StageB)

            global varMC1

            varMC1 = IntVar();

            global medianC
            
            medianC = Checkbutton(variable=varMC1)
            medianC.configure(relief = FLAT, text = "Show Medians",bg=tcgabackground, activebackground = tcgabackground)
            medianC_window = canvas.create_window(400, 830, anchor=NW,tag = 'all', window=medianC)

            global pDataL

            pDataL = Label(text = "Expression Graphs:",bg = tcgabackground)
            pDataL_window = canvas.create_window(30,890, anchor=NW ,tag = 'all',window=pDataL)

            global pDataB

            pDataB = Button(text = "Expression Graph", anchor = W, command = patient_data_graph)
            pDataB.configure(width = 30, relief = FLAT)
            pDataB_window = canvas.create_window(170, 930, anchor=NW,tag = 'all', window=pDataB)



##            global FRL,FRB,ORGAN,option1FR,IncExFR,option2FR
##
##            FRL = Label(text = "Full Report:",bg = tcgabackground)
##            FRL_window = canvas.create_window(30,970, anchor=NW ,tag = 'all',window=FRL)
##
##            FRB = Button(text = "Full Report", anchor = W, command = full_report_tcga)
##            FRB.configure(width = 13, relief = FLAT)
##            FRB_window = canvas.create_window(270, 1003, anchor=NW,tag = 'all', window=FRB)
##
##            ORGAN = StringVar(root)
##            ORGAN.set("Organ:") # initial value
##            
##            option1FR = OptionMenu(root,ORGAN, "Lung")
##            option1FR.configure(width = 8,height = 1, activebackground = "#33B5E5", relief = FLAT)
##            option1FR_window = canvas.create_window(170, 1000, anchor=NW,tag = 'all',  window=option1FR)
##
##            IncExFR = StringVar(root)
##            IncExFR.set("Include Outliers") # initial value
##            
##            option2FR = OptionMenu(root,IncExFR,"Include Outliers",'Exclude Outliers')
##            option2FR.configure(width = 15,height = 1, activebackground = "#33B5E5", relief = FLAT)
##            option2FR_window = canvas.create_window(378, 1000, anchor=NW,tag = 'all',  window=option2FR)
##
##            


            canvas.config(bg=tcgabackground)
            TCGAT.configure(bg = "#33B5E5")
            NCBIT.configure(bg = "#edefef")
            NCBIC.configure(bg = "#edefef")
            
        def NCBITUMOR():
            canvas.delete('all')
            canvas.delete('check')


            
            Button1 = Button(text = "Search", command = ncbiT_search, anchor = W)
            Button1.configure(width = 10, activebackground = "#33B5E5", relief = FLAT)
            button1_window = canvas.create_window(700, 10, anchor=NW,tag = 'all', window=Button1)

            global geneSearch
            geneSearch = Entry()
            geneSearch.insert(0,"Search Prob Number")
            geneSearch .configure(width = 25, relief = FLAT)
            geneSearch_window = canvas.create_window(525, 15, anchor=NW,tag = 'all',  window=geneSearch)

            global var
            var = StringVar(root)
            var.set("Data Source:") # initial value
            
            option = OptionMenu(root,var, "Small Cell Lung Cancer", "Lung Adenocarcinoma","Lung Squamous Carcinoma")
            option.configure(width = 50, activebackground = "#33B5E5", relief = FLAT)
            option_window = canvas.create_window(10, 10, anchor=NW,tag = 'all',  window=option)

            global out
            
            out = StringVar(root)
            out.set("Include Outliers") # initial value

            optionout = OptionMenu(root, out, "Include Outliers", "Exclude Outliers")
            optionout.configure(width = 15, activebackground = "#33B5E5", relief = FLAT)
            optionout_window = canvas.create_window(370, 10, anchor=NW,tag = 'all', window=optionout)

            #Normal Data
            
            global Average
            Average = Entry()
            Average.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Average_window = canvas.create_window(250, 75, anchor=NW,tag = 'all', window=Average)
            
            global STDEV
            STDEV = Entry()
            STDEV.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            STDEV_window = canvas.create_window(250, 95, anchor=NW,tag = 'all', window=STDEV)
            
            global IQR
            IQR = Entry()
            IQR.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            IQR_window = canvas.create_window(250, 115, anchor=NW,tag = 'all', window=IQR)
            
            global Range
            Range = Entry()
            Range.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Range_window = canvas.create_window(250, 135, anchor=NW,tag = 'all', window=Range)
            
            global Minimum
            Minimum = Entry()
            Minimum .configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Minimum_window = canvas.create_window(250, 155, anchor=NW,tag = 'all', window=Minimum)
            
            global FQM
            FQM = Entry()
            FQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            FQM_window = canvas.create_window(250, 175, anchor=NW,tag = 'all', window=FQM)
            
            global Median
            Median = Entry()
            Median.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Median_window = canvas.create_window(250, 195, anchor=NW,tag = 'all', window=Median)
            
            global TQM
            TQM = Entry()
            TQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            TQM_window = canvas.create_window(250, 215, anchor=NW,tag = 'all', window=TQM)

            global Maximum
            Maximum = Entry()
            Maximum.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Maximum_window = canvas.create_window(250, 235, anchor=NW,tag = 'all', window=Maximum)
            #Normal Labels

            standard = Label(text = "Standard Statistics:",bg = ncbitbackground)
            stan_window = canvas.create_window(30,55, anchor=NW ,tag = 'all',window=standard)

            ave = Label(text = "Average:",bg = ncbitbackground)
            ave_window = canvas.create_window(147,75, anchor=NW ,tag = 'all',window=ave)

            stdev = Label(text = "Standard Deviation:",bg = ncbitbackground)
            stdev_window = canvas.create_window(120,95, anchor=NW ,tag = 'all',window=stdev)

            iqr = Label(text = "Interquartile Range:",bg = ncbitbackground)
            iqr_window = canvas.create_window(120,115, anchor=NW ,tag = 'all',window=iqr)

            ran = Label(text = "Range:",bg = ncbitbackground)
            ran_window = canvas.create_window(150,135, anchor=NW ,tag = 'all',window=ran)

            minimum = Label(text = "Minimum:",bg = ncbitbackground)
            minimum_window = canvas.create_window(140,155, anchor=NW ,tag = 'all',window=minimum)

            fqm = Label(text = "First Quartile Median:",bg = ncbitbackground)
            fqm_window = canvas.create_window(115,175, anchor=NW ,tag = 'all',window=fqm)

            med = Label(text = "Median:",bg = ncbitbackground)
            med_window = canvas.create_window(145,195, anchor=NW ,tag = 'all',window=med)

            tqm = Label(text = "Third Quartile Median:",bg = ncbitbackground)
            tqm_window = canvas.create_window(115,215, anchor=NW ,tag = 'all',window=tqm)

            maximum = Label(text = "Maximum:",bg = ncbitbackground)
            maximum_window = canvas.create_window(137,235, anchor=NW ,tag = 'all',window=maximum)

            #correlation
            global corre
            
            corre = Label(text = "Correlations:",bg = ncbitbackground)
            corre_window = canvas.create_window(30,275, anchor=NW ,tag = 'all',window=corre)

            global YNREG
            global reglineYN
            
            YNREG = IntVar()

            reglineYN = Checkbutton(variable=YNREG)
            reglineYN.configure(relief = FLAT,text="Graph Regression Line", bg=ncbitbackground, activebackground = ncbitbackground)
            reglineYN_window = canvas.create_window(250, 360, anchor=NW,tag = 'all', window=reglineYN)

            global rsq

            rsq = Label(text = "R Squared:",bg = ncbitbackground)
            rsq_window = canvas.create_window(175,320, anchor=NW ,tag = 'all',window=rsq)

            global tt

            tt = Label(text = "t-test:",bg = ncbitbackground)
            tt_window = canvas.create_window(185,340, anchor=NW ,tag = 'all',window=tt)
            
            global geneCSearch
            geneCSearch = Entry()
            geneCSearch.insert(0,"Correlate to Prob")
            geneCSearch .configure(width = 24)
            geneCSearch_window = canvas.create_window(100, 295, anchor=NW, window=geneCSearch)

            global RSQ
            RSQ = Entry()
            RSQ.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            RSQ_window = canvas.create_window(250, 320, anchor=NW,tag = 'all', window=RSQ)

            global TT
            TT = Entry()
            TT.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            TT_window = canvas.create_window(250, 340, anchor=NW,tag = 'all', window=TT)

            global ButtonC

            ButtonC = Button(text = "Correlate to Gene", command = gene_correlate_ncbiT, anchor = W)
            ButtonC.configure(width = 21, activebackground = "#33B5E5", relief = FLAT)
            buttonC_window = canvas.create_window(250, 290, anchor=NW,tag = 'all', window=ButtonC)
            #heatmaps

            global heatmapSearch
            heatmapSearch = Entry()
            heatmapSearch.insert(0,'Heat Map Probes')
            heatmapSearch .configure(width = 40)
            heatmapSearch_window = canvas.create_window(500, 295, anchor=NW, window=heatmapSearch)

            global ButtonHeat

            ButtonHeat = Button(text = "Start Heat Map", command = ncbi_heatmap, anchor = W)
            ButtonHeat.configure(width = 21, activebackground = "#33B5E5", relief = FLAT)
            ButtonHeat_window = canvas.create_window(500, 315, anchor=NW,tag = 'all', window=ButtonHeat)

            #Expression Graph

            global pDataL

            pDataL = Label(text = "Expression Graphs:",bg = ncbitbackground)
            pDataL_window = canvas.create_window(30,395, anchor=NW ,tag = 'all',window=pDataL)

            global pDataB

            pDataB = Button(text = "Expression Graph", anchor = W, command = patient_data_graph)
            pDataB.configure(width = 30, relief = FLAT)
            pDataB_window = canvas.create_window(170,415, anchor=NW,tag = 'all', window=pDataB)



            

            canvas.config(bg=ncbitbackground)
            TCGAT.configure(bg = "#edefef")
            NCBIT.configure(bg = "#33B5E5")
            NCBIC.configure(bg = "#edefef")
        def NCBICELL():
            canvas.delete('all')
            canvas.delete('check')


            global button1_window
            
            Button1 = Button(text = "Search Cell Line", command = cellline_search_ncbi, anchor = W)
            Button1.configure(width = 22, activebackground = "#33B5E5",relief = FLAT)
            button1_window = canvas.create_window(175, 10, anchor=NW,tag = 'all', window=Button1)
            

            global CellLineE
            CellLineE = Entry()
            CellLineE.insert(0,"Cell Line")
            CellLineE .configure(width = 25, relief = FLAT)
            CellLineE_window = canvas.create_window(10, 12, anchor=NW,tag = 'all',  window=CellLineE)

            global SearchProbE
            SearchProbE = Entry()
            SearchProbE.insert(0,"Prob Number")
            SearchProbE .configure(width = 25, relief = FLAT)
            SearchProbE_window = canvas.create_window(155, 305, anchor=NW,tag = 'all',  window=SearchProbE)

            global button2_window
            
            Button2 = Button(text = "Search Prob", command = search_prob, anchor = W)
            Button2.configure(width = 16, activebackground = "#33B5E5",relief = FLAT)
            button2_window = canvas.create_window(315, 300, anchor=NW,tag = 'all', window=Button2)


            global out
            
            out = StringVar(root)
            out.set("Include Outliers") # initial value

            optionout = OptionMenu(root, out, "Include Outliers", "Exclude Outliers")
            optionout.configure(width = 15, activebackground = "#33B5E5", relief = FLAT)
            optionout_window = canvas.create_window(10, 300, anchor=NW,tag = 'all', window=optionout)

            #Normal Data
            
            global Average
            Average = Entry()
            Average.configure(width =25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Average_window = canvas.create_window(250, 375, anchor=NW,tag = 'all', window=Average)
            
            global STDEV
            STDEV = Entry()
            STDEV.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            STDEV_window = canvas.create_window(250, 395, anchor=NW,tag = 'all', window=STDEV)
            
            global IQR
            IQR = Entry()
            IQR.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            IQR_window = canvas.create_window(250, 415, anchor=NW,tag = 'all', window=IQR)
            
            global Range
            Range = Entry()
            Range.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Range_window = canvas.create_window(250, 435, anchor=NW,tag = 'all', window=Range)
            
            global Minimum
            Minimum = Entry()
            Minimum .configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Minimum_window = canvas.create_window(250, 455, anchor=NW,tag = 'all', window=Minimum)
            
            global FQM
            FQM = Entry()
            FQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            FQM_window = canvas.create_window(250, 475, anchor=NW,tag = 'all', window=FQM)
            
            global Median
            Median = Entry()
            Median.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Median_window = canvas.create_window(250, 495, anchor=NW,tag = 'all', window=Median)
            
            global TQM
            TQM = Entry()
            TQM.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            TQM_window = canvas.create_window(250, 515, anchor=NW,tag = 'all', window=TQM)

            global Maximum
            Maximum = Entry()
            Maximum.configure(width = 25, relief = FLAT, state = 'readonly',readonlybackground="White")
            Maximum_window = canvas.create_window(250, 535, anchor=NW,tag = 'all', window=Maximum)
            #Normal Labels

            global standard
            
            standard = Label(text = "Standard Statistics:",bg = ncbicbackground)
            stan_window = canvas.create_window(30,355, anchor=NW ,tag = 'all',window=standard)

            global ave

            ave = Label(text = "Average:",bg = ncbicbackground)
            ave_window = canvas.create_window(147,375, anchor=NW ,tag = 'all',window=ave)

            global stdev

            stdev = Label(text = "Standard Deviation:",bg = ncbicbackground)
            stdev_window = canvas.create_window(120,395, anchor=NW ,tag = 'all',window=stdev)

            global iqr

            iqr = Label(text = "Interquartile Range:",bg = ncbicbackground)
            iqr_window = canvas.create_window(120,415, anchor=NW ,tag = 'all',window=iqr)

            global ran

            ran = Label(text = "Range:",bg = ncbicbackground)
            ran_window = canvas.create_window(150,435, anchor=NW ,tag = 'all',window=ran)

            global minimum

            minimum = Label(text = "Minimum:",bg = ncbicbackground)
            minimum_window = canvas.create_window(140,455, anchor=NW ,tag = 'all',window=minimum)

            global fqm

            fqm = Label(text = "First Quartile Median:",bg = ncbicbackground)
            fqm_window = canvas.create_window(115,475, anchor=NW ,tag = 'all',window=fqm)

            global med

            med = Label(text = "Median:",bg = ncbicbackground)
            med_window = canvas.create_window(145,495, anchor=NW ,tag = 'all',window=med)

            global tqm

            tqm = Label(text = "Third Quartile Median:",bg = ncbicbackground)
            tqm_window = canvas.create_window(115,515, anchor=NW ,tag = 'all',window=tqm)

            global maximum

            maximum = Label(text = "Maximum:",bg = ncbicbackground)
            maximum_window = canvas.create_window(137,535, anchor=NW ,tag = 'all',window=maximum)

            
            
            

            global pDataL

            pDataL = Label(text = "Expression Graphs:",bg = ncbicbackground)
            pDataL_window = canvas.create_window(30,595, anchor=NW ,tag = 'all',window=pDataL)

            global pDataB

            pDataB = Button(text = "Expression Graph", anchor = W, command = cellline_data_graph)
            pDataB.configure(width = 30, relief = FLAT)
            pDataB_window = canvas.create_window(170,615, anchor=NW,tag = 'all', window=pDataB)
            #heatmaps
            global heatpDataL

            heatpDataL = Label(text = "Heat Map Graph:",bg = ncbicbackground)
            heatpDataL_window = canvas.create_window(30,650, anchor=NW ,tag = 'all',window=heatpDataL)

            global heatmapSearch
            heatmapSearch = Entry()
            heatmapSearch.insert(0,'Heat Map Prob Number')
            heatmapSearch .configure(width = 22)
            heatmapSearch_window = canvas.create_window(170, 673, anchor=NW, window=heatmapSearch)

            global ButtonHeat

            ButtonHeat = Button(text = "Start Heat Map", command = ncbic_heatmap, anchor = W)
            ButtonHeat.configure(width = 21, activebackground = "#33B5E5", relief = FLAT)
            ButtonHeat_window = canvas.create_window(310, 670, anchor=NW,tag = 'all', window=ButtonHeat)
            

            canvas.config(bg=ncbicbackground)
            TCGAT.configure(bg = "#edefef")
            NCBIT.configure(bg = "#edefef")
            NCBIC.configure(bg = "#33B5E5")

        
        #Choosing Database
        choosebuttons=Canvas(bg="GREY",width=200,height=50)

        TCGAT = Button(text = "TCGA Tumor",command = TCGATUMOR ,anchor = W)
        TCGAT.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        TCGAT_window = choosebuttons.create_window(10, 10, anchor=NW, window=TCGAT)

        NCBIT = Button(text = "NCBI Tumor",command = NCBITUMOR, anchor = W)
        NCBIT.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        NCBIT_window = choosebuttons.create_window(100, 10, anchor=NW, window=NCBIT)

        NCBIC = Button(text = "NCBI Cell Line",command = NCBICELL, anchor = W)
        NCBIC.configure(width = 11, activebackground = "#33B5E5", relief = FLAT)
        NCBIC_window = choosebuttons.create_window(190, 10, anchor=NW, window=NCBIC)


        
        choosebuttons.pack(side=TOP,fill=BOTH)




        
        



root = Tk()

root.geometry("1200x900")
app = Window(root)
root.mainloop()
