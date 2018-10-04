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

import correlatetcga
import Kaplan_Meier_Graphing_1_10_lines
import Kaplan_Meier_Graphing_expected
import kaplan_change
import bestkaplansalltwoline
import bestkaplanstwoline
import stage
import tcga_genesearch
import heatmapforexpression



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
