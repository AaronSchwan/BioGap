
import csv
import os
from tkinter import *



def kaplan_graph_change():

    t = Toplevel()
    t.wm_title("Kaplan Meier Graph")
    t.wm_geometry("460x500")
    t.resizable(width=False, height=False)
    kap=Canvas(bg="#e6ffff",width=200,height=200)
    kap.pack(t,expand=True)
     

