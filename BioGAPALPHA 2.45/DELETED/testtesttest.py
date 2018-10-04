
from statistics import mean
from statistics import median
import numpy as np
import csv
import os
from tkinter import *
import urllib.request
import urllib.parse
import re
import time
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM147894"
resp = urllib.request.urlopen(url)
respData = resp.read()

listall = respData.split()

listall = str(listall)
start = listall.index("b'nowrap>Scan")
end = listall.index("b'nowrap>Organism</td>")

xx=start;

while(xx<=end):
    listall.pop(xx)
    xx+=1;

print(listall)
