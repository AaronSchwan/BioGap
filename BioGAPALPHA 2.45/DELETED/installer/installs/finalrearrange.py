import os
import csv
from PIL import Image
import numpy as np

with open('BioGap\cell_lines.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\Celllines\cell_lines.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\cell_lines.csv")

with open('BioGap\lungadeno.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\_NCBI_\lungadeno.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\lungadeno.csv")
with open('BioGap\lungsquamous.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\_NCBI_\lungsquamous.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\lungsquamous.csv")
with open('BioGap\smallcell.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\_NCBI_\smallcell.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\smallcell.csv")


with open('BioGap\PDATALUAD.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\Data Lung Adenocarcinoma (TCGA, Provisional)\PDATALUAD.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\PDATALUAD.csv")
with open('BioGap\PDATALUSC.csv','r') as file:
    data = file.read();
    file.close();
with open('BioGap\Data Sources\Data Lung Squamous Cell Carcinoma (TCGA, Provisional)\PDATALUSC.csv','w') as file:
    file.write(data);
    file.close();
os.remove("BioGap\PDATALUSC.csv")

with open(r'BioGap\_titleimage.txt','r') as file:
    dataT = file.read().split('), (')
    file.close()
x=0;
data = [];
while(x<(len(dataT)-1)):
    p = dataT[x].replace('(','').split(',')
    p1 = (int(p[0]),int(p[1]),int(p[2]))
    data.insert(x,p1)

    
    x+=1;
    
img = Image.new("RGB", (587, 140), "white")

x = 0;

cordy = 0;

while(x<len(data)-1):
    
    cordx=0;

    while(cordx<587):
        try:
            img.putpixel((cordx, cordy),data[x])
            
            cordx+=1;
            x+=1;
        except:
            cordx+=1;
            x+=1;
            pass;
        
        

    cordy+=1;


img.save(r'BioGap\Program Setup\Pictures\Title.png')

#os.remove(r'BioGap\finalrearrange.py')
