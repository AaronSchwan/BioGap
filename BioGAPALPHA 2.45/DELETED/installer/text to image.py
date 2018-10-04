from PIL import Image
import numpy as np
import csv

with open('titleimage.txt','r') as file:
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


img.save('Title.png')

