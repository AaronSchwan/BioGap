from PIL import Image
import csv

im = Image.open('Title.png')

pixels = list(im.convert ("RGB").getdata())
width, height = im.size
print(width,height)
##pixels = [pixels[i * width:(i + 1) * width] for i in range(height)]

with open('titleimage.txt','w') as file:
    file.write(str(pixels).replace('[','').replace(']',''))
    file.close()
