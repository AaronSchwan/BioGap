import csv
import os

mypath = r'installs'

from os import walk

files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    files.extend(filenames)
    break
print(files)
with open('installs.py','w') as mainfile:
    mainfile.write("import os\n");
    mainfile.write("import csv\n");
    mainfile.write('''
import os

newpath = r'BioGap' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

newpath = r'BioGap\Data Sources' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Data Sources\_NCBI_' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Data Sources\Celllines' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Data Sources\Data Lung Adenocarcinoma (TCGA, Provisional)' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Data Sources\Data Lung Squamous Cell Carcinoma (TCGA, Provisional)' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

newpath = r'BioGap\Program Setup' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\Colors' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\KaplanMeier1-10' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\KaplanMeierAuto' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\Pictures' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\Standerds' 
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = r'BioGap\Program Setup\Titles' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

''')
    x=0;
    while(x<len(files)):
        
        with open(os.path.join('installs',files[x]),'r') as file:
            data = file.read();
            file.close();
        
        mainfile.write("\nwith open(os.path.join('BioGap','");
        mainfile.write(files[x]);
        mainfile.write("'),");
        mainfile.write('"w") as file:\n');
        mainfile.write("      file.write(str('''");
        mainfile.write(str(data));
        mainfile.write("'''))");
        mainfile.write('\n      file.close()');
   
        x+=1;
   
    mainfile.close();


mypath = r'installstxt'



files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    files.extend(filenames)
    break
print(files)
with open('installs.py','a') as mainfile:
    mainfile.write("\n");


    x=0;
    while(x<len(files)):
        try:
            with open(os.path.join('installstxt',files[x]),'r') as file:
                data = file.read();
                file.close();
            
            mainfile.write("\nwith open(os.path.join('BioGap','");
            mainfile.write(files[x]);
            mainfile.write("'),");
            mainfile.write('"w") as file:\n');
            mainfile.write("      file.write(str('''");
            mainfile.write(str(data));
            mainfile.write("'''))");
            mainfile.write('\n      file.close()');
        except:
            pass;
        x+=1;

    mainfile.close();

mypath = r'installscsv'



files = []
for (dirpath, dirnames, filenames) in walk(mypath):
    files.extend(filenames)
    break
print(files)
with open('installs.py','a') as mainfile:
    mainfile.write("\n");


    x=0;
    while(x<len(files)):
        try:
            with open(os.path.join('installscsv',files[x]),'r') as file:
                data = file.read();
                file.close();
            
            mainfile.write("\nwith open(os.path.join(r'BioGap','");
            mainfile.write(files[x]);
            mainfile.write("'),");
            mainfile.write('"w") as file:\n');
            mainfile.write("      file.write(str('''");
            mainfile.write(str(data));
            mainfile.write("'''))");
            mainfile.write('\n      file.close()');
        except:
            pass;
        x+=1;

    mainfile.close();
    
with open('installs.py','a') as mainfile:
    mainfile.write("\n");
    
    mainfile.write("\nimport BioGap.finalrearrange")

    mainfile.write('\nos.remove("installs.py")');
    mainfile.close();   
print(len(files))

