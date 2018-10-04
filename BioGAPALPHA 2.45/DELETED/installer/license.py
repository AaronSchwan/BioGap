import csv


password = [49,112,51,100,102,54 ,55,56,57,49,115,116,50,103,117,108]


with open("testfile.txt",'r') as file:
    data = list(file.read())
    file.close();

asciilist = list(map(ord,data))

passwordlist = [password[x%16] for x in range(len(data))]

spot = [(x%5)+1 for x in range(len(data))]

passmod = [(x+1)%passwordlist[x]+1 for x in range(len(data))]

powertop = [(passwordlist[x]*(x+1))+1 for x in range(len(data))]

powerbottom = [((x+1)*(spot[x])+passmod[x])+1 for x in range(len(data))]

power = [int(powertop[x]/powerbottom[x]) for x in range(len(data))]

end = [asciilist[x]*power[x]*sum(password) for x in range(len(data))]

with open("licensedecode.py",'r') as file:
    lice = file.read()
    file.close();


    
with open("installer.py",'w') as file:
    file.write('data = ')
    file.write(str(end))
    file.write('\n')
    file.write(str(lice))
    file.close()

with open("text.txt",'w') as file:
    file.write(str(end))
    file.close();

