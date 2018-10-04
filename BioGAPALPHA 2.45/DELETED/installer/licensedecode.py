import csv
import time

password = [49,112,51,100,102,54 ,55,56,57,49,115,116,50,103,117,108]

starttime =time.time()
##with open("testfilereturn.txt",'r') as file:
##    data = file.read().split(',')
##    file.close();


passwordlist = [password[x%16] for x in range(len(data))]

spot = [(x%5)+1 for x in range(len(data))]

passmod = [(x+1)%passwordlist[x]+1 for x in range(len(data))]

powertop = [(passwordlist[x]*(x+1))+1 for x in range(len(data))]

powerbottom = [((x+1)*(spot[x])+passmod[x])+1 for x in range(len(data))]

power = [int(powertop[x]/powerbottom[x]) for x in range(len(data))]

powerup = [int(power[x]*sum(password)) for x in range(len(data))]


end = [int(int(data[x])/powerup[x]) for x in range(len(data))]

a = [chr(end[x]) for x in range(len(data))];

characters = ''.join(a)




with open("testfilereturnfinal.txt",'w') as file:
    file.write(str(characters))
    file.close();

print(time.time()-starttime)
