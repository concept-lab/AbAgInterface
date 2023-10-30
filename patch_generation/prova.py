
from utils import *

# pqrFile = open("../structures/exposed/7mhy/7mhy_antibody_MN.pqr",'r')
exposedFile = open("../structures/exposed/7mhy/7mhy_antibody_MN_exposed.txt",'r')
triangFileName = "../structures/exposed/7mhy/7mhy_antibody_MN.off"

anumb_exposed = set()
anumb_vertExposed = set()

for line in exposedFile.readlines():
    anumb_exposed.add(int(line))

triangLines = loadTriang(triangFileName)
for line in triangLines:
    anumb_vertExposed.add(int(line.split()[6]))

print("exposed set = ",len(anumb_exposed))
print("closer to vert set = ",len(anumb_vertExposed))
print("intersection = ", len(anumb_vertExposed & anumb_exposed))

outFile=open("provaVertMap.txt",'w')

for line in triangLines:
    indAtomClose = int(line.split()[6]) #referred to index in antigen.xyzr
    flag = 0
    if indAtomClose in anumb_exposed:
        flag = 3
    outFile.write(str(flag)+'\n')
outFile.close()