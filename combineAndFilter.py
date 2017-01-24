import sys
import os

def filter(myFile):
    myLines = []

    with open(myFile) as f:
        next(f)
        for line in f:
            myLine = line.split('\t')
            myLines.append(str(myLine[0]) + '\t' + str(myLine[2]) + '\t' + str(myLine[4]) + '\t' + str(myLine[6]) + '\t' + str(myLine[8]) + '\t' + str(myLine[10]) + '\t' + str(myLine[12]) + '\t' + str(myLine[14]) + '\t' + str(myLine[16]) + '\t' + str(myLine[18]))
    f.close()

    return myLines

def combine():
    myFile = sys.argv[1]
    filtered = filter(myFile)

    newFile = open('combined2.txt','a+')
    
    if os.path.getsize('combined2.txt') > 0:
        for i in range(0,len(filtered)):
            newFile.write(str(filtered[i].strip()) + '\n')
    else:
        newFile.write('gene\tpve50\tpge50\tn_gamma50\tpve025\tpge025\tn_gamma025\tpve975\tpge975\tn_gamma975\n')
        for i in range(0,len(filtered)):
            newFile.write(str(filtered[i].strip()) + '\n')

    newFile.close()

combine()

