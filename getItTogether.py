import sys
import os

def addPop(myFile):
    myPop = myFile[0:3]
    fileName = myPop + '_temp.txt'
    newFile = open(fileName, 'w')

    f = open(myFile, 'r+')
    with open(myFile) as f:     	
	for line in f:
            if line.startswith('gene\t'):
                next(f)
            x = line.split('\t')
	    y = str(x[-1].strip('\n'))
	    x = x[0:9]
            x.append(y)
	    x.append(myPop)
            clean = str(x).strip(']').strip('[').replace(',','\t').replace("'",'').strip()
	    newFile.write(clean)
	    newFile.write('\n')

    f.close()	
    
    return fileName


def mergeFiles(a, b, file):
    tempA = addPop(a)
    tempB = addPop(b)
    
    file = file +  '.txt'	

    f = open(file, 'w')	
    
    f.write('gene\tpve50\tpge50\tn_gamma50\tpve025\tpge025\tn_gamma025\tpve975\tpge975\tn_gamma975\tPopulation\n')	

    with open(tempA) as tempA:
	next(tempA)
	for line in tempA:
            f.write(line)

    with open(tempB) as tempB: 
	next(tempB)
	for line in tempB:
	    f.write(line)
	
    #os.system('rm *_temp.txt')

    tempA.close()
    tempB.close()
    f.close()	
    
def getItTogether():
    fileName = sys.argv[1]
    pop1 = sys.argv[2]
    pop2 = sys.argv[3]

    mergeFiles(pop1, pop2, fileName)

getItTogether()    		
    		 
