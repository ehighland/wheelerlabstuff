def x():

    CHB_list = []
    JPT_list = []	

    CHB = open('CHB_Expression.txt', 'r')
    line = CHB.readline().strip('\n')	
    x = line.split('\t')
    for i in x[1:]:
	i = str(i)
	i += '-CHB'
        CHB_list.append(i)
    	

    JPT = open('JPT_Expression.txt', 'r')
    line = JPT.readline().strip('\n')
    x = line.split('\t')
    for j in x[1:]:
	j = str(j)
        j += '-JPT'
        JPT_list.append(j)
    	    
    newFile = open('AsianPop_Combined.txt', 'w')
    newFile.write('id\t' + str(CHB_list).strip('[').strip(']').replace("'", '').replace(',','\t') + '\t' \
    + str(JPT_list).strip('[').strip(']').replace("'",'').replace(',','\t'))
    newFile.write('\n')

    CHB_lines = CHB.readlines()
    JPT_lines = JPT.readlines()

    for i in CHB_lines:
        var = i[:13]
        if var not in JPT_lines:
            CHB_lines.remove(i)
        else:
            pass

    for j in JPT_lines:
        var = j[:13]
        if var not in CHB_lines:
            JPT_lines.remove(j)
        else:
	    pass 

    #print(totalLen)
    
    if len(CHB_lines) > len(JPT_lines):
        totalLen = len(CHB_lines)
    else:
	totalLen = len(JPT_lines)

    for i in range(1,totalLen):
        try:
    	    CHBVar = CHB_lines[i].strip('\n')
	except IndexError:
	    pass        
        try: 
            JPTVar = JPT_lines[i].strip('\n')
	except IndexError:
	    pass	
        line = CHBVar + '\t' +  JPTVar + '\n'	
        newFile.write(line)

    newFile.close()
    CHB.close()
    JPT.close()

x()
