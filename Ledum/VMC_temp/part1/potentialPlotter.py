import matplotlib.pyplot 	as plt
import numpy 				as np


R = []
V = []


#fileName = '../build-part1-Desktop_Qt_5_7_0_clang_64bit-Release/VMCData.dat'
fileNameList = ['../build-part1-Desktop_Qt_5_5_0_GCC_64bit-Release/VMCData05-15.dat',
                '../build-part1-Desktop_Qt_5_5_0_GCC_64bit-Release/VMCData15-55.dat',
                '../build-part1-Desktop_Qt_5_5_0_GCC_64bit-Release/VMCData25-35.dat',
                '../build-part1-Desktop_Qt_5_5_0_GCC_64bit-Release/VMCData35-45.dat']

for i in xrange(len(fileNameList)) :
    fileName = fileNameList[i]
    with open(fileName, 'r') as inFile :
        j=0
        for line in inFile :
            if (i == 0) or (i != 0 and j > 5):
		        line = line.split()
		        R.append(float(line[0]))
		        V.append(float(line[1]))
            j += 1


plt.figure()
plt.plot(R,V,'r.-')
plt.xlabel('R')
plt.ylabel('V')
plt.title('Potential as a function of intermolecular distance for H2 molecule.')
plt.show()


with open('VMC_H2.dat', 'w') as outFile :
    for i in xrange(len(R)) :
        outFile.write("%.15f %.15f\n" % (R[i], V[i]))

