import funcs
import sys
import math
import matplotlib.pyplot as plt
import random
import numpy
import os
import scipy.spatial	


#----------------------------------
THRESHOLD=0.001
MAXITERATIONS=30
EMMAXITERATIONS=30
EMTHRESHOLD=0.1
TRAININGDATALIMIT=0.75
CURRENTNOOFCLUSTERS=8
#----------------------------------

def belongsToClusture(dataPoint,means):
	distances=[]
	for i in range(len(means)):
		distances.append(scipy.spatial.distance.euclidean(dataPoint,means[i]))
	return distances.index(min(distances))

def emCheckConvergence(oldl,newl,interationNo):
	for clsi in range(len(oldl)):
		if(math.fabs(oldl[clsi]-newl[clsi])<EMTHRESHOLD):
			print ""
			print "NOTE:DIFFERENCE CROSSED THRESHOLD FOR CONVERGENCE -- EM method"
			return True
		if(interationNo>EMMAXITERATIONS):
			print ""
			print "NOTE:MAX ITERATIONS REACHED FOR CONVERGENCE -- EM method"
			return True
	return False

#to update
def checkConvergence(old,new,interationNo):
	if(interationNo>MAXITERATIONS):
		print "NOTE:MAX ITERATIONS REACHED FOR CONVERGENCE"
		return True
	return False

#problem of two means assigning same value -- solved
def kMeansClusture(k,data,datainclusters):

	curmeans=numpy.zeros((CURRENTNOOFCLUSTERS,VECTORDIMENSION,1),dtype=numpy.float)
	initialMeanIndex=[]

	i=0
	while(i<k):
		randindex=random.randrange(int(len(data)*TRAININGDATALIMIT))
		#checking for mean index repetitions
		try:
			initialMeanIndex.index(randindex)
			continue
		except:
			curmeans[i]=data[randindex]
			initialMeanIndex.append(randindex)
			datainclusters.append([])
			i+=1

	interationNo=0

	while(True):
		interationNo+=1

		for i in range(k):
			datainclusters[i]=[]

		for i in range(int(len(data)*TRAININGDATALIMIT)):
			datainclusters[belongsToClusture(data[i],curmeans)].append(data[i])

		#print float(len(datainclusters[0]))/(len(data)*0.75)+float(len(datainclusters[1]))/(len(data)*0.75)+float(len(datainclusters[2]))/(len(data)*0.75)+float(len(datainclusters[3]))/(len(data)*0.75)+float(len(datainclusters[4]))/(len(data)*0.75)

		newmeans=numpy.zeros((CURRENTNOOFCLUSTERS,VECTORDIMENSION,1),dtype=numpy.float)
		for i in range(len(datainclusters)):
			for j in range(0,len(datainclusters[i])):
				newmeans[i]=datainclusters[i][j]+newmeans[i];
			
			#no data point classified in cluster -- hopefully never executes
			if(len(datainclusters[i])==0):
				print "error iteration no :",interationNo
				print "error cluster index :",i
				print "error cluster data array :",datainclusters[i]
				print "no data point in a clusture"
				print "TERMINATING"
				exit()
				
			newmeans[i]=newmeans[i]/len(datainclusters[i])

		if(checkConvergence(curmeans,newmeans,interationNo)):
			return newmeans[:]

		curmeans=newmeans[:]

def calculateL(data,mean,covariance,piik):
	tempL=[]
	
	covINV=[]
	covINVdet=[]

	for clsi in range(len(data)):
		covINV.append([])
		covINVdet.append([])
		for clusi in range(len(mean[clsi])):
			covINV[clsi].append(numpy.linalg.inv(covariance[clsi][clusi]))
			covINVdet[clsi].append(numpy.linalg.det(covINV[clsi][clusi]))


	for clsi in range(len(data)):
		tempL.append(0)
		for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):

			#bs='\b'*1000
			#print bs,
			#print "\bCALCULATE-L CLASS :",(clsi+1),"PROGRESS :",round(100*float(datai+1)/float(len(data[clsi])),2),"%",

			tempSum=0.0
			for clusi in range(len(mean[clsi])):
				temp=funcs.subM(data[clsi][datai],mean[clsi][clusi])
				temp=numpy.dot(numpy.dot(funcs.transpose(temp),covINV[clsi][clusi]),temp)
				N=covINVdet[clsi][clusi]*math.exp(-0.5*temp[0][0])/math.sqrt(2*3.14)
				tempSum+=piik[clsi][clusi]*N
			tempSum=math.log(tempSum)
			tempL[clsi]+=tempSum
		#print ""
	return tempL

def NormalDist(data,mean,covariance):
	temp=funcs.subM(data,mean)
	temp=numpy.dot(numpy.dot(funcs.transpose(temp),numpy.linalg.inv(covariance)),temp)
	N=math.exp(-0.5*temp[0][0])/((2*3.14)*math.sqrt(numpy.linalg.det(covariance)))
	return N

#------------------------------------------------------------------------------------------MAIN-CODE

data=[]
datainclusters=[]

minx=0
miny=0
maxx=0
maxy=0

#reading files and storing data in data[]
DATAVALUESEPARATOR=" "
noofdir=len(sys.argv)-1
PRINTSTEPSFLAG=int(sys.argv[len(sys.argv)-1])
PRINTGRAPHSFLAG=0
VECTORDIMENSION=0

nooffiles=len(sys.argv)-1

PRINTSTEPSFLAG=int(sys.argv[len(sys.argv)-1])

#for each file
for filei in range(1,nooffiles):
	with open(sys.argv[filei]) as curFile:
		data.append([])
		datainclusters.append([])

		content = curFile.read().splitlines()
		lineno=len(content);
		#for each line
		for linei in range(lineno):
			value=content[linei].split(DATAVALUESEPARATOR);
			data[filei-1].append([])
			
			#for each data value in data point
			for valuei in range(len(value)):
				try:
					data[filei-1][linei].append([float(value[valuei])])
					if(filei==1 and linei==0):
							VECTORDIMENSION+=1
				except ValueError:
					continue
			#endoffor
			if(float(value[0])<minx):
				minx=float(value[0])
			if(float(value[0])>maxx):
				maxx=float(value[0])
			if(float(value[1])<miny):
				miny=float(value[1])
			if(float(value[1])>maxy):
				maxy=float(value[1])
		#endoffor
	#endofwith
#endoffor

data=numpy.array(data)
mean=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS,VECTORDIMENSION,1),dtype=numpy.float)
covariance=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS,VECTORDIMENSION,VECTORDIMENSION),dtype=numpy.float)
pik=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS),dtype=numpy.float)
oldL=numpy.zeros((noofdir-1),dtype=numpy.float)
newL=numpy.zeros((noofdir-1),dtype=numpy.float)

#no of classes
CLSNO=len(data)

for i in range(CLSNO):
	mean[i]=kMeansClusture(CURRENTNOOFCLUSTERS,data[i],datainclusters[i])

#print mean
print "CALCULTING INTITAL PIK"
for clsi in range(CLSNO):
	for clusi in range(CURRENTNOOFCLUSTERS):
		noofdata=len(datainclusters[clsi][clusi])
		for datai in range(noofdata):
			tempSub=datainclusters[clsi][clusi][datai]-mean[clsi][clusi]
			covariance[clsi][clusi]=numpy.dot(tempSub,numpy.transpose(tempSub))+covariance[clsi][clusi]
			
		covariance[clsi][clusi]=covariance[clsi][clusi]/noofdata
		pik[clsi][clusi]=float(len(datainclusters[clsi][clusi]))/float(len(data[clsi])*TRAININGDATALIMIT)

'''
print "CALCULATING L"
newL=calculateL(data,mean,covariance,pik)
'''

interationNo=0

clusterColors=[["#ffcdd2","#e57373","#f44336","#b71c1c","#f06292","#e91e63","#c2185b","#880e4f"],["#b2dfdb","#4db6ac","#009688","#00796b","#004d40","#a5d6a7","#66bb6a","#1b5e20"],["#303f9f","#1976d2","#0288d1","#0097a7","#64b5f6"]]

gammak=[]

for clsi in range(CLSNO):
	gammak.append([])
	for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
		gammak[clsi].append([])
		for clusi in range(CURRENTNOOFCLUSTERS):
			gammak[clsi][datai].append([0])

gammak=numpy.array(gammak)

while(True):
	interationNo+=1

	print "ITERATION :",interationNo

	#E STEP
	#oldL=newL[:]

	covINV=[]
	covINVdet=[]

	for clsi in range(len(data)):
		covINV.append([])
		covINVdet.append([])
		for clusi in range(CURRENTNOOFCLUSTERS):
			covINV[clsi].append(numpy.linalg.inv(covariance[clsi][clusi]))
			covINVdet[clsi].append(numpy.linalg.det(covINV[clsi][clusi]))

	print "E-STEP CALCULATING GAMMA-K"
	
	for clsi in range(CLSNO):
		for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):

			bs='\b'*1000
			print bs,
			print "\bGAMMA-K CLASS :",(clsi+1),"PROGRESS :",round(100*float(datai+1)/float(len(data[clsi])),2),"%",
			
			totalgamma=0;
			for clusi in range(CURRENTNOOFCLUSTERS):
				temp=data[clsi][datai]-mean[clsi][clusi]
				temp=numpy.dot(numpy.dot(numpy.transpose(temp),covINV[clsi][clusi]),temp)
				gammak[clsi][datai][clusi]=[pik[clsi][clusi]*covINVdet[clsi][clusi]*math.exp(-0.5*temp[0][0])/math.sqrt(2*3.14)]
				totalgamma+=gammak[clsi][datai][clusi][0]
			
			if(len(gammak[clsi][datai])!=CURRENTNOOFCLUSTERS):
				print "error"
			
			gammak[clsi][datai]=gammak[clsi][datai]/totalgamma

	#M STEP
	NK=[]

	for clsi in range(CLSNO):
		NK.append([])
		for clusi in range(CURRENTNOOFCLUSTERS):
			NK[clsi].append([])
			NK[clsi][clusi]=0
			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
				NK[clsi][clusi]+=gammak[clsi][datai][clusi][0]

	print "M-STEP MAXIMIZING DATA"

	for clsi in range(CLSNO):
		for clusi in range(CURRENTNOOFCLUSTERS):
			pik[clsi][clusi]=NK[clsi][clusi]/len(data[clsi])
			mean[clsi][clusi]=numpy.zeros((len(data[0][0]),1),dtype=numpy.float)

			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
#				bs='\b'*1000
#				print bs,
#				print "\bMAXIMIZING MEAN CLASS :",(clsi+1),"PROGRESS :",round(100*float(datai+1)/float(len(data[clsi])),2),"%",
				
				mean[clsi][clusi]=mean[clsi][clusi]+numpy.array(data[clsi][datai])*(gammak[clsi][datai][clusi][0]/NK[clsi][clusi])

#temp			mean[clsi][clusi]=funcs.divByConstM(mean[clsi][clusi],NK[clsi][clusi])

			covariance[clsi][clusi]=numpy.zeros((len(data[0][0]),len(data[0][0])), dtype=numpy.float)

			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
#				bs='\b'*1000
#				print bs,
#				print "\bMAXIMIZING COVARIANCE CLASS :",(clsi+1),"PROGRESS :",round(100*float(datai+1)/float(len(data[clsi])),2),"%",
				tempSub=data[clsi][datai]-mean[clsi][clusi]
				covariance[clsi][clusi]=covariance[clsi][clusi]+numpy.dot(tempSub,numpy.transpose(tempSub))*(gammak[clsi][datai][clusi][0]/NK[clsi][clusi])
#temp			covariance[clsi][clusi]=funcs.divByConstM(covariance[clsi][clusi],NK[clsi][clusi])
	print covariance
	#print "CALCULATING L"
	#newL=calculateL(data,mean,covariance,pik)

#as data is not going to converge -- REAL IMAGE DATA
	if(emCheckConvergence([1,1,1],[100,100,100],interationNo)):
	 	break

	# datainclusters=[]

	# if(PRINTSTEPSFLAG==1):
	# 	fig=plt.figure()
	# 	subplot=fig.add_subplot(111)
	# 	subplot.set_title(str(interationNo))

	# 	for clsi in range(CLSNO):
	# 		datainclusters.append([])
	# 		for clusi in range(CURRENTNOOFCLUSTERS):
	# 			datainclusters[clsi].append([])
	# 		for datai in range(len(data[clsi])):
	# 			clas=funcs.classify(data[clsi][datai],mean[clsi],covariance[clsi],pik[clsi])
	# 			subplot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=clusterColors[clsi][clas],marker="o")
	# 			datainclusters[clsi][clas].append(data[clsi][datai])

	# 	plt.axis([minx,maxx,miny,maxy])
	# 	fig.savefig(str(interationNo)+".png")


fig=plt.figure()
subplot=fig.add_subplot(111)

color=["#F49292","#A9F5AF","#A1A4FF","#F3F3AF"]
color1=["#E21818","#17E81F","#252BDF","#D2C81D"]
xinc=(maxx-minx)/50
yinc=(maxy-miny)/50

#abcd
gx=[]
for clsi in range(CLSNO):
	gx.append(0)

vecx=[]
vecy=[]
vecz=[]

mingx=0
maxgx=0

yi=0
y=miny
while(y<maxy):
	
	vecy.append(y)
	vecz.append([])

	bs='\b'*1000
	print bs,
	print "\bPLOTTING GRAPH :",round((y-miny)*100/(maxy-miny),2),"%",

	x=minx
	while(x<maxx):
		
		for clsi in range(CLSNO):
			temp=0.0
			for clusi in range(len(mean[clsi])):
				temp+=pik[clsi][clusi]*NormalDist([[x],[y]],mean[clsi][clusi],covariance[clsi][clusi])

			totalN=0
			for xx in range(CLSNO):
				totalN+=len(data[xx])

			temp*=float(len(data[clsi]))/float(totalN)

			#hopefully this does not execute
			if(temp>1):
				print "error probability :",temp
				print "probability more that 1 error"
				print "TERMINATING"
				exit()

			temp=math.log(temp)
			gx[clsi]=temp
		
		mgx=max(gx)

		vecz[yi].append(mgx)

		if(mingx>mgx):
			mingx=mgx
		if(maxgx<mgx):
			maxgx=mgx

		if(y==miny):
			vecx.append(x)
			if(x==minx):
				mingx=mgx
				maxgx=mingx

		clas=gx.index(mgx)
		subplot.plot(x,y,color=color[clas],marker="o",markeredgecolor=color[clas])
		x+=xinc
	y+=yinc
	yi+=1

print ""

plt.axis([minx,maxx,miny,maxy])

contourFig=plt.figure()
contourSubPlot=contourFig.add_subplot(111)
#CS = contourSubPlot.contour(vecx,vecy,vecz,zorder=2,level=funcs.generateGP(mingx,maxgx))
CS = contourSubPlot.contour(vecx,vecy,vecz,zorder=2,level=funcs.generateGP(mingx,maxgx))
contourSubPlot.clabel(CS, inline=1, fontsize=10)	

for clsi in range(CLSNO):
	for datai in range(int(TRAININGDATALIMIT*len(data[clsi]))):
		contourSubPlot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color[clsi],marker="o",markeredgecolor=color[clsi],zorder=1)
		subplot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color1[clsi],marker="o")

plt.axis([minx,maxx,miny,maxy])

print "SAVING CONTOUR GRAPH"
contourFig.savefig("contour.png")
print "SAVING BOUNDARIES GRAPH"
fig.savefig("boundaries.png")