import funcs

import sys
import math
import matplotlib.pyplot as plt
import random

THRESHOLD=0.001
MAXITERATIONS=50
EMMAXITERATIONS=50
EMTHRESHOLD=0.1
TRAININGDATALIMIT=0.75
CURRENTNOOFCLUSTERS=10

def belongsToClusture(dataPoint,means):
	distances=[]
	for i in range(len(means)):
		distances.append(math.sqrt(math.pow(dataPoint[0][0]-means[i][0][0],2)+math.pow(dataPoint[1][0]-means[i][1][0],2)))
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
	if(math.sqrt(math.pow(new[0][0][0]-old[0][0][0],2)+math.pow(new[0][1][0]-old[0][1][0],2))<THRESHOLD and math.sqrt(math.pow(new[1][0][0]-old[1][0][0],2)+math.pow(new[1][1][0]-old[1][1][0],2))<THRESHOLD):
		print "NOTE:DIFFERENCE CROSSED THRESHOLD FOR CONVERGENCE"
		return True
	if(interationNo>MAXITERATIONS):
		print "NOTE:MAX ITERATIONS REACHED FOR CONVERGENCE"
		return True
	return False

#problem of two means assigning same value
def kMeansClusture(k,data,minx,maxx,miny,maxy,datainclusters):
	#holding means
	curmeans=[]
	#holding data of different clusters
	#datainclusters=[]

	initialMeanIndex=[]

	i=0
	while(i<k):
		randindex=random.randrange(int(len(data)*TRAININGDATALIMIT))
		#checking for mean index repetitions
		try:
			initialMeanIndex.index(randindex)
			continue
		except:
			curmeans.append([[data[randindex][0][0]],[data[randindex][1][0]]])
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

		newmeans=[]
		for i in range(len(datainclusters)):
			newmeans.append([[0],[0]])
			for j in range(0,len(datainclusters[i])):
				newmeans[i]=funcs.addM([[datainclusters[i][j][0][0]],[datainclusters[i][j][1][0]]],newmeans[i]);
			
			#no data point classified in cluster -- hopefully never executes
			if(len(datainclusters[i])==0):
				print "error iteration no :",interationNo
				print "error cluster index :",i
				print "error cluster data array :",datainclusters[i]
				print "no data point in a clusture"
				print "TERMINATING"
				exit()
				
			newmeans[i]=funcs.divByConstM(newmeans[i],len(datainclusters[i]))

		if(checkConvergence(curmeans,newmeans,interationNo)):
			return newmeans[:]

		curmeans=newmeans[:]

def calculateL(data,mean,covariance,piik):
	tempL=[]
	for clsi in range(len(data)):
		tempL.append(0)
		for datai in range(len(data[clsi])):
			tempSum=0.0
			for clusi in range(len(mean[clsi])):
				temp=funcs.subM(data[clsi][datai],mean[clsi][clusi])
				temp=funcs.mulM(funcs.mulM(funcs.transpose(temp),funcs.inverse(covariance[clsi][clusi])),temp)
				N=funcs.det(funcs.inverse(covariance[clsi][clusi]))*math.exp(-0.5*temp[0][0])/math.sqrt(2*3.14)
				tempSum+=piik[clsi][clusi]*N
			tempSum=math.log(tempSum)
			tempL[clsi]+=tempSum
	return tempL

def NormalDist(data,mean,covariance):
	temp=funcs.subM(data,mean)
	temp=funcs.mulM(funcs.mulM(funcs.transpose(temp),funcs.inverse(covariance)),temp)
	N=math.exp(-0.5*temp[0][0])/math.sqrt(2*3.14)*math.sqrt(funcs.det(covariance))
	return N

#universal matrices
data=[]
mean=[]
datainclusters=[]
covariance=[]
pik=[]

oldL=[]
newL=[]

minx=0
miny=0
maxx=0
maxy=0
#reading files and storing data in data[]
DATAVALUESEPARATOR=" "
nooffiles=len(sys.argv)-1

PRINTSTEPSFLAG=int(sys.argv[len(sys.argv)-1])

#for each file
for filei in range(1,nooffiles):
	with open(sys.argv[filei]) as curFile:
		data.append([])
		mean.append([])
		datainclusters.append([])
		covariance.append([])
		pik.append([])
		newL.append(0.0)
		oldL.append(0.0)

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
				except ValueError:
					continue
			#endoffor
			#checking for ranges of data
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

#no of classes
clsno=len(data)

for i in range(clsno):
	mean[i]=kMeansClusture(CURRENTNOOFCLUSTERS,data[i],minx,maxx,miny,maxy,datainclusters[i])

#print mean

for clsi in range(clsno):
	for clusi in range(CURRENTNOOFCLUSTERS):
		covariance[clsi].append([])
		covariance[clsi][clusi]=[[0,0],[0,0]]
		for datai in range(len(datainclusters[clsi][clusi])):
			covariance[clsi][clusi]=funcs.addM(funcs.mulM(funcs.subM(datainclusters[clsi][clusi][datai],mean[clsi][clusi]),funcs.transpose(funcs.subM(datainclusters[clsi][clusi][datai],mean[clsi][clusi]))),covariance[clsi][clusi])
		covariance[clsi][clusi]=funcs.divByConstM(covariance[clsi][clusi],len(datainclusters[clsi][clusi]))
		pik[clsi].append([])
		pik[clsi][clusi]=float(len(datainclusters[clsi][clusi]))/float(len(data[clsi]))

print "PIK CALCULATED"
newL=calculateL(data,mean,covariance,pik)
print "L CALCULATED"

interationNo=0

clusterColors=[["#d32f2f","#c2185b","#7b1fa2","#512da8","#e57373"],["#00796b","#388e3c","#689f38","#afb42b","#81c784"],["#303f9f","#1976d2","#0288d1","#0097a7","#64b5f6"]]

while(True):
	interationNo+=1

	bs='\b'*1000
	print bs,
	print "\bITERATION :",interationNo,

	#E STEP
	oldL=newL[:]

	gammak=[]

	for clsi in range(clsno):
		gammak.append([])
		for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
			gammak[clsi].append([])
			totalgamma=0;
			for clusi in range(CURRENTNOOFCLUSTERS):
				gammak[clsi][datai].append([])
				temp=funcs.subM(data[clsi][datai],mean[clsi][clusi])
				temp=funcs.mulM(funcs.mulM(funcs.transpose(temp),funcs.inverse(covariance[clsi][clusi])),temp)
				gammak[clsi][datai][clusi]=[funcs.det(funcs.inverse(covariance[clsi][clusi]))*math.exp(-0.5*temp[0][0])/math.sqrt(2*3.14)]
				gammak[clsi][datai][clusi][0]*=pik[clsi][clusi]
				totalgamma+=gammak[clsi][datai][clusi][0]
			if(len(gammak[clsi][datai])!=CURRENTNOOFCLUSTERS):
				print "error"
			gammak[clsi][datai]=funcs.divByConstM(gammak[clsi][datai],totalgamma)

	#M STEP
	NK=[]

	for clsi in range(clsno):
		NK.append([])
		for clusi in range(CURRENTNOOFCLUSTERS):
			NK[clsi].append([])
			NK[clsi][clusi]=0
			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
				NK[clsi][clusi]+=gammak[clsi][datai][clusi][0]

	for clsi in range(clsno):
		for clusi in range(CURRENTNOOFCLUSTERS):
			pik[clsi][clusi]=NK[clsi][clusi]/len(data[clsi])
			mean[clsi][clusi]=[[0],[0]]
			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
				mean[clsi][clusi]=funcs.addM(mean[clsi][clusi],funcs.mulByConstM(data[clsi][datai],gammak[clsi][datai][clusi][0]))
			mean[clsi][clusi]=funcs.divByConstM(mean[clsi][clusi],NK[clsi][clusi])
			covariance[clsi][clusi]=[[0,0],[0,0]]
			for datai in range(int(len(data[clsi])*TRAININGDATALIMIT)):
				covariance[clsi][clusi]=funcs.addM(funcs.mulByConstM(funcs.mulM(funcs.subM(data[clsi][datai],mean[clsi][clusi]),funcs.transpose(funcs.subM(data[clsi][datai],mean[clsi][clusi]))),gammak[clsi][datai][clusi][0]),covariance[clsi][clusi])
			covariance[clsi][clusi]=funcs.divByConstM(covariance[clsi][clusi],NK[clsi][clusi])

	newL=calculateL(data,mean,covariance,pik)

#	print newL

	if(emCheckConvergence(oldL,newL,interationNo)):
		break

	datainclusters=[]

	if(PRINTSTEPSFLAG==1):
		fig=plt.figure()
		subplot=fig.add_subplot(111)
		subplot.set_title(str(interationNo))

		for clsi in range(clsno):
			datainclusters.append([])
			for clusi in range(CURRENTNOOFCLUSTERS):
				datainclusters[clsi].append([])
			for datai in range(len(data[clsi])):
				clas=funcs.classify(data[clsi][datai],mean[clsi],covariance[clsi],pik[clsi])
				subplot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=clusterColors[clsi][clas],marker="o")
				datainclusters[clsi][clas].append(data[clsi][datai])

		plt.axis([minx,maxx,miny,maxy])
		fig.savefig(str(interationNo)+".png")

fig=plt.figure()
subplot=fig.add_subplot(111)

color=["#F49292","#A9F5AF","#A1A4FF","#F3F3AF"]
color1=["#E21818","#17E81F","#252BDF","#D2C81D"]
xinc=(maxx-minx)/50
yinc=(maxy-miny)/50

#abcd
gx=[]
for clsi in range(clsno):
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
		
		for clsi in range(clsno):
			temp=0.0
			for clusi in range(len(mean[clsi])):
				temp+=pik[clsi][clusi]*NormalDist([[x],[y]],mean[clsi][clusi],covariance[clsi][clusi])

			totalN=0
			for xx in range(clsno):
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
CS = contourSubPlot.contour(vecx,vecy,vecz,zorder=2)
contourSubPlot.clabel(CS, inline=1, fontsize=10)	

for clsi in range(clsno):
	for datai in range(int(TRAININGDATALIMIT*len(data[clsi]))):
		contourSubPlot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color[clsi],marker="o",markeredgecolor=color[clsi],zorder=1)
		subplot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color1[clsi],marker="o")

plt.axis([minx,maxx,miny,maxy])

print "SAVING CONTOUR GRAPH"
contourFig.savefig("contour.png")
print "SAVING BOUNDARIES GRAPH"
fig.savefig("boundaries.png")
