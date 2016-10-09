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
MAXITERATIONS=10
EMMAXITERATIONS=10
EMTHRESHOLD=0.1
TRAININGDATALIMIT=0.75
CURRENTNOOFCLUSTERS=5
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
	temp=numpy.dot(numpy.dot(numpy.transpose(temp),numpy.linalg.inv(covariance)),temp)
	if(temp[0][0]<0):
		print "error less than zero"
	N=math.exp(-0.5*temp[0][0])/((2*3.14)*math.sqrt(numpy.linalg.det(covariance)))
	return N

#hardcoded
startIndexClass=[7020,8316,7884]

def confusionM(data,mean,covariance,pik):
	clsno=len(data)
	returnM=[]
	p=[]

	gx=[]
	for clsi in range(clsno):
		gx.append(0)
		p.append(0)

	for clsi in range(clsno):
		dtno=len(data[clsi])
		
		for cli in range(clsno):
			p[cli]=0

		datai=startIndexClass[clsi]
		while(datai<dtno):
			classified=[0,0,0]
			for iiii in range(36):
				for cli in range(clsno):
					temp=0.0
					for clusi in range(len(mean[cli])):
						temp+=pik[cli][clusi]*NormalDist(data[clsi][datai],mean[cli][clusi],covariance[cli][clusi])
						#print NormalDist(data[clsi][datai],mean[cli][clusi],covariance[cli][clusi])
					totalN=0
					for xx in range(clsno):
						totalN+=len(data[xx])

					temp*=float(len(data[cli]))/float(totalN)

					#hopefully this does not execute
					# if(temp>1):
					# 	print pik[cli]
					# 	print "error probability :",temp
					# 	print "probability more that 1 error"
					# 	print "TERMINATING"
					# 	exit()

					temp=math.log(temp)
					gx[cli]=temp
			
				mgx=gx.index(max(gx))
				classified[mgx]+=1
				datai+=1;
			p[classified.index(max(classified))]+=1
		returnM.append(p[:])
	return returnM

def accuracy(confmatrix,ntest):
	dmatrix=len(confmatrix)
	accu=0.0
	for i in range (dmatrix):
		accu=accu+confmatrix[i][i]
	accu=accu/ntest
	return accu

def per_accuracy(confmatrix,ntest):
	p_accu=accuracy(confmatrix,ntest)
	return p_accu*100

def precision(confmatrix):
	pre=[]
	dmatrix=len(confmatrix)
	for i in range (dmatrix):
		tp=0.0		
		for j in range (dmatrix):
			tp=tp+confmatrix[j][i]
		tc=confmatrix[i][i]
		if (tp==0):
			pre.append(-1)
		else:	
			pre.append(tc/tp)
	return pre

def mean_precision(confmatrix):
	pre=precision(confmatrix)
	temp=0.0
	lenn=len(pre)
	for i in range(lenn):
		if(pre[i]==-1):
			lenn-=1
		else:	
			temp=temp+pre[i]
	return (temp/lenn)

def recall(confmatrix):
	rec=[]
	dmatrix=len(confmatrix)
	for i in range(dmatrix):
		n=0.0
		for j in range(dmatrix):
			n=n+confmatrix[i][j]
		tc=confmatrix[i][i]
		rec.append(tc/n)
	return rec

def mean_recall(confmatrix):
	rec=recall(confmatrix)
	temp=0.0
	for i in range (len(rec)):
		temp=temp+rec[i]
	return (temp/(len(rec)))

def f_measure(confmatrix):
	f=[]
	pre=precision(confmatrix)
	rec=recall(confmatrix)
	temp=0.0
	for i in range (len(confmatrix)):
		if (pre[i]==-1):
			f.append(-1)
		else:
			temp=pre[i]*rec[i]*2/(pre[i]+rec[i])
			f.append(temp)

	return f

def mean_f_measure(confmatrix):
	f=f_measure(confmatrix)
	temp=0.0
	lenn=len(f)
	for i in range (lenn):
		if(f[i]==-1):
			lenn-=1
		else:
			temp=temp+f[i]
	return (temp/(lenn))

#------------------------------------------------------------------------------------------MAIN-CODE

data=[]
datainclusters=[]

#reading files and storing data in data[]
DATAVALUESEPARATOR=" "
noofdir=len(sys.argv)-1
PRINTSTEPSFLAG=int(sys.argv[len(sys.argv)-1])
PRINTGRAPHSFLAG=0
VECTORDIMENSION=0

nooffiles=len(sys.argv)-1

PRINTSTEPSFLAG=int(sys.argv[len(sys.argv)-1])

#for each file
for diri in range(1,noofdir):
	fileList=os.listdir(sys.argv[diri])
	data.append([])
	datainclusters.append([])
	lineiOffset=0
	for filei in range(len(fileList)):
		
		bs='\b'*1000
		print bs,
		print "\bDIRECTORY :",diri,"READ :",round(100*float(filei+1)/float(len(fileList)),2),"%",

		with open(str(sys.argv[diri])+fileList[filei]) as curFile:
			content = curFile.read().splitlines()
			lineno=len(content);
			for linei in range(lineno):
				value=content[linei].split(DATAVALUESEPARATOR);
				data[diri-1].append([])
				for valuei in range(len(value)):
					try:
						data[diri-1][lineiOffset+linei].append([float(value[valuei])])	
						if(diri==1 and filei==0 and linei==0):
							VECTORDIMENSION+=1
					except ValueError:
						continue
		lineiOffset+=lineno
	print ""

data=numpy.array(data)
mean=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS,VECTORDIMENSION,1),dtype=numpy.float)
covariance=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS,VECTORDIMENSION,VECTORDIMENSION),dtype=numpy.float)
pik=numpy.zeros((noofdir-1,CURRENTNOOFCLUSTERS),dtype=numpy.float)
oldL=numpy.zeros((noofdir-1),dtype=numpy.float)
newL=numpy.zeros((noofdir-1),dtype=numpy.float)

#no of classes
CLSNO=len(data)
totaltestdatano=0

for i in range(CLSNO):
	mean[i]=kMeansClusture(CURRENTNOOFCLUSTERS,data[i],datainclusters[i])
	totaltestdatano+=(1-TRAININGDATALIMIT)*len(data[i])

totaltestdatano/=36

#print mean
print "CALCULTING INTITAL PIK"
for clsi in range(CLSNO):
	for clusi in range(CURRENTNOOFCLUSTERS):
		noofdata=len(datainclusters[clsi][clusi])
		for datai in range(noofdata):
			tempSub=datainclusters[clsi][clusi][datai]-mean[clsi][clusi]
			covariance[clsi][clusi]=numpy.dot(tempSub,numpy.transpose(tempSub)/noofdata)+covariance[clsi][clusi]
			
		#covariance[clsi][clusi]=numpy.array(covariance[clsi][clusi])/noofdata
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

	#print covINVdet

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
				#covariance[clsi][clusi]=numpy.array(covariance[clsi][clusi])/NK[clsi][clusi]

	#print "CALCULATING L"
	#newL=calculateL(data,mean,covariance,pik)
	#print covariance
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

confM=confusionM(data,mean,covariance,pik)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print confM
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)

# fig=plt.figure()
# subplot=fig.add_subplot(111)

# color=["#F49292","#A9F5AF","#A1A4FF","#F3F3AF"]
# color1=["#E21818","#17E81F","#252BDF","#D2C81D"]
# xinc=(maxx-minx)/50
# yinc=(maxy-miny)/50

# #abcd
# gx=[]
# for clsi in range(CLSNO):
# 	gx.append(0)

# vecx=[]
# vecy=[]
# vecz=[]

# mingx=0
# maxgx=0

# yi=0
# y=miny
# while(y<maxy):
	
# 	vecy.append(y)
# 	vecz.append([])

# 	bs='\b'*1000
# 	print bs,
# 	print "\bPLOTTING GRAPH :",round((y-miny)*100/(maxy-miny),2),"%",

# 	x=minx
# 	while(x<maxx):
		
# 		for clsi in range(CLSNO):
# 			temp=0.0
# 			for clusi in range(len(mean[clsi])):
# 				temp+=pik[clsi][clusi]*NormalDist([[x],[y]],mean[clsi][clusi],covariance[clsi][clusi])

# 			totalN=0
# 			for xx in range(CLSNO):
# 				totalN+=len(data[xx])

# 			temp*=float(len(data[clsi]))/float(totalN)

# 			#hopefully this does not execute
# 			if(temp>1):
# 				print "error probability :",temp
# 				print "probability more that 1 error"
# 				print "TERMINATING"
# 				exit()

# 			temp=math.log(temp)
# 			gx[clsi]=temp
		
# 		mgx=max(gx)

# 		vecz[yi].append(mgx)

# 		if(mingx>mgx):
# 			mingx=mgx
# 		if(maxgx<mgx):
# 			maxgx=mgx

# 		if(y==miny):
# 			vecx.append(x)
# 			if(x==minx):
# 				mingx=mgx
# 				maxgx=mingx

# 		clas=gx.index(mgx)
# 		subplot.plot(x,y,color=color[clas],marker="o",markeredgecolor=color[clas])
# 		x+=xinc
# 	y+=yinc
# 	yi+=1

# print ""

# plt.axis([minx,maxx,miny,maxy])

# contourFig=plt.figure()
# contourSubPlot=contourFig.add_subplot(111)
# #CS = contourSubPlot.contour(vecx,vecy,vecz,zorder=2,level=funcs.generateGP(mingx,maxgx))
# CS = contourSubPlot.contour(vecx,vecy,vecz,zorder=2)
# contourSubPlot.clabel(CS, inline=1, fontsize=10)	

# for clsi in range(CLSNO):
# 	for datai in range(int(TRAININGDATALIMIT*len(data[clsi]))):
# 		contourSubPlot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color[clsi],marker="o",markeredgecolor=color[clsi],zorder=1)
# 		subplot.plot(data[clsi][datai][0][0],data[clsi][datai][1][0],color=color1[clsi],marker="o")

# plt.axis([minx,maxx,miny,maxy])

# print "SAVING CONTOUR GRAPH"
# contourFig.savefig("contour.png")
# print "SAVING BOUNDARIES GRAPH"
# fig.savefig("boundaries.png")