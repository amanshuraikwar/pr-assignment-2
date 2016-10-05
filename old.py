import sys
import math
import matplotlib.pyplot as plt

def generateGP(minv,maxv):
	print maxv,minv
	rv=math.exp(math.log(float(maxv)/float(minv))/float(9.0))
	gp=[]
	for i in range(10):
		gp.append(minv*math.pow(rv,i))
	return gp	

#function for addimg matrices
def addM(m1,m2):
	rno=len(m1)
	cno=len(m1[0])
	returnM=[]
	for r in range(rno):
		returnM.append([])
		for c in range(cno):
			s=m1[r][c]+m2[r][c]
			returnM[r].append(s)
		#endoffor
	#endoffor	
	return returnM

#function for subtracting matrices
def subM(m1,m2):
	rno=len(m1)
	cno=len(m1[0])
	returnM=[]
	for r in range(rno):
		returnM.append([])
		for c in range(cno):
			s=m1[r][c]-m2[r][c]
			returnM[r].append(s)
		#endoffor
	#endoffor	
	return returnM

#function for doing a transpose of matrix
def transpose(m):
	rno=len(m)
	cno=len(m[0])
	returnM=[]
	for c in range(cno):
		returnM.append([])
		for r in range(rno):
			returnM[c].append(m[r][c])
		#endoffor
	#endoffor
	return returnM

#function for multiplying two matrices
def mulM(m1,m2):
	if(len(m1[0])!=len(m2)):
		return False
	#endofif
	returnM=[]
	rno1=len(m1)
	cno1=len(m1[0])
	rno2=len(m2)
	cno2=len(m2[0])
	for r1 in range(rno1):
		returnM.append([])
		for c2 in range(cno2):
			s=0
			for c1 in range(cno1):
				s+=m1[r1][c1]*m2[c1][c2]
			#endoffor
			returnM[r1].append(s)
		#endoffor
	#endoffor	
	return returnM	

#function for dividing the while matrix by a constant
def divByConstM(m,const):
	rno=len(m)
	cno=len(m[0])
	returnM=[]
	for r in range(rno):
		returnM.append([])
		for c in range(cno):
			returnM[r].append(m[r][c]/const)
		#endoffor
	#endoffor	
	return returnM;		

#returns identity matrix for given order as argument
def I(x):
	returnM=[]
	for r in range(x):
		returnM.append([])
		for c in range(x):
			if(r==c):
				returnM[r].append(1)
			#endofif
			else:
				returnM[r].append(0)
			#endofelse
		#endoffor
	#endoffor		
	return returnM

#converts a normal matrix to identity matrix
def convertToDiagonal(m):
	o=len(m)
	returnM=[]
	for r in range(o):
		returnM.append([])
		for c in range(o):
			if(r==c):
				returnM[r].append(m[r][c])
			#endofif
			else:
				returnM[r].append(0)
			#endofelse
		#endoffor
	#endoffor		
	return returnM

#to update, currently only for 2X2 matrix
def inverse(m):
	returnM=[[0,0],[0,0]]
	const=m[0][0]*m[1][1]-m[0][1]*m[1][0]
	returnM[1][1]=m[0][0]
	returnM[0][0]=m[1][1]

	returnM[0][1]=-m[0][1]
	returnM[1][0]=-m[1][0]

	returnM=divByConstM(returnM,const)

	return returnM

def det(m):
	return m[0][0]*m[1][1]-m[1][0]*m[0][1]

#discriminating function
def g(x,co,me,p):
	temp=subM(x,me)

	returnM=mulM(mulM(transpose(temp),inverse(co)),temp)

	return (-0.5*returnM[0][0]-0.5*math.log(det(co))-math.log(2*3.14)+math.log(p))

#function that classifies the data point
def classify(data,me,cov):
	p=[]
	length=len(me)
	#checking for every class
	for cl in range(length):
		
		p.append(g(data,cov[cl],me[cl],classP[cl]))
	#endoffor

	return p.index(max(p))

#classifies a large data set
def runClassifier(data,hitRate,me,cov,subplot):
	print "runing classifier"
	clsno=len(data)
	#for data of every class
	for cl in range(clsno):
		dtpntno=len(data[cl])
		trainingdataboundary=int(math.floor(0.75*dtpntno))	
		#dtpnti=int(math.floor(0.75*length1))
		#for every data point of that class
		for dtpnti in range(dtpntno):
			clas=classify(data[cl][dtpnti],me,cov)
			if(dtpnti<trainingdataboundary):
				subplot.plot(data[cl][dtpnti][0][0],data[cl][dtpnti][1][0],color=color1[clas],marker='o')
			else:
				if(cl==clas):
					hitRate[cl]+=1
				#enfofif	
		#endoffor
	#endoffor

#classifies a large data set
def plotTrainingData(data,me,cov,subplot):
	print "plotting training data"
	clsno=len(data)
	#for data of every class
	for cl in range(clsno):
		dtpntno=len(data[cl])	
		#dtpnti=int(math.floor(0.75*length1))
		#for every data point of that class
		for dtpnti in range(int(math.floor(0.75*dtpntno))):
			clas=classify(data[cl][dtpnti],me,cov)
			subplot.plot(data[cl][dtpnti][0][0],data[cl][dtpnti][1][0],color=color[clas],marker='o',markeredgecolor=color[clas],zorder=1)
			#enfofif	
		#endoffor
	#endoffor

#this will only work for upto 3 classes
def plotForEveryPair(data,me,cov,title):
	print "plotting every pair"
	fig=plt.figure();
	fig.suptitle(title)
	clsno=len(data)
	spi=221
	for clsi in range(clsno):
		for cli in range (clsi+1,clsno):
			print "plotting for ",str(clsi)," and ",str(cli)
			currentclasses=[clsi,cli]
			ax=fig.add_subplot(spi)
			ax.set_title(str(clsi)+" and "+str(cli))
			spi+=1
			y=miny
			while (y < maxy):
				x=minx
				while (x < maxx):
					p=[g([[x],[y]],cov[clsi],me[clsi],classP[clsi]),g([[x],[y]],cov[cli],me[cli],classP[cli])]
					clas=currentclasses[p.index(max(p))]
					ax.plot(x,y,color=color[clas],marker="o",markeredgecolor=color[clas])
					x+=increamentorx
				y+=increamentory
			for ccli in [clsi,cli]:
				dtpntno=len(data[ccli])	
				#for every data point of that class
				for dtpnti in range(int(math.floor(0.75*dtpntno))):
					clas=classify(data[ccli][dtpnti],me,cov)
					ax.plot(data[ccli][dtpnti][0][0],data[ccli][dtpnti][1][0],color=color1[clas],marker='o')
	fig.savefig(str(sys.argv[len(sys.argv)-1])+"/"+title+".png")

def plotBoundries(subplot,me,cov,subplot1):
	print "plotting boundaries"
	vecx=[]
	vecy=[]
	vecz=[]
	yi=0
	y=miny
	mingx=g([[minx],[miny]],cov[0],me[0],classP[0])
	maxgx=mingx
	plt.figure()
	#print minx,miny,maxx,maxy
	while (y < maxy):
		vecz.append([])
		x=minx
		vecy.append(y)
		while (x < maxx):
			if(y==miny):
				vecx.append(x)
	
			clas=classify([[x],[y]],me,cov)
			gx=g([[x],[y]],cov[clas],me[clas],classP[clas])
		

			if(gx<mingx):
				mingx=gx
			if(gx>maxgx):
				maxgx=gx

			vecz[yi].append(gx)
			subplot.plot(x,y,color=color[clas],marker="o",markeredgecolor=color[clas])
			x+=increamentorx
		#print "came out"	
		y+=increamentory
		yi+=1
	CS = subplot1.contour(vecx,vecy,vecz,zorder=2,level=generateGP(mingx,maxgx))
	subplot1.clabel(CS, inline=1, fontsize=10)	

def confusionM(data,me,co):
	clsno=len(data)
	returnM=[]
	p=[]
	for clsi in range(clsno):
		p.append(0)

	for clsi in range(clsno):
		dtno=len(data[clsi])
		for cli in range(clsno):
			p[cli]=0
		for datai in range (int(math.floor(0.75*dtno)),dtno):
			k=classify(data[clsi][datai],me,co)
			p[k]+=1
		#bakchodi
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

#universal matrices
data=[]
mean=[]
covariance=[]
dtdim=0
classP=[]

fig=plt.figure()
fig1=plt.figure()
fig2=plt.figure()
fig3=plt.figure()
ax6=fig2.add_subplot(221)
ax6.set_title(" 1 a i ")
ax7=fig2.add_subplot(222)
ax7.set_title(" 1 a ii ")
ax8=fig2.add_subplot(223)
ax8.set_title(" 1 b ")
ax9=fig3.add_subplot(221)
ax9.set_title(" 2 a ")
ax10=fig3.add_subplot(222)
ax10.set_title(" 2 b ")
ax11=fig3.add_subplot(223)
ax11.set_title(" 2 c ")

ax1=fig.add_subplot(221)
ax1.set_title(" 1 a i ")
ax12=fig.add_subplot(222)
ax12.set_title(" 1 a ii ")
ax2=fig.add_subplot(223)
ax2.set_title(" 1 b ")
ax3=fig1.add_subplot(221)
ax3.set_title(" 2 a ")
ax4=fig1.add_subplot(222)
ax4.set_title(" 2 b ")
ax5=fig1.add_subplot(223)
ax5.set_title(" 2 c ")

#ranges of data
minx=0
maxx=0
miny=0
maxy=0

color1=["#E21818","#17E81F","#252BDF","#D2C81D"]
color=["#F49292","#A9F5AF","#A1A4FF","#F3F3AF"]

#reading files and storing data in data[]
datavalueseparator=" "
nooffiles=len(sys.argv)-1
#for each file
for filei in range(1,nooffiles):
	with open(sys.argv[filei]) as curFile:
		data.append([])
		mean.append([])
		covariance.append([])
		content = curFile.read().splitlines()
		lineno=len(content);
		#for each line
		for linei in range(lineno):
			value=content[linei].split(datavalueseparator);
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
increamentorx=(maxx-minx)/50
increamentory=(maxy-miny)/50

dtdim=len(data[0][0])

totaldatano=0
for i in range(clsno):
	totaldatano+=len(data[i])
for i in range(clsno):
	classP.append(float(len(data[i]))/float(totaldatano))	

#print classP

#calculating mean
#for each class
for clsi in range(clsno):
	tempSum=data[clsi][0]
	#how much data to be taken as training	
	dtpntno=int(math.floor(0.75*len(data[clsi])))
	for dtpnti in range(1,dtpntno):
		tempSum=addM(data[clsi][dtpnti],tempSum)
	#endoffor
	mean[clsi]=divByConstM(tempSum,dtpntno)
#endoffor
#print mean
#findnig covariance matrix
#for each class
for clsi in range(clsno):
	covariance[clsi]=mulM(subM(data[clsi][0],mean[clsi]),transpose(subM(data[clsi][0],mean[clsi])))
	#how much data to be taken as training
	dtpntno=int(math.floor(0.75*len(data[clsi])))
	#for each data point
	for dtpnti in range(1,dtpntno):
		covariance[clsi]=addM(mulM(subM(data[clsi][dtpnti],mean[clsi]),transpose(subM(data[clsi][dtpnti],mean[clsi]))),covariance[clsi])
	#endoffor
	covariance[clsi]=divByConstM(covariance[clsi],dtpntno)	
#endoffor
#print covariance
hitRate=[]
#initialising hitRate
for i in range(clsno):
	hitRate.append(0)
#endoffor

totaltestdatano=0
for clsi in range(clsno):
	totaltestdatano+=int(math.ceil(0.25*len(data[clsi])))
print "-------------------------------------------------------------------"

#part 1 (a) (i)
covariance1ai=[]
covariance1ai.append([])
covariance1ai[0]=covariance[0]
#average of each class
for clsi in range(1,clsno):
	covariance1ai.append([])
	covariance1ai[0]=addM(covariance[clsi],covariance1ai[0])
#endoffor
covariance1ai[0]=divByConstM(covariance1ai[0],clsno)
#making cov for each class same
for clsi in range(1,clsno):
	covariance1ai[clsi]=covariance1ai[0]
#endoffor
#-------------------------------------------------------------------
plotTrainingData(data,mean,covariance1ai,ax6)
plotBoundries(ax1,mean,covariance1ai,ax6)
plotForEveryPair(data,mean,covariance1ai,"1 (a) (i) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance1ai,ax1)
print "1 (a) (i)"
confM=confusionM(data,mean,covariance1ai)
print "confusion matrix : ",str(confM)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"

#part i (a) (ii)
mean1aii=[]
covariance1aii=[]
mean1aii.append([])
mean1aii[0]=divByConstM(mean[0],1/float(len(data[0])))
lendata=len(data[0])
#average of each class
for clsi in range(1,clsno):
	mean1aii.append([])
	mean1aii[0]=addM(divByConstM(mean[clsi],1/float(len(data[clsi]))),mean1aii[0])
	lendata+=len(data[clsi])
#endoffor
mean1aii[0]=divByConstM(mean1aii[0],lendata)
for clsi in range(1,clsno):
	mean1aii[clsi]=mean1aii[0]
#endoffor
covariance1aii.append([])
covariance1aii[0]=subM(I(dtdim),I(dtdim))
for clsi in range(clsno):
	#how much data to be taken as training
	dtpntno=int(math.floor(0.75*len(data[clsi])))
	#for each data point
	for dtpnti in range(dtpntno):
		covariance1aii[0]=addM(mulM(subM(data[clsi][dtpnti],mean1aii[0]),transpose(subM(data[clsi][dtpnti],mean1aii[0]))),covariance1aii[0])
	#endoffor
covariance1aii[0]=divByConstM(covariance1aii[0],lendata)
for clsi in range(1,clsno):
	covariance1aii.append([])
	covariance1aii[clsi]=covariance1aii[0]
#endoffor
#initialising hitRate
for i in range(clsno):
	hitRate[i]=0
#endoffor
#-------------------------------------------------------------------
plotBoundries(ax12,mean,covariance1aii,ax7)
plotTrainingData(data,mean,covariance1aii,ax7)
plotForEveryPair(data,mean,covariance1aii,"1 (a) (ii) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance1aii,ax12)
print "1 (a) (ii)"
confM=confusionM(data,mean,covariance1aii)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "confusion matrix : ",str(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"

#part 1 (b)
#initialising hitRate
for i in range(clsno):
	hitRate[i]=0
#endoffor
#-------------------------------------------------------------------
plotBoundries(ax2,mean,covariance,ax8)
plotTrainingData(data,mean,covariance,ax8)
plotForEveryPair(data,mean,covariance,"1 (b) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance,ax2)
print "1 (b) "
confM=confusionM(data,mean,covariance)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "confusion matrix : ",str(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"

#part 2 (a)
covariance2a=[]
covariance2a.append([])
covariance2a[0]=convertToDiagonal(covariance[0])
#average of each class
for clsi in range(1,clsno):
	covariance2a.append([])
	covariance2a[0]=addM(convertToDiagonal(covariance[clsi]),covariance2a[0])
#endoffor
covariance2a[0]=divByConstM(covariance2a[0],clsno)
covAvg=0
#taking average of diagonal values
for r in range(dtdim):
	covAvg+=covariance2a[0][r][r]
#endoffor
covAvg/=dtdim
#setting diagonals elements to avg
for r in range(dtdim):
	covariance2a[0][r][r]=covAvg	
#endoffor
for clsi in range(1,clsno):
	covariance2a[clsi]=covariance2a[0]
#eendoffor
#initialising hitRate
for i in range(clsno):
	hitRate[i]=0
#endoffor
#-------------------------------------------------------------------
plotBoundries(ax3,mean,covariance2a,ax9)
plotTrainingData(data,mean,covariance2a,ax9)
plotForEveryPair(data,mean,covariance2a,"2 (a) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance2a,ax3)
print "2 (a) "
confM=confusionM(data,mean,covariance2a)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "confusion matrix : ",str(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"
#part 2 (b)
covariance2b=[]
covariance2b.append([])
covariance2b[0]=convertToDiagonal(covariance[0])
#taking average of each class
for clsi in range(1,clsno):
	covariance2b.append([])
	covariance2b[0]=addM(convertToDiagonal(covariance[clsi]),covariance2b[0])
#endoffor
covariance2b[0]=divByConstM(covariance2b[0],clsno)
#setting each class to same cov
for clsi in range(1,clsno):
	covariance2b[clsi]=covariance2b[0]
#initialising hitRate
for i in range(clsno):
	hitRate[i]=0
#endoffor
#-------------------------------------------------------------------
plotBoundries(ax4,mean,covariance2b,ax10)
plotTrainingData(data,mean,covariance2b,ax10)
plotForEveryPair(data,mean,covariance2b,"2 (b) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance2b,ax4)
print "2 (b) "
confM=confusionM(data,mean,covariance2b)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "confusion matrix : ",str(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"
#part 2 (c)
covariance2c=covariance
#coverting cov to diagonal for each class
for clsi in range(clsno):
	covariance2c[clsi]=convertToDiagonal(covariance2c[clsi])
#endoffor
#initialising hitRate
for i in range(clsno):
	hitRate[i]=0
#endoffor
#-------------------------------------------------------------------
plotBoundries(ax5,mean,covariance2c,ax11)
plotTrainingData(data,mean,covariance2c,ax11)
plotForEveryPair(data,mean,covariance2c,"2 (c) ")
#-------------------------------------------------------------------
runClassifier(data,hitRate,mean,covariance2c,ax5)
print "2 (c) "
confM=confusionM(data,mean,covariance2c)
acc=accuracy(confM,totaltestdatano)
prec=precision(confM)
meprec=mean_precision(confM)
rec=recall(confM)
merec=mean_recall(confM)
fmes=f_measure(confM)
mefmes=mean_f_measure(confM)
print "confusion matrix : ",str(confM)
print "accuracy : ",str(acc)
print "precision : ",str(prec)
print "mean precision : ",str(meprec)
print "recall : ",str(rec)
print "mean recall : ",str(merec)
print "f measure : ",str(fmes)
print "mean fmeasure",str(mefmes)
print "-------------------------------------------------------------------"
#-------------------------------------------------------------------
plt.axis([minx,maxx,miny,maxy])

fig.savefig(str(sys.argv[len(sys.argv)-1])+"/1.png")
fig1.savefig(str(sys.argv[len(sys.argv)-1])+"/2.png")
fig2.savefig(str(sys.argv[len(sys.argv)-1])+"/3.png")
fig3.savefig(str(sys.argv[len(sys.argv)-1])+"/4.png")

#plt.show()
#-------------------------------------------------------------------
