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

#function for dividing the while matrix by a constant
def mulByConstM(m,const):
	rno=len(m)
	cno=len(m[0])
	returnM=[]
	for r in range(rno):
		returnM.append([])
		for c in range(cno):
			returnM[r].append(m[r][c]*const)
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
def classify(data,me,cov,pik):
	p=[]
	length=len(me)
	#checking for every class
	for cl in range(length):
		
		p.append(g(data,cov[cl],me[cl],pik[cl]))
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

		for datai in range (int(math.floor(0.75*dtno)),dtno):
			for cli in range(clsno):
				temp=0.0
				for clusi in range(len(mean[cli])):
					temp+=pik[cli][clusi]*NormalDist([[x],[y]],mean[cli][clusi],covariance[cli][clusi])

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
		
			mgx=gx.index(max(gx))
			p[mgx]+=1

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