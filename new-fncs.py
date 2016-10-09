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