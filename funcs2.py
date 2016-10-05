import math

THRESHOLD=0.1

def belongsToClusture(dataPoint,means):
	distaces=[]
	for i in range(len(means)):
		distaces.append(math.sqrt(math.pow(dataPoint[0][0]-means[i][0],2)+math.pow(dataPoint[1][0]-means[i][1],2)))
	return indexof(min(distaces))

def checkConvergence(old,new):
	if(math.sqrt(math.pow(new[0]-old[0],2)+math.pow(new[1]-old[1],2))<THRESHOLD):
		return true
	else:
		return false

def kMeansClusture(k,data,minx,maxx,miny,maxy):
	#holding means
	curmeans=[]
	#holding data of different clustures
	datainclusture=[]

	for i in range(k):
		mx=random.uniform(minx,maxx)
		my=random.uniform(miny,maxy)
		curmeans.append([mx,my])
		datainclusture.apend([])

	while(true):
		for i in range(len(data)):
			datainclusture[belongsToClusture(data[i],curmeans)].append(data[i])

		newmeans=[]
		for i in range(len(datainclusture)):
			newmeans.append[datainclusture[i][0][0][0],datainclusture[i][0][1][0]]
			for j in range(1,len(datainclusture[i])):
				newmeans[i]=[newmeans[i][0]+datainclusture[i][j][0][0],newmeans[i][1]+datainclusture[i][j][1][0]];
			newmeans[i][0]=newmeans[i][0]/len(datainclusture[i])
			newmeans[i][1]=newmeans[i][1]/len(datainclusture[i])

		if(checkConvergence(curmeans,newmeans)):
			break

		curmeans=newmeans[:]