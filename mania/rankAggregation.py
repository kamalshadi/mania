# rank aggregation used by group mania
import numpy as num 

def KendalMatrix(L):
	# L is the list of lists, each element being a confidence of the object
	m=len(L[0])
	W=num.zeros((m,m))
	for i in range(m-1):
		for j in range(i+1,m):
			for lis in L:
				if lis[i]>lis[j]:
					W[i,j]=W[i,j]+1
				elif lis[i]<lis[j]:
					W[j,i]=W[j,i]+1
				else:
					pass
	return W
	
def KendalMatrix(L):
	# L is the list of lists, each element being a confidence of the object
	m=len(L[0])
	W=num.zeros((m,m))
	for i in range(m-1):
		for j in range(i+1,m):
			for lis in L:
				if lis[i]>lis[j]:
					W[i,j]=W[i,j]+1
				elif lis[i]<lis[j]:
					W[j,i]=W[j,i]+1
				else:
					pass
	return W
	
def KendalMatrix2(L):
	# L is the list of lists, each element being a confidence of the object
	m=len(L[0])
	W=num.zeros((m,m))
	for i in range(m-1):
		for j in range(i+1,m):
			for lis in L:
				if lis.index(i)<lis.index(j):
					W[i,j]=W[i,j]+1
				elif lis.index(i)>lis.index(j):
					W[j,i]=W[j,i]+1
				else:
					pass
	return W
	
	
def agg(W,L):
	l=len(L)
	if l==2:
		if W[L[0],L[1]]>W[L[1],L[0]]:
			return L
		else:
			return [L[1],L[0]]
	elif l<2:
		return L
	else:
		piv=num.random.randint(0,l)
		Ll=[]
		Lr=[]
		for i in range(l):
			if i==piv:
				continue
			if W[L[i],L[piv]]>W[L[piv],L[i]]:
				Ll.append(L[i])
			else:
				Lr.append(L[i])
		return agg(W,Ll)+[L[piv]]+agg(W,Lr)
	
	
						
	
