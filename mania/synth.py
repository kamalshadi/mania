import numpy as num
import pickle as pk
from synth_tools import *
from tqdm import tqdm
#~ import pylab as pl
#Note all matrices and vectors must be provided as numpy arrays

def randomGraph(P,dim=None,Sym=True):
	# p is a numpy array
	# Create a random Adjacency Matrix based on probability matrix p
	if dim is None:
		p=P
		l1,l2=num.shape(p)
		if l1!=l2:
			print 'Error: Input matrix not square'
			return None
		l=l1
		A=num.zeros((l,l),dtype='int')
		for i in range(l):
			for j in range(i+1,l):
				r1=p[i,j]
				r2=p[j,i]
				x1=num.random.rand()
				x2=num.random.rand()
				if x1<=r1:
					A[i,j]=1
				if x2<=r2:
					A[j,i]=1
	else:
		l=dim
		A=num.zeros((l,l),dtype='int')
		for i in range(l):
			for j in range(i+1,l):
				if Sym:
					r=P
					x=num.random.rand()
					if x<=r:
						A[i,j]=1
						A[j,i]=1
				else:
					r1=P
					r2=P
					x1=num.random.rand()
					x2=num.random.rand()
					if x1<=r1:
						A[i,j]=1
					if x2<=r2:
						A[j,i]=1
	return A	
		
def kamalModel(a1,a2,A):
	l1,l2=num.shape(A)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	B=num.zeros((l,l))
	for i in range(l):
		for j in range(l):
			if i==j:
				continue
			else:
				if A[i,j]==0:
					P=maxent(a1)
					B[i,j]=P
				else:
					P=maxent(a2)
					B[i,j]=1-P
	return B
	
	
	
def sim(A1,A2):
	# Jaccard similarity on edges
	l1,l2=num.shape(A1)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	l1,l2=num.shape(A1)
	if l1!=l or l2!=l:
		print 'Error: Matrices must be of the same size'
		return None
	P=A1+A2
	P=P-num.diag(num.diag(P))
	P[P>0]=1
	S=A1*A2
	S=S-num.diag(num.diag(S))
	return float(sum(sum(S)))/sum(sum(P))
	

def density(A):
	# returns the density of edges in A
	l1,l2=num.shape(A)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	return float(total_edges)/(l*(l-1))
	
def miss_detection(A1,A2,norm=0):
	# A1 is the ground truth
	l1,l2=num.shape(A1)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	l1,l2=num.shape(A2)
	if l1!=l or l2!=l:
		print 'Error: Matrices must be of the same shape'
		return None
	np=sum(sum(A1-num.diag(num.diag(A1))))
	B=A1-A2
	B[B<0]=0
	nm=sum(sum(B-num.diag(num.diag(B))))
	if norm==0:
		return float(nm)/np
	else:
		return float(nm)*density(A1)/np
	
def numberOfEdges(A):
	return sum(sum(A-num.diag(num.diag(A))))
	
	
def false_alarm(A1,A2,norm=0):
	# A1 is the ground truth
	l1,l2=num.shape(A1)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	l1,l2=num.shape(A2)
	if l1!=l or l2!=l:
		print 'Error: Matrices must be of the same shape'
		return None
	A=1-A1
	na=sum(sum(A-num.diag(num.diag(A))))
	B=A2-A1
	B[B<0]=0
	nm=sum(sum(B-num.diag(num.diag(B))))
	if norm==0:
		return float(nm)/na
	else:
		return float(nm)*density(A)/na
		
def err(A1,A2):
	return false_alarm(A1,A2,norm=1)+miss_detection(A1,A2,norm=1)

def AR(A):
	# Compute Assymetry ratio
	B=A-A.transpose()
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	eu=float(sum(sum(abs(B))))/2
	if total_edges==0:
		return 1.0
	return eu/total_edges

def NAR(A,mode=0):
	# Computes the normalized asymmetry ratio
	l1,l2=num.shape(A)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	eta=1.0/(l*(l-1))
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	AR_rand=1-(float(total_edges)/(l*(l-1)))
	if AR_rand==0:
		return float('inf')
	if mode==1:
		return (AR(A)+eta)/AR_rand
	return AR(A)/AR_rand

def FDR(A,A1):
	Discovered=sum(sum(A1-num.diag(num.diag(A1))))
	B=A1-A
	B[B<0.5]=0
	e=sum(sum(B-num.diag(num.diag(B))))
	if Discovered==0:
		return 0
	return float(e)/Discovered
	
def FNR(A,A1):
	V=1-A1
	nD=sum(sum(V-num.diag(num.diag(V))))
	B=A-A1
	B[B<0.5]=0
	e=sum(sum(B-num.diag(num.diag(B))))
	if nD==0:
		return 0
	return float(e)/nD
	
def myFlat(A):
	l1,l2=num.shape(A)
	l=min(l1,l2)
	D=[0.0]*(l1*l2-l)
	k=0
	for i in range(l1):
		for j in range(l2):
			if i==j:
				continue
			D[k]=A[i,j]
			k=k+1
	return D
	
def U(W,th,sym = 'None'):
	l1,l2=num.shape(W)
	B=num.zeros((l1,l2))
	for i in range(l1):
		for j in range(l2):
			if i==j:
				continue
			else:
				if sym == 'None':
					if W[i,j]>th:
						B[i,j] = 1
				elif sym == "abs":
					# symmetrizing by absolute distance
					if W[i,j]>th and W[j,i]>th:
						B[i,j] = 1
						B[j,i] = 1
					elif (W[i,j]>th and W[j,i]<th) or (W[i,j]<th and W[j,i]>th):
						if (W[i,j]+W[j,i]) > (2*th):
							B[i,j] = 1
							B[j,i] = 1
					else:
						pass
				elif sym == "norm":
					# symmetrizing by absolute distance
					if W[i,j]>th and W[j,i]>th:
						B[i,j] = 1
						B[j,i] = 1
					elif (W[i,j]>th and W[j,i]<th):
						if (W[i,j]-th)/(1-th) > (th-W[j,i])/th:
							B[i,j] = 1
							B[j,i] = 1
					elif (W[i,j]<th and W[j,i]>th):
						if (W[j,i]-th)/(1-th) > (th-W[i,j])/th:
							B[i,j] = 1
							B[j,i] = 1
					else:
						pass
				else:
					raise MANIA_ERROR('sym')
						
	return B
					
	

def NARdetect(W,r=0.001,fg=True, sym = 'None'):
	#r is the resolution
	l1,l2=num.shape(W)
	l=int(num.ceil(1/r))
	t=num.linspace(0,1,l)
	B=num.zeros((l1,l2))
	his=float('inf')
	for w in t:
		C=U(W,w)
		if fg:
			mn=NAR(C)
		else:
			mn=AR(C)
		if mn<his:
			B[:]=C
			his=mn
			th = w
	B = U(W, th, sym = sym)
	return B

	
	
def MLdetect(a1,a2,W,sym = 'None'):
	c1=C(a1)
	c2=C(a2)
	th=(math.log(c2/c1)+a2)/(a1+a2)
	return U(W,th, sym = sym)
	
def bayesDetect(W,eg,a1,a2, sym = 'None'):
	c1=C(a1)*(1-eg)
	c2=C(a2)*(eg)
	th=(math.log(c2/c1)+a2)/(a1+a2)
	return U(W,th, sym = sym)
	
	
def add_edge(A,sim=True):
	l1,l2=num.shape(A)
	fg=False
	for i in range(l1):
		for j in range(i,l2):
			if A[i,j]==0:
				A[i,j]=1
				if sim:
					A[j,i]=1
				fg=True
				break
		if fg:
			break
	return A
				
def remove_edge(A,sim=True):
	l1,l2=num.shape(A)
	fg=False
	for i in range(l1):
		for j in range(i,l2):
			if A[i,j]==1:
				A[i,j]=0
				if sim:
					A[j,i]=0
				fg=True
				break
		if fg:
			break
		
	return A
				
def make_GT(eg,N,sim=True):
	A=randomGraph(eg,dim=N,Sym=True)
	e=density(A)
	if e<eg:
		while e<eg:
			A=add_edge(A,sim)
			e=density(A)
	else:
		while e>eg:
			A=remove_edge(A,sim)
			e=density(A)
	return A
	
def egsim():
	T = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
	fname = "egsim"
	mu1=0.1
	mu2=0.1
	a2=finda(mu2,False)
	a1=a2
	l = 10
	egL=num.linspace(0.1,0.9,l)
	M = 1000
	j1=[0.0]*l
	j2=[0.0]*l
	j3=[0.0]*l
	Data = {}
	ll = l*len(T)
	q = 0
	for eg in egL:
		for tt in T:
			jmp1 = [0.0]*M
			jmp2 = [0.0]*M
			jmp3 = [0.0]*M
			jmp4 = [0.0]*M
			print str(q+1)+"/"+str(ll)+":eg="+str(eg)+"; T="+str(tt)
			q = q + 1
			A=make_GT(eg,50)
			for k in tqdm(range(M),total = M):
				B = kamalModel(a1,a2,A)
				Ah1 = NARdetect(B,0.001,True, sym = "norm")
				Ah2 = NARdetect(B,0.001,True)
				Ah3 = U(B,tt,sym = "norm")
				Ah4 = U(B,tt)
				jmp1[k] = sim(A,Ah1)
				jmp2[k] = sim(A,Ah2)
				jmp3[k] = sim(A,Ah3)
				jmp4[k] = sim(A,Ah4)
			Data[(eg,tt)] = (jmp1, jmp2, jmp3, jmp4)
	with open(fname+'.pk','w') as f:
		pk.dump(Data,f)
		
def usim():
	T = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
	fname = "usim"
	mu2=0.1
	eg = 0.1
	a2=finda(mu2,False)
	muL=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45]
	a1L=[finda(xx,False) for xx in muL]
	l = len(muL)
	M = 100
	Data = {}
	ll = l*len(T)
	q = 0
	for a1_ind, a1 in enumerate(a1L):
		for tt in T:
			jmp1 = [0.0]*M
			jmp2 = [0.0]*M
			jmp3 = [0.0]*M
			jmp4 = [0.0]*M
			print str(q+1)+"/"+str(ll)+":mu="+str(muL[a1_ind])+"; T="+str(tt)
			q = q + 1
			A=make_GT(eg,50)
			for k in tqdm(range(M),total = M):
				B = kamalModel(a1,a2,A)
				Ah1 = NARdetect(B,0.001,True, sym = "norm")
				Ah2 = NARdetect(B,0.001,True)
				Ah3 = U(B,tt,sym = "norm")
				Ah4 = U(B,tt)
				jmp1[k] = sim(A,Ah1)
				jmp2[k] = sim(A,Ah2)
				jmp3[k] = sim(A,Ah3)
				jmp4[k] = sim(A,Ah4)
			Data[(muL[a1_ind],tt)] = (jmp1, jmp2, jmp3, jmp4)
	with open(fname+'.pk','w') as f:
		pk.dump(Data,f)
			
def fixT():
	T=[0.3,0.4,0.5,.6,.7]
	D={}
	M=100
	for t in T:
		j1 = [0.0]*M
		j2 = [0.0]*M
		j3 = [0.0]*M
		j4 = [0.0]*M
		print 'threshold:'+str(t)
		for i in tqdm(range(M),total = M):
			eg=num.random.rand()
			u1=0
			u2=0
			u1=0.45*num.random.rand()
			if u1<0.01:
				u1=u1+0.01
			u2=0.45*num.random.rand()
			if u2<0.01:
				u2=u2+0.01
			a1=finda(u1,False)
			a2=finda(u2,False)
			A=make_GT(eg,50)
			B=kamalModel(a1,a2,A)
			Ah1 = NARdetect(B,r=0.001,fg=True,sym = 'norm')
			Ah2 = NARdetect(B,r=0.001,fg=True)
			Ah3 = U(B,t,sym = 'norm')
			Ah4 = U(B,t)
			j1[i] = sim(A,Ah1)
			j2[i] = sim(A,Ah2)
			j3[i] = sim(A,Ah3)
			j4[i] = sim(A,Ah4)
		D[t] = (j1,j2,j3,j4)
	with open('fixT.pk','w') as f:
		pk.dump(D,f)
		
def synth_probabilistic_anatomy(N,eg,mu1,mu2):
	A=make_GT(eg,N)
	a1=finda(mu1,False)
	a2=finda(mu2,False)
	B = kamalModel(a1,a2,A)
	return B
	
		
#~ if __name__ == '__main__':
	#~ print synth_probabilistic_anatomy(10,.1,.1,.1)
