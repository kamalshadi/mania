import pylab as pl
import numpy as num
import math
from scipy.optimize import broyden1


def PDF(a,nBins=100):
	pdf,bins=num.histogram(a,nBins,density=True)
	x= bins[0:len(pdf)]
	re=x[1]-x[0]
	print num.trapz(pdf,x)
	return [[o for o in pdf],x]
	
	
def C(a):
	if a==0:
		return 1.0
	return a/(math.exp(a)-1)

def mua(u):
	return lambda a:(a*u+1)*(math.exp(a)-1)-a*math.exp(a)
	
def finda(u,lookup=True):
	if lookup:
		import pickle as pk
		with open('au_lookup') as f:
			a,ul=pk.load(f)
		for j,uu in enumerate(ul):
			if uu>u:
				break
		if j>0:
			return a[j-1]
		else:
			return a[0]
			
	f=mua(u)
	if u==0.5:
		return 0.0
	if u<.5:
		return broyden1(f,-10)
	return broyden1(f,10)
	
def mu(a):
	if a==0:
		return 0.5
	return math.exp(a)/(math.exp(a)-1)-(1.0/a)
	
def fpdf(a):
	if a==0:
		return lambda x:0.5
	c=a/(math.exp(a)-1)
	return lambda x:c*math.exp(a*x)
	
	
def var(a):
	if a==0:
		return 1.0/12
	c=C(a)
	u=mu(a)
	b=1.0/a
	return c*(math.exp(a)*(b-2*b**2+2*b**3)-2*b**3) - u**2
	
	
def maxent(a):
	# a sample from truncated exponential distribution
	w=num.random.rand()
	c=C(a)
	return (math.log((a*w/c)+1))/a
	
def ME_fdr(eg,a,th):
	z=C(a)
	ea=math.exp(a)
	x=z*(ea-math.exp(a*th))/a
	y=z*(ea*math.exp(-a*th)-1)/a
	return x/(eg*y+(1-eg)*x)

def ME_fnr(eg,a,th):
	z=C(a)
	ea=math.exp(a)
	x=z*ea*(1-math.exp(-a*th))/a
	y=z*(math.exp(a*th)-1)/a
	return x/(eg*x+(1-eg)*y)
	
def FDR_th(eg,a,fdr):
	z=(1-fdr*(1-eg))/fdr
	ea=math.exp(a)
	zea=z*ea
	D=math.sqrt((zea+1)**2-4*zea)
	g1=(zea+D)/(2*z)
	g2=(zea-D)/(2*z)
	fg1=True
	fg2=True
	try:
		th1=math.log(g1)/a
	except ValueError:
		fg1=False
	try:
		th2=math.log(g2)/a
	except ValueError:
		fg2=False
	if fg1 and fg2:
		if 0<=th1<=1:
			return th1
		elif  0<=th2<=1:
			return th2
		else:
			print 'Error:All value of theta out of bound'
			return 1.0
	elif fg1:
		if 0<=fg1<=1:
			return th1
		else:
			print 'Error:the_1 out of bound and the_2 non-computable'
			return None
	elif fg2:
		if 0<=fg2<=1:
			return th2
		else:
			print 'Error:the_2 out of bound and the_1 non-computable'
			return None
	else:
		print 'Error: theta non-computable'
		return None
			
	
def kamal_md(a,th):
	z=C(a)
	ea=math.exp(a)
	return z*ea*(1-math.exp(-a*th))/a
	
def kamal_fa(a,th):
	z=C(a)
	ea=math.exp(a)
	return z*(ea-math.exp(a*th))/a
	
def FDRalg(pvalues,q):
	# q is the maximum tolerated FDR
	pv=sorted(pvalues)
	e=-1
	pr=-1
	for i,w in enumerate(pv):
		if w<=(i*q/len(pv)):
			e=i
			pr=pv[i]
	return pr
		
def ME_pvalue(a,t):
	z=C(a)	
	ea=math.exp(a)
	return z*(ea-math.exp(a*t))/a
	

		


