###################Tools and metrics##################

'''The functions in this python module compute 
algebric metrics from adjacency matrices used in MANIA'''

import numpy as num

def AR(A):
	# Computes Assymetry ratio
	B=A-A.transpose()
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	eu=float(sum(sum(abs(B))))/2
	if total_edges==0:
		return 1.0
	return eu/total_edges

def NAR(A):
	# Computes the normalized asymmetry ratio
	l1,l2=num.shape(A)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	AR_rand=1-(float(total_edges)/(l*(l-1)))
	if AR_rand==0:
		return float('inf')
	return AR(A)/AR_rand

def sim(A1,A2,mode='jac'):
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
	if mode=='jac':
		P=A1+A2
		P=P-num.diag(num.diag(P))
		P[P>0]=1
		S=A1*A2
		S=S-num.diag(num.diag(S))
		return float(sum(sum(S)))/sum(sum(P))
	else:
		S=A1*A2
		S=S-num.diag(num.diag(S))
		A1=A1-num.diag(num.diag(A1))
		A2=A2-num.diag(num.diag(A2))
		denom=min(sum(sum(A1)),sum(sum(A2)))
		return float(sum(sum(S)))/denom

def density(A):
	# returns the density of edges in A
	l1,l2=num.shape(A)
	if l1!=l2:
		print 'Error: Input matrix not square'
		return None
	l=l1
	total_edges=sum(sum(A-num.diag(num.diag(A))))
	return float(total_edges)/(l*(l-1))
