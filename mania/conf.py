# This module compute confidence for mania binary network
import numpy as num
import pickle as pk

def conf(fp):
	if fp[-1] != '/':
		fp = fp + '/'
	# load mania network
	fn = fp + 'MANIA/mania_binary.net'
	with open(fn) as f:
		net = pk.load(f)
	# load mania meta data
	fn = fp + 'MANIA/mania_binary.meta'
	with open(fn) as f:
		meta = pk.load(f)
	
	rhos=meta['density']
	
	l1,l2=net.shape
	N=l1

	fn = fp + 'MANIA/raw.res'
	with open(fn) as f:
		DD = pk.load(f)

	curval=1
	D={}
	for d in DD:
		den=d[2]
		net=d[0]
		for i in range(N):
			for j in range(N):
				if net[i,j]==0 or i==j:
					continue
				try:
					D[(i,j)]
				except KeyError:
					D[(i,j)]=curval
		if den<rhos:
			curval=(rhos-den)/rhos
		else:
			curval=(rhos-den)/(1-rhos)
	C = num.zeros((N,N))
	for i in range(N):
		for j in range(N):
			try:
				C[i,j]= D[(i,j)]
			except KeyError:
				C[i,j] = None
				
	# save results in numpy array
	fn = fp + 'MANIA/mania_binary.conf'
	with open(fn,'w') as f:
		DD = pk.dump(C,f)
