# optimization engine of MANIA

import pylab as pl
import pickle as pk

def localmin(ar,tol=-1,nei=3):
	minInd=[]
	if tol<0:
		minn=min(ar)
		for i,w in enumerate(ar):
			if w==minn:
				minInd=minInd+[i]
		return minInd
	l=len(ar)
	lastMin=-nei-1
	for i in range(l):
		if ar[i] > float(tol)/100:
			continue
		minF=True
		for j in range(1,nei+1):
			if i-j<0:
				pri=0
			else:
				pri=i-j
			if i+j >= l:
				aft=l-1
			else:
				aft=i+j
			if ar[i] > ar[pri] or ar[i] > ar[aft]:
				minF=False
				break
		if minF:
			if i-lastMin <= nei:
				lastMin=i
				if len(minInd)>0:
					minInd[-1]=i
				else:
					minInd.append(i)
			else:
				lastMin=i
				minInd.append(i)
	his = float('inf')
	#~ for q,w in enumerate(ar):
		#~ if w<.04 or w > .5:
			#~ continue
		#~ if w < his:
			#~ his = w
			#~ o = q-1
	#~ minInd.append(o)
	#~ print minInd
	return minInd
	
def finalNet(fp, v = 0):
	if v == 0:
		path_to_out = fp + '/MANIA/'
		fName = path_to_out + 'raw.res'
	else:
		path_to_out = fp
		fName = path_to_out + 'agg.res'
		
	with open(fName) as f:
		D=pk.load(f)
	den=[w[2] for w in D]
	nar=[w[1] for w in D]
	T=[w[-1] for w in D]
	
	#~ pl.plot(den,nar,'-*')
	#~ pl.xlabel('density',fontsize = 20)
	#~ pl.ylabel('Normalized asymmetry ratio', fontsize = 20)
	# avoiding corner cases
	l = len(den)
	for i in range(l):
		if den[i]<.05 or den[i]>.95:
			nar[i] = float('Inf') 
	Ind=localmin(nar)
	#~ print set(Ind)
	#~ print '-----'
	#~ for j in Ind:
		#~ pl.plot(den[j],nar[j],'rs')
	#~ pl.show()
	#~ raw_input('----->')
	his= Ind[-1]
	

	if his is None:
		met=None
		net=None
		den=None
		Ts=None
		Vs=None
	else:
		met=D[his][1]
		net=D[his][0]
		den=D[his][2]
		Ts=D[his][-1]
		Vs=1
	if v==0:
		fName = path_to_out + 'mania_binary.net'
	else:
		fName = path_to_out + 'agg.net'
	with open(fName,'w') as f:
		D=pk.dump(net,f)
		
	if v==0:
		fName = path_to_out + 'mania_binary.meta'
	else:
		fName = path_to_out + 'agg.meta'
	with open(fName,'w') as f:
		tmp = {}
		tmp['density'] = den
		tmp['threshold'] = Ts
		tmp['normalized_asymmetry_ratio'] = met
		D=pk.dump(tmp,f)
	
	
	
	
