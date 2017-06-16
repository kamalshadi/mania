# This fils hosts the main functions to generate binary anatomical network from
# probtrackx output

from utils import *
from optimization_tools import finalNet
import numpy as num
import pickle as pk
from tqdm import tqdm
import os


def parseSubj(filepath, nroi):
	# this function parses the 'matrix_seeds_to_all_targets' file from probtrackx output
	W = {}
	for i in range(nroi):
		q = str(i+1)
		ad = '/'+filepath.strip('/')+'/probtrackx/' + q + '.s2t'
		with open(ad) as f:
			W1=[]
			for line in f:
				tmp=line.strip().split(' ')
				vox=[int(float(xx)) for xx in tmp if len(xx)>0]
				if len(vox)!=nroi:
					print 'Voxel rejected for not having all ROI data'
					continue
				W1.append(vox)
				W[q]=num.array(W1)
	return W

def vox_roi(W,th):
	K={}
	for idk in W.keys():
		tup=num.shape(W[idk])
		tmp=num.zeros(tup)
		tmp[W[idk] > th]=1
		K[idk]=tmp
	return K

def roi_roi(W):
	com=W.keys()[0].split('.')[-1]
	lis1=[int(xx.split('.')[0]) for xx in W.keys()]
	lis=sorted(lis1)
	l=len(lis)
	net=num.zeros((l,l))
	for i,idk1 in enumerate(lis):
		idk=str(idk1)
		K=W[idk]
		tmp=num.sum(K,0)
		net[i,:]=[1 if xx > 0 else 0 for xx in tmp]
	return net
		
def brainMap(W,fName,ns):
	nroi=len(W)
	B=num.zeros((nroi,nroi))
	tpe=nroi*nroi-nroi	# total possible edges
	res=[()]*(ns-1)
	q = 0
	for T in tqdm(reversed(range(1,ns)), total = ns-1):
			K=vox_roi(W,T)
			BM=roi_roi(K)
			den=density(BM)
			met=NAR(BM)
			res[q]=(BM,met,den,T)
			q = q + 1

	with open(fName,'w') as f:
		pk.dump(res,f)
		
def mania_on_subject(fp, ns = 5000):
	print 'loading data for '+ fp
	W=parseSubj(fp, 18)
	if fp[-1] != '/':
		fp = fp + '/'
	path_to_out = fp + '/MANIA/'
	fName = path_to_out + 'raw.res'
	if not os.path.isdir(path_to_out):
		os.mkdir(path_to_out)
	print 'Scanning the threshold space for ' + fp
	brainMap(W,fName,ns)
	print 'Running the optimization for ' + fp
	finalNet(fp)
	
	
	
#~ if __name__ == '__main__':
	#~ mania_on_subject('/Users/kshadi/Desktop/MANIA_pack/S1', ns = 5000)
	
