import numpy as np
from joblib import Parallel,delayed
import time

npts = 20000
nbins = 10
lo = 0
hi = 1

x1 = np.random.random(npts)
x2 = np.random.random(npts)

x2chunks = []
nprocesses = 8
for i in range(nprocesses):
    x2chunks.append(x2.copy())


################################################################################
# NORMAL
################################################################################
def dist(a,b):
    dists = []
    for ival in a:
        dx = ival-b
        dx2 = dx*dx
        dist = np.sqrt(dx2)
        dists += dist.tolist()
    hdists = np.histogram(dists,bins=nbins,range=(lo,hi))[0]
    return hdists
################################################################################

start = time.time()
# Normal distance calculations
ndists = dist(x1,x2)
print("NORMAL run time:   ",time.time()-start)

#print(ndists)
#print(len(ndists))


################################################################################
# Parallel
################################################################################
def worker(i):
    #x2index = i%nprocesses
    #print(x2index)

    dists = []
    ival = x1[i]
    #dx = ival-x2chunks[x2index]
    dx = ival-x2
    dx2 = dx*dx
    dist = np.sqrt(dx2)
    dists += dist.tolist()
    #return dists
    hdists = np.histogram(dists,bins=nbins,range=(lo,hi))[0]
    return hdists
################################################################################

start = time.time()
# Parallel distance calculations
pdists = Parallel(n_jobs=nprocesses)(delayed(worker)(i) for i in range(len(x1)))
print("PARALLEL run time: ",time.time()-start)

#print(pdists)
#print(len(pdists))

totpdists = np.zeros(nbins,dtype=int)
for p in pdists:
	totpdists += p

#print(len(totpdists))
print(totpdists)
print(ndists)

print(sum(ndists))
print(sum(totpdists))
