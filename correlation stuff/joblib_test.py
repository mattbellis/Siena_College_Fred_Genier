import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import seaborn as sn
import time
from joblib import Parallel,delayed
sn.set_style("ticks",sn.axes_style({'axes.grid': True}))

def distance(a,b):
    d = np.sqrt( (b[0]-a[0])**2 + (b[1]-a[1])**2 + (b[2]-a[2])**2 )
    return d

random = np.loadtxt('10k_weighted_random.dat')
data = np.loadtxt('10k_weighted_north_cmass.dat')

dt=data.transpose()
rt=random.transpose()

from spherical_to_cartesian import spherical_to_cartesian as s2c

#import astropy
#from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
cd = cosmo.comoving_distance(dt[2])
cdist = cd.value * 0.7

rcd = cosmo.comoving_distance(rt[2])
rcdist = rcd.value * 0.7

d_conv=s2c(dt)

r_conv=s2c(rt)

ndivs = 4

rxpts = np.linspace(np.round(min(r_conv[0])),np.round(max(r_conv[0])),ndivs)
rypts = np.linspace(np.round(min(r_conv[1])),np.round(max(r_conv[1])),ndivs)
rzpts = np.linspace(np.round(min(r_conv[2])),np.round(max(r_conv[2])),ndivs)

dxpts = np.linspace(np.round(min(d_conv[0])),np.round(max(d_conv[0])),ndivs)
dypts = np.linspace(np.round(min(d_conv[1])),np.round(max(d_conv[1])),ndivs)
dzpts = np.linspace(np.round(min(d_conv[2])),np.round(max(d_conv[2])),ndivs)

ddivs_master = {}

xwidth=dxpts[1]-dxpts[0]
ywidth=dypts[1]-dypts[0]
zwidth=dzpts[1]-dzpts[0]

x,y,z=d_conv[0],d_conv[1],d_conv[2]

ddivs = []

for a,i in enumerate(dxpts[range(ndivs)]):

    index_x0 = x>i
    index_x1 = x<= i+xwidth

    for b,j in enumerate(dypts[range(ndivs)]):

        index_y0 = y>j
        index_y1 = y<= j+ywidth
        
        for c,k in enumerate(dzpts[range(ndivs)]):
                  
            index_z0 = z>k
            index_z1 = z<= k+zwidth
            
            index = index_x0*index_x1 * index_y0*index_y1 * index_z0*index_z1
            
            xsub = x[index]
            ysub = y[index]
            zsub = z[index]
            
            ddivs.append([xsub,ysub,zsub])

            key = "%02d%02d%02d" % (a,b,c)
            ddivs_master[key] = [xsub,ysub,zsub]
            
rdivs_master = {}

xwidth=rxpts[1]-rxpts[0]
ywidth=rypts[1]-rypts[0]
zwidth=rzpts[1]-rzpts[0]

x,y,z=r_conv[0],r_conv[1],r_conv[2]

rdivs = []

for a,i in enumerate(rxpts[range(ndivs)]):

    index_x0 = x>i
    index_x1 = x<= i+xwidth
    
    #print(i,ndivs)

    for b,j in enumerate(rypts[range(ndivs)]):

        index_y0 = y>j
        index_y1 = y<= j+ywidth
        
        for c,k in enumerate(rzpts[range(ndivs)]):
                  
            index_z0 = z>k
            index_z1 = z<= k+zwidth
            
            index = index_x0*index_x1 * index_y0*index_y1 * index_z0*index_z1
            
            xsub = x[index]
            ysub = y[index]
            zsub = z[index]
            
            #rdivs.append([xsub,ysub,zsub])

            key = "%02d%02d%02d" % (a,b,c)
            rdivs_master[key] = [xsub,ysub,zsub]

def worker1(dpt):
    if __name__ == '__main__':
        pos1=d_home_subdiv[dpt]
        pos2 = None
        if nni==i and nnj==j and nnk==k:
            #print("here")
            pos2=d_nn_subdiv[dpt+1:].transpose()
        else:
            #print("there")
            pos2=d_nn_subdiv[:].transpose()


        # DD
        d_dist=distance(pos1,pos2)
        
        # DR
        dr_dist=distance(pos1,r_nn_subdiv[:].transpose())
        
        return d_dist,dr_dist
        
        
def worker2(rpt):
    if __name__ == '__main__':
        pos1=r_home_subdiv[rpt]
        pos2 = None
        if nni==i and nnj==j and nnk==k:
            pos2=r_nn_subdiv[rpt+1:].transpose()
        else:
            pos2=r_nn_subdiv[:].transpose()
            
            # RR
        r_dist=distance(pos1,pos2)
        
        # DR FOR ALL OF THEM, but don't double count!
        if not (nni==i and nnj==j and nnk==k):
            dr_dist=distance(pos1,d_nn_subdiv[:].transpose())
            
            return r_dist,dr_dist
            
        return r_dist
            
    

######################################################################################
# nearest neighbor
######################################################################################

nsubs = ndivs

nbins = 50

hist_dd = np.zeros(nbins)
hist_dr = np.zeros(nbins)
hist_rr = np.zeros(nbins)

start = time.time()

rsep = None
for i in range(nsubs):
    
    print("i: ",i)
    
    for j in range(nsubs):
        
        for k in range(nsubs):
            
            home = "%02d%02d%02d" % (i,j,k)
            #print(home)

            d_home_subdiv = np.array(ddivs_master[home]).transpose()
            r_home_subdiv = np.array(rdivs_master[home]).transpose()

            
            for nni in range(i,i+2):
                
                for nnj in range(j,j+2):
                    
                    for nnk in range(k,k+2):
                        
                        if nni<nsubs and nnj<nsubs and nnk<nsubs:
                            nn = "%02d%02d%02d" % (nni,nnj,nnk)
                            #print("\t",nn)
                            
                            d_nn_subdiv = np.array(ddivs_master[nn]).transpose()
                            r_nn_subdiv = np.array(rdivs_master[nn]).transpose()
                            
                            #print("dlen: ",len(d_nn_subdiv))
                            #print("rlen: ",len(r_nn_subdiv))

                            #hdd = []
                            #hdr = []
                            #hrr = []
                            
                            results=Parallel(n_jobs=4,verbose=11)(delayed(worker1)(dpt)for dpt in range(0,len(d_home_subdiv)))
                            res_t = np.array(results).transpose() 

                            #h_dr = np.histogram(res_t[1].transpose(),bins=nbins,range=(0,200))
                            #hist_dr += h_dr[0]
                            
                                
                            #h_dd = np.histogram(res_t[0].transpose(),bins=nbins,range=(0,200))
                            #hist_dd += h_dd[0]
                            
                            results=Parallel(n_jobs=4,verbose=11)(delayed(worker2)(rpt)for rpt in range(0,len(r_home_subdiv)))
                            res_t = np.array(results).transpose()
                            
                            # Histogram after getting the numbers
                            #h_rr = np.histogram(res_t[0].transpose(),bins=nbins,range=(0,200))
                            #hist_rr += h_rr[0]

                            #h_dr = np.histogram(res_t[1].transpose(),bins=nbins,range=(0,200))
                            #hist_dr += h_dr[0]

                            
                            

                                
end=time.time()
print("Done! Took %.2f seconds" % (end-start))

#plt.hist(hist_dd)
#plt.show()