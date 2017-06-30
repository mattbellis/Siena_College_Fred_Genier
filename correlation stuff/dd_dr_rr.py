def distance(a,b):
    import numpy as np
    d = np.sqrt( (b[0]-a[0])**2 + (b[1]-a[1])**2 + (b[2]-a[2])**2 )
    return d

def dd_dr_rr(_data_,_random_,nsubs,nbins=50):
    
    import numpy as np 
    import time
    from joblib import Parallel, delayed
    import multiprocessing
    
    num_cores = multiprocessing.cpu_count()
    
    def dd_dr(d_home,d_nn,r_nn,ind):
        
        pos1=d_home[ind]
        pos2 = None
                                    
        if nni==i and nnj==j and nnk==k:
                                        
            pos2=d_nn[ind+1:].transpose()
        else:
                                        
            pos2=d_nn[:].transpose()
                                    
        d_dist=distance(pos1,pos2).tolist()
    
        dr_dist=distance(pos1,r_nn[:].transpose()).tolist()
        
        return d_dist,dr_dist
    
    def rr_dr(r_home,d_nn,r_nn,ind):
        pos1=r_home[ind]
        pos2 = None
        
        if nni==i and nnj==j and nnk==k:
            
            pos2=r_nn[ind+1:].transpose()
        else:
            
            pos2=r_nn[:].transpose()
            
        # RR
        r_dist=distance(pos1,pos2)
    
        if not (nni==i and nnj==j and nnk==k):
            
            dr_dist=distance(pos1,d_nn[:].transpose())
            
        return r_dist,dr_dist

                                    
    
   # nsubs = len(_data_)-1
    
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
    
                d_home_subdiv = np.array( _data_[home] ).transpose()
                r_home_subdiv = np.array( _random_[home] ).transpose()
                
                for nni in range(i,i+2):
                    
                    for nnj in range(j,j+2):
                        
                        for nnk in range(k,k+2):
                            
                            if nni<nsubs and nnj<nsubs and nnk<nsubs:
                                
                                nn = "%02d%02d%02d" % (nni,nnj,nnk)
                                
                                d_nn_subdiv = np.array( _data_ [nn]).transpose()
                                r_nn_subdiv = np.array( _random_ [nn]).transpose()
                                
                                hdd=[]
                                hdr=[]
                                hrr=[]
                                
                                subdiv_range =  range(len(d_home_subdiv))
                            
                                results = Parallel(n_jobs=num_cores)(delayed(dd_dr)(d_home_subdiv,d_nn_subdiv,r_nn_subdiv,dpt) for dpt in subdiv_range)
                                
                                
                                
                                if len(results)>1:
                                
                                    hdd,hdr = results[0],results[1]
                                    
                                    h_dd = np.histogram(hdd,bins=nbins,range=(0,200))
                                    hist_dd += h_dd[0]
                                
                                    h_dr = np.histogram(hdr,bins=nbins,range=(0,200))
                                    hist_dr += h_dr[0]
                                
                                    rsep = h_dd[1][0:-1]
                                    print(results,type(results))
                                
                                subdiv_range =  range(len(r_home_subdiv))
                                results = Parallel(n_jobs=num_cores)(delayed(rr_dr)(r_home_subdiv,d_nn_subdiv,r_nn_subdiv,rpt) for rpt in subdiv_range)
                                
                                if len(results)>1:
                                    hrr,hdr=results[0],results[1]
    
                                    # Histogram after getting the numbers
                                    h_rr = np.histogram(hrr,bins=nbins,range=(0,200))
                                    hist_rr += h_rr[0]
    
                                    h_dr = np.histogram(hdr,bins=nbins,range=(0,200))
                                    hist_dr += h_dr[0]
    
    
                                    
    end=time.time()
    print("Done! Took %.2f seconds" % (end-start))
    
    return hist_dd,hist_dr,hist_rr,rsep

