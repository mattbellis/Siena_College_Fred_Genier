def get_collisions(infile):
    
    print "Reading in data..."
    
    start = time.time()
    
    with h5py.File(infile,'r') as hf:

        npions = hf.get('npions').value
        nkaons = hf.get('nkaons').value
        nmuons = hf.get('nmuons').value
        nelecs = hf.get('nelectrons').value
        nphots = hf.get('nphotons').value

        pions = hf.get('particles/pions').value
        kaons = hf.get('particles/kaons').value
        muons = hf.get('particles/muons').value
        elecs = hf.get('particles/electrons').value
        phots = hf.get('particles/photons').value

        #print "Read in the all the data in %f seconds" % (time.time()-start)

        nlist = [npions,nkaons,nmuons,nelecs,nphots]
        plist = [pion,kaon,muon,elec,phot]
        elist = [pions,kaons,muons,electrons,photons]

        for i in range(0,len(nlist)):

            n_particles = nlist[i]
            particles = plist[i]
            outlist = elist[i]
            n_start = 0

            for j in n_particles:

                n_end = n_start + j
                #print n_start,n_end

                event_particles = particles[n_start:n_end]
                