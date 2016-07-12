################################################################################

def get_collisions_from_filename(infilename,verbose=False):

    infile = None
    if zipfile.is_zipfile(infilename) is True:
        z = zipfile.ZipFile(infilename,'r')
        infile = z.open(z.namelist()[0],'r')
    else:
        infile = open(infilename)

    collisions = get_collisions(infile,verbose)

    return collisions

################################################################################

def get_collisions(infile, filetype=None, returntype=None):
    
    # File type can be:
    #    text (default)
    #    text zipped 
    #
    #    npy
    #    npz
    #    npz_compressed
    #
    #    hdf5
    #    hdf5 compressed
    #
    # Return type can be:
    #    list of lists (default)
    #    dictionary of lists/arrays
    #    dictionary of dictionaries 
    #    rec arrary (???????? only given time)
    
    ################################################################################
    
    if filetype == None or filetype == "txt":
        
        collisions = []

        not_at_end = True
        collision_count = 0
        new_collision = True
        while ( not_at_end ):
            
            # Read in one collision
            
            line = infile.readline()

            if collision_count%1000==0 and verbose:
                print "collision count: ",collision_count

            if line=="":
                not_at_end = False

            if line.find("Event")>=0:
                new_collision = True

            if new_collision==True:

                # Read in the pion info for this collision.
                pions = []
                line = infile.readline()
                npions = int(line)
                for i in xrange(npions):
                    line = infile.readline()
                    vals = line.split()
                    e = float(vals[0])
                    px = float(vals[1])
                    py = float(vals[2])
                    pz = float(vals[3])
                    q = int(vals[4])
                    sigpi = float(vals[5])
                    sigka = float(vals[6])
                    likpi = float(vals[7])
                    likka = float(vals[8])
                    nphopi = int(vals[9])
                    nphoka = int(vals[10])
                    depthmu = float(vals[11])
                    cluster_energy = float(vals[12])
                    pions.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

                # Read in the kaon info for this collision.
                kaons = []
                line = infile.readline()
                nkaons = int(line)
                for i in xrange(nkaons):
                    line = infile.readline()
                    vals = line.split()
                    e = float(vals[0])
                    px = float(vals[1])
                    py = float(vals[2])
                    pz = float(vals[3])
                    q = int(vals[4])
                    sigpi = float(vals[5])
                    sigka = float(vals[6])
                    likpi = float(vals[7])
                    likka = float(vals[8])
                    nphopi = int(vals[9])
                    nphoka = int(vals[10])
                    depthmu = float(vals[11])
                    cluster_energy = float(vals[12])
                    kaons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

                # Read in the muon info for this collision.
                muons = []
                line = infile.readline()
                nmuons = int(line)
                for i in xrange(nmuons):
                    line = infile.readline()
                    vals = line.split()
                    e = float(vals[0])
                    px = float(vals[1])
                    py = float(vals[2])
                    pz = float(vals[3])
                    q = int(vals[4])
                    sigpi = float(vals[5])
                    sigka = float(vals[6])
                    likpi = float(vals[7])
                    likka = float(vals[8])
                    nphopi = int(vals[9])
                    nphoka = int(vals[10])
                    depthmu = float(vals[11])
                    cluster_energy = float(vals[12])
                    muons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])

                # Read in the electron info for this collision.
                electrons = []
                line = infile.readline()
                nelectrons = int(line)
                for i in xrange(nelectrons):
                    line = infile.readline()
                    vals = line.split()
                    e = float(vals[0])
                    px = float(vals[1])
                    py = float(vals[2])
                    pz = float(vals[3])
                    q = int(vals[4])
                    sigpi = float(vals[5])
                    sigka = float(vals[6])
                    likpi = float(vals[7])
                    likka = float(vals[8])
                    nphopi = int(vals[9])
                    nphoka = int(vals[10])
                    depthmu = float(vals[11])
                    cluster_energy = float(vals[12])
                    electrons.append([e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy])


                # Read in the photon info for this collision.
                photons = []
                line = infile.readline()
                nphotons = int(line)
                for i in xrange(nphotons):
                    line = infile.readline()
                    vals = line.split()
                    e = float(vals[0])
                    px = float(vals[1])
                    py = float(vals[2])
                    pz = float(vals[3])
                    photons.append([e,px,py,pz])


                # Read in the information about the missing transverse energy (MET) in the collision.
                # This is really the x and y direction for the missing momentum.
                #line = infile.readline()
                #vals = line.split()
                #met_px = float(vals[0])
                #met_py = float(vals[1])

                new_collision = False
                collision_count += 1

                collisions.append([pions,kaons,muons,electrons,photons])

        return collisions
    
    ################################################################################
    
    elif filetype=="hdf5":

        with h5py.File(infile,'r') as hf:

            npions = hf.get('npions').value
            pions = hf.get('particles/pions').value

            nkaons = hf.get('nkaons').value
            kaons = hf.get('particles/kaons').value

            nmuons = hf.get('nmuons').value
            muons = hf.get('particles/muons').value

            nelectrons = hf.get('nelectrons').value
            electrons = hf.get('particles/electrons').value

            nphotons = hf.get('nphotons').value
            photons = hf.get('particles/photons').value

            h5_events = []
            pion = []
            kaon = []
            muon = []
            electron = []
            photon = []

            nlist = [npions,nkaons,nmuons,nelectrons,nphotons]
            plist = [pions,kaons,muons,electrons,photons]
            elist = [pion,kaon,muon,electron,photon]

            for i in range(0,len(nlist)):

                n_particles = nlist[i]
                particles = plist[i]
                outlist = elist[i]

                n_start = 0

                for j in n_particles:

                    n_end = n_start + j

                    event_particles = particles[n_start:n_end]

                    for k in event_particles:

                        ep = []
                        ep.append(k)

                    outlist.append(ep)

                    n_start = n_end

            for l in range(0,len(pion)):

                ev = []

                ev.append(pion[l])
                ev.append(kaon[l])
                ev.append(muon[l])
                ev.append(electron[l])
                ev.append(photon[l])

                h5_events.append(ev)
                
    elif filetype == "npy":
        
        npy = np.load(infile)
        
    elif filetype == "npz":
        
        npz = np.load(infile)['arr_0']
        
        
        
        
        
        
        
        
        
           