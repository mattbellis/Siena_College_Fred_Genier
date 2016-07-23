import zipfile

################################################################################

def is_zip_txt(infilename,verbose=False):
    

    infile = None
    if zipfile.is_zipfile(infilename) is True:
        print "Unzipping!!!!"
        z = zipfile.ZipFile(infilename,'r')
        infile = z.open(z.namelist()[0],'r')
    else:
        infile = open(infilename)

    return infile

################################################################################

def returntype(var):
    
    if returntype == None or returntype == "list":
        
        return var
    
    elif returntype == "dictionary": # dictionary of dictionaries
        
        d_events = []
        
        for i in var:
    
            pions,kaons,muons,electrons,photons = i

            pis = []
            kas = []
            mus = []
            els = []
            phs = []
            
            
            for pion in pions:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=pion
                pidi = {"e":e,"px":px,"py":py,"pz":pz,"q":q,"likpi":likpi,"likka":likka,"nphoka":nphoka,"depthmu":depthmu,"cluster_energy":cluster_energy}
                pis.append(pidi)

            for kaon in kaons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=kaon
                kadi = {"e":e,"px":px,"py":py,"pz":pz,"q":q,"likpi":likpi,"likka":likka,"nphoka":nphoka,"depthmu":depthmu,"cluster_energy":cluster_energy}
                kas.append(kadi)

            for muon in muons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=muon
                mudi = {"e":e,"px":px,"py":py,"pz":pz,"q":q,"likpi":likpi,"likka":likka,"nphoka":nphoka,"depthmu":depthmu,"cluster_energy":cluster_energy}
                mus.append(mudi)

            for electron in electrons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=electron
                eldi = {"e":e,"px":px,"py":py,"pz":pz,"q":q,"likpi":likpi,"likka":likka,"nphoka":nphoka,"depthmu":depthmu,"cluster_energy":cluster_energy}
                els.append(eldi)

            for photon in photons:
                e,px,py,pz=photon
                phdi = {"e":e,"px":px,"py":py,"pz":pz}
                phs.append(phdi)

            ev = {"pions":pis,"kaons":kas,"muons":mus,"electrons":els,"photons":phs}

            d_events.append(ev)
            
        return d_events
    
    elif returntype == "dictionary_list": # dictionary containing lists
        
        d_events = []
        
        for i in var:
    
            pions,kaons,muons,electrons,photons = i

            pis = []
            kas = []
            mus = []
            els = []
            phs = []
            
            
            for pion in pions:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=pion
                
                pis.append(pion)

            for kaon in kaons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=kaon
                
                kas.append(kaon)

            for muon in muons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=muon
                
                mus.append(muon)

            for electron in electrons:
                e,px,py,pz,q,sigpi,sigka,likpi,likka,nphopi,nphoka,depthmu,cluster_energy=electron
                
                els.append(electron)

            for photon in photons:
                e,px,py,pz=photon
                
                phs.append(photon)

            ev = {"pions":pis,"kaons":kas,"muons":mus,"electrons":els,"photons":phs}

            d_events.append(ev)
            
        return d_events
    
    else:
        
        print "ERROR: Invalid returntype. Valid returntypes include: 'list' (default), 'dictionary', and 'dictionary_list'."

################################################################################

def get_collisions(infilename, filetype=None, returntype=None):
    
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
        
        infile = is_zip_txt(infilename) # this may raise an issue -  needs testing
        
        collisions = []

        not_at_end = True
        collision_count = 0
        new_collision = True
        while ( not_at_end ):
            
            # Read in one collision
            
            line = infile.readline()
            
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
        
        if returntype=="lists" or returntype==None:
            return collisions
        
        #returntype(collisions)
    
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
                
        returntype(h5_events)
                
    elif filetype == "npy":
        
        npy = np.load(infile)
        
        returntype(npy)
        
    elif filetype == "npz":
        
        npz = np.load(infile)['arr_0']
        
        returntype(npz)
   
    else:
        
        print "ERROR: Invalid filetype. Valid filetypes include: 'txt' (default), 'hdf5', 'npy', and 'npz'."
        
        
        
        
        
        
        
        
           