def get_collisions(infile):
    
    with h5py.File('cleo_test.hdf5','r') as hf:

        npions = hf.get('npions')
        nkaons = hf.get('nkaons')
        nmuons = hf.get('nmuons')
        nelecs = hf.get('nelectrons')
        nphots = hf.get('nphotons')

        pions = hf.get('particles/pions')
        kaons = hf.get('particles/kaons')
        muons = hf.get('particles/muons')
        elecs = hf.get('particles/electrons')plt.hist(masses,bins=50)
        phots = hf.get('particles/photons')


        nlist = [npions,nkaons,nmuons,nelecs,nphots]
        plist = [pis,kas,mus,els,phs]


        for i in range(0,len(nlist)):

            n_particles = nlist[i]
            particles = plist[i]
            n_start = 0

            for j in n_particles:

                n_end = n_start + j
                event_particles = particles[n_start:n_end]

                for k in event_particles: