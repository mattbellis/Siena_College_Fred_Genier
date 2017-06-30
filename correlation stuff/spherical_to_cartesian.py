def spherical_to_cartesian(ar,alt_cmd=None):
    
    # Requires array with 4 columns of arbitrary length
    # If you don't have a comoving distance in your array, set alt_cmd to
    # an array of comoving distances of the same length as the inumpyut array
    
    import numpy
    
    newpts=[]
    
    if alt_cmd is None:
        
        for i in range(len(ar[0])):

            ra = numpy.deg2rad(ar[0,i])
            dec = numpy.deg2rad(ar[1,i])
            cmd=ar[3,i]
 
            ############################
    
            x=cmd*numpy.cos(ra)*numpy.cos(dec)
            y=cmd*numpy.sin(ra)*numpy.cos(dec)
            z=cmd*numpy.sin(dec)

            ############################
            xyz=[x,y,z]

            newpts.append(xyz)
            
    else:
            
        for i in range(len(ar[0])):
              
            ra = numpy.deg2rad(ar[0,i])
            dec = numpy.deg2rad(ar[1,i])
            cmd=alt_cmd[i]

            ############################
   
            x=cmd*numpy.cos(ra)*numpy.cos(dec)
            y=cmd*numpy.sin(ra)*numpy.cos(dec)
            z=cmd*numpy.sin(dec)
   
            ############################
    

            xyz=[x,y,z]
 
            newpts.append(xyz)
   

    return numpy.array(newpts).transpose()

