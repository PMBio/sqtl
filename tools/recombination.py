




def get_approx_map(chr,start,end):
    locs, rates = cl("/Users/leopold/data/projects/fourway/genotype/average-rec_%s.pickle"%chr) # chr is in roman form, e.g. X
    locs = SP.array(locs)
    I = SP.where(abs(locs - 0.5*(start+end)) < 0.5*(end-start))[0]
    cumr = SP.cumsum(rates[I])
    R = SP.zeros([len(I), len(I)])
    
    for i in range(len(I) - 1):
        for j in range(i+1, len(I)):
            R[i,j] = R[j,i] = cumr[j] - cumr[i]

    return locs[I], R
