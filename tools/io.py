import os
import pdb
import scipy as SP
#from rapid_qtl_yeast.tools.frequency import *
#from sqtl.model.stats_helper import *
from sqtl.tools.common import *

MEAN_GEN_REC_RATE = (66-30)/(5.*1.2E7)
OUT_DIR = "/Users/leopold/data/projects/sqtl/output"

def get_sift(chrm, loci):
    result = SP.zeros(len(loci))
    ifh = file('sift_predictions_complete.tab','r')
    header = ifh.next()
    scores = SP.array([SP.array(l.split('\t'))[[7,8,19]] for l in ifh])
    sift_locs = SP.array(scores[:,1], int)
    chr = SP.array([10*(ord(x[3]) - ord('0')) + ord(x[4]) - ord('0') for x in scores[:,0]], int)
    score = scores[:,2]
    Ic = SP.where(chr == chrm)[0]

    for i,l in enumerate(loci):
        for j in Ic:
            if l == sift_locs[j]:
                result[i] = 1. - float(score[j])

    result += 0.01
    return result/result.sum()


def get_prior(chrm, loci):
    Q = SP.ones(len(loci))
    Q += get_sift(chrm, loci)
    return Q/sum(Q)


""" chrm is in 1..16, loci are in base pairs"""
""" return """
def get_genetic_map(chrm, loci):
    ifh = file("%s/map/2way_v4.All.chr%d_0.5"%(DATA_DIR, chrm),"r")
    res = SP.zeros(100000)
    for l in ifh:
        d = l.strip().split()
        res[int(2.*float(d[0]) - 0.49)] = 12.*float(d[2])/1000. # median r per segregant per base pair per generation; location is in kb, indexes into slices of 500

    r = SP.zeros(len(loci)) # fill in recombination rates. r[i] = r between i-1 and i
    for i in range(len(r)-1):
        # check if the loci are within the same region
        if loci[i]/500 == loci[i+1]/500:
            r[i+1] = (loci[i+1] - loci[i])*res[loci[i]/500]
        else:
            r[i+1] = (loci[i+1] - 500*(loci[i+1]/500))*res[loci[i+1]/500] # bit from last slice beginning to loci[i+1]
            r[i+1] += (500*(loci[i]/500 + 1) - loci[i])*res[loci[i]/500] # bit from loci[i] to first slice end
            for j in range(loci[i]/500 + 1, loci[i+1]/500): # intermediate bits
                r[i+1] += 500*res[j]
    return r


def get_all_counts(sample):
    ifh = file("/Users/lp2/prog/projects/rapid_qtl_yeast/analysis/allele_freq/%s_count.tab"%sample, 'r')
    ifh.next() # get rid of header
    res = {}
    for l in ifh:
        chr, loc, n1, n2 = l.strip().split()
        if chr not in res: res[chr] = []
        res[chr].append((loc,n1,n2))
    for c in res: res[c] = SP.array(res[c], int)
    return res


def get_counts(sample, chrm, start, end):
    all_chr = get_all_counts(sample)[chrm]
    I = SP.where((all_chr[:,0] >= start) & (all_chr[:,0] <= end))[0]
    return all_chr[I,:]


def read_counts(file_name):
    res = {} # map of samples->chrms->(locs, seq, data)
    samples = []
    ifh = file(file_name, 'r')
    for line in ifh:
        if line[0] == "#": 
            if line[0:4] == "#Chr": # if header
                samples = [x.replace("_ref","") for x in line.strip().split()[4:-1:2]]
                res = dict((s,{}) for s in samples)
        else:
            d = line.strip().split()
            chrm, loc, ref_base, nonref_base = d[0:4]
            if chrm not in res.values()[0]:
                for sample in res: res[sample][chrm] = {'L':[],'seq':[],'D':[]}
            for s in range(len(samples)):
                res[samples[s]][chrm]['L'].append(int(loc))
                res[samples[s]][chrm]['seq'].append([ref_base, nonref_base])
                res[samples[s]][chrm]['D'].append([int(d[4+2*s]), int(d[5+2*s])])
    return res



def read_posterior(file_name):
    res = {} # map of samples->chrms->(locs, seq, data, posterior_mean, posterior_var, badness, coverage)
    samples = []
    ifh = file(file_name, 'r')
    for line in ifh:
        if line[0] == "#": 
            if line[0:4] == "#Chr": # if header
                samples = [x.replace("_ref","") for x in line.strip().split()[4:-1:6]]
                res = dict((s,{}) for s in samples)
        else:
            d = line.strip().split()
            chrm, loc, ref_base, nonref_base = d[0:4]
            if chrm not in res.values()[0]:
                for sample in res: res[sample][chrm] = {'L':[],'seq':[],'D':[], 'ML':[], 'AF':[], 'SD':[], 'bad':[]}
            for s in range(len(samples)):
                res[samples[s]][chrm]['L'].append(int(loc))
                res[samples[s]][chrm]['seq'].append([ref_base, nonref_base])
                res[samples[s]][chrm]['D'].append([int(float(d[4+6*s])), int(float(d[4+6*s + 1]))])
                res[samples[s]][chrm]['ML'].append(float(d[4+6*s+2]))
                res[samples[s]][chrm]['AF'].append(float(d[4+6*s+3]))
                res[samples[s]][chrm]['SD'].append(float(d[4+6*s+4]))
                res[samples[s]][chrm]['bad'].append(eval(d[4+6*s+5]))
    return res


def write_posterior(out_file, header, means, counts):
    ofh = file(out_file, 'w')
    ofh.write(header)
    ofh.write("#Chrm\tLoc\tRef\tNonref")
    for s in sorted(means):
        ofh.write("\t%s_ref\t%s_nonref\t%s_ML-AF\t%s_posterior-AF\t%s_posterior-AF-SD\t%s_loc_bad"%(s,s,s,s,s,s))
    ofh.write("\n")
    example = counts.values()[0]
    for chrm in sorted(example): # ["ref|NC_001133|", "ref|NC_001138|"]: #
        for l in range(len(example[chrm]['L'])):
            ofh.write("%s\t%d\t%s\t%s"%(chrm, example[chrm]['L'][l], example[chrm]['seq'][l][0], example[chrm]['seq'][l][1]))
            for sample in sorted(means):
                i,m,v,p,b,loc,c = means[sample][chrm]
                ofh.write("\t%d\t%d\t%.3f\t%.3f\t%.3f\t%s"%(int(counts[sample][chrm]['D'][l][0]), int(counts[sample][chrm]['D'][l][1]), i[l], m[l], v[l]**0.5, str(b[l])))
            ofh.write("\n")
    ofh.close()


def read_diffs():
    res = {} # map of sample-pair->chrm->(locs, counts, mean(diff), sd(diff))
    sample_pairs = []
    ifh = file(file_name, 'r')
    for line in ifh:
        if line[0] == "#": 
            if line[0:4] == "#Chr": # if header
                samples = [x.replace("_ref","") for x in line.strip().split()[4:-1:6]]
                res = dict((s,{}) for s in samples)
        else:
            d = line.strip().split()
            chrm, loc, ref_base, nonref_base = d[0:4]
            if chrm not in res.values()[0]:
                for sample in res: res[sample][chrm] = {'L':[],'seq':[],'D':[], 'ML':[], 'AF':[], 'SD':[], 'bad':[]}
            for s in range(len(samples)):
                res[samples[s]][chrm]['L'].append(int(loc))
                res[samples[s]][chrm]['seq'].append([ref_base, nonref_base])
                res[samples[s]][chrm]['D'].append([int(float(d[4+6*s])), int(float(d[4+6*s + 1]))])
                res[samples[s]][chrm]['ML'].append(float(d[4+6*s+2]))
                res[samples[s]][chrm]['AF'].append(float(d[4+6*s+3]))
                res[samples[s]][chrm]['SD'].append(float(d[4+6*s+4]))
                res[samples[s]][chrm]['bad'].append(eval(d[4+6*s+5]))
    return res


def write_diffs():
    pass


def read_qtls(filename):
    ifh = file(filename, 'r')
    qtls = []
    for l in ifh:
        if l[0] == "#": continue # skip comments and header
        qtls.append(l.strip("\n").split("\t"))
    return qtls


def write_qtls(out_file, qtls, header):
    ofh = file(out_file, "w")
    ofh.write(header)
    ofh.write("#Set\tChrm\tPeak\tStart\tCentre_start\tCentre_end\tEnd\tLength\tAF_peak\tSD_peak\tnumSD_peak\tCentre_genes\n")
    for q, (chrm, peak, d, s, sds, start, end, c_start, c_end, genes) in enumerate(qtls):
        ofh.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.1f\t%s\n"%("QTL-%d"%(q+1), chrm, peak, start, c_start, c_end, end, end-start, d, s, sds, ",".join(genes)))
    ofh.close()



def read_gff(file_name=""):
    pass


def read_pileup(pileup_file, out_file=None, base_qual_cutoff=20, map_qual_cutoff=20):
    if not os.path.exists(pileup_file):
        LOG.debug("No pileup file - %s"%pileup_file)
        return
    allele_index = {'A':0, 'C':1, 'G':2, 'T':3, '*':4}
    D, Q, MQ, L, refseq = [], [], [], [], [] # alleles, quality, mapping quality, loci, reference at loci
    pfh = file(pileup_file, 'r')
    index = 0

    for line in pfh: # read through pileup
        d = line.strip().split()
        if d[2] == "*": continue
        L.append([d[0], int(d[1])])  # init new site
        refseq.append(d[2].upper())
        D.append([0]*len(allele_index))
        Q.append([[] for i in range(len(allele_index))])
        MQ.append([[] for i in range(len(allele_index))])

        data, base_qual = d[8:10]         # init data values
        i, j = 0, 0

        while i < len(data): # for the mapped reads, parse alignment
            if data[i] in '+-': # insertion or deletion definition
                length = 0 # parse length from ensuing numbers
                while data[i + 1] in '0123456789' and i + 1 < len(data):
                    length *= 10
                    length += int(data[i + 1])
                    i += 1
                i += length
            elif data[i] == '$': pass # read beginning, end, and deletion markers markers
            elif data[i] == '^': i += 1
            else: # actual base
                base = data[i].upper()
                if data[i] in ',.': base = d[2].upper() # take reference as the base if ., in alignment
                if ord(base_qual[j]) - 33 > base_qual_cutoff:
                    D[index][allele_index[base]] += 1
                    Q[index][allele_index[base]].append(ord(base_qual[j]) - 33)
                j += 1 # only move qual for actual bases
            i += 1 # move on in pileup for this locus

        if sum(D[index]) == 0: 
            for x in [D,Q,MQ,L,refseq]: x.pop()     # if no bases left, don't add the base
        else: index += 1

    pfh.close()
    result = {'D': SP.array(D), 'Q':Q, 'L':L, 'refseq':refseq}  # create and return result hash
    if out_file is not None: cdm(result, out_file)
    return result



def get_peak_data(chrm, loc, sams=['WA-NA_Initial_R2_F12_T0', 'WA-NA_Heat_R2_F12_T4'], peak_window=1.5E4, gen=12):
    rec_rate=30/1.2E7 + (gen - 1)*MEAN_GEN_REC_RATE
    ds = [get_posterior(sam, chr=chrm, type='uniform', rec=0.9, filter=0.1) for sam in sams]

    combined_loci, I = combine_loci([ds[0][1], ds[1][1]])
    combined_loci = SP.array(combined_loci)
    data = [ds[i][3][1][I[i]] for i in range(len(sams))]
    Iq = SP.where(abs(combined_loci - loc) < peak_window)[0]
    dat = data[1][Iq]
    r_mean = SP.zeros(Iq.shape[0] - 1)
    for i in range(len(r_mean)): r_mean[i] = rec_rate*(combined_loci[Iq[i+1]] - combined_loci[Iq[i]])
    reverse = False
    f0 = ds[0][0][0][I[0]][Iq]
    f1 = ds[1][0][0][I[1]][Iq]
    if f0[len(f0)/2] > f1[len(f1)/2]: 
        reverse = True
        f0 = 1 - f0
        f1 = 1 - f1
    for i in range(len(Iq)): dat[i] = [dat[i, reverse],dat[i, 1 - reverse]]

    rparams = get_gamma_params(r_mean, 0.00001*SP.ones(len(r_mean)))

    return dat,  rparams, f0, get_prior(chrm, combined_loci[Iq]), combined_loci[Iq], f1
