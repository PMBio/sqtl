import scipy as SP
import time as time
import pdb
from sqtl.tools.common import *
from sqtl.tools.io import read_counts, write_posterior
import time

STR_LOG_SMOOTH = "smooth - Sample %s, chrm %s. n_rounds=%d, n_sites=%d, rec_rate=%.2e, rec_cutoff=%.2f (length_cuoff=%d), outlier_cutoff=%.2f, nearby_snp_cutoff=%d bp, median_coverage=%d, coverage_cap=%d include_bad=%s"
STR_HEADER_SMOOTH = """#sQTL version %s. smooth samples
# Samples=%s
# number of inference rounds=%d
# recombination rate=%.2e (events/bp)
# recombination rate cutoff=%.2f
# nearby SNP cutoff=%d bp (min distance between nearby sites)
# coverage_cap=%.1f
# 
"""


""" Filter list of locations that are nearby, but not too close
@param locs length-L list of integer site locations
@param l site
@param length_cutoff maximum distance from l
@param nearby_snp_cutoff minimum distance between nearby snps
@param bad length-L list of booleans of whether the sites are of good quality
@return list I of indexes such that locs[I[i+1]] - locs[I[i]] >= nearby_snp_cutoff and abs(locs[I[i]] - locs[l]) <= length_cutoff and not bad[I[i]]
"""
def filter_nearby(locs, l, length_cutoff, nearby_snp_cutoff, bad):
    # 1. find all loci that are no more than cutoff away
    start, end = l, l
    while start >= 0 and abs(locs[l] - locs[start]) < length_cutoff: start -= 1
    while end < len(locs) and abs(locs[l] - locs[end]) < length_cutoff: end += 1

    # 2. Filter to retain only good ones that are
    last = start
    while last < len(locs) and bad[last]: last += 1 # start off with at least one good site
    res = []    
    for i in range(last+1, end):
        if ((locs[i] - locs[last]) > nearby_snp_cutoff) and (~bad[i]):
            last = i
            res.append(i)

    return res


def smooth(data, loc, chrm_locs={}, sample="", chrm="", n_rounds=2, rec_rate=80./(12000000), rec_cutoff=0.9, outlier_cutoff=0.2, nearby_snp_cutoff=200, max_num_median_coverage=2., fixed_weight=0.05, include_bad=False):
    length_cutoff, coverage, coverage_cap = -SP.log(rec_cutoff)/rec_rate, data.sum(axis=1), SP.median(data.sum(axis=1))*max_num_median_coverage
    LOG.info(STR_LOG_SMOOTH%(sample, chrm, n_rounds, len(data), rec_rate, rec_cutoff, length_cutoff, outlier_cutoff, nearby_snp_cutoff, SP.median(coverage), coverage_cap, include_bad))
    res, fixed, bad = SP.zeros([data.shape[0], 2]), (data.prod(axis=1) == 0), SP.zeros(data.shape[0], bool) # only one allele observed, and outliers compared to mean
    
    for n in range(n_rounds): # for a number of rounds to hone bad loci and outliers
        for i in range(len(data)): # smooth each locus
            res[i,:] = 0
            I = chrm_locs[i] if i in chrm_locs else filter_nearby(loc, i, length_cutoff, nearby_snp_cutoff, bad)
            chrm_locs[i] = I
            for j in I: # add together information across all strain SNPs that are not too close
                d = data[j]
                d = d*(min(coverage_cap, d.sum())/(d.sum() + EPS)) # rescale observations so that no more than total of coverage_cap alleles are observed
                d = data[j]*(fixed_weight**(fixed[j])) + 1 # if allele fixed, divide coverage by 20
                res[i] += SP.exp(-rec_rate*abs(loc[i] - loc[j]))*d   # if not high, just add all total mapped bases
            
        init_mean = data[:,0]/(data.sum(axis=1) + EPS)
        post_mean = res[:,0]/(res.sum(axis=1) + EPS)
        post_var = (res.prod(axis=1) + EPS)/((res.sum(axis=1)**2 + EPS)*(res.sum(axis=1) + 1))
        bad = (abs(init_mean - post_mean) > outlier_cutoff) | (data.sum(axis=1) == 0)
    return (init_mean, post_mean, post_var, res, bad, loc, coverage), chrm_locs


""" Perform genome-wide smoothing of samples in the count_file, and output them to out_file
@param out_file path to output file
@param count_file path to input file of allele counts
@param args arguments for the smooth() function
@effects creates out_file with reference and nonreference counts, as well as maximum likelihood allele frequency, posterior allele frequency, posterior standard deviation, and bad locus for each sample"""
def smooth_samples(out_file, count_file, **args):
    LOG.info("smooth_samples - out_file=%s count_file=%s args=%s"%(out_file, count_file, str(args)))
    counts = read_counts(count_file)
    
    # 1. Calculate smoothed values
    means, i, total_chrms = {}, 1, len(counts)*len(counts.values()[0])
    for chrm in counts.values()[0]: # ["ref|NC_001133|", "ref|NC_001138|"]: #
        chrm_locs = {}
        for sample in counts:
            if sample not in means: means[sample] = {}
            means[sample][chrm], chrm_locs = smooth(SP.array(counts[sample][chrm]['D']), SP.array(counts[sample][chrm]['L']), chrm_locs, sample, chrm + " (%d of %d total)"%(i, total_chrms), **args)
            i += 1

    # 2. Output
    header = STR_HEADER_SMOOTH%(SQTL_VERSION, str(counts.keys()), args['n_rounds'], args['rec_rate'], args['rec_cutoff'], args['nearby_snp_cutoff'], args['max_num_median_coverage'])
    write_posterior(out_file, header, means, counts)
    LOG.info("smooth_samples - done")
