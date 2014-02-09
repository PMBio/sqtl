from sqtl.tools.common import *
from sqtl.tools.io import *

LOG_CALL_REGIONS = "call_qtl_regions out_file=%s, posterior_file=%s, sample_pairs=%s, af_lenient=%.2f, sd_lenient=%.1f, af_stringent=%.2f, sd_stringent=%.1f, length_cutoff=%d"
HEADER_CALL_REGIONS = """#sQTL version %s - call QTL regions output
#==========================================
# out_file=%s
# posterior_file=%s
# sample_pairs=%s
# af_lenient=%.2f
# sd_lenient=%.1f
# af_stringent=%.2f
# sd_stringent=%.1f
# length_cutoff=%d
#
"""

""" Calculate allele frequency differences and their confidences for given chromosome
@param posteriors map of sample->chrm->{data, allele frequency, locations, sds, bad}; see sqtl.tools.io read_posterior
@param chrm chromosome
@param sample_pairs list of selected_sample,unselected_sample pairs (see call_qtl_regions)
@return locs, delta, sds - length L arrays of site locations, average allele frequency difference, and average standard deviation of the difference """
def calc_diffs(posteriors, chrm, sample_pairs):
    L = len(posteriors.values()[0][chrm]['D']) # number of sites in chromosome
    locs, delta, vars, n_samples = SP.array(posteriors.values()[0][chrm]['L']), SP.zeros(L), SP.zeros(L), SP.zeros(L)

    for selected_sample, unselected_sample in sample_pairs:
        if selected_sample not in posteriors: # test preconditions
            LOG.error("Selected sample %s not found in file %s"%(selected_sample, posterior_file))
            return
        if unselected_sample not in posteriors:
            LOG.error("Unselected sample %s not found in file %s"%(unselected_sample, posterior_file))
            return

        d1, d0 = posteriors[selected_sample][chrm], posteriors[unselected_sample][chrm] # selected and unselected data
        I = SP.where(~(SP.array(d1['bad'], bool) | SP.array(d0['bad'], bool)))[0] # index of good sites
        delta[I] += (SP.array(d1['AF'], float) - SP.array(d0['AF'], float))[I] # store difference
        vars[I] += (SP.array(d1['SD'], float)**2 + SP.array(d0['SD'], float)**2)[I] # and variance for averaging
        n_samples[I] += 1 # Also remember how many sample pairs gave info at each site to average properly

    I = SP.where(n_samples == 0)[0] # bad sites
    delta[I] = vars[I] = SP.nan # are assigned SP.nan
    return locs, delta/(n_samples + 1E-6), (vars/(n_samples + 1e-6))**0.5 # return sites and averaged means and SDs


""" Call QTLs for a sequencing set
@param out_file path to output file
@param posterior_file path to input file of posterior allele frequencies
@param sample_pairs list of selected_sample,unselected_sample pairs
@require selected_sample in posterior_file sample list
@require unselected_sample in posterior_file sample list
@param af_lenient minimum allele frequency change to start a putative QTL region
@param sd_lenient minimum SD of allele frequency change to start a putative QTL region
@param af_stringent minimum allele frequency change in a putative QTL region to call it a QTL
@param sd_stringent minimum SD of allele frequency change in a putative QTL region to call it a QTL
@param length_cutoff minimum length of a putative QTL region to be called a QTL
@param peak_cutoff "relaxation" from maximum allele frequency change at the QTL peak that is considered a candidate region for the causal gene
"""
def call_qtl_regions(out_file, posterior_file, sample_pairs, af_lenient=0.1, sd_lenient=3, af_stringent=0.15, sd_stringent=5, length_cutoff=1000, peak_cutoff=0.03):
    LOG.info(LOG_CALL_REGIONS%(out_file, posterior_file, str(sample_pairs), af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff))
    posteriors = read_posterior(posterior_file)
    sample_pairs = [x.split(",") for x in sample_pairs]
    qtls = []

    # First, call QTLs
    for chrm in sorted(posteriors.values()[0]): # for each chromosome
        locs, delta, sd = calc_diffs(posteriors, chrm, sample_pairs) # calculate allele frequency differences between all sample pairs
        lenient = ((abs(delta) > af_lenient) & (abs(delta)/sd > sd_lenient)) | SP.isnan(delta) # and whether the differences satisfy the stringent
        stringent = ((abs(delta) > af_stringent) & (abs(delta)/sd > sd_stringent)) | SP.isnan(delta) # and lenient criteria for large change.
        
        for i in range(len(locs)): # scan loci left to right
            if lenient[i] and (not SP.isnan(delta[i])): # if change at the locus leniently significant by all measures
                is_qtl = False # not a QTL yet, but start scanning
                start, end, peak = i,i,i # start, end, strongest signal of peak region, region we're confident in having the gene in
                while (end < len(locs) - 1) and lenient[end]: # while can extend peak to right, do so
                    end += 1
                    if (not SP.isnan(delta[end])) and stringent[end]:
                        is_qtl = True # actual QTL if the stretch of leniently significant loci has at least one properly significant change
                        if abs(delta[end]) > abs(delta[peak]):   peak = end # store location of strongest peak

                if is_qtl and (locs[end] - locs[start] > length_cutoff): # if strong enough change exists, and the region is long enough,
                    centre = SP.where((locs >= locs[start]) & (locs <= locs[end]) & (abs(delta[peak]) <= abs(delta) + peak_cutoff))[0] # all loci with change within peak_cutoff
                    genes = [""]
                    qtls.append((chrm, locs[peak], delta[peak], sd[peak], abs(delta[peak])/(sd[peak] + 1e-6), locs[start], locs[end], locs[centre.min()], locs[centre.max()], genes)) # store QTL
                elif is_qtl: LOG.debug("Skipping QTL chrm %s %d (af change=%.2f) as too short (%d)"%(chrm, locs[peak], delta[peak], locs[end] - locs[start]))
                
                lenient[start:end+1] = stringent[start:end+1] = False # remove the stretch just observed from future consideration

    # After calling the QTLs, output them all to the desired file
    write_qtls(out_file, qtls, HEADER_CALL_REGIONS%(SQTL_VERSION, out_file, posterior_file, str(sample_pairs), af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff))
