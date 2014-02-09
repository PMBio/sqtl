import sys
import os
from sqtl.tools.common import *
from sqtl.tools.io import *

BWA = "bwa"
BWA = "~/bin/bwamem"
SAMTOOLS = "samtools"

LOG_COUNT = "create_counts - out_file=%s, # pileup_files=%d, base quality cutoff=%d, mapping quality cutoff=%d "
HEADER_COUNT = """#sQTL version %s - create_counts output
#======================================
# output file=%s
# pileup files=%s
# base quality cutoff=%d
# mapping quality cutoff=%d
#
"""


def find_second_readfile(file_i, all_files):
    read_file_base = all_files[file_i].replace("_R1_", "_").replace("_R2_", "_")

    for f, filename in enumerate(all_files):
        if (f != file_i) and (read_file_base == filename.replace("_R1_", "").replace("_R2_", "_")):
            return f, filename

    return (file_i, "")

""" Map one set of read files to a reference, and merge mappings to a single .bam file
@param out_file output file path
@param ref_file reference sequence (fasta) file path
@param read_files list of read file paths (can be zipped, as long as bwa can handle them)
@param bwamem_executable command used to map; default=bwa
@param bwamem_params additional parameters to bwa mem; default=""
@param samtools_executable command used to merge mappings; default=samtools
@param sort_mem maximum memory for sorting mappings; default=512Mb
@effects creates bwa index of reference, if it does not exist"""
def map_reads(out_file, ref_file, read_files, bwamem_executable=BWA, bwamem_params="", samtools_executable=SAMTOOLS, sort_mem=512000000):
    LOG.info("map_reads - outfile=%s, ref_file=%s, %d read files. bwamem_executable=%s, bwamem_params=%s, sort_mem=%d, samtools_executable=%s"%(out_file, ref_file, len(read_files), bwamem_executable, bwamem_params, sort_mem, samtools_executable)) 
    readfile_done = [False for f in read_files] # array to keep track of whether read file has been processed
    tmp_outfiles = [] # sorted .bam files created for each set of separately mapped read files
    if not os.path.exists("%s.bwt"%ref_file):
        LOG.debug("map_reads - no reference index found; creating")
        os.system("%s index %s"%(bwamem_executable, ref_file)) # if reference is not indexed, do so

    # Go through list of read files and map to the reference
    done_files = 0
    for f in range(len(read_files)):
        if readfile_done[f]: continue # if file already processed, skip
        tmp_outfiles.append(read_files[f] + "_tmp_out_bam")
        second_read_file_i, second_read_file = find_second_readfile(f, read_files) # this is "" if there is no second read
        map_cmd = "%s mem %s %s %s %s | %s view -bSu -F 0x04 - | %s sort -m %d - %s"%(bwamem_executable, bwamem_params, ref_file, read_files[f], second_read_file, samtools_executable, samtools_executable, sort_mem, tmp_outfiles[-1])
        LOG.debug("map_reads - mapping and sorting, command=%s"%(map_cmd))
        os.system(map_cmd)
        done_files += 1
        readfile_done[f] = readfile_done[second_read_file_i] = True # keep track of which files processed to avoid double work. If no second read file, second_read_file_i = f

    # Merge individual sorted .bam files into one output, and remove intermediates
    if done_files > 1:
        merge_cmd = "%s merge %s "%(samtools_executable, out_file) + " ".join([x + '.bam' for x in tmp_outfiles]) # <samtools-exec> merge <outfile> <infile1> <infile2> ...
        LOG.debug("map_reads - Merging: %s"%merge_cmd)
        os.system(merge_cmd)
    else:
        os.system("mv %s %s"%(tmp_outfiles[0], out_file))
    LOG.debug("map_reads - Removing intermediates")
    #for f in tmp_outfiles: os.system("rm %s.bam"%f)
    LOG.debug("map_reads - Done.")


""" Create pileup of a .bam file
@param out_file output file path
@param bam_file input file path
@param ref_file reference file (fasta) path
@param samtools_executable command used to merge mappings; default=samtools
@param samtools_pileup_params additional parameters to pileup; default=" -vcs -N 2 -d 100 "
@param site_file path to file with locations of segregating sites (format one site per line, tab-delimited, chrm and location). Default: None
@param delete_bam whether to delete bam files after creating pileups. Default: False.
"""
def create_pileup(out_file, bam_file, ref_file, samtools_executable="samtools", samtools_pileup_params=" -vcs -N 2 -d 100 ", site_file=None, delete_bam=False):
    LOG.info("create_pileup - out_file=%s bam_file=%s samtools_executable=%s pileup_params=%s site_file=%s"%(out_file, bam_file, samtools_executable, samtools_pileup_params, str(site_file)))
    pileup_cmd = "%s pileup %s -f %s %s > %s"%(samtools_executable, samtools_pileup_params, ref_file, bam_file, out_file)
    if site_file is not None:
        pileup_cmd = "%s pileup %s -l %s -f %s %s > %s"%(samtools_executable, samtools_pileup_params, site_file, ref_file, bam_file, out_file)
    LOG.debug("create_pileup - command: %s"%pileup_cmd)
    os.system(pileup_cmd)
    if delete_bam: os.system("rm %s"%bam_file)
    LOG.debug("create_pileup - Done.")



""" Combine multiple pileup files; take the union of site locations
@param pileups list of pileup hashes
@return locs (list of chrm, site pairs), data (NxLx2 array of sample site allele counts (ref, nonref)), refseq (Lx2 list of ref and nonref bases at the locs)"""
def combine_pileups(pileups):
    # 0.0 Create the union of all sites across pileups
    all_locs = set([])
    for pileup in pileups: all_locs = all_locs | set(map(tuple, pileup['L']))
    # 0.1 Create indexes for each site of where it is in each pileup array; -1 if not present
    I = [dict((l,-1) for l in all_locs) for pileup in pileups] # default not present
    for p, pileup in enumerate(pileups): # fill in observed sites with their index
        for l,loc in enumerate(pileup['L']): I[p][tuple(loc)] = l

    # 1. Go through all locs, and store the counts from pileups
    locs, refseq, data = sorted(all_locs), [], SP.zeros([len(pileups), len(all_locs), 2]) # final return variables
    allele_index = {'A':0, 'C':1, 'G':2, 'T':3, '*':4}

    for l, loc in enumerate(locs): # for each site
        seq = None
        for p, pileup in enumerate(pileups): # for each sample
            if I[p][loc] >= 0: # if site observed in sample
                idx = I[p][loc]
                d,q,seq = pileup['D'][idx], pileup['Q'][idx], pileup['refseq'][idx]
                data[p,l,0] = d[allele_index[seq]]
                d[allele_index[seq]] = 0
                data[p,l,1] = d.max()
        refseq.append(tuple([seq, "ACGT*"[d.argmax()]]))
        
    return locs, data, refseq


""" Create counts from a set of pileup files. First, find the union of all sites, then output the allele counts at that site for all samples in one line.
@param out_file path to output file
@param pileup_files list of paths to pileup files
@param base_qual_cutoff minimum base quality to be included in output
@param map_qual_cutoff minimum mapping quality to be included in output"""
def create_counts(out_file, pileup_files, base_qual_cutoff=20, map_qual_cutoff=20):
    LOG.info(LOG_COUNT%(out_file, len(pileup_files), base_qual_cutoff, map_qual_cutoff))
    samples, pileups = [], []
    ofh = file(out_file, 'w')

    # 0. Create header, read pileups and sample names
    for p,pileup_file in enumerate(pileup_files):
        samples.append(pileup_file.split("/")[-1].split(".")[0])
        pileups.append(read_pileup(pileup_file, base_qual_cutoff=base_qual_cutoff, map_qual_cutoff=map_qual_cutoff))
        
    ofh.write(HEADER_COUNT%(SQTL_VERSION, out_file, str(pileup_files), base_qual_cutoff, map_qual_cutoff))
    ofh.write("#Chrm\tLoc\tRef\tNonref")
    for p,pileup_file in enumerate(pileup_files):
        ofh.write("\t%s_ref\t%s_nonref"%(samples[p], samples[p]))
    ofh.write("\n")

    # 1. combine pileups, output at the union of all loci across samples
    locs, counts, ref = combine_pileups(pileups)
    for l in range(len(locs)):
        ofh.write("%s\t%d"%(locs[l]))
        ofh.write("\t%s\t%s"%(ref[l]))
        for s in range(len(samples)):
            ofh.write("\t%d\t%d"%(tuple(counts[s][l])))
        ofh.write("\n")
    ofh.close()
    LOG.debug("create_counts - Done.")
