#!/usr/bin/python
import sys
sys.path.append("../")
from optparse import OptionParser
from sqtl.tools.map import *
from sqtl.tools.call import *
from sqtl.model.smooth import *
from sqtl.plot.genome import *
from sqtl.plot.qtl import *
# possible commands to include: demo, finemap, call
# usage: \tdemo\t\tShow example commands to run on the included test dataset
# \t\tfinemap\t\tFinemap QTLs
# \t\tcall\t\tRun entire pipeline from reads to allele frequency estimates, QTLs, and plots


USAGE_STR = """Program: sqtl (selection quantitative trait locus mapping)
Version: %s
Contact: Leopold Parts <leopold.parts@gmail.com>

Usage: sqtl <command> [options]

Command:
\t\tmap\t\tMap fastq files to reference
\t\tcalc_afs\tCreate pileups at segregating sites, count reference and non-reference alleles, and infer allele frequencies
\t\t    pileup\t    Create pileups at segregating sites
\t\t    count\t    Count reference and non-reference alleles
\t\t    infer\t    Infer posterior allele frequencies at segregating sites
\t\tcall_regions\tCall regions with large change between pairs of samples

\t\tplot_genome\tPlot genome-wide allele frequencies
\t\tplot_chrms\tPlot allele frequencies separately for each chromosome
\t\tplot_qtls\tPlot allele frequencies at QTLs
"""%SQTL_VERSION

USAGE_MAP = """Map reads from fastq to a reference genome to produce a sorted merged .bam file.
Usage:
sqtl map [options] out_file reference_file read_file1 [read_file2 [read_file3 ...]]

Options:
        -b STR    (path to) BWA 0.7.5+ executable [bwa]
        -p STR    BWA mem parameters [""]
        -s STR    (path to) samtools 0.1.12 executable [samtools]
        -m INT    maximum memory for merging bam files [512000000]
"""

USAGE_PILEUP = """Create read pileups at segregating sites in the cross from alignment .bam file.
Usage:
sqtl pileup [options] out_file reference_file bam_file 

Options: 
        -s STR    (path to) samtools 0.1.12 executable [samtools]
        -p STR    samtools pileup parameters [" -vcs -N 2 -d 100 "]
        -d        delete the BAM file [False]
        -f STR    (path to) tab-delimited file of segregating site locations (chrm location) [""]
"""

USAGE_COUNT = """Count reference and non-reference alleles from given pileup files.  
Usage:
sqtl count [options] out_file pileup_file1 [pileup_file2 [pileup_file3 ...]]

Options:
        -b INT    base quality cutoff [20]
        -q INT    mapping quality cutoff [20]
"""

USAGE_INFER = """Infer posterior allele frequencies from observed allele counts.
Usage:
sqtl infer [options] out_file count_file

Options: 
        -r INT    number of rounds of inference to iterate between allele frequency calls and bad site filtering [2]
        -a FLOAT  per-base recombination rate [80/12M]
        -c FLOAT  recombination rate cutoff in window [0.9]
        -o FLOAT  outlier cutoff; difference between max likelihood allele frequency estimate and the average to be considered ok [0.2]
        -n INT    nearby SNP cutoff; minimum distance between nearby segregating sites. If set, in-between sites assumed to come from the same fragment are skipped [200]
        -m FLOAT  max number of median coverage. Limit the maximum coverage at any one site to be this much times the median genome-wide coverage [2.]
        -w FLOAT  fixed weight; weight relative to 1.0 that sites that appear fixed in the population have. These could be selected, or false SNP calls. [0.05]
"""

USAGE_AFS = """Calculate pileups, count alleles, and infer posterior allele frequencies from observed allele counts.
Usage:
sqtl calc_afs [options] out_file ref_file bam_file1 [bam_file2 [bam_file3 ...]]

Options: 
        -s STR    (path to) samtools 0.1.12 executable [samtools]
        -p STR    samtools pileup parameters [" -vcs -N 2 -d 100 "]
        -f STR    (path to) tab-delimited file of segregating site locations (chrm location) [""]
        -b INT    base quality cutoff [20]
        -q INT    mapping quality cutoff [20]
        -r INT    number of rounds of inference to iterate between allele frequency calls and bad site filtering [2]
        -a FLOAT  per-base recombination rate [80/12M]
        -c FLOAT  recombination rate cutoff in window [0.9]
        -o FLOAT  outlier cutoff; difference between max likelihood allele frequency estimate and the average to be considered ok [0.2]
        -n INT    nearby SNP cutoff; minimum distance between nearby segregating sites. If set, in-between sites assumed to come from the same fragment are skipped [200]
        -m FLOAT  max number of median coverage. Limit the maximum coverage at any one site to be this much times the median genome-wide coverage [2.]
        -w FLOAT  fixed weight; weight relative to 1.0 that sites that appear fixed in the population have. These could be selected, or false SNP calls. [0.05]
"""


USAGE_REGION = """Call QTL regions that are on average changing in the same direction in allele frequency in the contrasted pairs.
Usage:
sqtl call_regions [options] out_file posterior_file selected_sample1,unselected_sample1 [selected_sample2,unselected_sample2 [selected_sample3,unselected_sample3 [...]]]

Options: 
        -l FLOAT  lenient absolute change in allele frequency cutoff [0.1]
        -s FLOAT  strict absolute change in allele frequency cutoff [0.2]
        -p FLOAT  peak height (maximum allowed allele frequency "distance" from the peak) [0.03]

        --sd_lenient   FLOAT lenient change z-score cutoff   [3.0]
        --sd_stringent FLOAT strict  change z-score cutoff   [6.0]
        --min_length   INT   minimum bp length of QTL region [0]
"""


USAGE_PLOT_GENOME = """Plot entire genome allele frequency for given samples. The extension of out_file determines the file type
Usage:
sqtl plot_genome [options] out_file data_file sample1 [sample2 [sample3 [...]]]

Options:
        -d INT    DPI (dots per inch); resolution [300]
        -s STR    plot with screen parameters instead of publication plot parameters [False]
"""


USAGE_PLOT_CHRM = """Plot individual chromosome allele frequencies for given samples.
Usage:
sqtl plot_chromosomes [options] out_dir data_file sample1 [sample2 [sample3 [...]]]

Options:
        -e STR    extension [.png]
        -d INT    DPI (dots per inch); resolution [300]
        -m STR    plot dots for maximum likelihood estimates (raw data) [False]
        -s STR    plot with screen parameters instead of publication plot parameters [False]
"""

USAGE_PLOT_QTLS = """Plot QTL allele frequencies and combined signal for given samples.
Usage:
sqtl plot_qtls [options] out_dir qtl_file posterior_file sample1 [sample2 [sample3 [...]]]

Options:
        -e STR    extension [.png]
        -d INT    DPI (dots per inch); resolution [300]
        -m STR    plot dots for maximum likelihood estimates (raw data) [False]
        -s STR    plot with screen parameters instead of publication plot parameters [False]
"""


"""Map sequencing reads to the reference for multiple sets of reads """ 
def sqtl_map():
    parser = OptionParser(usage=USAGE_MAP)
    parser.add_option("-b", "--bwa", type="string", dest="bwamem_executable", default="bwa")
    parser.add_option("-p", "--bwa_params", type="string", dest="bwamem_params", default="")
    parser.add_option("-s", "--samtools", type="string", dest="samtools_executable", default="samtools")
    parser.add_option("-m", "--memory", type="int", dest="sort_mem", default=512000000)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_MAP
        return
    map_reads(out_file=args[1], ref_file=args[2], read_files=args[3:], **(eval(str(opts))))


""" Pile mapped reads up at segregating sites for one mapping file """
def sqtl_pileup():
    parser = OptionParser(usage=USAGE_PILEUP)
    parser.add_option("-s", "--samtools", type="string", dest="samtools_executable", default="samtools")
    parser.add_option("-p", "--pileup_params", type="string", dest="samtools_pileup_params", default=" -vcs -N 2 -d 100 ")
    parser.add_option("-d", "--delete_bamfile", dest="delete_bam", default=False)
    parser.add_option("-f", "--site_file", type="string", dest="site_file", default=None)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_PILEUP
        return
    create_pileup(out_file=args[1], ref_file=args[2], bam_file=args[3], **(eval(str(opts))))


""" Create allele counts for multiple pileup files """
def sqtl_count():
    parser = OptionParser(usage=USAGE_COUNT)
    parser.add_option("-b", "--basequal", type="int", dest="base_qual_cutoff", default=20)
    parser.add_option("-q", "--mapqual", type="int", dest="map_qual_cutoff", default=20)
    opts, args = parser.parse_args()
    if len(args) < 3:
        print USAGE_COUNT
        return
    create_counts(out_file=args[1], pileup_files=args[2:], **(eval(str(opts))))


""" Infer allele frequency across the genome without QTL calling """
def sqtl_infer():
    parser = OptionParser(usage=USAGE_INFER)
    parser.add_option("-r", "--rounds", type="int", dest="n_rounds", default=2)
    parser.add_option("-a", "--rec_rate", type="float", dest="rec_rate", default=80./12000000)
    parser.add_option("-c", "--rec_cutoff", type="float", dest="rec_cutoff", default=0.9)
    parser.add_option("-o", "--outlier_cutoff", type="float", dest="outlier_cutoff", default=0.2)
    parser.add_option("-n", "--nearby_snp_cutoff", type="int", dest="nearby_snp_cutoff", default=200)
    parser.add_option("-m", "--max_median_coverage", type="float", dest="max_num_median_coverage", default=2.)
    parser.add_option("-w", "--fixed_site_weight", type="float", dest="fixed_weight", default=0.05)
    opts, args = parser.parse_args()
    if len(args) < 3:
        print USAGE_INFER
        return    
    smooth_samples(out_file=args[1], count_file=args[2], **(eval(str(opts))))


""" Create pileups, counts and infer allele frequency across the genome without QTL calling """
def sqtl_calc_afs():
    parser = OptionParser(usage=USAGE_AFS)
    parser.add_option("-s", "--samtools", type="string", dest="samtools", default="samtools")
    parser.add_option("-p", "--pileup_params", type="string", dest="pileup_params", default=" -vcs -N 2 -d 100 ")
    parser.add_option("-f", "--site_file", type="string", dest="site_file", default=None)
    parser.add_option("-b", "--basequal", type="int", dest="basequal", default=20)
    parser.add_option("-q", "--mapqual", type="int", dest="mapqual", default=20)
    parser.add_option("-r", "--rounds", type="int", dest="rounds", default=2)
    parser.add_option("-a", "--rec_rate", type="float", dest="rec_rate", default=80./12000000)
    parser.add_option("-c", "--rec_cutoff", type="float", dest="rec_cutoff", default=0.9)
    parser.add_option("-o", "--outlier_cutoff", type="float", dest="outlier_cutoff", default=0.2)
    parser.add_option("-n", "--nearby_snp_cutoff", type="int", dest="snp_cutoff", default=200)
    parser.add_option("-m", "--max_median_coverage", type="float", dest="max_num_median_coverage", default=2.)
    parser.add_option("-w", "--fixed_site_weight", type="float", dest="fixed_weight", default=0.05)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_AFS
        return

    pileup_files = []
    for bam_file in args[3:]:
        pileup_files.append(bam_file + ".pileup")
        create_pileup(out_file=pileup_files[-1], ref_file=args[2], bam_file=bam_file, samtools_executable=opts.samtools, samtools_pileup_params=opts.pileup_params, site_file=opts.site_file, delete_bam=False)
    count_file = args[1] + ".counts"
    create_counts(out_file=count_file, pileup_files=pileup_files, base_qual_cutoff=opts.basequal, map_qual_cutoff=opts.mapqual)
    smooth_samples(out_file=args[1], count_file=count_file, n_rounds=opts.rounds, rec_rate=opts.rec_rate, rec_cutoff=opts.rec_cutoff, outlier_cutoff=opts.outlier_cutoff, nearby_snp_cutoff=opts.snp_cutoff, max_num_median_coverage=opts.max_num_median_coverage, fixed_weight=opts.fixed_weight)
    os.system("rm %s"%count_file)

    
""" Call regions that are under selection """ 
def sqtl_call_regions():
    parser = OptionParser(usage=USAGE_REGION)
    parser.add_option("-l", "--af_lenient", type="float", dest="af_lenient", default=0.1)
    parser.add_option("-s", "--af_stringent", type="float", dest="af_stringent", default=0.2)
    parser.add_option("", "--sd_lenient", type="float", dest="sd_lenient", default=3)
    parser.add_option("", "--sd_stringent", type="float", dest="sd_stringent", default=6)
    parser.add_option("", "--min_length", type="float", dest="length_cutoff", default=0)
    parser.add_option("-p", "--peak_height", type="float", dest="peak_cutoff", default=0.03)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_REGION
        return    
    call_qtl_regions(out_file=args[1], posterior_file=args[2], sample_pairs=args[3:], **(eval(str(opts))))
    

""" Create genome-wide plots of requested samples, one image """ 
def sqtl_plot_genome():
    parser = OptionParser(usage=USAGE_PLOT_GENOME)
    parser.add_option("-d", "--dpi", type="int", dest="dpi", default=300)
    parser.add_option("-s", "--screen", dest="screen", default=False)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_PLOT_GENOME
        return
    plot_genome(out_file=args[1], data_file=args[2], samples=args[3:], **(eval(str(opts))))


""" Create genome-wide plots of requested samples, one image per chromosome """ 
def sqtl_plot_chromosomes():
    parser = OptionParser(usage=USAGE_PLOT_CHRM)
    parser.add_option("-e", "--extension", type="str", dest="extension", default="png")
    parser.add_option("-d", "--dpi", type="int", dest="dpi", default=300)
    parser.add_option("-m", "--skip_ml_estimates", dest="skip_plot_ml", default=False)
    parser.add_option("-s", "--screen", dest="screen", default=False)
    opts, args = parser.parse_args()
    if len(args) < 4:
        print USAGE_PLOT_CHRM
        return
    plot_chromosomes(out_dir=args[1], data_file=args[2], samples=args[3:], **(eval(str(opts))))


""" Plot QTLs from a file given QTL info, allele frequencies and samples """
def sqtl_plot_qtls():
    parser = OptionParser(usage=USAGE_PLOT_QTLS)
    parser.add_option("-e", "--extension", type="str", dest="extension", default="png")
    parser.add_option("-d", "--dpi", type="int", dest="dpi", default=300)
    parser.add_option("-m", "--skip_ml_estimates", dest="skip_plot_ml", default=False)
    parser.add_option("-c", "--skip_confidence_intervals", dest="skip_plot_ci", default=False)
    parser.add_option("-s", "--screen", dest="screen", default=False)
    opts, args = parser.parse_args()
    if len(args) < 5:
        print USAGE_PLOT_QTLS
        return
    plot_qtls(out_dir=args[1], qtl_file=args[2], posterior_file=args[3], samples=args[4:], **(eval(str(opts))))


def sqtl_call():
    pass


def main():
    if len(sys.argv) < 2:
        print USAGE_STR
        return
    
    if sys.argv[1] == "map": sqtl_map()
    elif sys.argv[1] == "pileup": sqtl_pileup()
    elif sys.argv[1] == "count": sqtl_count()
    elif sys.argv[1] == "infer": sqtl_infer()
    elif sys.argv[1] == "calc_afs": sqtl_calc_afs()
    elif sys.argv[1] == "call_regions": sqtl_call_regions()
    elif sys.argv[1] == "plot_genome": sqtl_plot_genome()
    elif sys.argv[1] == "plot_chrms": sqtl_plot_chromosomes()
    elif sys.argv[1] == "plot_qtls": sqtl_plot_qtls()
    else:
	print USAGE_STR
	return


if __name__ == '__main__':
    main()
