import scipy as SP
import pylab as PL
import glob
import os
import sys
from sqtl.tools.io import *
from sqtl.tools.common import *

PLOT_PARAMS_SCREEN = {'text.fontsize':40, 'xtick.labelsize':24, 'ytick.labelsize':24, 'text.size':40, 'axes.titlesize':30, 'axes.labelsize':26, 'figure.figsize':(16,10), 'legend.fontsize':24}
PLOT_PARAMS_PAPER = {'text.fontsize':10,'xtick.labelsize':8, 'ytick.labelsize':8, 'text.size':10, 'axes.titlesize':10, 'axes.labelsize':10, 'figure.figsize':(7,6), 'legend.fontsize':10,'backend':'ps'}
PL.rcParams.update(PLOT_PARAMS_PAPER)

def plot_chromosomes(out_dir, data_file, samples, skip_plot_ml=False, extension="png", dpi=300, screen=False):
    if screen: PL.rcParams.update(PLOT_PARAMS_SCREEN)
    LOG.info("plot_chromosomes - out_dir=%s, data_file=%s, samples=%s, skip_plot_ml=%s, extension=%s, dpi=%d"%(out_dir, data_file, str(samples), str(skip_plot_ml), extension, dpi))
    colours = 'bgryckbgryck'
    if not os.path.exists(out_dir):
        LOG.debug("plot_genome - Output directory %s does not exist; creating."%out_dir)
        os.system("mkdir -p %s"%out_dir)

    data = read_posterior(data_file)
    if samples is None or len(samples) == 0: samples = data.keys()
    if len(samples) == 0: return

    for chrm in data.values()[0]:
        max_site = max(data[samples[0]][chrm]['L']) # length of chromosome
        PL.figure(None, [10.*max_site/1e6, 5]) # assume yeast, width is 10inches per megabase
        for s, sample in enumerate(samples):
            plot_chrm_af(data[sample][chrm]['L'], data[sample][chrm]['ML'], data[sample][chrm]['AF'], colours[s%6], chrm, " ".join(samples), plot_ml=not skip_plot_ml)
        PL.savefig("%s/%s.%s"%(out_dir, chrm, extension), dpi=dpi)


""" Plot allele frequencies for the entire genome in a single plot
@param out_file file to save output image to
@param data_file file containinig posterior allele frequencies
@param samples list of sample IDs to plot
@require s in data_file for each s in samples
@param dpi resolution parameter
@param screen if True, use screen plot settings, else publication (larger font etc)"""
def plot_genome(out_file, data_file, samples, dpi=300, screen=False):
    if screen: PL.rcParams.update(PLOT_PARAMS_SCREEN)
    LOG.info("plot_genome - out_file=%s, data_file=%s, samples=%s, dpi=%d"%(out_file, data_file, str(samples), dpi))
    colors = 'bgryckbgryck'

    data = read_posterior(data_file)
    if samples is None or len(samples) == 0: samples = data.keys()
    if len(samples) == 0: return

    PL.figure(None, [14, 4])
    right_end = 0 # rightmost plotted base pair
    for chrm in sorted(data.values()[0]): # for chromosomes in ascending order
        max_site = max(data[samples[0]][chrm]['L']) # length of chromosome
        for s, sample in enumerate(samples): # plot all samples
            I = SP.where(SP.array(data[sample][chrm]['SD']) < 0.3)[0] # at sites that have confident posteriors
            PL.plot(SP.array(data[sample][chrm]['L'])[I] + right_end, SP.array(data[sample][chrm]['AF'])[I], alpha=0.4, color=colors[s], lw=2) # offset by the end of last chromosome
        if right_end > 0: PL.plot([right_end, right_end], [0,1], 'k--', lw=0.4, alpha=0.2) # plot separators between chromosomes
        right_end = right_end + max(data[sample][chrm]['L']) # update rightmost end
    PL.plot([0,right_end], [0.5,0.5], 'k--', alpha=0.3)
    PL.xlim(0,right_end)
    xrange = SP.arange(0,right_end, 1000000)
    PL.xticks(xrange, ["%d"%(int(x/1000000)) for x in xrange])
    PL.xlabel("Genome (Mb)"), PL.ylabel("Reference allele frequency")
    PL.savefig(out_file, dpi=dpi)
    

""" Plot allele frequency for one chromosome. Assume all plotting things are taken care of. Plot dots for ML AF estimates, solid line for smoothed versions.
@param locs length-L array of chromosomal coordinates
@param init length-L array of max likelihood estimates of locus allele frequencies (a/(a+b))"""
def plot_chrm_af(locs, init, post, color, chrm, title, plot_ml=True):
    if plot_ml: PL.plot(locs, init, ".", markersize=1, alpha=1, color=color)
    PL.plot(locs, post, alpha=1, color=color, lw=4)
    PL.plot([min(locs), max(locs)], [0,0], 'k-')
    PL.plot([min(locs), max(locs)], [1,1], 'k-')
    PL.ylim(-0.1,1.1)
    PL.yticks(SP.arange(0,1.01,0.2))
    PL.ylabel("Reference allele frequency")
    PL.xlim(min(locs), max(locs))
    PL.xlabel(chrm)
    PL.title(title)
    

def plot_posterior_ci(locs, mean, sd, color, alpha_multiplier=0.1, rm=True):
    x_ci = SP.array(list(locs) + list(locs)[::-1])
    y_ci = SP.array(list(mean) + list(mean)[::-1])
    if rm: y_ci = 1. - y_ci
    sds = SP.array(list(sd) + list(-sd)[::-1])
    PL.fill(x_ci, y_ci + sds, color, alpha=alpha_multiplier)
    PL.fill(x_ci, y_ci + 2*sds, color, alpha=2*alpha_multiplier) 


def plot_posterior_diff(locs, means, sds, color, plot_ci=True, alpha_multiplier=0.1, diffline=0.23, rm=True):
    overlap = set(locs[0]) & set(locs[1])
    I = [[i for i in range(len(locs[j])) if locs[j][i] in overlap] for j in range(2)]
    d = means[0][I[0]] - means[1][I[1]]
    l = locs[0][I[0]]
    sd = 2*(sds[0][I[0]]**2 + sds[1][I[1]]**2)**0.5
    if rm: d = -d
    PL.plot(l, d, color + '-', lw=4)
    PL.plot([min(l), max(l)], [0,0], 'k-')
    if diffline is not None:
        PL.plot([min(l), max(l)], [diffline,diffline], 'k--')
        PL.plot([min(l), max(l)], [-diffline,-diffline], 'k--')
    if plot_ci:
        plot_posterior_ci(l, d, sd, color, alpha_multiplier)


def plot_qtl_old(set_name, sample_high, sample_low, chrm, peakloc, diff, sd, startloc, endloc, add=False, show=True, colors="rb", rec_cutoff=None, diff_ylim=None, diffline=None, title=True, filter=False, separate_plots=True, rm=True):
    data = [get_sample_seq_data(sample, rec_cutoff)[chrm] for sample in sample_high, sample_low]
    I = [SP.where((data[i][4] >= startloc) & (data[i][4] <= endloc))[0] for i in range(2)]
    locs = [data[i][4][I[i]] for i in range(2)]
    init, mean, post, bad = [data[i][0][I[i]] for i in range(2)], [data[i][1][I[i]] for i in range(2)], [data[i][2][I[i]] for i in range(2)], [data[i][3][I[i]] for i in range(2)]
    vars = [p.prod(axis=1)/((p.sum(axis=1)**2+0.1)*(p.sum(axis=1) + 1.1)) for p in post]
    chrnames = get_seq_chr_names()
    
    if not add: PL.figure(figsize=(11, 4 + 2.5*separate_plots))
    if separate_plots:
        PL.subplot(211)
        PL.plot([peakloc, peakloc], [0,1], "k--") # highlight the peak
    for i in range(2):
        if len(I[i]) < 2: continue
        J = range(len(locs[i]))
        if filter: J = SP.where(~bad[i])[0]
        plot_chrm_af(locs[i][J], init[i][J], mean[i][J], colors[i], rm=rm)
        plot_posterior_ci(locs[i][J], mean[i][J], 2.*(vars[i][J]**0.5), colors[i], rm=rm)
    PL.xlim(startloc, endloc)
    if title and (not add):  PL.title("%s QTL at chrm %s %d-%d (peak at %d, change=%.2f; %d sd)"%(set_name, chrnames[chrm], startloc, endloc, peakloc, diff, sd))
    PL.ylabel("RM allele frequency")
    if separate_plots:
        PL.subplot(212)
        plot_posterior_diff(locs, mean, [v**0.5 for v in vars], colors[0], diffline=diffline, rm=rm)
        PL.plot([peakloc, peakloc], [-2,2], "k--") # highlight the peak
        if diff_ylim is None: diff_ylim = (-0.5,0.5)
        PL.ylabel("Difference in BY allele frequency")
    else:
        overlap = set(locs[0]) & set(locs[1])
        I = [[i for i in range(len(locs[j])) if locs[j][i] in overlap] for j in range(2)]
        d = mean[0][I[0]] - mean[1][I[1]]
        if rm: d = -d
        PL.plot(locs[0][I[0]], d, "k-", lw=2, alpha=0.7)
        if peakloc is None:
            peakloc = locs[0][I[0]][SP.argmax(abs(d))]
        PL.plot([peakloc, peakloc], [0,1], "k--") # highlight the peak
    PL.xlim(startloc, endloc)
    PL.ylim(*diff_ylim)
    PL.xlabel("Chr %s"%(chrnames[chrm]))
    if show:
        PL.show()


def plot_genes(chrm, start, end, y, h):
    all_genelocs = read_gene_locs()
    genelocs = []
    for (c, s, e) in all_genelocs.values():
        if (chrm == c) and (s > start) and (e < end):
            genelocs.append([s,e])
    for (x1,x2) in genelocs:
        PL.fill([x1,x2,x2,x1], [y+h, y+h, y-h, y-h], lw=0.5, fill=True, color="k", alpha=0.3)



def plot_paper_qtl1():
    peakloc = 215000
    window = 90000
    colors = 'rb' #['#33ff33','#006600']
    #plot_qtl("X1_sR1", "Sample4_32", "Sample4_44", "ref|NC_001136|", peakloc, 0.2, 2, peakloc-window, peakloc+window, add=False, show=False, colors=colors, rec_cutoff=0.9, diffline=None, title=False, filter=True, separate_plots=False, diff_ylim=(-0.1,1.1))
    plot_qtl("X1_sR2", "Sample4_93", "Sample4_81", "ref|NC_001136|", peakloc, 0.2, 2, peakloc-window, peakloc+window, add=True, show=False, colors=colors, rec_cutoff=0.9, diff_ylim=(-0.1,1.1), diffline=None, title=False, filter=True, separate_plots=False, rm=True)
    y,h = -0.05, 0.03
    plot_genes(4, peakloc-window, peakloc+window, y=y, h=h)
    x1,x2 = 213000, 216000
    PL.fill([x1,x2,x2,x1], [y+h, y+h, y-h, y-h], lw=0.5, fill=True, color="k")
    PL.plot([x2, x2+3000],[y+h, 0.05], 'k-')
    PL.text(x2+4000, y+h+0.07, "RGT2")#, fontsize=16)
    PL.savefig("RGT2.svg")
    PL.show()


def plot_paper_qtl2():
    peakloc = 474000
    window = 100000
    colors = 'rb' #['#33ff33','#006600']
    #plot_qtl("X3_sR1", "Sample4_58", "Sample4_70", "ref|NC_001139|", peakloc, 0.2, 2, peakloc-window, peakloc+window, add=False, show=False, colors="rb", rec_cutoff=0.9, diffline=None, title=False, filter=True)
    plot_qtl("X3_sR2", "Sample4_11", "Sample4_23", "ref|NC_001139|", peakloc, 0.2, 2, peakloc-window, peakloc+window, add=False, show=False, colors=colors, rec_cutoff=0.9, diff_ylim=(-0.1,1.1), diffline=None, title=False, filter=True, separate_plots=False)
    y,h = -0.05, 0.03
    plot_genes(4, peakloc-window, peakloc+window, y=y, h=h)
    x1,x2 = 469000, 473000
    PL.fill([x1,x2,x2,x1], [y+h, y+h, y-h, y-h], lw=0.5, fill=True, color="k")
    PL.plot([x2, x2+3000],[y+h, 0.05], 'k-')
    PL.text(x2+4000, y+h+0.07, "PDR1", fontsize=16)
    PL.savefig("PDR1.svg")
    PL.show()
    
