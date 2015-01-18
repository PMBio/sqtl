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

def plot_qtls(out_dir, qtl_file, posterior_file, samples, skip_plot_ml=False, skip_plot_ci=False, genefile="", extension=".png", dpi=300, screen=False):
    posteriors = read_posterior(posterior_file)
    for qtlset, chrm, peak, start, c_start, c_end, end, length, delta, sd, sds, genes in read_qtls(qtl_file):
        plot_qtl(posteriors, samples, chrm, int(peak), int(start), int(end), int(c_start), int(c_end), out_dir, skip_plot_ml=skip_plot_ml, skip_plot_ci=skip_plot_ci, extension=extension, dpi=dpi, screen=screen, genefile=genefile)

        
def plot_singlepanel_qtl(posteriors, samples, chrm, peak, start, end, c_start, c_end, out_dir, screen=False, extension="png", dpi=300, skip_plot_ml=False, skip_plot_ci=False):
    if screen: PL.rcParams.update(PLOT_PARAMS_SCREEN)
    LOG.info("plot_qtl - out_dir=%s, chrm=%s, start-csctart-peak-cend-end=%d %d %d %d %d, samples=%s, skip_plot_ml=%s, skip_plot_ci=%s, extension=%s, dpi=%d"%(out_dir, chrm, start, c_start, peak, c_end, end, str(samples), str(skip_plot_ml), str(skip_plot_ci), extension, dpi))
    colours = 'bgryck'
    if not os.path.exists(out_dir):
        LOG.debug("plot_qtl - Output directory %s does not exist; creating."%out_dir)
        os.system("mkdir -p %s"%out_dir)

    PL.figure(None, [8,3.5])
    for s, sample in enumerate(samples):
        col = colours[s%6]
        if not skip_plot_ml: plot_qtl_ml(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['ML'], col)
        plot_qtl_af(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['AF'], col)
        if not skip_plot_ci: plot_qtl_ci(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['AF'], posteriors[sample][chrm]['SD'], col)
    #PL.plot([start, end], [0,0], 'k-')
    #PL.plot([start, end], [1,1], 'k-')
    PL.plot([start, end], [0.5,0.5], 'k--', alpha=0.5)
    PL.plot([peak, peak], [0,1], 'r--', alpha=0.5)
    PL.plot([c_start, c_start], [0,1], 'k--', alpha=0.5)
    PL.plot([c_end, c_end], [0,1], 'k--', alpha=0.5)
    PL.ylim(0,1)
    PL.yticks(SP.arange(0,1.01,0.2))
    PL.ylabel("Reference allele frequency")
    PL.xlim(start, end)
    PL.xlabel(chrm)
    PL.title("")
    PL.savefig("%s/%s_%d-%d.%s"%(out_dir, chrm, start, end, extension), dpi=dpi)

    

""" Plot allele frequency for one chromosome. Assume all plotting things are taken care of. Plot dots for ML AF estimates, solid line for smoothed versions.
@param locs length-L array of chromosomal coordinates
@param init length-L array of max likelihood estimates of locus allele frequencies (a/(a+b))"""
def plot_qtl_af(locs, post, color):
    PL.plot(locs, post, alpha=1, color=color, lw=4)

def plot_qtl_ml(locs, init, color):
    PL.plot(locs, init, ".", markersize=1, alpha=1, color=color)

def plot_qtl_ci(locs, post, sd, color, alpha_multiplier=0.06):
    sd = SP.array(sd)
    x_ci = SP.array(list(locs) + list(locs)[::-1])
    y_ci = SP.array(list(post) + list(post)[::-1])
    sds = SP.array(list(sd) + list(-sd)[::-1])

    PL.fill(x_ci, [max(0,y) for y in y_ci + 2*sds], color, alpha=alpha_multiplier)
    PL.fill(x_ci, [max(0,y) for y in y_ci + 4*sds], color, alpha=2*alpha_multiplier) 


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


def plot_qtl2(set_name, sample_high, sample_low, chrm, peakloc, diff, sd, startloc, endloc, add=False, show=True, colors="rb", rec_cutoff=None, diff_ylim=None, diffline=None, title=True, filter=False, separate_plots=True, rm=True):
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


def plot_qtl(posteriors, samples, chrm, peak, start, end, c_start, c_end, out_dir, screen=False, extension="png", dpi=300, skip_plot_ml=False, skip_plot_ci=False, genefile=None):
    if screen: PL.rcParams.update(PLOT_PARAMS_SCREEN)
    LOG.info("plot_twopanel_qtl - out_dir=%s, chrm=%s, start-csctart-peak-cend-end=%d %d %d %d %d, samples=%s, skip_plot_ml=%s, skip_plot_ci=%s, extension=%s, dpi=%d"%(out_dir, chrm, start, c_start, peak, c_end, end, str(samples), str(skip_plot_ml), str(skip_plot_ci), extension, dpi))
    colours = 'bgryck'
    if not os.path.exists(out_dir):
        LOG.debug("plot_qtl - Output directory %s does not exist; creating."%out_dir)
        os.system("mkdir -p %s"%out_dir)

    PL.figure(None, [8,7])
    for w,window in enumerate([200000, 20000]): # two window sizes - 200kb, and ~16 genes (or 20kb, if genes not given)
        PL.subplot(2,1,w+1)
        for s, sample in enumerate(samples):
            col = colours[s%6]
            if not skip_plot_ml: plot_qtl_ml(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['ML'], col)
            plot_qtl_af(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['AF'], col)
            if not skip_plot_ci: plot_qtl_ci(posteriors[sample][chrm]['L'], posteriors[sample][chrm]['AF'], posteriors[sample][chrm]['SD'], col)
        PL.plot([-1e5,1e7], [0.5,0.5], 'k--', alpha=0.5)
        PL.plot([peak, peak], [0,1], 'r--', alpha=0.5)
        PL.plot([c_start, c_start], [0,1], 'k--', alpha=0.5)
        PL.plot([c_end, c_end], [0,1], 'k--', alpha=0.5)
        PL.ylim(0,1)
        PL.yticks(SP.arange(0,1.01,0.2))
        PL.ylabel("Reference allele frequency")
        if (w == 1) and (genefile is not None):
            PL.xlim(*plot_gene_boxes(genefile, chrm, peak, n_genes=16))
            PL.ylim(-0.2,1)
        else: PL.xlim(peak-window/2, peak+window/2)
        PL.title(chrm)
    PL.savefig("%s/%s_%d-%d.%s"%(out_dir, chrm, start, end, extension), dpi=dpi)


def plot_gene_boxes(genefile, chrm, peak, n_genes=16, h=0.06, y=-0.1):
    genes = read_genes(genefile)
    if chrm not in genes:
        LOG.error("QTL location chromosome %s not in the list of genes in file %s. Resorting to 20kb window."%(chrm, genefile))
        return peak-10000, peak+10000
    loc, name = genes[chrm]
    
    mids = loc.mean(axis=1)
    I = SP.argsort(abs(mids - peak))
    for j,i in enumerate(sorted(I[0:n_genes])):
        x1, x2 = loc[i]
        PL.fill([x1,x2,x2,x1], [y+h, y+h, y-h, y-h], lw=0.5, fill=True, color="k", alpha=0.3)
        PL.text(mids[i], y - 0.02 + 0.03*((-1)**j), name[i], horizontalalignment='center', size=7, style="italic")
    win_start, win_end = min(loc[I[0:n_genes]].min(), peak-5000), max(loc[I[0:n_genes]].max(), peak+5000) # at least 10kb window
    PL.plot([win_start, win_end], [0,0], 'k-', alpha=0.1)
    return win_start, win_end


def read_genes(genefile):
    res = {}
    for l in file(genefile, 'r'):
        if l[0] == "#": continue
        chrm, start, stop, name = l.strip().split("\t")
        if chrm not in res: res[chrm] = [],[]
        res[chrm][0].append([int(start), int(stop)])
        res[chrm][1].append(name)
    for chrm in res: res[chrm] = SP.array(res[chrm][0]), res[chrm][1]
    return res


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

