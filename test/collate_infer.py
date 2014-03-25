from leo.common import *
import os
import glob
import pylab as PL
import scipy as SP
import pdb

DATA_DIR = "/home/morphology/shared/lparts/data/projects/sqtl"
if os.popen("hostname").next().strip() == "can1.local":
    DATA_DIR = "/Users/leopold/data/projects/sqtl"

def collate(r_set="truth"):
    infl_sets = list(set([x.split("_")[-1] for x in glob.glob("%s/sim-3/seed-*/r-%s_*.pickle"%(DATA_DIR,r_set))])) # seed-66_r-fixed_config-config_rincrease-4.00.pickle .. => rincrease-4.00.pickle, ...
    coverages = sorted(list(set([int(x.split("_")[0][4:]) for x in cl(glob.glob("%s/sim-3/seed-*/r-%s_*.pickle"%(DATA_DIR,r_set))[0])]))) # 'cov-900_af-0.70'... => 100,...,900
    afs = map(lambda x:"%.2f"%x, [0.7,0.85,0.99])
    models = ['sQTL','Smooth','ML','MP']

    for infl_set in infl_sets:
        files = glob.glob("%s/sim-3/seed-*/r-%s_config-config_%s"%(DATA_DIR,r_set,infl_set))
        print "%s/sim-3/seed-*_r-%s_config-%s"%(DATA_DIR,r_set,infl_set)
        print "%s - %d files"%(infl_set, len(files))
        result = {}
        for af in afs:
            result[af] = {}
            for m in models:
                result[af][m] = {'F':SP.zeros([len(coverages),len(files)]), 'X':SP.zeros([len(coverages),len(files)]), 'L':SP.zeros([len(coverages),len(files)])}

        for i,f in enumerate(files):
            d = cl(f)
            for conf in d: #"cov-100_af-0.99"
                coverage = conf.split("_")[0].split("-")[1]
                af = conf.split("_")[1].split("-")[1]
                for m in d[conf][0]: # For each model
                    result[af][m]['F'][coverages.index(int(coverage)), i] = d[conf][0][m] # store difference in F
                    result[af][m]['X'][coverages.index(int(coverage)), i] = d[conf][1][m] # and difference in X
                    result[af][m]['L'][coverages.index(int(coverage)), i] = d[conf][2][m] # and difference in L

        cdm(result, "%s/output-2013/sim3-results_r-%s_%s"%(DATA_DIR,r_set,infl_set))

def plot_all_collated(r_set):
    files = glob.glob("%s/output-2013/sim3-results_r-%s_*.pickle"%(DATA_DIR,r_set))
    sets = list(set([f.split("_")[-1] for f in files]))
    for s in sets: plot_collated(r_set, s)


def plot_collated(r_set="truth", infl_set="varinfl-0.25", subplots=True, save=False):
    d = cl("%s/output-2013/sim3-results_r-%s_%s"%(DATA_DIR,r_set, infl_set))
    coverages = SP.array(range(20,200,20) + range(200,1001,100)) #range(200,500,50) + range(500,1001,100))
    if r_set == "truth": coverages = SP.array(range(20,200,20) + range(200,500,50) + range(500,1001,100))
    afs = map(lambda x:"%.2f"%x, [0.7,0.85,0.99])
    models = ['sQTL','Smooth','ML','MP']
    p = 0
    colors = 'bgry'
    if subplots: PL.figure(figsize=(14,10))
    for feature in 'FX':
        for af in afs:
            if subplots: PL.subplot(2,3,p+1)
            else: PL.figure()
            p += 1
            lines = []
            
            for i,model in enumerate(models):
                I = SP.where(d[af][model][feature].var(axis=0) > 1e-10)[0]
                err = d[af][model][feature][:,I].var(axis=1)**0.5
                lines.append(PL.plot(coverages + 2*i,SP.median(d[af][model][feature][:,I],axis=1), "-o", linewidth=3, markersize=9, color=colors[i])[0])
                PL.errorbar(coverages + 2*i, SP.median(d[af][model][feature][:,I],axis=1), yerr=err, fmt="-o", linewidth=1, markersize=9,color=colors[i])
            PL.xticks(coverages)
            #PL.xlim(min(coverages),max(coverages))
            PL.title("%s %s - %s"%(infl_set, feature, af))
            PL.xlim(15,220)

            if feature == "X": PL.ylim(0,8)
            if p == 1:  PL.legend(lines, models)
            if save: PL.savefig("/Users/leopold/doc/write/manuscripts/2011_X_sQTL/figures/figure2013-3_2%s.pdf"%("ABCDEF"[p-1:p]))
    PL.show()
