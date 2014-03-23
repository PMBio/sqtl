import os
import scipy as SP
import scipy.stats as ST
import pylab as PL
import pdb
from sqtl.model.basic import *
from sqtl.model.infer import *
from sqtl.tools.common import *
N_eff = 1E7
RND_SEED = 10 # was 11

def sample_f0(R, f_init=0.5, f_init_loc=0, unif_f0=False):
    L = R.shape[0]
    F0 = SP.zeros(L)
    n = SP.zeros([L, 2])
    n[f_init_loc,0] = f_init*N_eff
    n[f_init_loc,1] = (1-f_init)*N_eff
    F0[f_init_loc] = f_init
    for l in range(f_init_loc + 1, n.shape[0]):
        F0[l] = 0.5 + 0.5*(F0[l-1] - 0.5)*(1 + SP.exp(-2*R[l-1]))
    for l in range(0,f_init_loc)[::-1]:
        F0[l] = 0.5 + 0.5*(F0[l+1] - 0.5)*(1 + SP.exp(-2*R[l]))
    if unif_f0: F0[:] = f_init
    return F0

    for i in range(f_init_loc + 1, n.shape[0]):
        n[i] = n[i-1]
        rec = SP.random.poisson(R[i-1]*N_eff)
        n0 = sum(SP.rand(rec) <= F0[i-1])
        n[i][1] += n0
        n[i][1] -= (rec - n0)
        n[i][0] -= n0
        n[i][0] += (rec - n0)
        F0[i] = n[i][0]/N_eff

    for i in range(0, f_init_loc)[::-1]:
        n[i] = n[i+1]
        rec = SP.random.poisson(R[i+1]*N_eff)
        n0 = sum(SP.rand(rec) <= F0[i+1])
        n[i][1] += n0
        n[i][1] -= (rec - n0)
        n[i][0] -= n0
        n[i][0] += (rec - n0)
        F0[i] = n[i][0]/N_eff

    if unif_f0: F0[:] = f_init
    return F0
        

# return location and allele frequency at qtl
def sample_qtl(F0):
    q = len(F0)/3 # cheat
    f = (F0[q] + 3)/4.
    return q,f

# return 
def sample_r(L, r_mean, r_var):
    R = SP.zeros(L)
    a = (r_mean**2)/r_var
    b = r_mean/r_var # mean r = 0.1, var = 0.03 (sd=0.055)
    for i in range(L): 
        #R[i] = -SP.log(SP.rand(min(a,1))).sum()/b
        #R[i] = SP.randn()*(r_var**0.5) + r_mean
        R[i] = SP.random.gamma(a,1./b)
        #R[i] = 0.5*(1-SP.exp(-2*R[i])) # P(odd number of recombinations)
        #if R[i] < 0: R[i] = 1E-5
    debug = False
    if debug:
        print R
        import pylab as PL
        PL.hist(R, bins=40), PL.show()
        pdb.set_trace()
    return R


def sample_f1(F0, R, q, fq):
    F1 = SP.zeros(F0.shape)
    F1[q] = fq
    d = F1[q] - F0[q]
    r = 0

    for i in range(q + 1, R.shape[0]):
        r += R[i]
        rho = 0.5*(1. - SP.exp(-2.*r))
        F1[i] = F0[i] + 0.5*d*((F0[i] - rho)/F0[q] + (1. - F0[i] - rho)/(1 - F0[q]))

    r = 0
    for i in range(0,q)[::-1]:
        r += R[i]
        rho = 0.5*(1. - SP.exp(-2.*r))
        F1[i] = F0[i] + 0.5*d*((F0[i] - rho)/F0[q] + (1. - F0[i] - rho)/(1 - F0[q]))

    return F1
        


def sample_D(F1, mean_coverage, sd_coverage):
    D = SP.zeros([F1.shape[0], 2])
    coverage = SP.array(SP.randn(len(D))*sd_coverage + mean_coverage, int)
    for i in range(D.shape[0]):
        if coverage[i] <= 0: coverage[i] = 1
        D[i][0] = ST.binom.rvs(coverage[i], F1[i])
        D[i][1] = coverage[i] - D[i][0]
    return D

def calc_rmat(rho):
    L = len(rho)
    mat = SP.zeros([L,L])
    cr = SP.cumsum(rho)
    for l1 in range(L):
        for l2 in range(l1, L):
            mat[l1,l2] = 0.5*(1-SP.exp(-2*(cr[l2] - cr[l1])))
    return mat


def generate_data(L=51, mean_coverage=50, sd_coverage=3, r_mean=0.1, r_var=0.001, fqs=None, qs=None):
    if fqs is None: fqs = [0.85]
    if qs is None: qs = L/2

    sim = cSim()
    sim.mean_coverage = mean_coverage
    sim.sd_coverage = sd_coverage
    sim.R = sample_r(L, r_mean, r_var)
    sim.Rmat = calc_rmat(sim.R)
    F1 = sample_f0(sim.R, 0.5, 0, unif_f0=True)
    sim.Qloc, sim.fq = qs[0], fqs[0]

    for i, fq in enumerate(fqs):
        sim.F0 = F1
        F1 = sample_f1(sim.F0, sim.R, qs[i], fqs[i])

    sim.F1 = F1
    sim.Qloc, sim.fq = qs[-1], fqs[-1]

    sim.D = sample_D(sim.F1, sim.mean_coverage, sim.sd_coverage)
    sim.Q = SP.ones(len(sim.F1))
    sim.Rparam = SP.zeros(sim.D.shape)
    sim.Rparam[:,0] = r_mean
    sim.Rparam[:,1] = r_var
    sim.r_mean = r_mean
    sim.r_var = r_var

    return sim


def get_paper_sim_data(L=250, qtl_candidates=None, total_rec=2.5, var_r=2E-4, pX=None, norho=True):
    sim = None
    mean_r = total_rec/L
    qtl_loci = map(int, [L*0.4, L*0.6, L*0.75])
    qtl_afs = [0.6,0.5,0.85]
    if qtl_candidates is None: qtl_candidates = qtl_loci
    coverage = 100
    filename = "sim_%d_%s_%.2f_%.3e_%s_%s.pickle"%(L, str(qtl_candidates), total_rec, var_r, str(pX), str(norho))
    if os.path.exists(filename): sim = cl(filename)
    else:
        sim = generate_data(L=L, r_mean=mean_r, r_var=var_r, mean_coverage=coverage,fqs=qtl_afs, qs=qtl_loci)
        cdm(sim, filename)
    if pX is None: pX = sim.Q
    model = cSQtlModel(D=sim.D, rho_basemean=mean_r, rho_basevar=var_r, F0=sim.F0, pX=pX)#, qtl_loci=qtl_candidates)
    model_norho = None
    if norho:
        model_norho = cSQtlModel(sim.D, mean_r, 0, sim.F0, sim.Q, qtl_loci=qtl_candidates)
    return sim, model_norho, model



def get_paper_comparison_data(L=160, qtl_candidates=None, qtl_afs=None, total_rec=0.5, r_inflation=1, rvar_inflation=1, coverage=100, dataonly=True, norho=True):
    mean_r = total_rec/L
    var_r = 2E-6*total_rec

    sim = generate_data(L=L, r_mean=mean_r, r_var=var_r, mean_coverage=coverage,fqs=qtl_afs, qs=qtl_candidates)

    models = {}
    models['sQTL'] = cSQtlModel(sim.D, mean_r*r_inflation, var_r*rvar_inflation, sim.F0, sim.Q, qtl_loci=qtl_candidates)
    models['ML'] = cMLModel(sim.D)
    if dataonly:
        models['dataonly'] = cDataOnlyModel(sim.D, var_r*rvar_inflation)
    if norho:
        models['norho'] = cSQtlModel(sim.D, mean_r, 0, sim.F0, sim.Q, qtl_loci=qtl_candidates)
    return sim, models



def get_paper_X_data(L=160, qtl_af=0.85, total_rec=0.5, coverage=100):
    mean_r = total_rec/L
    var_r = 1E-6
    qtl_candidates = range(L/2-5, L/2+5)
    sim = generate_data(L=L, r_mean=mean_r, r_var=var_r, mean_coverage=coverage,fqs=[qtl_af], qs=[L/2])

    models = {}
    models['sQTL'] = cSQtlModel(sim.D, mean_r, var_r, sim.F0, sim.Q, qtl_loci=qtl_candidates)
    models['dataonly'] = cDataOnlyModel(sim.D, var_r)
    return sim, models


class cSim:
    def __init__(self): pass

