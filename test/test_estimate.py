import sys
import os
import scipy as SP
from generate import *
from run_infer import read_config
from sqtl.model.infer import *
from sqtl.plot.freq import *


def get_simulation(seed, config_file, coverage, qtl_af):
    LOG.debug("get_simulation - seed=%d, coverage=%d, qtl_af=%.2f"%(seed, coverage, qtl_af))
    L, qtl_loci, r_mean, r_mean_infl, r_var, r_var_infl, small_qtl_afs = read_config(config_file)
    SP.random.seed(seed)
    sim = generate_data(L=L, r_mean=r_mean, r_var=r_var, mean_coverage=coverage, sd_coverage=0.1, fqs=small_qtl_afs+[qtl_af], qs=qtl_loci)
    LOG.debug("get_simulation - Done.")
    return sim

def get_simple_estimate(net, f_mean, f_var, l1, l2, qtl_l2):
    if (not qtl_l2) and SP.isnan(net.alpha[qtl_l2, 0, l1, l2]):
        calc_qtlside_estimates(net, l1, l2)
    est_mean = net.alpha[qtl_l2, 0, l1, l2]*f_mean +  net.alpha[qtl_l2, 1, l1, l2]# l1 - locus we are looking at now, estimating frequency from l2
    est_var = f_var*net.alphavar[qtl_l2,0,l1,l2] + f_mean*f_mean*net.alphavar[qtl_l2,1,l1,l2] + f_mean*net.alphavar[qtl_l2,2,l1,l2] + net.alphavar[qtl_l2, 3,l1,l2] #(net.alpha[qtl_l1, 0, l2, l1]**2)


def test_estimates():
    sim = get_simulation(11, "config.txt", 100000, 0.95)
    net = cSQtlModel(D=sim.D, rho_basemean=sim.r_mean, rho_basevar=sim.r_var, F0=sim.F0, calc_all_afs=False)
    l, q = 10, 188
    
    # Four combinations
    # 1. predict from QTL
    plot_prediction_fromlocus(net, sim.F1, q, q)
    # 2. predict from non-QTL locus
    plot_prediction_fromlocus(net, sim.F1, l, q)
    # 3. predict for QTL
    plot_prediction_forlocus(net, sim.F1, q, q)
    # 4. predict for non-QTL locus
    plot_prediction_forlocus(net, sim.F1, l, q)

if __name__ == '__main__':
    test_estimates()
