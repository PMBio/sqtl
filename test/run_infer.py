import sys
import os
import scipy as SP
from leo.common import *
from generate import *
from sqtl.model.infer import *
from sqtl.model.infer_local import *
from sqtl.model.model_dataonly import *
from sqtl.plot.freq import *
DATA_DIR = '/home/morphology/shared/lparts/data/projects/sqtl'
if os.popen("hostname").next().strip() == 'can1.local': DATA_DIR = "/Users/leopold/data/projects/sqtl"


def read_config(config_file):
    ifh = file(config_file, 'r')
    L = int(ifh.next().strip())
    qtl_loci = map(int, ifh.next().strip().split(","))
    r_mean = float(ifh.next().strip())
    r_mean_infl = float(ifh.next().strip())
    r_var = float(ifh.next().strip())
    r_var_infl = float(ifh.next().strip())
    small_qtl_afs = map(float, ifh.next().strip().split(","))
    return L, qtl_loci, r_mean, r_mean_infl, r_var, r_var_infl, small_qtl_afs


def get_simulation(seed, config_file, coverage, qtl_af):
    LOG.debug("get_simulation - seed=%d, coverage=%d, qtl_af=%.2f"%(seed, coverage, qtl_af))
    L, qtl_loci, r_mean, r_mean_infl, r_var, r_var_infl, small_qtl_afs = read_config(config_file)
    SP.random.seed(seed)
    sim = generate_data(L=L, r_mean=r_mean, r_var=r_var, mean_coverage=coverage, sd_coverage=0.1, fqs=small_qtl_afs+[qtl_af], qs=qtl_loci)
    sim.r_mean_infl = r_mean_infl
    sim.r_var_infl = r_var_infl
    LOG.debug("get_simulation - Done.")
    return sim


def get_models(sim, rho_mode="fixed", debug=False):
    models, qtl_loci = {}, None
    if debug: qtl_loci = [1,176,187]
    if rho_mode == "fixed":
        #pdb.set_trace()
        models['MP'] = cSQtlModelMP(D=sim.D, rho_basemean=sim.r_mean*sim.r_mean_infl, rho_basevar=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False)
        models['sQTL'] = cSQtlModel(D=sim.D, rho_basemean=sim.r_mean*sim.r_mean_infl, rho_basevar=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False)
    elif rho_mode == "truth":
        models['MP'] = cSQtlModelMP(D=sim.D, rho_basemean=sim.r_mean*sim.r_mean_infl, rho_basevar=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False, rho=sim.R)
        models['sQTL'] = cSQtlModel(D=sim.D, rho_basemean=sim.r_mean*sim.r_mean_infl, rho_basevar=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False, rho=sim.R)
    models['Smooth'] = cDataOnlyModel(D=sim.D, rho_basemean=sim.r_mean*sim.r_mean_infl, rho_basevar=sim.r_var*sim.r_var_infl, F0=sim.F0)
    models['ML'] = cMLModel(sim.D, sim.F0)
    for m in models:
        if m != "sQTL":
            models[m].infer()
    return models


def calculate_stats(sim, models):
    result_f = {}
    result_q = {}
    result_l = {}

    for m in models:
        qhat = models[m].X.E1.argmax()
        if m in ("sQTL", "MP"):
            result_f[m] = (abs(models[m].F.E1[qhat] - sim.F1)).mean()
        else:
            result_f[m] = (abs(models[m].F.E1 - sim.F1)).mean()
        s,e = min(qhat, sim.Qloc), max(qhat, sim.Qloc)
        result_q[m] = abs(sim.R[s:e].sum()*100.) # distance in cM
        result_l[m] = qhat - sim.Qloc
    return result_f, result_q, result_l


def run_simulation(seed, config_file, debug=False, rho_mode="fixed"):
    coverages = range(20,200,20) + range(200,1001,100) #range(200,500,50) + range(500,1001,100)
    qtl_afs = [0.7,0.85,0.99]
    if debug: coverages, qtl_afs = [80], [0.99] # [90000],[0.99]
    result = {}
    
    for coverage in coverages:
        for qtl_af in qtl_afs:
            sim = get_simulation(seed, config_file, coverage, qtl_af)
            models = get_models(sim, rho_mode=rho_mode, debug=debug)
            model_result = calculate_stats(sim, models)
            if debug:
                print seed, model_result
            #if False:
                #if model_result[1]['sQTL'] < 10: return
                #return
                import pylab as PL
                L = len(sim.F0)
                PL.figure()
                PL.plot(range(L),sim.F0,'k')
                PL.plot(range(L),sim.F1,'k')
                lX,lX2 = models['sQTL'].X.lnX, models['MP'].X.lnX
                L4 = [180,188]
                PL.plot(range(L),(lX -lX.min())/(abs(lX.min())),'g')
                PL.plot(range(L),(lX2 -lX2.min())/(abs(lX2.min())),'g--')
                for l in L4:
                    PL.plot(range(L),models['sQTL'].F.E1[l],'b')
                    PL.plot(range(L), models['MP'].F.E1[l], 'b--')
                PL.plot(range(L), sim.D[:,0]/(sim.D.sum(axis=1)), 'r.',markersize=5)
                PL.plot(range(L), models['Smooth'].F.E1, 'k--')
                PL.show()
                net = models['MP']
                qq = 180
                qq = 188
                PL.plot(net.message[0,qq,:,0], 'b-o')
                PL.plot(net.message[1,qq,:,0], 'g-o')
                PL.plot(net.F.E1[qq,:], 'r-o')
                PL.plot(sim.F1,'k')
                PL.plot(net.mu_d, 'y.')
                PL.show()                
                #net = models['MP']
                #PL.plot(net.mu_d, net.alpha[True,0,:,qq]*sim.F1[qq] + net.alpha[True,1,:,qq], ".")
                #PL.plot([0.5,1],[0.5,1],"-"),PL.show()
                #plot_prediction_fromlocus(net, sim.F1, 10, 188)
                #plot_prediction_forlocus(net, sim.F1, 10, 188)
                #a = SP.array([calc_posterior_estimate(net, qq, qq) for qq in range(250)])
                #PL.plot(a[:,0]/a.sum(axis=1))
                #PL.plot(range(250), net.mu_d, '.')
                #PL.plot(range(250), sim.F1)
                #PL.show()
                pdb.set_trace()
                pass
            result["cov-%d_af-%.2f"%(coverage,qtl_af)] = model_result
    if not os.path.exists("%s/sim/seed-%d"%(DATA_DIR, seed)): os.system("mkdir -p %s/sim/seed-%d"%(DATA_DIR, seed))
    cdm(result, "%s/sim/seed-%d/r-%s_config-%s.pickle"%(DATA_DIR, seed, rho_mode, config_file.split("/")[-1][0:-4]))


def main():
    seed, config_file = int(sys.argv[1]), sys.argv[2]    
    run_simulation(seed, config_file, len(sys.argv)>3)


if __name__ == '__main__':
    main()
