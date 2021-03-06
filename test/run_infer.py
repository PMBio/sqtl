import sys
import os
import scipy as SP
from leo.common import *
from generate import *
from sqtl.model.infer import *
from sqtl.model.smooth import *
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


def get_models(sim):
    smooth_F1 = smooth(sim.D, range(len(sim.D)), n_rounds=2, rec_rate=sim.r_mean*sim.r_mean_infl, nearby_snp_cutoff=0)[0][1]
    sqtl0 = cSQtlModel(D=sim.D, r_init_mean=sim.r_mean*sim.r_mean_infl, r_init_var=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False, init_from_smoothed=False)
    sqtl1 = cSQtlModel(D=sim.D, r_init_mean=sim.r_mean*sim.r_mean_infl, r_init_var=sim.r_var*sim.r_var_infl, F0=sim.F0, calc_all_afs=False, F1_init=smooth_F1)
    sqtl0.infer()
    sqtl1.infer()
    sqtl_q0 = sqtl0.X.E1.argmax()
    sqtl_q1 = sqtl1.X.E1.argmax()
    
    result_qtl = {'sQTL0':sqtl_q0, 'sQTL1':sqtl_q1, 'ML':abs(sqtl0.mu_d - sim.F0).argmax(), 'Smooth':abs(smooth_F1 - sim.F0).argmax()}
    result_f = {'sQTL0': abs(sqtl0.F.E1[sqtl_q0] - sim.F1).mean(), 'sQTL1': abs(sqtl1.F.E1[sqtl_q1] - sim.F1).mean(), 'ML':abs(sqtl0.mu_d - sim.F1).mean(), 'Smooth':abs(smooth_F1 - sim.F1).mean()}
    result_r, result_l = {}, {}
    for model in result_qtl:
        qhat = result_qtl[model]
        s,e = min(qhat, sim.Qloc), max(qhat, sim.Qloc)
        result_r[model] = abs(sim.R[s:e].sum()*100.) # distance in cM
        result_l[model] = qhat - sim.Qloc
    return result_f, result_r, result_l



def run_simulation(seed, config_file, debug=False, rho_mode="fixed"):
    coverages = range(20,200,20) + range(200,1001,100) #range(200,500,50) + range(500,1001,100)
    qtl_afs = [0.7,0.85,0.99]
    if debug: coverages, qtl_afs = [60], [0.85]#[80], [0.99] # [90000],[0.99]
    result = {}
    
    for coverage in coverages:
        for qtl_af in qtl_afs:
            sim = get_simulation(seed, config_file, coverage, qtl_af)
            model_result = get_models(sim)
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
                    PL.plot(range(L),models['sQTL'].F.E1[l],'b--')
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
    if not os.path.exists("%s/sim-3/seed-%d"%(DATA_DIR, seed)): os.system("mkdir -p %s/sim-3/seed-%d"%(DATA_DIR, seed))
    cdm(result, "%s/sim-3/seed-%d/r-%s_config-%s.pickle"%(DATA_DIR, seed, rho_mode, config_file.split("/")[-1][0:-4]))


def main():
    seed, config_file = int(sys.argv[1]), sys.argv[2]    
    run_simulation(seed, config_file, len(sys.argv)>3)


if __name__ == '__main__':
    main()
