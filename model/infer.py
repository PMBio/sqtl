import scipy as SP
import scipy.stats as ST
import time as time
import pdb
from sqtl.model.stats import *
from sqtl.tools.common import *

PREDICT_FROM_RIGHT, PREDICT_FROM_LEFT = [0,1]
PREDICT_TOWARDS_QTL, PREDICT_FROM_QTL = [0,1]


class cFNode:
    def __init__(self, D):
        self.L = len(D)
        self.E1, self.var = SP.zeros([self.L, self.L])*SP.nan, SP.zeros([self.L, self.L])*SP.nan
        self.E1[:] = (D[:,0] + EPS)/(D.sum(axis=1) + 2*EPS) # EPS < 1e-10, included just to avoid division by 0 in mean/variance computation
        self.var[:] = D.prod(axis=1)/((D.sum(axis=1)**2 + EPS)*(D.sum(axis=1) + 1))
        self.message = SP.zeros([2,self.L,self.L,2]) # for message passing to neighbours; E(f_l | left side data) and right side.

    def get_neighbour_estimate(self, net, f_mean, f_var, l1, l2, qtl_l2, q):
        if (not qtl_l2) and SP.isnan(net.alpha[qtl_l2, 0, l1, l2]):
            self.calc_qtlside_estimate(net, l1, l2)
        est_mean = net.alpha[qtl_l2, 0, l1, l2]*f_mean +  net.alpha[qtl_l2, 1, l1, l2]# l1 - locus we are looking at now, estimating frequency from l2
        est_var = f_var*net.alphavar[qtl_l2,0,l1,l2] + f_mean*f_mean*net.alphavar[qtl_l2,1,l1,l2] + f_mean*net.alphavar[qtl_l2,2,l1,l2] + net.alphavar[qtl_l2, 3,l1,l2] 
        #if (not qtl_l2) and (f_mean - net.F0[l2] > est_mean - net.F0[l1]): pdb.set_trace() 
        return est_mean, est_var


    def calc_qtlside_estimate(self, net, l1, l2):
        net.R.gm[l1,l2] = get_gmean(net.R.rm[l1,l2], net.R.rvar[l1,l2], net.F0[l1] + net.F0[l2] - 2*net.F0[l1]*net.F0[l2], n_bins=200)
        net.R.gvar[l1,l2] = get_gvar(net.R.rm[l1,l2], net.R.rvar[l1,l2], net.F0[l1] + net.F0[l2] - 2*net.F0[l1]*net.F0[l2], net.R.gm[l1,l2],n_bins=200)
        a, b = 2.*net.F0[l1]*(1 - net.F0[l1]), 0
        c,d  = -2.*net.F0[l1]*(1 - net.F0[l1])*net.F0[l2], net.F0[l1]
        net.alpha[0,0,l1,l2] = a*net.R.gm[l1,l2] + b # alpha[0, :, l1, l2] = coefficients for estimating f if QTL is towards l2 compared to l1 from l1
        net.alpha[0,1,l1,l2] = c*net.R.gm[l1,l2] + d
        net.alphavar[0,0,l1,l2] = (a**2)*(net.R.gvar[l1,l2] + net.R.gm[l1,l2]**2) + b**2
        net.alphavar[0,1,l1,l2] = (a**2)*net.R.gvar[l1,l2]
        net.alphavar[0,2,l1,l2] = 2*a*c*net.R.gvar[l1,l2]
        net.alphavar[0,3,l1,l2] = net.R.gvar[l1,l2]*(c**2)


    def update(self, net, debug=False):
        t0 = time.time()
        LOG.debug("Calculating allele frequencies")
        for q in range(net.L):
            for l in range(net.L):
                # first, left to right
                m1,t1,m2,t2 = net.mu_d[l], net.tau_d[l], 0, 0
                if l > 0: # if can predict from left, do so
                    m2,t2 = self.get_neighbour_estimate(net, self.message[PREDICT_FROM_LEFT, q, l-1, 0], 1./self.message[PREDICT_FROM_LEFT, q, l-1, 1], l, l-1, q < l, q)
                    t2 = 1./t2
                self.message[PREDICT_FROM_LEFT, q, l] = (m1*t1 + m2*t2)/(t1+t2), t1+t2
                
                # then, right to left from right hand side
                l2 = net.L - l - 1
                m1,t1,m2,t2 = net.mu_d[l2], net.tau_d[l2], 0, 0
                if l2 < net.L - 1: # if can predict form right, do so
                    m2,t2 = self.get_neighbour_estimate(net, self.message[PREDICT_FROM_RIGHT, q, l2+1, 0], 1./self.message[PREDICT_FROM_RIGHT, q, l2+1, 1], l2, l2+1, q > l2, q)
                    t2 = 1./t2
                self.message[PREDICT_FROM_RIGHT, q, l2]  = (m1*t1 + m2*t2)/(t1+t2), t1+t2

            for l in range(net.L):
                # once all messages left and right computed, form the actual prediction by combining estimates from left, right, and observed data. This can be shortened at expense of readability.
                m1,t1 = net.mu_d[l], net.tau_d[l]
                m2,t2 = self.message[PREDICT_FROM_LEFT, q, l]
                m3,t3 = self.message[PREDICT_FROM_RIGHT, q, l]
                self.E1[q,l] = (m1*t1+m2*t2+m3*t3)/(t1+t2+t3)
                self.var[q,l] = 1./(t1+t2+t3)


class cXNode:
    def __init__(self, pX):
        self.pX = pX
        self.lnX = SP.log(pX)
        self.E1 = SP.zeros(pX.shape)
        self.update_moments()

    def update_moments(self):
        self.lnX = self.lnX - self.lnX.max()
        self.E1 = SP.exp(self.lnX)/(SP.exp(self.lnX).sum()) 

    def update(self, net, debug=False):
        LOG.debug("Calculating log-likelihoods of QTL location")
        for q in range(net.L): # correct?
            tau = 1./(net.F.var[q] + 1./net.tau_d)
            self.lnX[q] = (-0.5*tau*(net.F.E1[q] - net.mu_d)**2 + 0.5*SP.log(tau)).sum()
        self.update_moments()


class cRNode:
    def __init__(self, locus_positions, r_init_mean, r_init_var, init_from_smoothed, F1_init, F0):
        self.L = len(locus_positions)
        self.dist = SP.zeros([self.L, self.L])
        for l1 in range(self.L):
            for l2 in range(l1, self.L):
                self.dist[l1,l2] = self.dist[l2,l1] = abs(locus_positions[l1] - locus_positions[l2])
        self.rm, self.rvar = SP.zeros([self.L,self.L]), SP.zeros([self.L,self.L])
        self.rhom, self.rhovar = SP.zeros([self.L,self.L]), SP.zeros([self.L,self.L])
        if init_from_smoothed and (F1_init is not None):
            self.update(F0, F1_init, abs(F1_init - F0).argmax())
        else:
            self.update_estimates(r_init_mean, r_init_var)


    """ Calculate required moments for genetic map """
    def update_estimates(self, r_basemean, r_basevar):
        LOG.debug("Calculating rho expectations - basemean=%.2e"%(r_basemean))
        self.r_init, self.r_basemean, self.r_basevar = r_basemean, r_basemean, r_basevar  # base total rates and variance. Assume Poisson variance (= mean)
        b = self.r_basemean/(self.r_basevar + 1e-20) # Gamma distribution parameters
        a = b*self.r_basemean
        self.rm = self.dist*self.r_basemean # raw event rate between loci
        self.rvar = self.dist*self.r_basevar
        self.rhom = 0.5*(1 - (b/(b + 2*self.dist))**a) # p(g1 = g2) for each pair of loci
        self.rhovar = 0.25*((b/(b + 4*self.dist))**a - (b/(b + 2*self.dist))**(2*a))
        self.gm, self.gvar = SP.zeros(self.rhom.shape)*SP.nan, SP.ones(self.rhom.shape)*SP.nan
        
        if (self.rhom < 0).any() or (self.rhovar < 0).any():
            LOG.error("average recombination estimates negative. %s %s"%(str(self.rhom), str(self.rhovar)))
            pdb.set_trace()
            pass


    def update(self, F0, F1, q=None, debug=False):
        rhos = SP.zeros(self.L-1)            
        for l in range(1,self.L): # calculate estimate of recombination rate for each pair of adjacent loci
            if q >= l: l1,l2 = l-1,l
            else: l1,l2 = l, l-1
            a = -0.5/(F0[l2]*(1 - F0[l2]))
            b = 0.5*(F0[l1])/F0[l2] + 0.5*(1.-F0[l1])/(1. - F0[l2])
            c = 0.5/(1 - F0[l2])
            d = 0.5*(F0[l1] - F0[l2])/(1 - F0[l2])
            rhos[l-1] = (F1[l1] - b*F1[l2] - d)/(a*F1[l2] + c)
        self.r_basemean = rhos.mean()
        if self.r_basemean < 0:
            if debug:  pdb.set_trace()
            else: self.r_basemean = self.r_init
        self.update_estimates(self.r_basemean, self.r_basemean/100.) # assume low variance in recombination rate - alternative is bad


class cSQtlModel:
    
    def __init__(self, D, r_init_mean, r_init_var, F0, F1_init=None, pX=None, locus_positions=None, qtl_loci=None, calc_all_afs=False, n_iterations=3, init_from_smoothed=True):
        # data
        self.D = D
        eps_tau = 0.5 # Prior of two prior observations of both alleles
        self.mu_d = 1.*(self.D[:,0]+eps_tau)/(self.D.sum(axis=1)+2*eps_tau) # mean and variance of data
        self.tau_d = 1./((self.D[:,0]+eps_tau)*(self.D[:,1]+eps_tau)/((self.D.sum(axis=1)+2*eps_tau)*(self.D.sum(axis=1)+2*eps_tau)*(self.D.sum(axis=1) + 1 + 2*eps_tau))) # precision of data

        # Locations of sites
        self.L = len(D)   
        self.locus_positions = locus_positions
        if locus_positions is None: self.locus_positions = SP.arange(self.L)
        self.qtl_loci = qtl_loci
        if qtl_loci is None: qtl_loci = range(self.L)

        # Initial parameters
        if pX is None: pX = 1.*SP.ones(self.L)/self.L
        self.F0 = F0
        self.n_iterations = n_iterations
        self.r_init_mean, self.r_init_var = r_init_mean, r_init_var
        self.init_from_smoothed = init_from_smoothed

        # Nodes
        self.F = cFNode(self.D)
        self.X = cXNode(pX)
        self.R = cRNode(self.locus_positions, r_init_mean, r_init_var, init_from_smoothed, F1_init, F0)
        

    # calculate allele frequency estimates at locus l1 given data at locus l2, and whether l1 is between qtl and l2, or l2 between qtl and l1 (assume qtl most extreme)
    def calc_estimates(self):
        LOG.debug("Calculating model coefficients")
        # Model constants
        a = SP.zeros([self.L, self.L])
        b = SP.zeros([self.L, self.L])
        c = SP.zeros([self.L, self.L])
        d = SP.zeros([self.L, self.L])
        self.alpha = SP.zeros([2, 2, self.L, self.L])*SP.nan
        self.alphavar = SP.zeros([2, 4, self.L, self.L])*SP.nan
        for l1 in range(self.L):
            a[l1,:] = -0.5/(self.F0*(1 - self.F0))
            b[l1,:] = 0.5*(self.F0[l1])/self.F0 + 0.5*(1.-self.F0[l1])/(1. - self.F0)
            c[l1,:] = 0.5/(1 - self.F0)
            d[l1,:] = 0.5*(self.F0[l1] - self.F0)/(1 - self.F0)

        # for alpha,alphavar - f_{1,l} = alpha[0,0,l,q]f_{1,q} + alpha[0,1,l,q] if l is towards QTL (l is "higher", or more diverged from f_0)
        # f_{1,l} = alpha[1,0,l,q]f_{1,q} + alpha[1,1,l,q] if q is towards QTL (l is "lower", or less diverged from f_0)
        self.alpha[1,0,:,:] = a*self.R.rhom + b # alpha[0, :, l1, l2] = coefficients for estimating f if QTL is towards l2 compared to l1 from l1
        self.alpha[1,1,:,:] = c*self.R.rhom + d
        self.alphavar[1,0,:,:] = (a**2)*(self.R.rhovar + self.R.rhom**2) + b**2
        self.alphavar[1,1,:,:] = (a**2)*self.R.rhovar
        self.alphavar[1,2,:,:] = 2*a*c*self.R.rhovar
        self.alphavar[1,3,:,:] = self.R.rhovar*(c**2)


    def infer(self, debug=False):
        for i in range(self.n_iterations):
            if i > 0:
                q = self.X.E1.argmax()
                self.R.update(F0=self.F0, F1=self.F.E1[q], q=q, debug=debug)
            self.calc_estimates() # estimates of F given 
            self.F.update(net=self, debug=debug)
            self.X.update(net=self, debug=debug)
