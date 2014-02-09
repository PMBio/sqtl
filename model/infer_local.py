import scipy as SP
import scipy.stats as ST
import time as time
import pdb
from sqtl.model.stats import *
from sqtl.tools.common import *

PREDICT_FROM_RIGHT, PREDICT_FROM_LEFT = [0,1]
MU, TAU = [0,1]

class cFNode:
    def __init__(self, D):
        self.D = D.copy()
        self.a = SP.zeros([D.shape[0], D.shape[0], D.shape[1]],float)
        self.a[:] = D
        self.E1 = (self.a[:,:,0] + EPS)/(self.a.sum(axis=2) + 2*EPS)
        self.var = self.a.prod(axis=2)/((self.a.sum(axis=2)**2 + EPS)*(self.a.sum(axis=2) + 1))


class cXNode:
    def __init__(self, pX):
        self.pX = pX
        self.lnX = SP.log(pX)
        self.E1 = SP.zeros(pX.shape)
        self.update_moments()

    def update_moments(self):
        self.lnX = self.lnX - self.lnX.max()
        self.E1 = SP.exp(self.lnX)/(SP.exp(self.lnX).sum()) 


class cSQtlModelMP:
    
    def __init__(self, D, rho_basemean, rho_basevar, F0, pX=None, locus_positions=None, qtl_loci=None, calc_all_afs=False, rho=None):
        # data
        self.D = D + 1
        self.L = len(D)        
        self.locus_positions = locus_positions
        if locus_positions is None: self.locus_positions = SP.arange(self.L)
        if pX is None: pX = 1.*SP.ones(self.L)/self.L
        self.F0 = F0
        self.rho_basemean, self.rho_basevar, self.rho = rho_basemean, rho_basevar, rho
        self.alpha = None #, self.f_est, self.tau = None, None, None
        if qtl_loci is None: qtl_loci = range(self.L)
        self.qtl_loci = qtl_loci
        self.calc_all_afs = calc_all_afs
        self.n_iterations = 2
        self.message = SP.zeros([2,self.L,self.L,2])

        # nodes
        self.F = cFNode(self.D)
        self.X = cXNode(pX)

        # precalculate convenience variables
        self.precalc_rho(rho)
        self.precalc_constants()
        self.precalc_estimates()
        self.precalc_tau()


    """ Calculate required moments for genetic map """
    def precalc_rho(self, rho=None):
        t0 = time.time()
        b = self.rho_basemean/(self.rho_basevar + 1e-20)
        a = b*self.rho_basemean
        LOG.debug("Init rho - basemean=%.5f basevar=%.5f a=%.5f b=%5f"%(self.rho_basemean, self.rho_basevar, a, b))        
        self.dist = SP.zeros([self.L, self.L])
        for l1 in range(self.L):
            for l2 in range(l1, self.L):
                self.dist[l1,l2] = self.dist[l2,l1] = abs(self.locus_positions[l1] - self.locus_positions[l2])
        self.rm = self.dist*self.rho_basemean
        self.rvar = self.dist*self.rho_basevar
        self.rhom = 0.5*(1 - (b/(b + 2*self.dist))**a)
        self.rhovar = 0.25*((b/(b + 4*self.dist))**a - (b/(b + 2*self.dist))**(2*a))
        if self.rho_basevar <= 0: self.rhovar[:,:] = 0
        self.gm, self.gvar = SP.zeros(self.dist.shape)*SP.nan, SP.ones(self.dist.shape)*SP.nan
        if (self.rhom < 0).any() or (self.rhovar[l1,l2] < 0).any():
            LOG.error("average recombination estimates negative. l1=%d l2=%d rm=%.1e rvar=%.1e gm=%.1e gvar=%.1e"%(l1, l2, self.rhom[l1,l2], self.rhovar[l1,l2], self.gm[l1,l2], self.gvar[l1,l2]))
            pdb.set_trace()
            pass
        if rho != None: # if fixed recombination rate with no variation
            cr = SP.cumsum(rho)
            for l1 in range(self.L):
                for l2 in range(l1, self.L):
                    self.rhom[l1,l2] = 0.5*(1-SP.exp(-2*(cr[l2] - cr[l1])))
                    self.rhovar[l1,l2] = 0

        LOG.debug("Done init rho, time=%.2f"%(time.time() - t0))


    """ Initialise inference parameters. """
    def precalc_constants(self):
        LOG.debug("Calculating constants")
        t0 = time.time()
        a = SP.zeros([self.L, self.L])
        b = SP.zeros([self.L, self.L])
        c = SP.zeros([self.L, self.L])
        d = SP.zeros([self.L, self.L])
        self.alpha = SP.zeros([2, 2, self.L, self.L])*SP.nan
        self.alphavar = SP.zeros([2, 4, self.L, self.L])*SP.nan
        self.message[:,:,:,:] = 0

        for l1 in range(self.L):
            a[l1,:] = -0.5/(self.F0*(1 - self.F0))
            b[l1,:] = 0.5*(self.F0[l1])/self.F0 + 0.5*(1.-self.F0[l1])/(1. - self.F0)
            c[l1,:] = 0.5/(1 - self.F0)
            d[l1,:] = 0.5*(self.F0[l1] - self.F0)/(1 - self.F0)

        # for alpha,alphavar - f_{1,l} = alpha[0,0,l,q]f_{1,q} + alpha[0,1,l,q] if l is towards QTL (l is "higher", or more diverged from f_0)
        # f_{1,l} = alpha[1,0,l,q]f_{1,q} + alpha[1,1,l,q] if q is towards QTL (l is "lower", or less diverged from f_0)
        self.alpha[1,0,:,:] = a*self.rhom + b # alpha[0, :, l1, l2] = coefficients for estimating f if QTL is towards l2 compared to l1 from l1
        self.alpha[1,1,:,:] = c*self.rhom + d
        self.alphavar[1,0,:,:] = (a**2)*(self.rhovar + self.rhom**2) + b**2
        self.alphavar[1,1,:,:] = (a**2)*self.rhovar
        self.alphavar[1,2,:,:] = 2*a*c*self.rhovar
        self.alphavar[1,3,:,:] = self.rhovar*(c**2)
        LOG.debug("Done calculate constants, t=%.2f"%(time.time() - t0))
        return


    def precalc_tau(self):
        LOG.debug("Calculating tau")
        t0 = time.time()
        eps_tau = 0.5 # Prior of two prior observations of both alleles
        self.mu_d = 1.*(self.D[:,0]+eps_tau)/(self.D.sum(axis=1)+2*eps_tau) # mean of data
        self.tau_d = 1./((self.D[:,0]+eps_tau)*(self.D[:,1]+eps_tau)/((self.D.sum(axis=1)+2*eps_tau)*(self.D.sum(axis=1)+2*eps_tau)*(self.D.sum(axis=1) + 1 + 2*eps_tau))) # precision of data
        tau_fq = 1./self.f_est[True, :, :, 1] # precision of estimates
        self.tau_tilde = SP.zeros([self.L, self.L]) # l1 from data, l2 from q
        for l1 in range(self.L):
            self.tau_tilde[l1,:] = self.tau_d[l1]*tau_fq[l1,:]/(self.tau_d[l1] + tau_fq[l1,:]) # t = t1t2/(t1+t2) where t1 is from data, t2 is from prediction from q
        LOG.debug("Done calculating tau, t=%.2f"%(time.time() - t0))
        

    # calculate allele frequency estimates at locus l1 given data at locus l2, and whether l1 is between qtl and l2, or l2 between qtl and l1 (assume qtl most extreme)
    def precalc_estimates(self):
        LOG.debug("Calculating allele frequency estimates")
        t0 = time.time()
        self.f_est = SP.zeros([2, self.L, self.L, 2])*SP.nan
        for l2 in range(self.L):
            f_mean = 1.*(self.D[l2,0]+1)/(self.D[l2].sum() + 2)
            f_var = 1.*(self.D[l2,0]+1)*(self.D[l2,1]+1)/((self.D[l2].sum()+2)*(self.D[l2].sum()+2)*(self.D[l2].sum()+3))
            for l1 in range(self.L):
                for qtl_l2 in [True]: # the others require numerical integration, and will only be filled in on demand
                    self.f_est[qtl_l2, l1, l2] = get_simple_estimate(self, f_mean, f_var, l1, l2, qtl_l2)
        LOG.debug("Done calculating allele frequency estimates, t=%.2f"%(time.time() - t0))


    def recalc_rho(self):
        L = self.X.E1.shape[0]
        gamma_b = self.rho_basemean/(self.rho_basevar + 1e-20)
        gamma_a = gamma_b*self.rho_basemean
        allrhos, bhat = [],[]
        q = self.X.lnX.argmax()
        locs = SP.array(SP.arange(0,L, (L-1)/10), int)
        locs = range(L/6) + range(5*L/6, L)
        for i,l1 in enumerate(locs):#range(L): # average of five spots
            rho, l2 = -1, -1
            if (q > l1) and (i+1 < len(locs)) and (q > locs[i+1]):
                l2 = locs[i+1]
            elif (q < l1) and (i > 0) and (q < locs[i-1]):
                l2 = locs[i-1]
            if l2 < 0: continue
            a = -0.5/(self.F0[l2]*(1 - self.F0[l2]))
            b = 0.5*(self.F0[l1])/self.F0[l2] + 0.5*(1.-self.F0[l1])/(1. - self.F0[l2])
            c = 0.5/(1 - self.F0[l2])
            d = 0.5*(self.F0[l1] - self.F0[l2])/(1 - self.F0[l2])
            rho = (self.F.E1[q,l1] - b*self.F.E1[q,l2] - d)/(a*self.F.E1[q,l2] + c)
            allrhos.append(-0.5*SP.log(1 - 2*rho))
            if (rho > 0):
                delta = abs(l1 - l2)
                k = SP.exp(SP.log(1 - 2*rho)/gamma_a)
                bhat.append(2*delta*k/(1-k))
                #print self.rhom[l1,l2], rho, gamma_a/bhat[-1], bhat[-1]
        gamma_b = SP.exp(ST.nanmean(SP.log(bhat)))
        self.rho_basemean = gamma_a/gamma_b
        self.rho_basemean = SP.median(allrhos)
        if self.rho_basemean < 0: pdb.set_trace()
        #self.rho_basevar = gamma_a/(gamma_b**2)
        LOG.debug("Update rho - basemean=%.5f basevar=%.5f a=%.5f b=%5f"%(self.rho_basemean, self.rho_basevar, gamma_a, gamma_b))        
        self.precalc_rho()
        self.precalc_constants()
        self.precalc_estimates()
        self.precalc_tau()


    def infer(self, debug=False):
        for i in range(self.n_iterations+10):            
            if i > 0:
                self.recalc_rho()
            self.update_LX(debug=debug)
            

    def update_LX(self, debug=False):
        t0 = time.time()
        LOG.debug("Calculating allele frequencies")
        for q in range(self.L):
            for l in range(self.L):
                # first, left to right
                m1,t1,m2,t2 = self.mu_d[l], self.tau_d[l], 0, 0
                if l > 0:
                    m2,t2 = get_simple_estimate(self, self.message[PREDICT_FROM_LEFT, q, l-1, 0], 1./self.message[PREDICT_FROM_LEFT, q, l-1, 1], l, l-1, q < l)
                    t2 = 1./t2
                self.message[PREDICT_FROM_LEFT, q, l] = (m1*t1 + m2*t2)/(t1+t2), t1+t2
                # then, right to left from right hand side
                l2 = self.L - l - 1
                if q == 0 and l2 == 72:
                    pass#pdb.set_trace()
                m1,t1,m2,t2 = self.mu_d[l2], self.tau_d[l2], 0, 0
                if l2 < self.L - 1:
                    if l2 == 248 and q == 0: pass #pdb.set_trace()
                    m2,t2 = get_simple_estimate(self, self.message[PREDICT_FROM_RIGHT, q, l2+1, 0], 1./self.message[PREDICT_FROM_RIGHT, q, l2+1, 1], l2, l2+1, q > l2)
                    t2 = 1./t2
                #if q == 188: print l, m1, t1, m2, t2
                #if t1 > 1e8 or t2 > 1e8: pdb.set_trace()
                self.message[PREDICT_FROM_RIGHT, q, l2]  = (m1*t1 + m2*t2)/(t1+t2), t1+t2
        LOG.debug("Calculating log-likelihoods of QTL location")
        for q in range(self.L):
            for l in range(self.L):
                m1,t1 = self.message[PREDICT_FROM_LEFT, q, l]
                m2,t2 = self.message[PREDICT_FROM_RIGHT, q, l]
                m3,t3 = self.mu_d[l],self.tau_d[l]
                self.F.E1[q,l] = (m1*t1+m2*t2+m3*t3)/(t1+t2+t3)
                self.F.var[q,l] = 1./(t1+t2+t3)
            tau = 1./(self.F.var[q] + 1./self.tau_d)
            self.X.lnX[q] = (-0.5*tau*(self.F.E1[q] - self.mu_d)**2 + 0.5*SP.log(tau)).sum()
        self.X.update_moments()
        if debug:
            import pylab as PL
            for q in [188, self.X.lnX.argmax()]:
                PL.plot(self.F.E1[q,:])
            PL.plot(self.mu_d, 'k.')
            PL.plot(0.5+0.5*self.X.lnX/self.X.lnX.max(), 'r--')
            PL.show()
            pdb.set_trace()
            pass
        LOG.debug("Done sQTL(MP) inference, t=%.2f"%(time.time() - t0))
