import scipy as SP
import pdb

# functions and approximations 
C_LIMIT = 49
S_LIMIT = 1e-5
GAMMA = 0.577215664901532860606512090082
SMALL_MU = 1E-5
LARGE_TAU = 1E5


def incomplete_beta_int_array(x1a,x2a,a,b,n_steps=100):
    if type(x1a) is not type(SP.zeros(2)): return incomplete_beta_int(x1a,x2a,a,b,n_steps)

    res = SP.zeros(x1a.shape)
    for i in range(len(x1a)):
        res[i] = incomplete_beta_int(x1a[i], x2a[i], a,b,n_steps)
    return res

def incomplete_beta_int(x1,x2,a,b, n_steps=100):
    if x2 > 1: x2 = 1.
    if x1 <= 0: x1 = 0.
    return approx_int(lambda t:SP.log(t)*(a-1) + SP.log(1-t)*(b-1), x1, x2, 1.*(x2-x1)/n_steps, True)


# ln(gamma(x))
def loggamma(x):
    tmp = (x - 0.5)*SP.log(x + 4.5) - (x + 4.5)
    ser = 1.0 + 76.18009173/(x + 0) - 86.50532033/(x + 1) + 24.01409822/(x + 2) - 1.231739516/(x + 3) + 0.00120858003/(x + 4) - 0.00000536382/(x + 5);
    return tmp + SP.log(ser) + 0.5*SP.log(2*SP.pi);


def digamma(x):
    is_array = type(x) is type(SP.ones(1))
    if not is_array: x = SP.ones(1)*x
    result = SP.zeros(x.shape)
    result[SP.where(x <= 0)[0]] = 0
    I1 = SP.where((x > 0) & (x <= S_LIMIT))[0]
    I2 = SP.where(x > S_LIMIT)[0]
    result[I1] = -GAMMA - 1./x[I1] # small value

    for i in I2:
        xb = x[i]
        while xb < C_LIMIT: 
            result[i] -= 1./xb
            xb += 1
        x2 = xb**(-2.)
        result[i] += SP.log(xb) - 0.5/xb - x2*(1./12 - x2*(1./120 - x2/252))

    if is_array: return result
    return result[0]


def digamma_num(x):
    if x > 0 & x <= S_LIMIT: return -GAMMA - 1./x[I] # small value
    if x >= C_LIMIT: # big enough value - close enough to log(x) + correction terms, remaining error O(x^-8)
        x2 = x**(-2.)
        return SP.log(x) - 0.5/x2 - inv*(1./12 + x2*(1./120 - x2/252))
    return digamma(x + 1.) - 1./x # else use recursive relation


def logbetafun(a,b): return -loggamma(a+b) + loggamma(a) + loggamma(b)


def betafun(a,b):
    return SP.exp(logbetafun(a,b))


# beta distribution
def beta(x, a):
    return (x**(a[0]-1))*((1-x)**(a[1]-1))/betafun(a[0],a[1])


def normal(x, m, v):
    return SP.exp(-0.5*((x-m)**2)/v)/((2.*SP.pi*v)**0.5)


""" log(x1+x2+x3) given log(x1), log(x2), ..."""
def logsum(x):
    x = SP.array(x)
    maxlog = x.max()
    x -= maxlog
    return maxlog + SP.log(SP.exp(x).sum())
    

# approximation of (int_start^stop f(x) dx) according to Simpson's rule
# result_logscale implies f_logscale
def approx_int(f, start, stop, delta, f_logscale=False, result_logscale=False):
    result = 0.
    xrange = SP.arange(start + delta,stop,delta) # N+1 delimiters of the N smaller chunks
    xrangemid = SP.arange(start + 0.5*delta, stop, delta) # midpoints of the N chunks

    if f_logscale:
        if result_logscale:
            r1 = logsum(f(xrange)) + SP.log(delta/3.)
            r2 = logsum(f(xrangemid)) + SP.log(2.*delta/3.)
            r3 = SP.log(delta/6.) + logsum([f(start),f(stop)])
            result = logsum([r1,r2,r3])
        else:
            result += (SP.exp(f(xrange))*delta/3.).sum()
            result += (2.*SP.exp(f(xrangemid))*delta/3.).sum()
            result += delta*(SP.exp(f(start)) + SP.exp(f(stop)))/6.
    else:
        result += (f(xrange)*delta/3.).sum()
        result += (2.*f(xrangemid)*delta/3.).sum()
        result += delta*(f(start) + f(stop))/6.

    return result


# approximation of (int_start^stop f(x) g(x) I(filter(x)) dx) / int_start^stop f(x)I(filter(x)) dx according to Simpson's rule
def approx_mean(fn3d, start, stop, delta):
    result, result_z = 0., 0.
    xrange = SP.arange(start + delta,stop,delta) # N+1 delimiters of the N smaller chunks
    xrangemid = SP.arange(start + 0.5*delta, stop, delta) # midpoints of the N chunks
    fval = SP.array(fn3d(xrange), float)
    fvalm = SP.array(fn3d(xrangemid), float)
    Inan = SP.where(SP.isnan(fval[1]))
    Inanm = SP.where(SP.isnan(fvalm[1]))
    fval[2][Inan] = False
    fval[1][Inan] = 0.
    fvalm[2][Inanm] = False
    fvalm[1][Inanm] = 0.

    result += ((fval.prod(axis=0))*delta/3.).sum()
    result += (2.*(fvalm.prod(axis=0))*delta/3.).sum()
    result += delta*(SP.prod(fn3d(start)) + SP.prod(fn3d(stop)))/6.
    result_z += (fval[0]*fval[2]*delta/3.).sum()
    result_z += (2.*fvalm[0]*fvalm[2]*delta/3.).sum()
    result_z += delta*(fn3d(start)[0]*fn3d(start)[2] + fn3d(stop)[0]*fn3d(stop)[2])/6.

    return result #/result_z



# Monte Carlo integration in two dimensions
def twod_mc_int(f, startx, stopx, starty, stopy, n_samples=100):
    rangex = stopx - startx
    rangey = stopy - starty
    locs = SP.rand(n_samples, 2)
    return SP.array([f(l[0]*rangex + startx, l[1]*rangey + starty) for l in locs]).sum()*rangex*rangey/n_samples

# Monte Carlo integration in two dimensions
def twod_mc_beta_int(a1, a2, f, n_samples=100):
    m1 = a1[0]/a1.sum()
    m2 = a2[0]/a2.sum()
    v1 = a1[0]*a1[1]/(a1.sum()**2)/(a1.sum()+1)
    v2 = a1[0]*a2[1]/(a2.sum()**2)/(a2.sum()+1)
    startx = max(0, m1 - 3*(v1**0.5))
    rangex = min(6*(v1**0.5), 1. - startx)
    starty = max(0, m2 - 3*(v2**0.5))
    rangey = min(6*(v2**0.5), 1. - starty)
    locs = SP.rand(n_samples, 2)
    print startx, rangex, starty, rangey
    return SP.array([f(l[0]*rangex + startx, l[1]*rangey + starty) for l in locs]).sum()*rangex*rangey/n_samples

# Monte Carlo integration in two dimensions using importance sampling; assume beta variables 
def twod_mcis_beta_int(a1,a2, f, n_samples=100):
    m1 = a1[0]/a1.sum()
    m2 = a2[0]/a2.sum()
    v1 = a1[0]*a1[1]/(a1.sum()**2)/(a1.sum()+1)
    v2 = a1[0]*a2[1]/(a2.sum()**2)/(a2.sum()+1)
    locs = SP.randn(n_samples, 2)
    return SP.array([f(l[0]*rangex + startx, l[1]*rangey + starty) for l in locs]).sum()*rangex*rangey/n_samples

 
# parameter estimates for Beta(a,b) such that the resulting distribution has specified mean m and variance v
def get_beta_params(m,v, min_count=None):
    if v < 0: v = 1E-4 # return SP.array([1E-6, 1E-6])
    if m < 0: v, m = v/2., 1E-4
    if m > 1: v, m = v/2., 1-1E-4
        
    t = m/v - m*m/v - 1
    if min_count is not None: t = max(t, min_count)
    return SP.array([t*m, t*(1-m)])


def get_gamma_params(m,v):
    return SP.array([m*m/v, m/v])

def get_beta_coeff_from_triv_normal(coeff, debug=False):
    tau = -2*coeff[:,:,1].sum(axis=0)
    mu = coeff[:,:,0].sum(axis=0)/tau
    if debug: print -2*coeff[:,0,1], -2*coeff[:,0,1].sum(axis=0), coeff[:,0,0]/(-2*coeff[:,0,1]), coeff[:,0,0].sum()/tau
#    if debug: pdb.set_trace()
    I = SP.where(mu < 0)[0]
    if len(I) > 0:
        print "beta mu < 0: %d sites"%(len(I))
        mu[I] = SMALL_MU
    I = SP.where(mu > 1)[0]
    if len(I) > 0:
        print "beta mu > 1: %d sites"%(len(I))
        mu[I] = 1 - SMALL_MU
    return get_beta_params(mu, 1./tau)


def get_beta_coeff_from_normal(coeff):
    tau = -2*coeff[:,1]
    mu = coeff[:,0]/tau
    I = SP.where(mu < 0)[0]
    if len(I) > 0:
        print "beta mu < 0: %d sites"%(len(I))
        mu[I] = SMALL_MU
    I = SP.where(mu > 1)[0]
    if len(I) > 0:
        print "beta mu > 1: %d sites"%(len(I))
        mu[I] = 1 - SMALL_MU
    return get_beta_params(mu, 1./tau)

def get_gamma_coeff_from_normal(coeff):
    tau = -2*coeff[:,1]
    mu = coeff[:,0]/tau
    I = SP.where(mu < 0)[0]
    if len(I) > 0:
        print "gamma mu < 0: %d sites"%(len(I))
        mu[I] = SMALL_MU
        tau[I] = LARGE_TAU
    return get_gamma_params(mu, 1./tau)



#### Model-specific helpers

max_rho_val = 1E3 # infinity for our purposes, used to keep integrals within check
min_rho_f_var = 0.001 # zero for our purposes, used to keep estimates from collapsing

# Beta estimate of P(r | f0q, f0l, f1q, f1l) according to q being between l and QTL 
def get_free_side_beta_approximation(f1q, f0l, f0q, (rho_m, rho_var)):
    rho_a, rho_b = rho_m*rho_m/rho_var, rho_m/rho_var # distribution of rho
    rho_f_a = 2*f0q*(1-f0q)/(f0q - f1q) # rho = a*f + b
    rho_f_b = -2*f0l*f1q*(1-f0q)/(f0q-f1q) + f0q - f0l
    rho_f_mean = (rho_a/rho_b - rho_f_b)/rho_f_a
    rho_f_var = rho_a/(rho_b*rho_b*rho_f_a*rho_f_a)
    rho_f_var = max(min_rho_f_var, rho_f_var)
    return get_beta_params(rho_f_mean, rho_f_var)


# Beta estimate of P(r | f0q, f0l, f1q, f1l) according to l being between q and QTL
def get_qtl_side_beta_approximation(f1q, f0l, f0q, (rho_m, rho_var)):
    rho_a, rho_b = rho_m*rho_m/rho_var, rho_m/rho_var # gamma parameters for distribution of rho
    delta = 0.001 # integration precision
    rho_f_a = f0q + f0l - 2.*f0q*f0l # r = (af + b)/(f + c)
    rho_f_b = f0l*(2.*f1q*f0l - 2.*f1q - f0l + f0q)
    rho_f_c = -f0l
    # r = (a*f + b)/(f+c) (ignoring inferred r too large or too small) 
    def rhofun(x):
        r = (rho_f_a*x + rho_f_b)/(x + rho_f_c + 1E-10)
        return (abs(r) < max_rho_val)*r

    # log distribution of r - gamma, with parameters calculated via rhofun
    logrhodist = lambda x: rho_a*SP.log(rho_b) - loggamma(rho_a) -rho_b*rhofun(x) + (rho_a - 1)*SP.log(rhofun(x))

    points = sorted([0,1,-rho_f_b/rho_f_a, -rho_f_c]) # intervals of f according to corresponding r. 
    m0,m1,m2 = 0,0,0

    for i in range(len(points) - 1): # estimate 3 moments
        if rhofun(0.5*(points[i] + points[i+1])) > 0: # if r is positive, include interval in moment estimation
            m0 += approx_int(lambda x:SP.exp(logrhodist(x)), points[i], points[i+1], delta)
            m1 += approx_int(lambda x:x*SP.exp(logrhodist(x)), points[i], points[i+1], delta)
            m2 += approx_int(lambda x:x*x*SP.exp(logrhodist(x)), points[i], points[i+1], delta)

    rho_f_mean = m1/m0 # calculate mean and variance from moments, and corresponding beta parameters
    rho_f_var = m2/m0 - rho_f_mean**2
    rho_f_var = max(min_rho_f_var, rho_f_var)
    return get_beta_params(rho_f_mean, rho_f_var)



# Numerically calculate int_x1 int_x2 Beta(x1; a1, b1) Beta(x2; a2, b2) ln(Gamma(cx/(d-y) + e + g/(d-y); a, b)))
# Approach: first, approximate Gamma with Beta matching moments in x; then integrate out x by evaluating <ln P'(x)>
# This results in a function of y.
# Finally, use Simpson's rule on the resulting function to calculate the integral over x2
def doubleint_betas_ln_approxgamma(a1,b1, a2,b2, rm, rv, c, d, e, g, fn='ln'):
    a, b = rm*rm/rv, rm/rv # gamma parameters for distribution of rho
    afun = lambda x: c/(d-x)
    bfun = lambda x: e + g/(d-x)
    def approx_tfun(x):
        aa, bb = afun(x), bfun(x)
        return - 1 - a + aa*b + 2*b*bb - aa*bb*b*b/a - bb*bb*b*b/a
    def approx_exp(x):
        aa = approx_tfun(x)*(a-b*bfun(x))/(afun(x)*b)
        bb = -aa + approx_tfun(x)
#        print "approx", aa, bb
        if fn == 'ln':
            return beta(x, (a2, b2)), (digamma(aa) - digamma(aa + bb)), 1. - (1.-(aa > 0)*(bb > 0)*(aa/(aa+bb)*afun(x) + bfun(x) > 0))*(1000.)
        else:
            return beta(x, (a2, b2)), aa/(aa+bb), 1. - (1.-(aa > 0)*(bb > 0)*(aa/(aa+bb)*afun(x) + bfun(x) > 0))*(1000.)
    return approx_mean(approx_exp, 0, 1, 0.01)
    

def doubleint_betas_entropy_approxgamma(a1,b1, a2,b2, rm, rv, c, d, e, g):
    a, b = rm*rm/rv, rm/rv # gamma parameters for distribution of rho
    afun = lambda x: c/(d-x)
    bfun = lambda x: e + g/(d-x)
    def approx_tfun(x):
        aa, bb = afun(x), bfun(x)
        return - 1 - a + aa*b + 2*b*bb - aa*bb*b*b/a - bb*bb*b*b/a
    def approx_exp(x):
        aa = approx_tfun(x)*(a-b*bfun(x))/(afun(x)*b) # parameters for approximate distribution of f_l
        bb = -aa + approx_tfun(x)
#        print "approx", aa, bb
        filter = (aa > 0)*(bb > 0)*(aa/(aa+bb)*afun(x) + bfun(x) > 0) # beta params ok and r > 0
        vals = logbetafun(aa,bb) - (aa-1)*digamma(aa) - (bb-1)*digamma(bb) + (aa+bb-2)*digamma(aa+bb)
        return beta(x, (a2, b2)), vals*filter - 0*(1-filter), abs(x) < SP.inf
    return approx_mean(approx_exp, 0, 1, 0.01)
    


def test_stats():
    print twod_mc_int(lambda x,y: x+y, 0, 3, 0, 2, 100)


def integrate(f, p, a, b, n_bins=20):
    delta = 1.*(b-a)/n_bins
    x = SP.arange(a + delta/2, b + delta/2, delta)
    return sum(f(x)*p(x)*delta)

def pnorm(x,m,v):
    return SP.exp(logpnorm(x,m,v))

def logpnorm(x,m,v):
    return -0.5*(x-m)*(x-m)/v - 0.5*SP.log(2.*3.141592652*v)

def cumnorm(x,m,v):
    flip = False
    x0 = (x-m)/(v**0.5) 
    if x0 < 0: # need x0 > 0 - otherwise return 1 - Phi(-x0)
        x0 = -x0
        flip = True
    t = 1./(1+0.2316419*x0)
    b = [0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429]
    ts = [t, t*t, t*t*t, t*t*t*t, t*t*t*t*t]
    if flip: return pnorm(x0,0,1)*SP.dot(b,ts)
    else: return 1 - pnorm(x0,0,1)*SP.dot(b,ts)

#TODO : calculate gmean, gvar exactly
def get_gvar(rm,rv,a, gm, min_var=1E-4, n_bins=20):
    rsd = rv**0.5
    gmhat = 0.5*(1-SP.exp(-2.*rm))
    if rv < 1E-6 and abs(a/(rsd+1e-9)) > 7: 
        return max(min_var, ((a-gmhat)**(-2) - gm**2))
    n_sd = 5
    rsd = rv**0.5
    lower = rm - n_sd*rsd
    upper = rm + n_sd*rsd
    total_area = 1
    if lower < 0:
        lower = 0
        total_area = cumnorm(2*rm,rm,rv)
    if gm < 0: return 10000

#    if (a-lower)*(a-upper) < 0: return 1
    max_var = 1000
    return min(max_var, max(integrate(lambda x:((a-0.5*(1-SP.exp(-2.*x)))**(-2)), lambda x:pnorm(x,rm,rv), lower, upper, n_bins=n_bins)/total_area - gm**2, min_var))


def get_gmean(rm,rv,a, n_bins=20):
    gmhat = 0.5*(1-SP.exp(-2.*rm))
    if rv < 1E-6: return 1./(a-gmhat)
    n_sd = 5
    rsd = rv**0.5
    lower = rm - n_sd*rsd    
    upper = rm + n_sd*rsd
    total_area = 1
    if lower < 0:
        lower = 0
        total_area = cumnorm(2*rm,rm,rv)
    res = 1./(a-gmhat)    
    if (a - lower)*(a - upper) > 0: res = integrate(lambda x:1./(a-0.5*(1-SP.exp(-2.*x))), lambda x:pnorm(x,rm,rv), lower, upper, n_bins=n_bins)/total_area
    if res < 0:
        print rm,rv, a, lower, upper, rsd
        pdb.set_trace()
#    if res < -1: 
#        print rm, rv, a, lower, upper, 1./(a-rm)
#        pdb.set_trace()
    return res



#TODO : calculate gmean, gvar exactly
def get_gvar_old(rm,rv,a, gm, min_var=1E-4, n_bins=20):
    rsd = rv**0.5
    if rv < 1E-6:# and abs(a/rsd) > 7: 
        return max(min_var, ((a-rm)**(-2) - gm**2))
    n_sd = 5
    rsd = rv**0.5
    lower = rm - n_sd*rsd
    upper = rm + n_sd*rsd
    total_area = 1
    if lower < 0:
        lower = 0
        total_area = cumnorm(2*rm,rm,rv)
    if gm < 0: return 10000

#    if (a-lower)*(a-upper) < 0: return 1
    max_var = 1000
    return min(max_var, max(integrate(lambda x:((a-x)**(-2)), lambda x:pnorm(x,rm,rv), lower, upper, n_bins=n_bins)/total_area - gm**2, min_var))


def get_gmean_old(rm,rv,a, n_bins=20):
    if rv < 1E-6: return 1./(a-rm)
    n_sd = 5
    rsd = rv**0.5
    lower = rm - n_sd*rsd    
    upper = rm + n_sd*rsd
    total_area = 1
    if lower < 0:
        lower = 0
        total_area = cumnorm(2*rm,rm,rv)
    res = 1./(a-rm)
    if (a - lower)*(a - upper) > 0: res = integrate(lambda x:1./(a-x), lambda x:pnorm(x,rm,rv), lower, upper, n_bins=n_bins)/total_area
    if res < 0:
        print rm,rv, a, lower, upper, rsd
        pdb.set_trace()
#    if res < -1: 
#        print rm, rv, a, lower, upper, 1./(a-rm)
#        pdb.set_trace()
    return res
