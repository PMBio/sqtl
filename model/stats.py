import scipy as SP
import pdb

# functions and approximations 
C_LIMIT = 49
S_LIMIT = 1e-5
GAMMA = 0.577215664901532860606512090082
SMALL_MU = 1E-5
LARGE_TAU = 1E5
EPS = 1E-10


def partial_loglikelihood(net,x,q):
    log_res = 0 # keep track of answer in log scale to avoid numerical issues
    for l in range(net.L): # have to do in the loop, as x is a vector that's incompatible with other dimensions
        tt = 1./(1./net.tau_d[l] + 0.25*net.rhovar[l,q]*(x-net.F0[q])*(x-net.F0[q])/((net.F0[q]**2)*(1-net.F0[q])**2))
        log_res += 0.5*SP.log(tt)
        mu_e = (net.alpha[True,0,l,q]*x + net.alpha[True,1,l,q])
        log_res -= 0.5*tt*(net.mu_d[l] - mu_e)**2 # squared difference in the estimates from data and QTL, scaled by precision
        if SP.isnan(log_res).any(): pdb.set_trace()
        pass
    return log_res

def partial_loglikelihood_separate(net,x,q):
    log_res = SP.zeros(net.L) # keep track of answer in log scale to avoid numerical issues
    for l in range(net.L):
        tt = 1./(1./net.tau_d[l] + 0.25*net.rhovar[l,q]*(x-net.F0[q])*(x-net.F0[q])/((net.F0[q]**2)*(1-net.F0[q])**2))
        log_res[l] += 0.5*SP.log(tt)
        log_res[l] -= 0.5*tt*(net.mu_d[l] - net.alpha[True,0,l,q]*x - net.alpha[True,1,l,q])**2
        #if q in [188,198] and x in [0.92,0.95,0.98] and l == 189:
        #    print q, x, tt, log_res[l], "Fd=%.3f Fe=%.3f"%(net.mu_d[l], net.alpha[True,0,l,q]*x + net.alpha[True,1,l,q])
    return log_res


def calc_posterior_estimate(net, q, l):
    i, total_prec, last_prec = 0, 1E-6, SP.ones(2)*(1E-6) # distance from l, total precision of all estimates, last precision
    means, precs = [], []

    while i < net.L and (last_prec > 0.001*total_prec).any() and ((i < 10) or (last_prec < 0.1*precs[0]).all()): # for each offset from l1, while the last addition made a difference, and it does not get overly confident and cocky about things far away
        for j, direction in enumerate([[0], [1, -1]][i > 0]): # in both directions, but only if distance from QTL > 0
            l2 = l + direction*i 
            if l2 < 0 or l2 >= net.L: continue # skip ones outside limits
            if i > 1 and last_prec[(direction + 1)/2] <= 0.001*total_prec: continue # skip uninformative ones
            #if l == 175 and l2 == 184 and q == 183: pdb.set_trace()
            m,v = get_estimate(net, l, l2, q)  # Estimate mean and variance of the allele frequency of l1 from l2
            last_prec[(direction + 1)/2] = 1./v # Store statistics of estimates
            total_prec += 1./v
            means.append(m)
            if SP.isnan(m): pdb.set_trace()
            precs.append(1./v)
        i += 1
    return get_beta_params(SP.dot(means,precs)/total_prec, 1./total_prec) # Once all sites have made estimates, get parameters of the posterior 


# Calculate numbers necessary to estimate frequency at l1 given data at l2.
def calc_qtlside_estimates(net, l1, l2):
    net.gm[l1,l2] = get_gmean(net.rm[l1,l2], net.rvar[l1,l2], net.F0[l1] + net.F0[l2] - 2*net.F0[l1]*net.F0[l2], n_bins=200)
    net.gvar[l1,l2] = get_gvar(net.rm[l1,l2], net.rvar[l1,l2], net.F0[l1] + net.F0[l2] - 2*net.F0[l1]*net.F0[l2], net.gm[l1,l2],n_bins=200)
    a, b = 2.*net.F0[l1]*(1 - net.F0[l1]), 0
    c,d  = -2.*net.F0[l1]*(1 - net.F0[l1])*net.F0[l2], net.F0[l1]
    net.alpha[0,0,l1,l2] = a*net.gm[l1,l2] + b # alpha[0, :, l1, l2] = coefficients for estimating f if QTL is towards l2 compared to l1 from l1
    net.alpha[0,1,l1,l2] = c*net.gm[l1,l2] + d
    net.alphavar[0,0,l1,l2] = (a**2)*(net.gvar[l1,l2] + net.gm[l1,l2]**2) + b**2
    net.alphavar[0,1,l1,l2] = (a**2)*net.gvar[l1,l2]
    net.alphavar[0,2,l1,l2] = 2*a*c*net.gvar[l1,l2]
    net.alphavar[0,3,l1,l2] = net.gvar[l1,l2]*(c**2)
    

# estimate allele frequency posterior at locus l1 given allele frequency mean f_mean and variance f_var at locus l2. If qtl_l2, l2 is between l1 and qtl 
def get_simple_estimate(net, f_mean, f_var, l1, l2, qtl_l2):
    if (not qtl_l2) and SP.isnan(net.alpha[qtl_l2, 0, l1, l2]):
        calc_qtlside_estimates(net, l1, l2)
    est_mean = net.alpha[qtl_l2, 0, l1, l2]*f_mean +  net.alpha[qtl_l2, 1, l1, l2]# l1 - locus we are looking at now, estimating frequency from l2
    est_var = f_var*net.alphavar[qtl_l2,0,l1,l2] + f_mean*f_mean*net.alphavar[qtl_l2,1,l1,l2] + f_mean*net.alphavar[qtl_l2,2,l1,l2] + net.alphavar[qtl_l2, 3,l1,l2] #(net.alpha[qtl_l1, 0, l2, l1]**2)
    if (not qtl_l2) and (f_mean - net.F0[l2] > est_mean - net.F0[l1]): pdb.set_trace() 
    #if est_var < 1e-8: pdb.set_trace()
    #if SP.isnan(est_mean) or SP.isnan(est_var): pdb.set_trace()
    return est_mean, est_var


# estimate allele frequency posterior at locus l1 from data at locus l2, regardless of position of qtl locus lx
def get_estimate(net, l1, l2, lx):
    qtl_l1 = (lx <= l1 and lx <= l2 and l1 <= l2) or (lx >= l1 and lx >= l2 and l1 >= l2) # qtl towards l1, and not between l1 and l2
    qtl_l2 = (lx <= l1 and lx <= l2 and l1 >= l2) or (lx >= l1 and lx >= l2 and l1 <= l2) # qtl towards l2, and not between l1 and l2
    qtl_between = not(qtl_l1 or qtl_l2)
    if qtl_between:
        calc_qtlside_estimates(net, lx, l2) # If QTL between two loci, precalculate prediction of lx from l2. These are NaN by default.
    elif not qtl_l2: calc_qtlside_estimates(net, l2, l1) # If QTL not between two loci, but towards l1 that we are trying to predict, precalculate l2 from l1. These are NaN by default if QTL not towards l2.

    if qtl_between:  # if QTL between two loci, 
        f_mean = 1.*(net.D[l2,0] + 0.1)/(net.D[l2].sum() + 0.2)
        f_var = 1.*(net.D[l2,0] + 0.1)*(net.D[l2,1] + 0.1)/((net.D[l2].sum() + 0.2)*(net.D[l2].sum() + 0.2)*(net.D[l2].sum() + 2.2))
        #f_mean, f_var = net.f_est[True, lx, l2] # estimate the QTL allele frequency from observed locus l2
        fq_mean, fq_var = net.f_est[False, lx, l2] = get_simple_estimate(net, f_mean, f_var, lx, l2, False)
        return get_simple_estimate(net, fq_mean, fq_var, l1, lx, True) # and calculate estimate at desired locus l1 from the result.
    else: # return simple estimate
        if SP.isnan(net.f_est[qtl_l2, l1, l2,0]):
            f_mean = 1.*(net.D[l2,0] + 0.1)/(net.D[l2].sum() + 0.2)
            f_var = 1.*(net.D[l2,0] + 0.1)*(net.D[l2,1] + 0.1)/((net.D[l2].sum() + 0.2)*(net.D[l2].sum() + 0.2)*(net.D[l2].sum() + 2.2))
            net.f_est[qtl_l2, l1, l2] = get_simple_estimate(net, f_mean, f_var, l1, l2, qtl_l2)            
        return net.f_est[qtl_l2, l1, l2]


# estimate allele frequency posterior at locus l1 from data at locus l2, regardless of position of qtl locus lx
def get_estimate_old(net, l1, l2, lx):
    qtl_l2 = (lx <= l1 and lx <= l2 and l1 >= l2) or (lx >= l1 and lx >= l2 and l1 <= l2) # qtl towards l2

    if (l1 < lx and lx < l2) or (l2 < lx and lx < l1):  # if QTL between two loci, 
        f_mean, f_var = net.f_est[qtl_l2, lx, l2] # estimate the QTL allele frequency from observed locus
        return get_simple_estimate(net, f_mean, f_var, l1, lx, True)
    else: # return simple estimate
        return net.f_est[qtl_l2, l1, l2]


def get_gvar(rm,rv, a, gm, min_var=1E-4, n_bins=20):
    gamma_b = rm/(rv + 1e-20)
    gamma_a = gamma_b*rm
    rhohat = 0.5*(1-SP.exp(-2.*rm))
    rsd = rv**0.5
    #rv = 0.
    if rv < 1E-4 and abs(a/(rsd+1e-9)) > 7: 
        return max(min_var, ((a-rhohat)**(-2) - gm**2))
    n_sd = 5
    lower = max(0,rm - n_sd*rsd)
    upper = rm + n_sd*rsd
    lower, upper = rm - n_sd*rsd, rm + n_sd*rsd
    if gm < 0: return 10000
    max_var = 1000
    delta = (upper-lower)/n_bins
    rhofun = lambda x: 0.5*(1-SP.exp(-2.*x))
    #expfun = lambda x:gamma_dist(x,gamma_a, gamma_b)*((a - rhofun(x))**(-2))
    expfun = lambda x:norm_dist(x,rm, rv)*((a - rhofun(x))**(-2))
    total_area = region_cumnorm_dist(lower+delta, upper, rm, rsd) + region_cumnorm_dist(upper, upper + 5*rsd, rm, rsd)
    res = (approx_int(expfun, lower+delta, upper, delta) + approx_int(expfun, upper, upper+5*rsd, rsd/4.))/total_area - gm**2
    if res > max_var or res < min_var: pdb.set_trace()
    return min(max_var, max(res, min_var))


def get_gmean(rm,rv,a, n_bins=20):
    rhohat = 0.5*(1-SP.exp(-2.*rm))
    gamma_b = rm/(rv + 1e-20)
    gamma_a = gamma_b*rm
    if gamma_a < 0 or gamma_b < 0: pdb.set_trace()
    rhofun = lambda x: 0.5*(1-SP.exp(-2.*x))
    #expfun = lambda x:gamma_dist(x,gamma_a, gamma_b)/(a - rhofun(x))
    expfun = lambda x:norm_dist(x,rm, rv)/(a - rhofun(x))
    res = 1./(a-rhohat)
    if SP.isnan(res): pdb.set_trace()
    #rv = 0.
    if True or rv < 1E-4: return res
    rsd, n_sd = rv**0.5, 5
    #lower, upper = max(0, rm - n_sd*rsd), rm + n_sd*rsd
    lower, upper = rm - n_sd*rsd, rm + n_sd*rsd
    if (a - rhofun(lower))*(a - rhofun(upper)) > 0: # if the to-be-integrated crosses 0 in denominator and blows up        
        delta = (upper-lower)/n_bins
        total_area = region_cumnorm_dist(lower+delta, upper, rm, rsd) + region_cumnorm_dist(upper, upper + 5*rsd, rm, rsd)
        res = (approx_int(expfun, lower+delta, upper, delta) + approx_int(expfun, upper, upper + 5*rsd, rsd/4.))/total_area
    if res < 0 or SP.isnan(res):
        print rm,rv, a, lower, upper, rsd
        pdb.set_trace()
    return res


""" Calculate log(y1+y2+y3) given x1=log(y1), x2=log(y2), ...; using
log(y1+y2+y3) = log(exp(log(y1)) + exp(log(y2)) + ...) = m + log(exp(x1)/exp(m) + exp(x2)/exp(m) + ...) = m + log(exp(x1 - m) + exp(x2 - m) + ...)
"""
def logsum(x):
    x = SP.array(x)
    maxlog = x.max()
    return maxlog + SP.log(SP.exp(x - maxlog).sum())

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

def loggamma(x):
    tmp = (x - 0.5)*SP.log(x + 4.5) - (x + 4.5)
    ser = 1.0 + 76.18009173/(x + 0) - 86.50532033/(x + 1) + 24.01409822/(x + 2) - 1.231739516/(x + 3) + 0.00120858003/(x + 4) - 0.00000536382/(x + 5);
    return tmp + SP.log(ser) + 0.5*SP.log(2*SP.pi);

def log_gamma_dist(x,a,b):
    return a*SP.log(b) - loggamma(a) + (a-1)*SP.log(x) - b*x

def gamma_dist(x,a,b):
    return SP.exp(log_gamma_dist(x,a,b))

def norm_dist(x,m,s):
    return SP.exp(log_norm_dist(x,m,s))

def log_norm_dist(x,m,s):
    return -0.5*(SP.log(2*SP.pi) + SP.log(s) + ((x-m)**2)/s)

def region_cumnorm_dist(start, end, m, s):
    return cumnorm_dist(end,m,s) - cumnorm_dist(start,m,s)


def cumnorm_dist(x,m,v):
    flip = False
    s = abs(v)**0.5
    x0 = (x-m)/s 
    if x0 < 0: # need x0 > 0 - otherwise return 1 - Phi(-x0)
        x0 = -x0
        flip = True
    t = 1./(1+0.2316419*x0)
    b = [0.319381530, -0.356563782, 1.781477937, -1.821255978, 1.330274429]
    ts = [t, t*t, t*t*t, t*t*t*t, t*t*t*t*t]
    if flip: return norm_dist(x0,0,1)*SP.dot(b,ts)
    else: return 1 - norm_dist(x0,0,1)*SP.dot(b,ts)


# parameter estimates for Beta(a,b) such that the resulting distribution has specified mean m and variance v
def get_beta_params(m,v, min_count=None):
    if v < 0: v = 1E-4 # return SP.array([1E-6, 1E-6])
    if m < 0: v, m = v/2., 1E-4
    if m > 1: v, m = v/2., 1-1E-4
    t = m/v - m*m/v - 1
    if min_count is not None: t = max(t, min_count)
    return SP.array([t*m, t*(1-m)])

def log_beta_dist(x,a,b):
    res = (a-1)*SP.log(x) + (b-1)*SP.log(1-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
    if SP.isnan(res).any():
        print a,b,x
        pdb.set_trace()
    return (a-1)*SP.log(x) + (b-1)*SP.log(1-x) - loggamma(a) - loggamma(b) + loggamma(a+b)
