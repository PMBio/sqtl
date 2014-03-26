import scipy as SP


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


def region_cumnorm_dist(start, end, m, s):
    return cumnorm_dist(end,m,s) - cumnorm_dist(start,m,s)


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


def norm_dist(x,m,s):
    return SP.exp(log_norm_dist(x,m,s))

def log_norm_dist(x,m,s):
    return -0.5*(SP.log(2*SP.pi) + SP.log(s) + ((x-m)**2)/s)

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

