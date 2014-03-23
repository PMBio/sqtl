import pylab as PL
import scipy as SP

def plot_inferred(models, model_names, sim, step=0.02, name=None, title=None, X_logscale=False, qtl_loci=None):
    PL.figure()
    plot_sim(sim)
    colors = 'rc'
    tau_mean = 0
    f1_err = 0
    f10_err = 0
    qtl_mode = 0

    for i in range(len(models)): # AF estimates
        if len(models[i].F.E1.shape) == 1:
            PL.plot(models[i].F.E1, colors[i])
        else:
            if qtl_loci is None: qtl_loci = [models[i].F.E1.shape[0]/2]
            for q in qtl_loci:
                PL.plot(models[i].F.E1[q], colors[i])

    for i, m in enumerate(models): # QTL estimates
        if model_names[i] == 'Inferred F1':
            x = m.X.E1
            if X_logscale:
                x = (m.X.lnX - min(m.X.lnX))/(abs(m.X.lnX - min(m.X.lnX)).max())
            PL.plot(m.qtl_loci, x[m.qtl_loci], 'k')
            #PL.plot(models[i].Rho.E1*2 + 0.3, 'm--') # true rec rate
            #PL.plot(models[i].F.E1[1:] - models[i].F.E1[0:-1] + 0.2, 'y') # inferred change
            #tau_mean = models[i].Tau.E1
            f1_err = abs(m.F.E1[m.X.E1.argmax()] - sim.F1).mean()
            qtl_mode = m.X.E1[len(m.X.E1)/2]
        if model_names[i] == 'Smoothed F1':
            f10_err = abs(m.F.E1 - sim.F1).mean()

    PL.legend(['F0', 'True F1', 'Data'] + model_names)
    T = models[0].F.E1.shape[0]/50.
    PL.xticks([5*T,15*T,25*T,35*T,45*T], [-100,-50,0,50,100])
    PL.xlabel('Distance to QTL (cM)')
    PL.ylabel('A allele frequency')
    if title is None: title = ''
#    title = title + ' $\\langle \\tau \\rangle = %.3f$'%tau_mean
    title = title + ' $|\Delta f_1|=%.3f$ '%f1_err
    title = title + ' $|\Delta f^0_1|=%.3f$ '%f10_err
    title = title + ' $p(X)=%.3f$ '%qtl_mode
    PL.title(title)
    if name is not None: PL.savefig('plot_%s.pdf'%name)


def plot_single(af, sd=None, locs=None, color="b", linewidth=2, alpha_multiplier=0.1):
    if locs is None: locs = range(len(af))
    line = PL.plot(locs, af, color+'-', linewidth=linewidth)

    if sd != None:
        x_ci = SP.array(list(locs) + list(locs)[::-1])
        y_ci = SP.array(list(af) + list(af)[::-1])
        sds = SP.array(list(sd) + list(-sd)[::-1])
        PL.fill(x_ci, y_ci + sds, color, alpha=alpha_multiplier)
        PL.fill(x_ci, y_ci + 2*sds, color, alpha=2*alpha_multiplier) 

    return line


def plot_sim(sim):
    PL.plot(sim.F0, '--')
    PL.plot(sim.F1, '--')
    PL.plot(sim.D[:,0]/sim.D.sum(axis=1), 'r.', markersize=8)



def plot_prediction_forlocus(model, F1, l, q):
    L = len(model.F.E1)
    res = SP.zeros(L)
    for l2 in range(L):
        #if l != q and l2 == 230: pdb.set_trace()
        res[l2] = get_estimate(model, l, l2, q)[0]## estimate allele frequency posterior at locus l from data at locus l2, regardless of position of qtl locus q
    PL.plot(range(L), SP.ones(L)*F1[l], "k--")
    PL.plot(range(L), res, "b.")
    PL.show()


def plot_prediction_fromlocus(model, F1, l, q):
    L = len(model.F.E1)
    res = SP.zeros(L)
    for l2 in range(L):
        qtl_l = (q <= l and l <= l2) or (q >= l and l >= l2) # qtl towards l
        if (l <= q and q <= l2) or (l2 <= q and q <= l):  # if QTL between two loci, estimate from q
            res[l2] = get_simple_estimate(model, F1[q], 0, l2, q, True)[0]
        else:
            res[l2] = get_simple_estimate(model, F1[l], 0, l2, l, qtl_l)[0]
    PL.plot(range(L), F1, "k--")
    PL.plot(range(L), res, "b.")
    PL.show()
