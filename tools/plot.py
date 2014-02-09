import pylab as PL


def init_plot_params():
    plparams = {'backend': 'ps',
          'axes.labelsize': 26,
          'font.size': 40,
          'legend.fontsize': 34,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'text.fontsize':40,
          'text.size':40,
          'text.usetex': True
    }
    PL.rcParams.update(plparams)



# def plot_qtlmapping():
#         print "%d\t%.1f\t%.1f\t%.1f\t%.1f"%(q, SP.log(tt).sum(), -SP.log(tau), -(tt*(m**2 - mu**2)).sum(),k)
#         if q in cands:
#                 #if abs(q-50) < 4: pdb.set_trace()
#                 PL.plot(tt*(mu**2 - m**2))
#                 taus.append(tt)
#                 #mudiffs.append(m**2-mu**2)
#                 mudiffs.append(m)
#         PL.plot(net.qtl_loci, self.lnX[net.qtl_loci] - self.lnX[net.qtl_loci].max(), "r--")
#         PL.subplot(3,1,2)
#         for i,tt in enumerate(mudiffs):
#             PL.plot(range(len(net.D)), tt)
#             PL.plot([cands[i]],[net.F.E1[cands[i],cands[i]]], 'k.', markersize=10)
#         PL.plot(range(len(net.D)), net.D[:,0]/(net.D.sum(axis=1)), 'b.')
#         PL.ylim(0,2)
#         PL.subplot(3,1,3)
#         for tt in taus: PL.plot(range(len(net.D)), tt)
