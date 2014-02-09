import pylab as PL
import scipy as SP

def init_plot_params():
    plparams = {'backend': 'ps',
          'axes.labelsize': 18,
          'font.size': 26,
          'legend.fontsize': 24,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.fontsize':28,
          'text.size':28,
          'text.usetex': False
    }
    PL.rcParams.update(plparams)


def plot_toy():
    PL.figure()
    xs = SP.array([1,4,7,8,9,10,11])
    ys_0 = SP.array([0.58,0.55,0.54,0.51,0.5,0.54,0.57])
    ys_1 = SP.array([0.61,0.71,0.84,0.92,0.83,0.77,0.73])
    I_marker = [1,2,3]
#    PL.plot([0,11],[1,1],'k--', linewidth=1)
    l1 = PL.plot(xs,ys_1, 'b-', linewidth=6)
    l2 = PL.plot(xs,ys_0, 'r-', linewidth=6)
    PL.plot(xs[I_marker], ys_1[I_marker], 'bo', markersize=14)
    PL.plot(xs[I_marker], ys_0[I_marker], 'ro', markersize=14)
    PL.ylabel('Parent "A" allele frequency')
    PL.xlabel("Locus")
    PL.xticks(xs[I_marker], ['$l$','$l+1$','$X$'])
    PL.yticks([0.5,0.75, 1], [0.5, 0.75,1])
    PL.text(3.6, 0.73, '$f_{1,l}$')
    PL.text(3.6, 0.45, '$f_{0,l}$')
    PL.text(5.7, 0.85, '$f_{1,l+1}$')
    PL.text(6.25, 0.45, '$f_{0,l+1}$')
    PL.legend([l1,l2], ["$t = 1$", "$t = 0$"], loc=2)
    PL.xlim(1,11)
    PL.ylim(0.4, 1)
    PL.savefig('toy_af.pdf')
    PL.show()



def main():
    init_plot_params()
    plot_toy()


if __name__ == '__main__': 
    main()
