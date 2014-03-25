import os
import sys
from sqtl.test.collate_infer import *
N_JOBS = 40
JOB_START = 0 # 160 + 40 on viimane; j2rgmised peaks tulema 200+
CONFIG_FILE = "/home/morphology/shared/lparts/data/projects/sqtl/sim/config/config.txt"

def main():
    if sys.argv[1] == "run":
        i = 0
        for seed in range(JOB_START, JOB_START + N_JOBS):
            for f in glob.glob("%s"%(CONFIG_FILE.replace(".txt","_*.txt"))):
                outfilename = "%s/sim-3/seed-%d/r-%s_config-%s.pickle"%(DATA_DIR, seed, "fixed", f.split("/")[-1][0:-4])
                if not os.path.exists(outfilename):
                    i += 1
                    if i <= 1000:
                        os.system("submitjob python run_infer.py %d %s"%(seed, f))
#        print i, "not done from earlier"
    if sys.argv[1] == "collate":
        if len(sys.argv) < 3: 
            print "Need rho set as parameter - either fixed or truth"
            return 
        collate(sys.argv[2])
    if sys.argv[1] == "plot":
        if len(sys.argv) < 3: 
            print "Need rho set as parameter - either fixed or truth"
            return
        plot_all_collated(sys.argv[2])
    if sys.argv[1] == "create_configs":
        create_configs()


def create_configs():
    config_base = CONFIG_FILE.replace(".txt","")
    os.system("cp %s %s"%(CONFIG_FILE, CONFIG_FILE.replace(".txt","_meaninfl-0.00.txt")))
    for infl in [0.25, 0.5, 2,4]:
        ofh = file("%s_meaninfl-%.2f.txt"%(config_base, infl),'w')
        ofh.write("250\n100,150,188\n0.005\n%.2f\n0.00002\n1\n0.6,0.5\n"%infl)
        ofh.close()
        ofh = file("%s_varinfl-%.2f.txt"%(config_base, infl),'w')
        ofh.write("250\n100,150,188\n0.005\n1\n0.00002\n%.2f\n0.6,0.5\n"%infl)
        ofh.close()
        ofh = file("%s_rincrease-%.2f.txt"%(config_base, infl),'w')
        ofh.write("250\n100,150,188\n%.4f\n1\n0.00002\n1\n0.6,0.5\n"%(0.005*infl))
        ofh.close()

    
if __name__ == '__main__':
    main()
