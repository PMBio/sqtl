import logging
import cPickle

SQTL_VERSION = "0.1"

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] %(name)-5s: %(levelname)-8s %(message)s')
LOG = logging.getLogger("sQTL")
DATA_DIR = "/Users/leopold/data/projects/sqtl"
EPS = 1e-10

def cl(f): return cPickle.load(open(f,'rb'))
def cdm(o,f): cPickle.dump(o, open(f, 'wb'), -1)
