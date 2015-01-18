import logging
import cPickle
import datetime

SQTL_VERSION = "0.1"

logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] %(name)-5s: %(levelname)-8s %(message)s')
LOG = logging.getLogger("sQTL")
DATA_DIR = "/Users/leopold/data/projects/sqtl"
EPS = 1e-10

def cl(f): return cPickle.load(open(f,'rb'))
def cdm(o,f): cPickle.dump(o, open(f, 'wb'), -1)
def current_time_str(): return datetime.datetime.today().strftime("[%Y-%m-%d (%a) %H:%M]")
