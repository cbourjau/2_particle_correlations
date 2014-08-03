"""Given one file the background, signal and total yield is computed"""

from bg_sig_total_yield import calc_bg, calc_signal_and_total_yield

import sys
import logging
from rootpy.io import root_open

if len(sys.argv) != 2:
    sys.exit('Wrong number of arguments, mate!')

fn = sys.argv[1]


with root_open(fn, 'update') as f:
    try:
        for dir_name in f.processed.walk(maxdepth=0).next()[1]:
            logging.info("""Delete """ + 'processed/' + dir_name +"""  directory""")
            f.rm('processed/' + dir_name)
    except:
        pass

logging.info('Calculating background')
calc_bg(fn)

logging.info('Calculating signal and total yield')
calc_signal_and_total_yield(fn)

