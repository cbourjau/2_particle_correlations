"""Divide the total yields in the first file by the total yields in the second one"""

from rootpy.io import root_open
from rootpy import asrootpy

import sys
import logging
from os.path import dirname, abspath


if len(sys.argv) != 3:
    sys.exit("Wrong number of arguments, dude!")

with root_open(sys.argv[1], 'read') as f1, root_open(sys.argv[2], 'read') as f2:
    out_dir = dirname(abspath(sys.argv[1]))
    with root_open(out_dir + '/output_effs.root', 'recreate') as f_eff:
        f_eff.mkdir('processed')
        f_eff.mkdir('processed/eff_from_total_yield')
        f_eff.cd('processed/eff_from_total_yield')
        names = ['total_yield_class_' + str(i) for i in range(0,4)]
        for name in names:
            h = (asrootpy(f1.processed.total_yield.get(name)) /
                 asrootpy(f2.processed.total_yield.get(name)))
            h.SetNameTitle(name[12:], 'Total yield division class ' + name[-1])
            f_eff.Write(h.name)

logging.info('output_effs.root was saved in the directory of the first file')
