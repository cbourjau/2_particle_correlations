"""Start of the analysis"""

# ROOT interfers with argv -> import it later
# see: http://root.cern.ch/phpBB3/viewtopic.php?t=8501
import argparse
import pickle
import time
from os.path import dirname, abspath
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--data", help="List of data files")
parser.add_argument("--workers", help="Number of workers",
                    type=int, default=None)
parser.add_argument("--assocs", help="interval for associated particles;"
                    + "min max",
                    nargs=2, type=float, default=[1.0, 2.0])
parser.add_argument("--triggers", help="interval for trigger particles;"
                    + "min max",
                    nargs=2, type=float, default=[2.0, 4.0])
parser.add_argument("--datatype", help="real, MC or reconstructed",
                    default='real')
parser.add_argument("--cut", help="golden or tpc",
                    default='golden')
parser.add_argument("--single_track",
                    help="Enable single track correction",
                    action='store_true')
parser.add_argument("--threshold", help="Threshold in GeV",
                    type=float, default=0.0)
parser.add_argument("--use_chain",
                    help="Process chain in single track without PROOF",
                    action='store_true')
parser.add_argument("--events_have_one_mc_trigger",
                    help="Valid events have exactly one mc trigger",
                    action='store_true')
parser.add_argument("--events_have_recon_triggers",
                    help="Valid events have at least one recon trigger",
                    action='store_true')
parser.add_argument("--allow_0_trig_in_pool",
                    help="Allow events in the pool which have 0 triggers",
                    action='store_true')

args = parser.parse_args()
if args.data:
    fn = args.data
else:
    fn = '/home/christian/msc/analysis/pyproof/ppb.dat'
print 'Using data file: ', fn

if args.datatype not in ['MC', 'reconstructed',  'real']:
    assert(0)
workers = args.workers
with open('parameters.pkl', 'wb') as pkl:
    pickle.dump(args, pkl)

tmpargv = sys.argv[:]    # [:] for a copy, not reference
sys.argv = []

from ROOT import TProof, TFileCollection, TChain, gROOT
sys.argv = tmpargv

if args.use_chain:
    chain = TChain("tree")
    gROOT.LoadMacro('DebugClassesMultESA2013.C+')
    with open(fn, 'read') as f:
        for l in f.readlines():
            chain.AddFile(l[:-1])
    chain.Process('TPySelector', 'selector')

else:  # Use Proof
    if workers:
        proof = TProof.Open("workers=" + str(workers))
    else:
        proof = TProof.Open('')

    proof.Load('./DebugClassesMultESA2013.C+,./parameters.pkl')
    # add the current folder to the python path of root (or something)
    proof.Exec('TPython::Exec("%s");' %
            ("import sys; sys.path.insert(0,'"+dirname(abspath("selector.py"))+"')"))

    time.sleep(1)  # needed for GUI to settle

    # giving some name to collection to shut up warning
    fc = TFileCollection('analysis_files')
    fc.AddFromFile(fn)
    proof.Process(fc, 'TPySelector', 'selector')

    # make the single track efficiency:
    #proof.Process(fc, 'TPySelector', 'single_eff_selector')
