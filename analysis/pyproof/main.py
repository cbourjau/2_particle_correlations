"""Start of the analysis"""

# ROOT interfers with argv -> import it later
# see: http://root.cern.ch/phpBB3/viewtopic.php?t=8501
import argparse
import pickle
import time
from os.path import dirname, abspath
import os
import sys
import json
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--data", help="List of data files")
parser.add_argument("--workers", help="Number of workers",
                    type=int, default=None)
parser.add_argument("--assocs", help="interval for associated particles;"
                    + "min max", dest='assoc_inter',
                    nargs=2, type=float, default=[0.5, 1.0])
parser.add_argument("--triggers", help="interval for trigger particles;"
                    + "min max", dest='trigger_inter',
                    nargs=2, type=float, default=[1.0, 2.0])
parser.add_argument("--datatype", help="real, MC or reconstructed",
                    default='real', dest='data_type')
parser.add_argument("--cut", help="golden or tpc", dest='track_cut',
                    default='golden')
# parser.add_argument("--single_track",
#                     help="Enable single track correction",
#                     action='store_true')
parser.add_argument("--threshold", help="Threshold in GeV. Negative for soft case",
                    type=float, default=0.0, dest='pt_threshold')
parser.add_argument("--use_chain",
                    help="Process chain in single track without PROOF",
                    action='store_true')
parser.add_argument("--events_have_one_mc_trigger",
                    help="Valid events have exactly one mc trigger",
                    action='store_true')
parser.add_argument("--events_have_recon_triggers",
                    help="Valid events have at least one recon trigger",
                    action='store_true')
parser.add_argument("--allow_0_trig_in_pool", dest='allow_0_trig',
                    help="Allow events in the pool which have 0 triggers",
                    action='store_true')
parser.add_argument("--exclude_tpc_border",
                    help="Exclude MC tracks too close to the tpc border",
                    action='store_true')
parser.add_argument("--only_pos", help="consider only positive charged tracks",
                    action='store_const', const=1, dest='only_charge')
parser.add_argument("--only_neg", help="consider only positive charged tracks",
                    action='store_const', const=-1, dest='only_charge')
parser.add_argument("--append_to_dir", help="Append extra info to dir name",
                    type=str, default='')
parser.add_argument("--mc_and_recon_valid", help="Both, gen and recon have to be valid",
                    action='store_true')

args = parser.parse_args()
if args.data:
    fn = args.data
else:
    fn = '/home/christian/msc/analysis/pyproof/ppb.dat'
print 'Using data file: ', fn

# Create a new director and copy the needed files there
directory = 'assocs{}-{}_trigs{}-{}_thresh{}_cut_{}_type_{}{}'.format(
    args.assoc_inter[0], args.assoc_inter[1],
    args.trigger_inter[0], args.trigger_inter[1], args.pt_threshold,
    args.track_cut, args.data_type, args.append_to_dir)
if not os.path.exists(directory):
    os.makedirs(directory)
else:
    print 'directory exists already!'
    assert 0 
args.directory = directory

if args.data_type not in ['MC', 'reconstructed',  'real']:
    assert(0)

workers = args.workers
with open('parameters.pkl', 'wb') as pkl:
    pickle.dump(args, pkl)

args.gitversion = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'])[:-1]

# save data to file for humans
with open(directory+"/info.txt", "w") as text_file:
    text_file.write(json.dumps(vars(args)))

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
            ("import sys; sys.path.insert(0,'"
              +dirname(abspath("selector.py"))+"')"))

    time.sleep(1)  # needed for GUI to settle

    # giving some name to collection to shut up warning
    fc = TFileCollection('analysis_files')
    fc.AddFromFile(fn)
    proof.Process(fc, 'TPySelector', 'selector')

    # make the single track efficiency:
    #proof.Process(fc, 'TPySelector', 'single_eff_selector')
