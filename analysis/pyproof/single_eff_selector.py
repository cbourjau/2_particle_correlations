"""Create one TH3F histo with eta, pt and zvtx"""

### selector module (selector.py, name has to match as per in main.py)
from selector import try_except
import selector
from ROOT import gROOT, THnSparseF, TH1I

from array import array
import numpy as np
from math import pi

from rootpy.io import root_open

# load MyClass definition macro (append '+' to use ACLiC)
gROOT.LoadMacro('DebugClassesMultESA2013.C+')
#from ROOT import DeDxEvent


class Efficiency_selector(selector.MyPySelector):

    def __init__(self):
        self.zvtx_bins = (10, -10, 10)
        self.trigger = 1  # 1 = V0AND !!! Double check this!!!
        self.max_dist_vtx_xy = 2.4
        self.max_dist_vtx_z = 3.2
        self.pt_threshold = 0.0
        self.mc_and_recon_valid = False  # Do manually to avoid counting triggers
        self.events_have_recon_triggers = False
        self.events_have_one_mc_trigger = False
        self.exclude_tpc_border = False
        self.only_charge = False
        self.track_filter = 1  # 1:golden; 2:TPC-only

    def Begin(self):
        pass

    def SlaveTerminate(self):
        pass

    def Init(self, tree):
        pass

    def SlaveBegin(self, tree):
        # register counters wrt z index
        print self.__class__.__module__+": SlaveBegin"

        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins
        pt_edges = array('d', (list(np.arange(0.5, 1, 0.01))
                                    + list(np.arange(1, 4.1, 0.1))))
        self.tracks_counter_mc = THnSparseF('tracks_mc',
                    'Valid_tracks;#phi;#eta;z_vtx;pt',
                    4,
                    array('i', array('i', [36, 32, 10, len(pt_edges)-1])),
                            array('d', array('d', [0, -.8, -10, .5])),
                            array('d', array('d', [2*pi, .8, 10, 4])))
        self.tracks_counter_re = THnSparseF('tracks_re',
                    'Valid_tracks;#phi;#eta;z_vtx;pt',
                    4,
                    array('i', array('i', [36, 32, 10, len(pt_edges)-1])),
                            array('d', array('d', [0, -.8, -10, .5])),
                            array('d', array('d', [2*pi, .8, 10, 4])))
        self.tracks_counter_mc.SetBinEdges(3, pt_edges)
        self.tracks_counter_mc.Sumw2()
        self.tracks_counter_re.SetBinEdges(3, pt_edges)
        self.tracks_counter_re.Sumw2()

        self.GetOutputList().Add(self.tracks_counter_mc)
        self.GetOutputList().Add(self.tracks_counter_re)

        # Helper histograms to find bins
        self._bin_zvtx = TH1I('bin_zvtx', 'used to find z vtx',
                              zvtx_nbins, zvtx_low, zvtx_high)

    @try_except
    def Process(self, entry):
        self.fChain.GetEntry(entry)
        # Validate on the event branch level
        if self.data_type == 'MC':
            if not self.validate_event_mc():
                return 1
            zvtx = self.fChain.event.zvtxMC
        if self.data_type in ['real', 'reconstructed']:
            if not self.validate_event_rec():
                return 1
            zvtx = self.fChain.event.zvtx
        # validate wrt number of triggers        
        if not self.validate_event_wrt_n_triggers():
            return 1
        
        tracks_mc = self.get_tracks_generator('MC')
        tracks = self.get_tracks_generator('reconstructed')

        fill_mc = self.tracks_counter_mc.Fill
        fill_re = self.tracks_counter_re.Fill

        zvtx_mc = self.fChain.event.zvtxMC
        zvtx = self.fChain.event.zvtx

        for track in tracks_mc:
            if self._validate_mc_track(track):
                fill_mc(array('d', [track.phiMC, track.etaMC, zvtx_mc, track.ptMC]))

        for track in tracks:
            if self._validate_recon_track(track):
                fill_re(array('d', [track.phi, track.eta, zvtx, track.pt]))

        return 0


    def Terminate(self):
        with root_open('single_track_eff/eff.root', 'recreate') as f:
            # write original histograms:
            f.mkdir('raw')
            f.raw.cd()
            for l in self.GetOutputList():
                l.Write()
        print 'fOutput in Terminate', self.GetOutputList().ls()


if __name__ == '__main__':
    from ROOT import TProof, TFileCollection, TChain
    from os.path import dirname, abspath
    import os
    import time
    from sys import argv

    # giving some name to collection to shut up warning
    fc = TFileCollection('analysis_files')
    try:
        fn = argv[1]
    except IndexError:
        print "no MC file given"
        assert(0)
    fc.AddFromFile(fn)

    proof = TProof.Open('')

    proof.Load('./DebugClassesMultESA2013.C+,./parameters.pkl')
    # add the current folder to the python path of root (or something)
    proof.Exec('TPython::Exec("%s");' %
            ("import sys; sys.path.insert(0,'"+dirname(abspath("selector.py"))+"')"))

    time.sleep(1)  # needed for GUI to settle
    
    directory = 'single_track_eff'
    if not os.path.exists(directory):
        os.makedirs(directory)

    # check the phi vs. pt acceptance with threshold
    proof.Process(fc, 'TPySelector', 'single_eff_selector')
    # chain = TChain("tree")
    # gROOT.LoadMacro('DebugClassesMultESA2013.C+')
    # with open(fn, 'read') as f:
    #     for l in f.readlines():
    #         chain.AddFile(l[:-1])
    # chain.Process('TPySelector', 'single_eff_selector')
