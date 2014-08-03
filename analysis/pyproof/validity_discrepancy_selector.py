"""Create one TH3F histo with eta, pt and zvtx"""

### selector module (selector.py, name has to match as per in main.py)
from selector import try_except
import selector
from ROOT import gROOT, THnSparseF, TH1I

import numpy as np

from rootpy.io import root_open
from rootpy.plotting import Hist, Profile

# load MyClass definition macro (append '+' to use ACLiC)
gROOT.LoadMacro('DebugClassesMultESA2013.C+')
#from ROOT import DeDxEvent


class Validation_discrep_selector(selector.MyPySelector):

    def __init__(self):
        self.zvtx_bins = (10, -10, 10)
        self.trigger = 1  # 1 = V0AND !!! Double check this!!!
        self.max_dist_vtx_xy = 2.4
        self.max_dist_vtx_z = 3.2
        self.pt_threshold = -3.0
        self.assoc_inter = [0.5, 1.0]
        self.trigger_inter =[1.0,2.0]
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
        """overwrite original one to have all branches"""
        pass

    def SlaveBegin(self, tree):
        print self.__class__.__module__+": SlaveBegin"
        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins
        # Count events of the following categories: valid_mc_valid_rec,
        # valid_mc_invalid_rec, invalid_mc_valid_rec, both_invalid
        self.valid_counter = Hist(4, 0, 4, name='MC_rec_discrepancy',
                                  title='Discrepancy between event validation of MCrec and MCtruth at thresh {}'.format(str(self.pt_threshold)))
        self.valid_counter.Sumw2()
        self.GetOutputList().Add(self.valid_counter)

        self.mean_pts = []
        for i in range(0,4):
            self.mean_pts.append(Profile(10, 0, 100,
                                         name="mean_pt_case_{}".format(str(i)),
                                         title="mean pt for case {}".format(str(i))))
            self.mean_pts[-1].Sumw2()
            self.GetOutputList().Add(self.mean_pts[-1])


        # Helper histograms to find bins
        self._bin_zvtx = TH1I('bin_zvtx', 'used to find z vtx',
                              zvtx_nbins, zvtx_low, zvtx_high)

    @try_except
    def Process(self, entry):
        self.data_type = 'dont_use_me_in_this_analysis'
        self.fChain.GetEntry(entry)

        valid_mc = (self.validate_event_mc() and self.get_n_triggers_for_type('MC') > 1)
        valid_rec =(self.validate_event_rec()
                    and self.get_n_triggers_for_type('reconstructed') > 1)

        # calculate mean_pt from the generator tracks
        tracks = self.get_valid_tracks('MC')
        if len(tracks) > 0:
            mean_pt = np.mean([track.ptMC for track in tracks])
        else:
            return 1
        if valid_mc and valid_rec:
            self.valid_counter.fill(.5)  # first bin
            self.mean_pts[0].Fill(self.fChain.event.trackmultMC, mean_pt)
        elif valid_mc and not valid_rec:
            self.valid_counter.fill(1.5)  # second bin
            self.mean_pts[1].Fill(self.fChain.event.trackmultMC, mean_pt)
        elif not valid_mc and valid_rec:
            self.valid_counter.fill(2.5)  # third bin
            self.mean_pts[2].Fill(self.fChain.event.trackmultMC, mean_pt)
        else:
            self.valid_counter.fill(3.5)  # forth, both invalid
            self.mean_pts[3].Fill(self.fChain.event.trackmultMC, mean_pt)
        return 0


    def Terminate(self):
        with root_open('event_validation_discrepancy{}.root'.format(str(self.pt_threshold)), 'recreate') as f:
            # write original histograms:
            f.mkdir('raw')
            f.raw.cd()
            for l in self.GetOutputList():
                l.Write()
        print 'fOutput in Terminate', self.GetOutputList().ls()

    def get_valid_tracks(self, track_type):
        """return list of valid tracks"""
        valid_tracks = []
        if track_type == 'MC':
            tracks = self.fChain.trackmc
            n = self.fChain.trackmc.GetEntriesFast()
            for i in xrange(0, n):
                t = tracks.At(i)
                if self._validate_mc_track(t):
                    valid_tracks.append(t)
        elif track_type == 'reconstructed':
            tracks = self.fChain.track
            n = self.fChain.track.GetEntriesFast()
            for i in xrange(0, n):
                t = tracks.At(i)
                if self._validate_recon_track(t):
                    valid_tracks.append(t)
        elif track_type == 'real':
            tracks = self.fChain.track
            n = self.fChain.track.GetEntriesFast()
            for i in xrange(0, n):
                t = tracks.At(i)
                if self._validat_real_track(t):
                    valid_tracks.append(t)
        else:
            assert(0)
        return valid_tracks


if __name__ == '__main__':
    from ROOT import TProof, TFileCollection, TChain
    from os.path import dirname, abspath
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

    proof = TProof.Open("workers=2")

    proof.Load('./DebugClassesMultESA2013.C+,./parameters.pkl')
    # add the current folder to the python path of root (or something)
    proof.Exec('TPython::Exec("%s");' %
            ("import sys; sys.path.insert(0,'"+dirname(abspath("selector.py"))+"')"))

    time.sleep(1)  # needed for GUI to settle
    
    # check the phi vs. pt acceptance with threshold
    proof.Process(fc, 'TPySelector', 'validity_discrepancy_selector')
    chain = TChain("tree")
    gROOT.LoadMacro('DebugClassesMultESA2013.C+')
    with open(fn, 'read') as f:
        for l in f.readlines():
            chain.AddFile(l[:-1])
    # chain.Process('TPySelector', 'validity_discrepancy_selector')
