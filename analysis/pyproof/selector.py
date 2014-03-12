### selector module (selector.py, name has to match as per in main.py)
from ROOT import TPySelector, gROOT, TH1I, THnSparseF

from array import array
import numpy as np
from math import pi
import pickle
from rootpy.io import root_open


# load MyClass definition macro (append '+' to use ACLiC)
gROOT.LoadMacro('DebugClassesMultESA2013.C+')
#from ROOT import DeDxEvent


def try_except(fn):
    """decorator for extra debugging output"""
    def wrapped(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except:
            import traceback
            traceback.print_exc()
            assert(0)
    return wrapped


class MyPySelector(TPySelector):

    def __init__(self):
        print self.__class__.__module__+": init"
        with open('parameters.pkl', 'read') as pkl:
            args = pickle.load(pkl)
        print 'Settings: ', args
        self.assoc_inter = args.assocs
        self.trigger_inter = args.triggers
        self.data_type = args.datatype
        self.single_track_correction = args.single_track
        self.track_cut = args.cut
        self.allow_0_trig = args.allow_0_trig_in_pool
        self.events_have_one_mc_trigger = args.events_have_one_mc_trigger
        self.events_have_reconstructed_triggers = args.events_have_recon_triggers
        if self.track_cut == 'golden':
            self.fn_single_tracks = (
                '~/msc/results/efficiencies/single_track_eff_golden/single_eff.root')
            self.track_filter = 1  # 1:golden; 2:TPC-only
        if self.track_cut == 'tpc':
            self.fn_single_tracks = (
                '~/msc/results/efficiencies/single_track_eff_TPC/single_eff.root')
            self.track_filter = 2  # 1:golden; 2:TPC-only
        self.pool_size = 1000
        self.pt_threshold = args.threshold
        self.zvtx_bins = (10, -10, 10)
        self.cent_bins = (4, 0, 100.01)  # edges are just dummies!
        self.cent_edges = array('d', [0, 20, 40, 60, 100.01])
        self.ptmax_bins = (20, 0, 10)
        self.phi_bins = (36, (-pi/2), (3*pi/2))
        self.eta_bins = (32, -1.6, 1.6)
        self.trigger = 1  # 1 = V0AND
        self.max_dist_vtx_xy = 2.4
        self.max_dist_vtx_z = 3.2
        # for counter histograms only:
        #number of tracks in the given interval (# of assocs or triggers)
        self.n_tracks_nbins = 10
        # last bin is really large
        self.n_tracks_edges = array('d', list(range(0, 10)) + [100])

    def Begin(self):
        print 'Data type in Begin(): ', self.data_type

    def SlaveBegin(self, tree):
        print 'py: slave beginning'
        print 'data type in slave: ', self.data_type
        # register counters wrt z index
        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins
        cent_nbins, cent_low, cent_high = self.cent_bins  # edges are ignored
        ptmax_nbins, ptmax_low, ptmax_high = self.ptmax_bins

        self.valid_counter = TH1I('valid', 'Valid events',
                                  zvtx_nbins, zvtx_low, zvtx_high)
        self.trig_counter = self.make_counter_histogram('trig',
                                        'Total number of triggers')

        # self.assoc_counter = self.make_counter_histogram('assocs',
        #                             'Total number of assocs')

        # Dimensions: phi, eta, z_vtx, cent, ptmax
        self.signal = self.make_diff_histogram('signal', 'Signal distribution')
        self.background = self.make_diff_histogram('background',
                                                   'Background distribution')
        self.GetOutputList().Add(self.signal)
        self.GetOutputList().Add(self.background)
        self.GetOutputList().Add(self.valid_counter)
        self.GetOutputList().Add(self.trig_counter)
        # self.GetOutputList().Add(self.assoc_counter)

        if self.single_track_correction:
            # load single weight histogram
            self.single_track_weight = (
                root_open(self.fn_single_tracks, 'read').Get('single_track_eff'))

        # The pool for each z section
        self.pool = [[]]
        for i in range(0, 10):
            for c in range(0, 4):
                self.pool[-1].append([])
            self.pool.append([])
        # Helper histograms to find bins
        self._bin_zvtx = TH1I('bin_zvtx', 'used to find z vtx',
                              zvtx_nbins, zvtx_low, zvtx_high)
        self._bin_cent = TH1I('bin_cent', 'used to find eclasses',
                              cent_nbins, self.cent_edges)

        # Set the validation function
        # Done like this for performance
        if self.data_type == 'MC':
            self.validate_track = self._validate_mc_track
        elif self.data_type == 'reconstructed':
            self.validate_track = self._validate_recon_track
        elif self.data_type == 'real':
            self.validate_track = self._validate_real_track

    @try_except
    def Process(self, entry):
        self.fChain.GetEntry(entry)
        if not self.validate_event(entry):
            return 1
        if self.data_type == 'MC':
            zvtx = self.fChain.event.zvtxMC
            ptmax = self.fChain.event.ptmaxMC
        if self.data_type in ['real', 'reconstructed']:
            zvtx = self.fChain.event.zvtx
            ptmax = self.fChain.event.ptmax
        # there is no centMC!
        cent = self.fChain.event.cent

        # stop here if there are no triggers for sure
        if not self.allow_0_trig and ptmax < self.trigger_inter[0]:
            return 1

        trigs, assocs = self.get_triggers_assocs()
        if not self.allow_0_trig and len(trigs) == 0:
            return 1

        self.valid_counter.Fill(zvtx, 1)

        if self.added_to_pool(assocs):
            return 1

        trig_fill = self.trig_counter.Fill
        # assoc_fill = self.assoc_counter.Fill

        # fill trigger and assoc counters
        nr_trigs = len(trigs)
        for phi, eta, pt, w in trigs:
            a = array('d', [phi, eta, zvtx, cent, ptmax, nr_trigs])
            trig_fill(a)
        # for phi, eta, pt, w in assocs:
        #     a = array('d', [phi, eta, zvtx, cent, ptmax, nr_trigs])
        #     assoc_fill(a, 1 / w)  # w != 0

        s_fill = self.signal.Fill
        bg_fill = self.background.Fill
        #cdef double trig_phi, trig_eta, assoc_phi, assoc_eta, bg_phi, bg_eta
        #cdef int i
        bgs_it = self.it_background_tracks()
        for trig_phi, trig_eta, trig_pt, trig_w in trigs:
            for assoc_phi, assoc_eta, assoc_pt, assoc_w in assocs:
                if assoc_pt >= trig_pt:
                    continue
                a = array('d', [wrap(assoc_phi-trig_phi),
                                assoc_eta-trig_eta,
                                zvtx, cent, ptmax])
                s_fill(a, 1.0/assoc_w)
            # 10 bg tracks per assoc
            for i in xrange(len(assocs)*10):
                bg_phi, bg_eta, bg_pt, bg_w = bgs_it.next()
                if bg_pt >= trig_pt:
                    continue
                a = array('d', [wrap(bg_phi - trig_phi),
                                bg_eta - trig_eta,
                                zvtx, cent, ptmax])
                bg_fill(a, 1.0/bg_w)

        # replace tracks in pool
        self.replace_pool(assocs)

        return 0

    def SlaveTerminate(self):
        print 'py: slave terminating'

    def Terminate(self):
        """currently ignoring the pt threshold part"""
        with root_open('output_' + self.data_type + '.root',
                       'recreate') as f:
            # write original histograms:
            f.mkdir('raw')
            f.raw.cd()
            for l in self.GetOutputList():
                l.Write()
        print 'fOutput in Terminate', self.GetOutputList().ls()
        print 'Successfully wrote results to output_' + self.data_type + '.root'

    ### Analysis functions
    def make_diff_histogram(self, name, title,
                       axis_title=';#Delta #varphi;#Delta #eta;z_vtx;cent;pt_max'):
        """create a n dim histogram"""
        phi_nbins, phi_low, phi_high = self.phi_bins
        eta_nbins, eta_low, eta_high = self.eta_bins
        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins
        cent_nbins, cent_low, cent_high = self.cent_bins  # edges are ignored
        ptmax_nbins, ptmax_low, ptmax_high = self.ptmax_bins

        tmp = THnSparseF(name,
                         title + axis_title,
                         5,  # dimensions
                         array('i', [phi_nbins, eta_nbins, zvtx_nbins,
                                     cent_nbins, ptmax_nbins]),
                         array('d', [phi_low, eta_low,  zvtx_low,
                                     cent_low, ptmax_low]),
                         array('d', [phi_high, eta_high, zvtx_high,
                                     cent_high, ptmax_high]))
        # Set centrality bin edegs
        tmp.SetBinEdges(3, self.cent_edges)
        tmp.Sumw2()
        return tmp

    def make_counter_histogram(self, name, title,
                               axis_title=(';#Delta #varphi;#Delta #eta;'
                                           'z_vtx;cent;pt_max;n_trig;p_t')):
        """create a n dim histogram"""
        phi_nbins, phi_low, phi_high = 36, 0, 2*pi
        eta_nbins, eta_low, eta_high = 32, -0.8, 0.8
        # pt bins are from analysis note
        # pt_nbins = 22
        # pt_edges = array('d', (list(np.arange(0.5, 1, 0.1))
        #                        + list(np.arange(1, 4, 0.25))
        #                        + list(np.arange(4, 5, 0.5))
        #                        + list(range(5, 9, 1))))
        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins
        cent_nbins, fake_low, fake_high = self.cent_bins  # edges are ignored
        ptmax_nbins, ptmax_low, ptmax_high = self.ptmax_bins
        tmp = THnSparseF(name,
                         title + axis_title,
                         6,  # dimensions
                         array('i', [phi_nbins, eta_nbins,  zvtx_nbins,
                                     cent_nbins, ptmax_nbins, self.n_tracks_nbins]),
                         array('d', [phi_low, eta_low, zvtx_low,
                                     fake_low, ptmax_low, fake_low]),
                         array('d', [phi_high, eta_high, zvtx_high,
                                     fake_high, ptmax_high, fake_high]))
        # Set centrality bin edegs
        tmp.SetBinEdges(3, self.cent_edges)
        tmp.SetBinEdges(5, self.n_tracks_edges)
        # tmp.SetBinEdges(6, pt_edges)
        tmp.Sumw2()
        return tmp

    def get_triggers_assocs(self):
        """Ret ([[trig_phi, trig_eta, pt, w]], [[assoc_phi, assoc_eta, pt, w]])
        for valid tracks
        """
        lower_trig, upper_trig = self.trigger_inter
        lower_assoc, upper_assoc = self.assoc_inter
        triggers, assocs = [], []
        single_track_correction = self.single_track_correction
        if single_track_correction:
            weight_hist = self.single_track_weight
        if self.data_type == 'MC':
            tracks = self.get_tracks_generator('MC')
            for track in tracks:
                if not self.validate_track(track):
                    continue
                pt = track.ptMC
                if (lower_assoc < pt < upper_assoc):
                    zvtx = self.fChain.event.zvtxMC
                    if single_track_correction:
                        w = weight_hist.GetBinContent(
                            weight_hist.GetBin(
                                array('d',
                                      [track.phiMC, track.etaMC, zvtx, pt])))
                        if w < 0.01:  # dont know if int or float
                            continue
                    else:
                        w = 1
                    assocs.append([track.phiMC, track.etaMC, pt, w])
                if (lower_trig < pt < upper_trig):
                    w = 1
                    triggers.append([track.phiMC, track.etaMC, pt, w])
        elif self.data_type in ['real', 'reconstructed']:
            tracks = self.get_tracks_generator(self.data_type)
            for track in tracks:
                if not self.validate_track(track):
                    continue
                pt = track.pt
                if (lower_assoc < pt < upper_assoc):
                    zvtx = self.fChain.event.zvtx
                    if single_track_correction:
                        w = weight_hist.GetBinContent(
                            weight_hist.GetBin(
                                array('d', [track.phi, track.eta, zvtx, pt])))
                        if w < 0.01:  # dont know if int or float
                            continue
                    else:
                        w = 1
                    assocs.append([track.phi, track.eta, pt, w])
                if (lower_trig < pt < upper_trig):
                    w = 1
                    triggers.append([track.phi, track.eta, pt, w])
        return (triggers, assocs)

    def it_background_tracks(self):
        """Iterator for background tracks"""
        zvtx = self.get_zvtx_id()
        cent = self.get_cent_id()
        pool = self.pool[zvtx][cent]
        while True:
            np.random.shuffle(pool)
            for track in pool:
                yield track

    def validate_event(self, entry):
        """Validate event on event level"""
        if not (self.fChain.event.trig & self.trigger):
            return False
        if self.fChain.event.vtxstatus < 1:
            return False
        if self.data_type == 'MC':
            if not (-10 < self.fChain.event.zvtxMC < 10):
                return False
            if self.fChain.event.ptmaxMC < self.pt_threshold:
                return False
            if (self.events_have_reconstructed_triggers and
                    not self._has_recon_trigger(entry)):
                return False
            if (self.events_have_one_mc_trigger
                    and not self._has_exactly_one_mc_trigger(entry)):
                return False

        elif self.data_type in ['real', 'reconstructed']:
            if not (-10 < self.fChain.event.zvtx < 10):
                return False
            if self.fChain.event.ptmax < self.pt_threshold:
                return False
            if (self.data_type == 'reconstructed'
                    and self.events_have_one_mc_trigger
                    and not self._has_exactly_one_mc_trigger(entry)):
                return False
        return True

    def _has_recon_trigger(self, entry):
        """check if at least one reconstructed track is presenten"""
        tracks = self.get_tracks_generator('reconstructed')
        lower_trig, upper_trig = self.trigger_inter
        for track in tracks:
            if not self._validate_recon_track(track):
                continue
            pt = track.pt
            if (lower_trig < pt < upper_trig):
                return True
        else:
            return False

    def _has_exactly_one_mc_trigger(self, entry):
        """check if there is exactly on mc trigger"""
        tracks = self.get_tracks_generator('MC')
        lower_trig, upper_trig = self.trigger_inter
        triggers = 0
        for track in tracks:
            if not self._validate_mc_track(track):
                continue
            pt = track.ptMC
            if (lower_trig < pt < upper_trig):
                triggers += 1
                if triggers > 1:
                    return False
        if triggers == 1:
            return True
        else:
            return False

    def added_to_pool(self, assocs):
        """See if the fitting z-vtx pool is already filled"""
        zvtx = self.get_zvtx_id()
        cent = self.get_cent_id()
        if len(self.pool[zvtx][cent]) > self.pool_size:
            return False
        else:
            # [self.bg_debug.Fill(t[0], t[1]) for t in assocs]
            self.pool[zvtx][cent] += assocs
            return True

    def replace_pool(self, assocs):
        """replace tracks in pool with given assocs"""
        if len(assocs) > 0:
            zvtx = self.get_zvtx_id()
            cent = self.get_cent_id()
            self.pool[zvtx][cent][:len(assocs)] = assocs

    def get_zvtx_id(self):
        """Return z vertex as index [0, 9] but requires events with
        limited z section"""
        if self.data_type == 'MC':
            return self._bin_zvtx.FindBin(self.fChain.event.zvtxMC) - 1
        elif self.data_type in ['real', 'reconstructed']:
            return self._bin_zvtx.FindBin(self.fChain.event.zvtx) - 1

    def get_cent_id(self):
        "Return the event class (bin - 1) for the current event"""
        return self._bin_cent.FindBin(self.fChain.event.cent) - 1

    def validate_track(self):
        """over write this function in SlaveBegin()"""
        pass

    def _validate_mc_track(self, track):
        """return True if the track is valid"""
        if (track.qMC == 0 or (abs(track.etaMC) > 0.8)):
            return False
        return True

    def _validate_recon_track(self, track):
        if not (track.filter & self.track_filter):
            return False
        if not (track.primary & 1):
            return False
        if abs(track.dcaxy) > self.max_dist_vtx_xy:
            return False
        if abs(track.dcaz) > self.max_dist_vtx_z:
            return False
        return True

    def _validate_real_track(self, track):
        if not (track.filter & self.track_filter):
            return False
        if abs(track.dcaxy) > self.max_dist_vtx_xy:
            return False
        if abs(track.dcaz) > self.max_dist_vtx_z:
            return False
        return True

    def get_tracks_generator(self, track_type):
        if track_type == 'MC':
            tracks = self.fChain.trackmc
            n = self.fChain.trackmc.GetEntriesFast()
        elif track_type in ['real', 'reconstructed']:
            tracks = self.fChain.track
            n = self.fChain.track.GetEntriesFast()
        else:
            assert(0)
        return (tracks.At(i) for i in xrange(0, n))

    def Init(self, tree):
        # event branch
        self.fChain.SetBranchStatus('*', 0)
        self.fChain.SetBranchStatus('ptmax', 1)
        self.fChain.SetBranchStatus('cent', 1)
        self.fChain.SetBranchStatus('zvtx', 1)
        self.fChain.SetBranchStatus('trig', 1)
        self.fChain.SetBranchStatus('vtxstatus', 1)
        # track branch
        self.fChain.SetBranchStatus('track.phi', 1)
        self.fChain.SetBranchStatus('track.eta', 1)
        self.fChain.SetBranchStatus('track.pt', 1)
        self.fChain.SetBranchStatus('track.filter', 1)
        self.fChain.SetBranchStatus('track.primary', 1)
        self.fChain.SetBranchStatus('track.dcaxy', 1)
        self.fChain.SetBranchStatus('track.dcaz', 1)
        if self.data_type in ['MC', 'reconstructed']:
            self.fChain.SetBranchStatus('zvtxMC', 1)
            self.fChain.SetBranchStatus('ptmaxMC', 1)
            self.fChain.SetBranchStatus('trackmc.phiMC', 1)
            self.fChain.SetBranchStatus('trackmc.etaMC', 1)
            self.fChain.SetBranchStatus('trackmc.ptMC', 1)
            self.fChain.SetBranchStatus('trackmc.qMC', 1)


def wrap(x):
    """wrap [rad] angle to [-PI/2..3PI/2)"""
    pi = 3.141592653589793
    return (x + pi/2) % (2 * pi) - pi/2


#MyPySelector.fill_histograms = process.fill_histograms
#MyPySelector.get_triggers_assocs = process.get_triggers_assocs
#MyPySelector.it_background_tracks = process.it_background_tracks
