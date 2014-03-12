"""Create one TH3F histo with eta, pt and zvtx"""

### selector module (selector.py, name has to match as per in main.py)
from ROOT import TPySelector, gROOT, THnSparseF, TH1I

from array import array
import numpy as np
from math import pi
import pickle


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
        self.data_type = args.datatype

        # self.phi_bins = [36, (-pi/2), (3*pi/2)]
        # self.eta_bins = [32, -.8, .8]
        self.zvtx_bins = [10, -10, 10]
        # self.pt_nbins = 22
        self.pt_edges = array('d', (list(np.arange(0.5, 1, 0.1))
                                    + list(np.arange(1, 4, 0.25))
                                    + list(np.arange(4, 5, 0.5))
                                    + list(range(5, 9, 1))))

        self.track_filter = 1  # 1:golden; 2:TPC-only
        self.trigger = 1  # 1 = V0AND !!! Double check this!!!
        self.max_dist_vtx_xy = 2.4
        self.max_dist_vtx_z = 3.2
        # for counter histograms only:

    def Begin(self):
        print 'Data type in Begin(): ', self.data_type

    def SlaveBegin(self, tree):
        print 'py: slave beginning'
        print 'data type in slave: ', self.data_type
        # register counters wrt z index
        zvtx_nbins, zvtx_low, zvtx_high = self.zvtx_bins

        self.event_counter = TH1I('valid', 'Valid events',
                                  zvtx_nbins, zvtx_low, zvtx_high)
        self.tracks_counter = THnSparseF('tracks',
                                         'Valid_tracks;#phi;#eta;z_vtx;pt',
                                         4,
                                         array('i', array('i', [36, 32, 10, 22])),
                                         array('d', array('d', [0, -.8, -10, .5])),
                                         array('d', array('d', [2*pi, .8, 10, 8])))
        self.tracks_counter.SetBinEdges(3, self.pt_edges)
        self.tracks_counter.Sumw2()

        self.GetOutputList().Add(self.event_counter)
        self.GetOutputList().Add(self.tracks_counter)

        # Helper histograms to find bins
        self._bin_zvtx = TH1I('bin_zvtx', 'used to find z vtx',
                              zvtx_nbins, zvtx_low, zvtx_high)

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
        self.fChain.GetBranch('event').GetEntry(entry)
        if not self.validate_event():
            return 1
        if self.data_type in ['real', 'reconstructed']:
            self.fChain.GetBranch('track').GetEntry(entry)
            zvtx = self.fChain.event.zvtx
        elif self.data_type == 'MC':
            self.fChain.GetBranch('trackmc').GetEntry(entry)
            zvtx = self.fChain.event.zvtxMC
        else:
            assert(0)

        self.event_counter.Fill(zvtx, 1)

        tracks = self.get_tracks()
        trig_fill = self.tracks_counter.Fill
        # fill trigger and assoc counters
        [trig_fill(array('d', [phi, eta, zvtx, pt])) for phi, eta, pt in tracks]

        return 0

    def SlaveTerminate(self):
        print 'py: slave terminating'

    def Terminate(self):
        """currently ignoring the pt threshold part"""
        print 'py: terminating'
        print 'fOutput in Terminate', self.GetOutputList().ls()
        from rootpy.io import root_open
        with root_open('output_' + self.data_type + '.root',
                       'recreate') as f:
            # write original histograms:
            f.mkdir('raw')
            f.raw.cd()
            for l in self.GetOutputList():
                l.Write()

    ### Analysis functions

    def get_tracks(self):
        """Return [[phi, eta, pt]]
        for valid tracks
        """
        tracks = []
        if self.data_type == 'MC':
            for track in self.fChain.trackmc:
                if not self.validate_track(track):
                    continue
                pt = track.ptMC
                tracks.append([track.phiMC, track.etaMC, pt])

        elif self.data_type in ['real', 'reconstructed']:
            for track in self.fChain.track:
                if not self.validate_track(track):
                    continue
                pt = track.pt
                tracks.append([track.phi, track.eta, pt])
        return tracks

    def validate_event(self):
        """Validate event on event level"""
        if self.data_type == 'MC':
            if not (-10 < self.fChain.event.zvtxMC < 10):
                return False
        if self.data_type in ['real', 'reconstructed']:
            if not (-10 < self.fChain.event.zvtx < 10):
                return False
        if not (self.fChain.event.trig & self.trigger):
            return False
        if self.fChain.event.vtxstatus < 1:
            return False
        else:
            return True

    def get_zvtx(self):
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


def wrap(x):
    """wrap [rad] angle to [-PI/2..3PI/2)"""
    pi = 3.141592653589793
    return (x + pi/2) % (2 * pi) - pi/2






























