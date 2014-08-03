### selector module (selector.py, name has to match as per in main.py)
from ROOT import TProfile3D
import ROOT
from array import array
import numpy as np
from numpy import pi

import pickle
from rootpy.io import root_open
from rootpy.plotting import Hist3D, Hist1D
from root_numpy import fill_profile


# load MyClass definition macro (append '+' to use ACLiC)
ROOT.gROOT.LoadMacro('DebugClassesMultESA2013.C+')
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


class MyPySelector(ROOT.TPySelector):

    def __init__(self):
        print self.__class__.__module__+": init"
        with open('parameters.pkl', 'read') as pkl:
            args = pickle.load(pkl)
        print 'Settings: ', args
        for (k, v) in vars(args).items():
            setattr(self, k, v)

        if self.track_cut == 'golden':
            self.fn_single_tracks = (
                #'~/msc/results/efficiencies/single_track_eff_golden/single_eff.root')
                '~/msc/results/efficiencies/single_track_eff_golden_z_pt/eff_z_pt_only.root')
            self.track_filter = 1  # 1:golden; 2:TPC-only
            # effs around the tpc border 
            self.fn_eff_sector_neg = (
                '~/msc/cern_stay/results/tpc_border_eff_golden/sector_border_eff_neg.root')
            self.fn_eff_sector_pos = (
                '~/msc/cern_stay/results/tpc_border_eff_golden/sector_border_eff_pos.root')
        if self.track_cut == 'tpc':
            raise NotImplementedError
            self.fn_single_tracks = (
                '~/msc/results/efficiencies/single_track_eff_TPC/single_eff.root')
            self.track_filter = 2  # 1:golden; 2:TPC-only

        self.pool_size = 1000
        self.zvtx_bins = (10, -10, 10)
        self.cent_bins = (4, 0, 100.01)  # edges are just dummies!
        self.cent_edges = array('d', [0, 20, 40, 60, 100.01])
        self.phi_bins = (36, (-pi/2), (3*pi/2))
        self.eta_bins = (32, -1.6, 1.6)
        self.trigger = 1  # 1 = V0AND !!! Double check this!!!
        self.max_dist_vtx_xy = 2.4
        self.max_dist_vtx_z = 3.2

    def Begin(self):
        print 'Data type in Begin(): ', self.data_type

    def SlaveBegin(self, tree):
        # trigger counters (with and without weight
        self.trig_counter = Hist1D(self.cent_edges, name='trigger_counter',
                                   title='Trigger counter')
        self.trig_counter.sumw2()
        self.GetOutputList().Add(self.trig_counter)
        self.trig_counter_weighted = Hist1D(self.cent_edges,
                                            name='weighted_trigger_counter',
                                            title='Weighted trigger counter')
        self.trig_counter_weighted.sumw2()
        self.GetOutputList().Add(self.trig_counter_weighted)

        # count how many events have tracks for which no efficiency is available
        self.no_eff_avail_counter = Hist1D(1, 0, 2, name='no_eff_avail_counter',
                                           title='No efficiency counter')
        self.GetOutputList().Add(self.no_eff_avail_counter)

        self.signals = [] 
        self.signals_weighted = []
        self.backgrounds = []
        self.backgrounds_weighted = []
        self.eta1eta2 = []
        self.signal_pt_profiles = []
        self.background_pt_profiles = []
        for i in range(len(self.cent_edges) - 1):
            self.signals.append(
                Hist3D(*(self.phi_bins + self.eta_bins + self.zvtx_bins), 
                       name='signal'+str(i),
                       title='Signal distribution;phi;eta;zvtx'))
            self.signals[-1].sumw2()

            self.signals_weighted.append(
                Hist3D(*(self.phi_bins + self.eta_bins + self.zvtx_bins), 
                       name='weighted_signal'+str(i),
                       title='Weighted signal distribution;phi;eta;zvtx'))
            self.signals[-1].sumw2()

            self.backgrounds.append(
                Hist3D(*(self.phi_bins + self.eta_bins + self.zvtx_bins),
                       name='background'+str(i),
                       title='Background distribution;phi;eta;zvtx'))
            self.backgrounds[-1].sumw2()

            self.backgrounds_weighted.append(
                Hist3D(*(self.phi_bins + self.eta_bins + self.zvtx_bins),
                       name='weighted_background'+str(i),
                       title='Weighted background distribution;phi;eta;zvtx'))
            self.backgrounds[-1].sumw2()

            #############
            ###   Special profiles   ###
            #############
            self.eta1eta2.append(
                TProfile3D('eta1eta2_pt_prof'+str(i),
                           'eta1 eta2 pt profile ;eta_a;eta_t;pt',
                           32, -.8, .8, 32, -.8, .8, 10, -10, 10))

            self.signal_pt_profiles.append(
                TProfile3D('signal_pt_prof'+str(i), 'Signal pt profile ;phi;eta;pt',
                           *(self.phi_bins + self.eta_bins + self.zvtx_bins)))

            self.background_pt_profiles.append(
                TProfile3D('background_pt_prof'+str(i),'Background pt profile ;phi;eta;pt',
                           *(self.phi_bins + self.eta_bins + self.zvtx_bins)))


        [self.GetOutputList().Add(s) for s in self.signals]
        [self.GetOutputList().Add(s) for s in self.signals_weighted]
        [self.GetOutputList().Add(b) for b in self.backgrounds]
        [self.GetOutputList().Add(b) for b in self.backgrounds_weighted]
        [self.GetOutputList().Add(h) for h in self.eta1eta2]
        [self.GetOutputList().Add(h) for h in self.signal_pt_profiles]
        [self.GetOutputList().Add(h) for h in self.background_pt_profiles]

        # load single weight histogram
        self.single_track_weight = (
            root_open(self.fn_single_tracks, 'read').Get('single_track_eff'))
        # effs around the tpc border 
        self.eff_sector_neg = root_open(self.fn_eff_sector_neg, 'read').Get(
            'pt_phi_neg_counter')
        self.eff_sector_pos = root_open(self.fn_eff_sector_pos, 'read').Get(
            'pt_phi_pos_counter')

        # The pool for each z section
        self.pool = [[]]
        for i in range(0, 10):
            for c in range(0, 4):
                self.pool[-1].append([])
            self.pool.append([])
        # Helper histograms to find bins
        self._bin_zvtx = Hist1D(*self.zvtx_bins, name='bin_zvtx',
                                 title='used to find zvtx')
        self._bin_cent = Hist1D(self.cent_edges, name='bin_cent',
                                title='used to find eclasses')

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
        # Validate on the event branch level
        if self.data_type == 'MC' or self.mc_and_recon_valid:
            if not self.validate_event_mc():
                return 1
            zvtx = self.fChain.event.zvtxMC  # reconstructed zvtx is used if mc_and_rec!
        if self.data_type in ['real', 'reconstructed'] or self.mc_and_recon_valid:
            if not self.validate_event_rec():
                return 1
            zvtx = self.fChain.event.zvtx
        # validate wrt number of triggers        
        if not self.validate_event_wrt_n_triggers():
            return 1

        cent = self.fChain.event.cent

        trigs, assocs = self.get_triggers_assocs()

        # stop here if there are no triggers for sure, those event MUST not make it 
        # into the pool!!!
        if not self.allow_0_trig and len(trigs)==0:
            return 1

        # Good tips for numpy:
        # www.astro.washington.edu/users/vanderplas/Astr599/notebooks/11_EfficientNumpy
        trigs = np.asarray(trigs)
        assocs = np.asarray(assocs)
        """ Example data: (phi, eta, pt, weight)
        assocs = [[ 4.9218154   0.11242714  0.57215816  1.        ]
                  [ 2.26730394  0.60090125  0.55954134  1.        ]]
        trigs  = [[ 4.02582455 -0.11806787  1.0537529   1.        ]]"""

        # Put numpy arrays int the pool
        if self.added_to_pool(assocs):
            return 1

        # the trigger counter must also include events without associated tracks but with
        # a trigger!!! Events of the pool are not included and we can stop after 
        # increasing the counter
        self.trig_counter.fill(cent, len(trigs))
        self.trig_counter_weighted.fill(cent, np.sum(1.0 / trigs[:, 3]))
        if len(assocs) == 0:
            return 1

        cent_id = self.get_cent_id()

        # extend the all arrays to form all permutations
        # Meshgrid gives array of arrays, so flatten them as well
        trigs_phi_ex, assocs_phi_ex = (a.flatten() for a in
                                       np.meshgrid(trigs[:,0], assocs[:,0]))
        trigs_eta_ex, assocs_eta_ex = (a.flatten() for a in
                                       np.meshgrid(trigs[:,1], assocs[:,1]))
        trigs_pt_ex, assocs_pt_ex = (a.flatten() for a in
                                     np.meshgrid(trigs[:,2], assocs[:,2]))
        trigs_weight_ex, assocs_weight_ex = (a.flatten() for a in
                                             np.meshgrid(trigs[:,3], assocs[:,3]))

        dphis = assocs_phi_ex - trigs_phi_ex
        dphis = (dphis + pi/2) % (2 * pi) - pi/2
        detas = assocs_eta_ex - trigs_eta_ex
        weights = 1 / (assocs_weight_ex * trigs_weight_ex)

        # Filling: The array needs to be in shape (3, N)
        # make a new array [[dphis], [detas], [zvtxs]] and 
        # transpose to get [[dphi, deta, zvtx], ...]
        # Example: [[ 0.89599085  0.23049501 -2.4725287 ]
        #           [ 4.5246647   0.71896911 -2.4725287 ]]

        #######   Signal   ##################
        zvtxs = np.full(dphis.shape, zvtx, dtype=float)
        self.signals_weighted[cent_id].fill_array(np.array(
            [dphis, detas, zvtxs]).T, weights)
        self.signals[cent_id].fill_array(np.array(
            [dphis, detas, zvtxs]).T)

        ### Special histograms:
        fill_profile(self.signal_pt_profiles[cent_id], 
                     np.array([dphis, detas, zvtxs, assocs_pt_ex]).T)

        ### eta1 eta2
        fill_profile(self.eta1eta2[cent_id],
                     np.array([assocs_eta_ex, trigs_eta_ex, zvtxs, assocs_pt_ex]).T)

        #######   Background   ################
        bgs = self.get_background_tracks(10*len(assocs))

        trigs_phi_ex, bgs_phi_ex = (a.flatten() for a in
                                       np.meshgrid(trigs[:,0], bgs[:,0]))
        trigs_eta_ex, bgs_eta_ex = (a.flatten() for a in
                                       np.meshgrid(trigs[:,1], bgs[:,1]))
        trigs_pt_ex, bgs_pt_ex = (a.flatten() for a in
                                     np.meshgrid(trigs[:,2], bgs[:,2]))
        trigs_weight_ex, bgs_weight_ex = (a.flatten() for a in
                                             np.meshgrid(trigs[:,3], bgs[:,3]))

        dphis = bgs_phi_ex - trigs_phi_ex
        dphis = (dphis + pi/2) % (2 * pi) - pi/2
        detas = bgs_eta_ex - trigs_eta_ex
        weights = 1 / (bgs_weight_ex * trigs_weight_ex)
        zvtxs = np.full(dphis.shape, zvtx, dtype=float)

        self.backgrounds_weighted[cent_id].fill_array(np.array(
            [dphis, detas, zvtxs]).T, weights)
        self.backgrounds[cent_id].fill_array(np.array([dphis, detas, zvtxs]).T)
        fill_profile(self.background_pt_profiles[cent_id],
                     np.array([dphis, detas, zvtxs, bgs_pt_ex]).T)

        # replace tracks in pool
        self.replace_pool(assocs)

        return 0

    def SlaveTerminate(self):
        print 'py: slave terminating'

    def Terminate(self):
        """currently ignoring the pt threshold part"""
        with root_open(self.directory+'/output_' + self.data_type + '.root',
                       'recreate') as f:
            # write original histograms:
            f.mkdir('raw')
            f.raw.cd()
            for l in self.GetOutputList():
                l.Write()
        print 'fOutput in Terminate', self.GetOutputList().ls()
        print 'Successfully wrote results to output_' + self.data_type + '.root'

    ### Analysis functions

    def get_triggers_assocs(self):
        """Ret ([[trig_phi, trig_eta, pt, w]], [[assoc_phi, assoc_eta, pt, w]])
        for valid tracks
        """
        lower_trig, upper_trig = self.trigger_inter
        lower_assoc, upper_assoc = self.assoc_inter
        triggers, assocs = [], []
        if self.data_type == 'MC':
            tracks = self.get_tracks_generator('MC')
            for track in tracks:
                if not self.validate_track(track):
                    continue
                pt = track.ptMC
                if (lower_assoc < pt < upper_assoc):
                    # MC efficiency is always 1
                    w = 1
                    assocs.append([track.phiMC, track.etaMC, pt, w])
                elif (lower_trig < pt < upper_trig):
                    # MC efficiency is always 1
                    w = 1
                    triggers.append([track.phiMC, track.etaMC, pt, w])
        elif self.data_type in ['real', 'reconstructed']:
            tracks = self.get_tracks_generator(self.data_type)
            for track in tracks:
                if not self.validate_track(track):
                    continue
                pt = track.pt
                if (lower_assoc < pt < upper_assoc):
                    w = self.get_efficiency_for_track(track)
                    if not w:
                        self.no_eff_avail_counter.fill(1)
                        continue
                    assocs.append([track.phi, track.eta, pt, w])
                elif (lower_trig < pt <upper_trig):
                    w = self.get_efficiency_for_track(track)
                    if not w:
                        self.no_eff_avail_counter.fill(1)
                        continue
                    triggers.append([track.phi, track.eta, pt, w])
        return (triggers, assocs)

    def get_efficiency_for_track(self, track):
        """return the single track efficiency (weight) of a given track
        return False if no value is available"""
        w = self.single_track_weight.GetBinContent(
            self.single_track_weight.FindBin(self.fChain.event.zvtx ,track.pt))
        if w < 0.01:  # dont know if int or float
            return False
        return w

    def get_background_tracks(self, n):
        """Return n fitting background tracks"""
        zvtx = self.get_zvtx_id()
        cent = self.get_cent_id()
        pool = self.pool[zvtx][cent]
        np.random.shuffle(pool)
        return pool[:n]    


    def validate_event_mc(self):
        """Validate event on event level"""
        if not (self.fChain.event.trig & self.trigger):
            return False
        if self.fChain.event.vtxstatus < 1:
            return False

        if not (-10 < self.fChain.event.zvtxMC < 10):
            return False
        if self.pt_threshold < 0.0:
            # The soft case
            if self.fChain.event.ptmaxMC > (-1)*self.pt_threshold:
                return False
        else:
            # the hard case
            if self.fChain.event.ptmaxMC < self.pt_threshold:
                return False
        return True

    def validate_event_rec(self):
        """Validate event on event level. This is for MC_rec and real reconstructed"""
        if not (self.fChain.event.trig & self.trigger):
            return False
        if self.fChain.event.vtxstatus < 1:
            return False

        if not (-10 < self.fChain.event.zvtx < 10):
            return False
        if self.pt_threshold < 0.0:
            # The soft case
            if self.fChain.event.ptmax > (-1)*self.pt_threshold:
                return False
        else:
            # the hard case
            if self.fChain.event.ptmax < self.pt_threshold:
                return False
        return True

    def validate_event_wrt_n_triggers(self):
        """Validate the event with respect to optinally given selections regarding
        the event's number of triggers"""
        if self.events_have_recon_triggers:
            if self.get_n_triggers_for_type('reconstructed') < 1:
                return False
        if self.events_have_one_mc_trigger:
            if 1 != self.get_n_triggers_for_type('MC'):
                return False
        if self.mc_and_recon_valid:
            if ((self.get_n_triggers_for_type('MC') < 1) or
                (self.get_n_triggers_for_type('reconstructed') <1)):
                # on or the other does not have valid events
                return False
        return True


    def get_n_triggers_for_type(self, data_type):
        """Return the number of triggers for given type"""
        tracks = self.get_tracks_generator(data_type)
        lower_trig, upper_trig = self.trigger_inter
        triggers = 0
        if data_type == 'MC':
            validator = self._validate_mc_track 
            for track in tracks:
                if not validator(track):
                    continue
                pt = track.ptMC
                if (lower_trig < pt < upper_trig):
                    triggers += 1

        elif data_type == 'reconstructed':
            validator = self._validate_recon_track
            for track in tracks:
                if not validator(track):
                    continue
                pt = track.pt 
                if (lower_trig < pt < upper_trig):
                    triggers += 1
        return triggers


    def added_to_pool(self, assocs):
        """See if the fitting z-vtx pool is already filled"""
        zvtx = self.get_zvtx_id()
        cent = self.get_cent_id()
        if len(self.pool[zvtx][cent]) > self.pool_size:
            return False
        elif len(assocs) == 0:
            # cannot concatenate an empty array, but this event class is not full yet,
            # thus return True
            return True
        else:
            if not len(self.pool[zvtx][cent]):
                self.pool[zvtx][cent]  = assocs
            else:
                self.pool[zvtx][cent] = np.concatenate((self.pool[zvtx][cent], assocs),
                                                       axis=0)
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

    def is_track_at_tpc_border(self, track):
        """Return true if the track is close to the TPC border. "Too close is defined
        as an efficiency < 0.65 which is rather restrictive."""
        # find the bin (both histograms are identical)
        min_eff = 0.65
        if self.data_type == 'MC':
            bin = self.eff_sector_neg.find_bin(track.ptMC, (track.phiMC + pi/18) % (pi/9))
            q = track.qMC
        else:
            bin = self.eff_sector_neg.find_bin(track.pt, (track.phi + pi/18) % (pi/9))
            q = track.q
        if q >= 1:
            if min_eff > self.eff_sector_pos.get_bin_content(bin):
                return True
        elif q <= -1:
            if min_eff > self.eff_sector_neg.get_bin_content(bin):
                return True
        return False

    def validate_track(self):
        """over write this function in SlaveBegin()"""
        pass

    def _validate_mc_track(self, track):
        """return True if the track is valid"""
        if (track.qMC == 0 or (abs(track.etaMC) > 0.8)):
            return False
        if self.exclude_tpc_border and self.is_track_at_tpc_border(track):
            return False
        # Diff charges? if track and requirement have not same sign the mult is negative
        if self.only_charge and self.only_charge*track.qMC < 0:
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
        if self.exclude_tpc_border and self.is_track_at_tpc_border(track):
            return False
        # Diff charges? if track and requirement have not same sign the mult is negative
        if self.only_charge and self.only_charge*track.q < 0:
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
