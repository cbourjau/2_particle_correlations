cimport cython
from array import array
#from math import pi
import numpy as np


@cython.binding(True)
def fill_histograms(self, trigs, assocs, double zvtx, double cent, double ptmax):
    # cache bound methods:
    # ftp://root.cern.ch/root/doc/19PythonRuby.pdf
    s_fill = self.signal.Fill
    bg_fill = self.background.Fill
    cdef double trig_phi, trig_eta, assoc_phi, assoc_eta, bg_phi, bg_eta
    cdef int i
    bgs_it = self.it_background_tracks()
    for trig_phi, trig_eta, pt in trigs:
        for assoc_phi, assoc_eta, pt in assocs:
            s_fill(array('d', [wrap(assoc_phi-trig_phi),
                               assoc_eta-trig_eta,
                               zvtx,
                               cent,
                               ptmax]))
        # 10 bg tracks per assoc
        for i in xrange(len(assocs)*10):
            bg_phi, bg_eta, pt = bgs_it.next()
            a = array('d', [wrap(bg_phi - trig_phi),
                            bg_eta - trig_eta,
                            zvtx, cent, ptmax])
            bg_fill(a)

@cython.binding(True)
def get_triggers_assocs(self):
    """Return ([[trig_phi, trig_eta, pt]], [[assoc_phi, assoc_eta, pt]])
    for valid tracks
    """
    lower_trig, upper_trig = self.trigger_inter
    lower_assoc, upper_assoc = self.assoc_inter
    triggers, assocs = [], []
    if self.data_type == 'MC':
        for track in self.fChain.trackmc:
            if not self.validate_track(track):
                continue
            pt = track.ptMC
            if (lower_assoc < pt < upper_assoc):
                assocs.append([track.phiMC, track.etaMC, pt])
            elif (lower_trig < pt < upper_trig):
                triggers.append([track.phiMC, track.etaMC, pt])
    elif self.data_type in ['real', 'reconstructed']:
        for track in self.fChain.track:
            if not self.validate_track(track):
                continue
            pt = track.pt
            if (lower_assoc < pt < upper_assoc):
                assocs.append([track.phi, track.eta, pt])
            elif (lower_trig < pt < upper_trig):
                triggers.append([track.phi, track.eta, pt])
    return (triggers, assocs)

@cython.binding(True)
@cython.boundscheck(False)
def it_background_tracks(self):
    """Iterator for background tracks"""
    cdef int zvtx = self.get_zvtx()
    cdef int cent = self.get_cent_id()
    pool = self.pool[zvtx][cent]
    while True:
        np.random.shuffle(pool)
        for track in pool:
            yield track


cdef double wrap(double x):
    """wrap [rad] angle to [-PI/2..3PI/2)"""
    cdef double pi = 3.141592653589793
    return (x + pi/2) % (2 * pi) - pi/2










