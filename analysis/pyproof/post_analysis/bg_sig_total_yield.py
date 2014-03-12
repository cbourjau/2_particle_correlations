"""Perform post analysis on 'raw' data"""
from rootpy.io import root_open
from rootpy import asrootpy
from ROOT import Double

import logging
import re
from os.path import dirname
# Suppress "WARNING"-level notices from THn not in rootpy
logging.getLogger("rootpy").setLevel(logging.ERROR)


def calc_bg(fn):
    print 'computing backgrounds'
    with root_open(fn, 'update') as f:
        try:
            f.mkdir('processed')
            f.mkdir('processed/background')
        except ValueError:
            print 'removing old dir'
            f.rm('processed/background')
            f.mkdir('processed/background')

        # get signal per z section per eclass
        # NOTE: range(bin_bin, bin_max)
        # Bin number starts at 1 !!!!!!
        bg_raw = f.Get('raw/background')
        for sec_bin in range(1, 11):
            for eclass_bin in range(1, 5):
                # get one background per section
                bg_raw.GetAxis(2).SetRange(sec_bin, sec_bin)
                bg_raw.GetAxis(3).SetRange(eclass_bin, eclass_bin)
                # include pt_max over and under flow
                bg_raw.GetAxis(4).SetRange(0, 0)
                bg = bg_raw.Projection(1, 0)
                bg.SetNameTitle((bg_raw.GetName()
                                 + '_z_sec_' + str(sec_bin - 1)
                                 + '_class_' + str(eclass_bin - 1)),
                                bg_raw.GetTitle() + ' z section'
                                + str(sec_bin - 1))
                scale_background(bg)
                f.processed.background.cd()
                bg.Write()


def calc_signal_and_total_yield(path, correct=None, assoc_inter=None, cut=''):
    """call after calc_bg"""
    print 'computing signal'
    if correct:
        sfx = '_corrected'
    else:
        sfx = ''
    with root_open(path, 'update') as f:
        try:
            f.mkdir('processed')
            f.mkdir('processed/signal')
            f.mkdir('processed/total_yield' + sfx)
        except ValueError:
            print 'removing old dir'
            f.rm('processed/total_yield' + sfx)
            f.rm('processed/signal')
            f.mkdir('processed/signal')
            f.mkdir('processed/total_yield' + sfx)

        sig_raw = f.Get('raw/signal')
        trig_raw = f.Get('raw/trig')
        threshs = [1] + range(4, 12)  # max of 10 GeV: 20 bins
        if correct:
            eff = get_corr(assoc_inter, cut)
        for pt_thresh in threshs:  # bins of ptmax
            for eclass_bin in range(1, 5):
                total_yield = None
                for sec_bin in range(1, 11):
                    bg = f.processed.background.Get(
                        ('background'
                         + '_z_sec_' + str(sec_bin - 1)
                         + '_class_' + str(eclass_bin - 1)))
                    # set z section, class and ptmax:
                    sig_raw.GetAxis(2).SetRange(sec_bin, sec_bin)
                    sig_raw.GetAxis(3).SetRange(eclass_bin, eclass_bin)
                    sig_raw.GetAxis(4).SetRange(pt_thresh, 21)
                    sig = sig_raw.Projection(1, 0)  # y, x

                    sig.SetNameTitle((sig_raw.GetName()
                                      + '_z_sec_' + str(sec_bin - 1)
                                      + '_class_' + str(eclass_bin - 1)
                                      + '_thresh_' + str((pt_thresh-1)*0.5)),
                                     sig_raw.GetTitle() + ' z section'
                                     + str(sec_bin - 1))
                    f.processed.signal.cd()
                    sig.Write()

                    # compute yield
                    sig.Divide(bg)  # divide by bg of current z-section
                    if not total_yield:
                        total_yield = sig
                        # print ('calculating total yield for event class',
                        #       str(eclass_bin - 1))
                    else:
                        total_yield.Add(sig)

                # Get counters:
                trig_raw.GetAxis(3).SetRange(eclass_bin, eclass_bin)
                trig_raw.GetAxis(4).SetRange(pt_thresh, 21)

                trig = trig_raw.Projection(0)
                n_trig = trig.Integral()
                trig.Delete()
                total_yield.Scale(1.0/(n_trig
                                       * sig.GetXaxis().GetBinWidth(1)
                                       * sig.GetYaxis().GetBinWidth(1)))

                # compute total yield ber threshold and eclass
                total_yield.SetNameTitle(('total_yield'
                                         + '_class_' + str(eclass_bin - 1)
                                         + '_thresh_' + str((pt_thresh-1)*0.5)),
                                         ('Total associated yield per' +
                                          'trigger particle\n '
                                          + 'class: ' + str(eclass_bin - 1)
                                          + 'thresh: ' + str((pt_thresh-1)*0.5)))
                # apply corrections
                if correct:
                    total_yield.Scale(1.0 / eff)
                    f.processed.total_yield_corrected.cd()
                else:
                    f.processed.total_yield.cd()
                total_yield.Write()


def scale_background(hist):
    """scale background ridge to a little less than one due to
    finite bin size. Expects a even bin number in eta"""
    bins_phi = hist.GetXaxis().GetLast()  # probably 36 or so
    #bins_eta = hist.GetYaxis().GetLast()  # probably 32 or so
    hist.GetYaxis().SetRange(16, 17)
    # two bins in eta
    s = 1.0 / (hist.Integral() / float(bins_phi * 2))
    s *= (-0.5 * 0.8 * 0.1 + 1)  # 0.96 see analysis note
    # reset range
    hist.GetYaxis().SetRange(0, 0)
    hist.Scale(s)


def get_corr(assoc_inter, cut):
    if cut == 'golden':
        f_eff_st = root_open('/home/christian/msc/results/efficiencies/'
                             'single_track_eff_golden/single_eff.root',
                             'read')
    if cut == 'TPC':
        f_eff_st = root_open('/home/christian/msc/results/efficiencies/'
                             'single_track_eff_TPC/single_eff.root',
                             'read')
    steff = f_eff_st.single_track_eff
    h_p3 = asrootpy(steff.Projection(3))
    h_p3.Scale(1.0 / (36*32*10))
    err = Double()
    eff = (h_p3.IntegralAndError(h_p3.FindBin(assoc_inter[0]),
                                 h_p3.FindBin(assoc_inter[1]-0.001),
                                 err, 'width')
           / (assoc_inter[1] - assoc_inter[0]))
    print 'correction single value: ', eff, ' +- ', err
    return eff
