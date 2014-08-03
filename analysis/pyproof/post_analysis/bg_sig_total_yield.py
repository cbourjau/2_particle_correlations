"""Perform post analysis on 'raw' data"""
from rootpy.io import root_open
import logging
import re

# Suppress "WARNING"-level notices from THn not in rootpy
logging.getLogger("rootpy").setLevel(logging.ERROR)


def calc_bg(fn):
    print 'computing backgrounds'
    with root_open(fn, 'update') as f:
        try:
            f.mkdir('processed')
        except ValueError:
            pass
        folders = ['background', 'weighted_background']
        [f.rm('processed/'+folder) for folder in folders]
        [f.mkdir('processed/'+folder) for folder in folders]

        # get signal per z section per eclass
        # NOTE: range(bin_bin, bin_max)
        # Bin number starts at 1 !!!!!!
        for w in ['', 'weighted_']:
            bgs = [f.Get('raw/'+w+'background'+str(i)) for i in range(0, 4)]
            for eclass, bg in enumerate(bgs):
                logging.info('Calculating '+w+'background for class ' + str(eclass))
                for sec_bin in range(1, 11):
                    # get one background per section
                    bg.GetZaxis().SetRange(sec_bin, sec_bin)
                    bg_tmp = bg.Project3D('yx')  # yes, 'yx'...
                    bg_tmp.SetNameTitle((bg.GetName()[:-1]
                                         + '_z_sec_' + str(sec_bin - 1)
                                         + '_class_' + str(eclass)),
                                    bg.GetTitle() + ' z section' + str(sec_bin - 1))
                    scale_background(bg_tmp)
                    f.cd('processed/'+w+'background')
                    bg_tmp.Write()


def calc_signal_and_total_yield(path, correct=None, path_to_mc=''):
    """call after calc_bg"""
    print 'computing signal'
    logging.info('Computing Signals')
    if correct:
        sfx = '_corrected'
    else:
        sfx = ''
    with root_open(path, 'update') as f:
        try:
            f.mkdir('processed')
        except ValueError:
            pass

        folders = ['signal', 'total_yield', 'weighted_signal', 'weighted_total_yield',
                   'total_yield_corrected']
        [f.rm('processed/'+folder) for folder in folders]
        [f.mkdir('processed/'+folder) for folder in folders]
        for w in ['', 'weighted_']:
            signals = [f.Get('raw/'+w+'signal'+str(i)) for i in range(0,4)]
            triggers = f.Get('raw/'+w+'trigger_counter')
            for eclass, sig in enumerate(signals):
                logging.info('Calculating '+w+'signal for class ' + str(eclass))
                total_yield = None
                for sec_bin in range(1, 11):
                    # set z section
                    sig.GetZaxis().SetRange(sec_bin, sec_bin)
                    sig_tmp = sig.Project3D('yx')  # yes, 'yx'...

                    sig_tmp.SetNameTitle((sig.GetName()[:-1]
                                      + '_z_sec_' + str(sec_bin - 1)
                                      + '_class_' + str(eclass)),
                                     (sig.GetTitle() + ' z section'
                                     + str(sec_bin - 1)))
                    f.cd('processed/'+w+'signal')
                    sig_tmp.Write()
                    bg_tmp = f.Get('processed/'+w+'background/'+w+'background'
                                   + '_z_sec_' + str(sec_bin - 1)
                                   + '_class_' + str(eclass))
                    # compute yield
                    sig_tmp.Divide(bg_tmp)  # divide by bg of current z-section
                    if not total_yield:
                        total_yield = sig_tmp
                    else:
                        total_yield.Add(sig_tmp)

                # Get counters:
                n_trig = triggers.GetBinContent(eclass + 1)
                scale = 1.0/(n_trig 
                             * sig.GetXaxis().GetBinWidth(1)
                             * sig.GetYaxis().GetBinWidth(1))
                logging.info('ntrig class ' + str(eclass) + ': '+ str(n_trig))
                logging.info('scalling class ' + str(eclass) + ' by ' + str(scale))
                total_yield.Scale(scale)

                # compute total yield per eclass
                total_yield.SetNameTitle(('total_yield'
                                         + '_class_' + str(eclass)),
                                         ('Total associated yield per' +
                                          'trigger particle\n '
                                          + 'class: ' + str(eclass)))
                # apply corrections
                if correct:
                    apply_corr(total_yield, path_to_mc=path_to_mc)
                    f.cd('processed/total_yield'+sfx)
                else:
                    f.cd('processed/'+w+'total_yield')
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


def apply_corr(hist, path_to_mc):
    # if corr == 'golden':
    #     eff = 0.764984  # central, no thresh, small data set!!!!!!!!!!
    #     print 'Corrected for golden cuts'
    #     hist.Scale(1/eff)
    # elif corr == 'TPC':
    #     eff = 0.909211  # central, no thresh, small data set!!!!!!!!!!
    #     print 'Corrected for TPC-only cuts'
    #     hist.Scale(1/eff)
    # else:
    #     print 'No corrections applied'
    path_eff = path_to_mc.replace('output_MC.root', 'output_effs.root')
    with root_open(path_eff) as f:
        re_class = r'_class_(\d)_'
        re_thresh = r'_thresh_(\d+\.\d)'
        c = re.findall(re_class, hist.GetName())[0]
        t = re.findall(re_thresh, hist.GetName())[0]
        h_eff = f.Get('processed/eff_from_ty/eff_yield_class_' + c +
                      '_thresh_' + t)
        hist.Divide(h_eff.GetFunction('peak'))










