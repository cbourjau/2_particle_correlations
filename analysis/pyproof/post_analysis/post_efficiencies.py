"""Calculate the efficiency for each event class from total yield"""
from rootpy.io import root_open
from rootpy import asrootpy
import ROOT
from os.path import dirname, abspath
from fit_2d import fit_2d_efficiency
import numpy as np


def calc_effs(path_mc, path_recon, out_dir=None):
    """
    Save the calculated effs for each z-section and eventclass
    to file
    """
    with root_open(path_mc, 'read') as f_mc:
        with root_open(path_recon, 'read') as f_re:
            if not out_dir:
                out_dir = dirname(abspath(path_mc))
            with root_open(out_dir + '/output_effs.root', 'recreate') as f_effs:
                f_effs.mkdir('processed')
                f_effs.mkdir('processed/eff_from_ty')
                f_effs.mkdir('processed/eff_from_signal')

                print 'calculating effs from total yield'
                f_effs.processed.eff_from_ty.cd()
                for tpl_mc in f_mc.walk('processed/total_yield'):
                    names = tpl_mc[2]
                    names.sort()
                    for h_name in tpl_mc[2]:
                        ty_mc = f_mc.Get('processed/total_yield/' + h_name)
                        ty_re = f_re.Get('processed/total_yield/' + h_name)
                        eff = ty_re / ty_mc
                        fit_2d_efficiency(eff)
                        name = 'eff_' + ty_re.name[6:]
                        eff.SetNameTitle(name,
                                         'efficiency' + ty_re.name[6:])
                        f_effs.Write(name)
                        err = ROOT.Double()
                        x_bins, y_bins = (eff.GetXaxis().GetNbins(),
                                          eff.GetYaxis().GetNbins())
                        inte = eff.IntegralAndError(1, x_bins, 1, y_bins, err)
                        print (eff.name,
                               inte / (x_bins * y_bins), 'Error: ',
                               err / (x_bins * y_bins))


def get_2d_eff_from_signal(f_mc, f_re, eclass, thresh, sec):
    """Return 2d hist for given parameters"""
    #import pdb; pdb.set_trace()
    trig_mc = f_mc.raw.get('trig')
    trig_re = f_re.raw.get('trig')
    try:
        sig_mc = asrootpy(f_mc.processed.signal.get(
            'signal_z_sec_' + str(sec) +
            '_class_' + str(eclass) +
            '_thresh_'+str(thresh)))
        sig_re = asrootpy(f_re.processed.signal.get(
            'signal_z_sec_' + str(sec) +
            '_class_'+str(eclass) +
            '_thresh_'+str(thresh)))
    except:
        assert(0)
    # Getting the number of triggers for recon and MC
    trig_mc.GetAxis(2).SetRange(sec+1, sec+1)
    trig_mc.GetAxis(3).SetRange(eclass+1, eclass+1)
    trig_mc.GetAxis(4).SetRange(int(thresh*2), 21)
    trig_pro_mc = trig_mc.Projection(0)
    n_trig_mc = trig_pro_mc.Integral()
    trig_pro_mc.Delete()

    trig_re.GetAxis(2).SetRange(sec+1, sec+1)
    trig_re.GetAxis(3).SetRange(eclass+1, eclass+1)
    trig_re.GetAxis(4).SetRange(int(thresh*2), 21)
    trig_pro_re = trig_re.Projection(0)
    n_trig_re = trig_pro_re.Integral()
    trig_pro_re.Delete()

    eff = ((sig_re / sig_mc) * (n_trig_mc / n_trig_re))
    name = 'eff_' + sig_mc.name
    eff.SetNameTitle(
        name, 'efficiency ' + sig_mc.name)
    # eff.SetDirectory(0)
    trig_mc.Delete()
    trig_re.Delete()
    return eff


def get_2d_eff_from_bg(f_mc, f_re, eclass, sec):
    """Return 2d hist for given parameters"""
    try:
        bg_mc = f_mc.processed.background.get(
            'background_z_sec_' + str(sec) +
            '_class_' + str(eclass))
        bg_re = f_re.processed.background.get(
            'background_z_sec_' + str(sec) +
            '_class_'+str(eclass))
    except:
        assert(0)
    # normalize to on on the ridge along phi
    bg_mc.GetYaxis().SetRange(16, 17)
    bg_mc.Scale((36 * 2) / bg_mc.Integral())
    bg_mc.GetYaxis().SetRange(0, 0)

    bg_re.GetYaxis().SetRange(16, 17)
    bg_re.Scale((36 * 2) / bg_re.Integral())
    bg_re.GetYaxis().SetRange(0, 0)

    eff = bg_re / bg_mc
    name = 'eff_' + bg_mc.name
    eff.SetNameTitle(
        name, 'efficiency ' + bg_mc.name)
    #eff.SetDirectory(0)
    del bg_re
    del bg_mc
    return eff


















