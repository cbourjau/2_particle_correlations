from rootpy.io import root_open

def calc_single_track_effs(path):
    """Calc the 2D effs for trigger and assocs"""
        with root_open(path + 'output_MC.root', 'read') as f_mc:
        with root_open(path + 'output_reconstructed.root', 'read') as f_re:
            with root_open(path + 'output_effs.root', 'update') as f_effs:
                f_effs.mkdir('processed')
                f_effs.mkdir('processed/efficiencies')

                f_effs.processed.efficiencies.cd()
                # names are identical in mc and recon
                # this might be a problem!
                print 'calculating effs'
                for tpl_mc in f_mc.walk('processed/total_yield'):
                    names = tpl_mc[2]
                    names.sort()
                    for h_name in tpl_mc[2]:
                        ty_mc = f_mc.Get('processed/total_yield/' + h_name)
                        ty_re = f_re.Get('processed/total_yield/' + h_name)
                        eff = ty_re / ty_mc
                        name = 'eff_' + ty_re.name[6:]
                        eff.SetNameTitle(name,
                                         'efficiency' + ty_re.name[6:])
                        f_effs.Write(name)
                        err = r.Double()
                        x_bins, y_bins = (eff.GetXaxis().GetNbins(),
                                          eff.GetYaxis().GetNbins())
                        inte = eff.IntegralAndError(1, x_bins, 1, y_bins, err)
                        print (eff.name,
                               inte / (x_bins * y_bins), 'Error: ',
                               err / (x_bins * y_bins))
