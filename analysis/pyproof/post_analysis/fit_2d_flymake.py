from rootpy.io import root_open
from rootpy import asrootpy
from ROOT import TF2, TF1
from math import pi


def fit_2d_efficiency(h):
    """2D fit for near and away side of the given histogram.
    The near side peak is fitted with a 2D Gaussian, the away side with
    a deta independent Gaussian
    """
    peak = TF2("peak", ('[0]+[1]*exp(-0.5*((x-[2])/[3])**2)*' +   # xygaus
                        'exp(-0.5*((y-[4])/[5])**2)' +
                        '+[6]*exp(-0.5*((x-[7])/[8])**2)'),
               (-pi/2), (3*pi/2), -1.5, 1.5)
    hist = h.Clone()
    hist.Rebin2D()
    hist.Scale(1/4.0)
    h_px = hist.ProjectionX()
    h_px.Scale(1.0/16)
    zyam = h_px.GetBinContent(h_px.FindBin(1.3))
    far_max = h_px.GetBinContent(h_px.FindBin(pi))
    peak_max = hist.GetBinContent(hist.FindBin(0, 0))

    peak.SetParameter(0, zyam)
    peak.FixParameter(1, peak_max - zyam)
    peak.SetParameter(2, 0)
    peak.SetParameter('p3', 0.25)
    peak.SetParameter(4, 0)
    peak.SetParameter('p5', 0.3)
    peak.FixParameter(6, far_max - zyam)
    peak.FixParameter(7, pi)
    peak.SetParameter('p8', 0.3)

    hist.Fit(peak, "R")
    f = hist.GetFunction('peak')
    f.SetRange(hist.GetXaxis().GetXmin(), -1.6,
               hist.GetXaxis().GetXmax(), 1.6)
    print 'Reduced chi^2: ', f.GetChisquare()/f.GetNDF()
    h.GetListOfFunctions().Add(f)


def func_to_hist(f, h=None):
    """Replace bin contents of the hists with values of f at the bin centers"""
    if not h:
        from rootpy.plotting import Hist2D
        h =Hist2D(36, (-pi/2), (3*pi/2), 32, -1.6, 1.6)
    for bin in asrootpy(h).bins():
        x = bin.x.center
        y = bin.y.center
        bin.value = f.Eval(x, y)
        bin.error = 0
    return h


if __name__ == '__main__':
    f = root_open('/home/christian/msc/results/140120/pool_with_4GeV_threshold/' + 
                  'output_effs.root', 'read')
    eclass = 3
    t = 4.0
    h = asrootpy(f.processed.eff_from_ty.get('eff_yield_class_' + str(eclass) +
                                             '_thresh_' + str(t)))
    fit_2d_efficiency(h)
    h.Draw('surf1')










