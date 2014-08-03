"""extract the yield from a given 2D histogram"""

from ROOT import TCutG, TF1
from math import pi
from rootpy import asrootpy


def extract_yield(sub, exclude_peak=True):
    """Extract the yield from the given subtraction hist
    return: TH1F projection onto delta phi with the
            fit functions "cos" and "double_cos" in the histo"""

    peak_cut = TCutG("peak", 5)
    peak_cut.SetPoint(0, -pi/3, 0.8)
    peak_cut.SetPoint(1, pi/3, 0.8)
    peak_cut.SetPoint(2, pi/3, -0.8)
    peak_cut.SetPoint(3, -pi/3, -0.8)
    peak_cut.SetPoint(4, -pi/3, 0.8)

    #nearside
    sub.GetXaxis().SetRange(sub.GetXaxis().FindBin(-pi/3 + 0.01),
                            sub.GetXaxis().FindBin(pi/3 - 0.01))
    y_n_bins = len([type(x) for x in sub.y()])
    if exclude_peak:
        px_near = asrootpy(sub.ProjectionX("near_without_peak",
                                           0, -1, "[-peak] o"))
        px_near.Scale(1.0/(y_n_bins - (sub.get_yaxis().FindBin(0.8 - 0.01)
                                 - sub.get_yaxis().FindBin(-0.8 + 0.01) + 1)))
    else:
        px_near = asrootpy(sub.ProjectionX("near_without_peak", 0, -1, "o"))
        px_near.Scale(1.0/y_n_bins)

    # far side
    sub.GetXaxis().SetRange(sub.GetXaxis().FindBin(2*pi/3 + 0.01),
                            sub.GetXaxis().FindBin(4*pi/3 - 0.01))
    px_far = asrootpy(sub.ProjectionX("farside", 0, -1, "o"))
    px_far.Scale(1.0/y_n_bins)

    # remaining
    sub.GetXaxis().SetRange(1, sub.GetXaxis().FindBin(-pi/3 - 0.01))
    px_remain = asrootpy(sub.ProjectionX("reminder_1", 0, -1, "o"))

    sub.GetXaxis().SetRange(sub.GetXaxis().FindBin(pi/3 + 0.01),
                            sub.GetXaxis().FindBin(2*pi/3 - 0.01))
    px_remain += asrootpy(sub.ProjectionX("reminder_2", 0, -1, "o"))

    sub.GetXaxis().SetRange(sub.GetXaxis().FindBin(4*pi/3 + 0.01), 36)
    px_remain += asrootpy(sub.ProjectionX("reminder_3", 0, -1, "o"))
    px_remain.Scale(1.0/y_n_bins)

    total = px_near + px_far + px_remain
    total.SetStats(0)

    cos = TF1("cos", "[0] + 2*[1]*cos(2*x)")
    cos.SetParameter(0, 0.8)
    cos.SetParameter(1, 1)

    double_cos = TF1("double_cos", "[0] + 2*[1]*cos(2*x) + 2*[2]*cos(3*x)")
    double_cos.SetParameter(0, 0.8)
    double_cos.SetParameter(1, 1)

    total.Fit(cos)
    total.Fit(double_cos, "+")  # + adds this fit to the list of fits in histo
    
    return total, peak_cut















