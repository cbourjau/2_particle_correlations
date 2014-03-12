"""Scale the peripheral class in a way that the near side ridge is
flat (no dip nor peak)"""

from rootpy.io import root_open


def get_subtraction_scalling(fn, thresh):
    """Return scalling factor for "smooth" subtraction of
    peripheral and central events"""
    x_range = (9, 10)  # bin 9 and 10
    y_range_peak = (16, 17)
    y_range1 = (1, 8)
    y_range2 = (25, 32)

    with root_open(fn, 'read') as f:
        h_high = f.Get(('processed/total_yield/'
                        'total_yield_class_0_'
                        'thresh_' + str(thresh*0.5)))
        h_low = f.Get(('processed/total_yield/'
                       'total_yield_class_3_'
                       'thresh_' + str(thresh*0.5)))

        h_high.GetXaxis().SetRange(*x_range)
        h_high.GetYaxis().SetRange(*y_range1)
        non_peak = h_high.Integral()
        h_high.GetYaxis().SetRange(*y_range2)
        non_peak += h_high.Integral()
        non_peak /= (2*2*8)  # divide by nr. of bins
        h_high.GetYaxis().SetRange(*y_range2)
        h_high.GetYaxis().SetRange(*y_range_peak)
        peak = h_high.Integral() / (2*2)
        cent_peak_hight = peak - non_peak

        h_low.GetXaxis().SetRange(*x_range)
        h_low.GetYaxis().SetRange(*y_range1)
        non_peak = h_low.Integral()
        h_low.GetYaxis().SetRange(*y_range2)
        non_peak += h_low.Integral()
        non_peak /= (2*2*8)  # divide by nr. of bins
        h_low.GetYaxis().SetRange(*y_range2)
        h_low.GetYaxis().SetRange(*y_range_peak)
        peak = h_low.Integral() / (2*2)
        peri_peak_hight = peak - non_peak

    return cent_peak_hight / peri_peak_hight

# fn = '/home/christian/msc/results/131116/golden_1-2_2-4/output_real.root'
# get_subtraction_scalling(fn, 0)
