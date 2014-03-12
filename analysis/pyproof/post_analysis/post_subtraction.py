"""Subtract the low multiplicity class from the high one"""

from rootpy.io import root_open, DoesNotExist
from rootpy.plotting import Hist
from scale_subtraction import get_subtraction_scalling


def high_minus_low(fn, scalling=None, correct=None):
    h_scale = Hist(21, -0.25, 10.25, name="scalling_threshold",
                   title="Threshold vs scalling factor;threshold;scalling")
    with root_open(fn, 'update') as f:
        sfx = ''
        h_name_extra = ''
        if scalling:
            sfx += '_with_scalling'
        if correct:
            sfx += '_corrected'
            h_name_extra = '_corrected'
        try:
            f.mkdir('processed')
            f.mkdir('processed/subtracted' + sfx)
        except:
            f.rm('processed/subtracted' + sfx)
            f.mkdir('processed/subtracted' + sfx)
        # loop through all thresholds and catch exceptions
        for thresh in range(0, 21):
            try:
                h_high = f.Get(('processed/total_yield' + h_name_extra + '/'
                                'total_yield_class_0_'
                                'thresh_' + str(thresh*0.5)))
                h_low = f.Get(('processed/total_yield' + h_name_extra + '/'
                               'total_yield_class_3_'
                               'thresh_' + str(thresh*0.5)))
            except DoesNotExist:
                continue
            if scalling:
                c = get_subtraction_scalling(fn, thresh)
                h_scale.Fill(thresh*0.5, c)
            else:
                c = 1
            sub = h_high - h_low * c
            sub.SetNameTitle('Subtraction_thresh_' + str(thresh*0.5),
                             ('Central minus peripheral at threshold '
                              + str(thresh*0.5) + ' GeV'))
            if scalling:
                sub.title += ' peripheral scaled by ' + str(c)
            f.cd('processed/subtracted' + sfx)
            sub.Write()
        if scalling:
            h_scale.Write()
















