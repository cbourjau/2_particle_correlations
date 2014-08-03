"""Run the post analysis"""
from bg_sig_total_yield import calc_bg, calc_signal_and_total_yield
from post_subtraction import high_minus_low
from post_efficiencies import calc_effs

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--real', type=str, help='Path to real data')
parser.add_argument('--mc', type=str, help='Path to MC data')
parser.add_argument('--reconstructed', type=str,
                    help='Path to reconstructed data')
parser.add_argument('--eff_out_dir',  default=None,
                    help='Path to safe efficiencies to')
parser.add_argument('--recalc_effs', action='store_true',
                    help='Recalculate efficiencies from MC and reconstructed')
args = parser.parse_args()


path_real = args.real
path_mc = args.mc
path_recon = args.reconstructed

if not (path_real or path_mc or path_recon):
    assert(0)
if path_recon:
    print 'calc total yield for reconstructed'
    calc_bg(path_recon)
    calc_signal_and_total_yield(path_recon)

if path_mc:
    print 'calc total yield for MC'
    calc_bg(path_mc)
    calc_signal_and_total_yield(path_mc)

if path_mc and path_recon:
    print 'calc 2D efficiencies from total yield and signal distributions'
    calc_effs(path_mc, path_recon)

if path_real:
    print 'calc UNcorrected real'
    calc_bg(path_real)
    calc_signal_and_total_yield(path_real)

# print 'calc corrected real'
# fn = 'output_real.root'
# #calc_bg(path_real)
# #calc_signal_and_total_yield(path_real, correct=True, path_to_mc=path_mc)

# print 'calculate uncorrected and unscaled high minus low histos'
# fn = 'output_real.root'
# high_minus_low(path_real)

# print 'calculate corrected and unscaled high minus low histos'
# fn = 'output_real.root'
# #high_minus_low(path_real, correct=True)

# print 'calculate high minus low histos SCALED'
# fn = 'output_real.root'
# high_minus_low(path_real, scalling=True)

# print 'calculate high minus low histos SCALED and corrected'
# fn = 'output_real.root'
# #high_minus_low(path_real, scalling=True, correct=True)

# #print 'extract yield over entire phi'

# #print 'extract yield only from near side'




















