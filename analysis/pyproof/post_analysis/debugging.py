from rootpy.io import root_open
from rootpy import asrootpy
from rootpy.plotting import root2matplotlib as rplt
from rootpy.plotting import Canvas

import sys
sys.path.append('/home/christian/msc/analysis/pyproof/post_analysis/')
#import ROOT
#from bg_sig_total_yield import calc_bg, calc_signal_and_total_yield
#from post_subtraction import high_minus_low
from post_efficiencies import calc_effs, get_2d_eff_from_signal, get_2d_eff_from_bg
from fit_2d import fit_2d_efficiency, func_to_hist
from combi_plot import mega_plot

#%load_ext autoreload
#%autoreload 2


# rcParams dict
import matplotlib as mpl
import matplotlib.pyplot as plt
colors = [(0.4, 0.7607843137254902, 0.6470588235294118), (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
          (0.5529411764705883, 0.6274509803921569, 0.796078431372549), (0.9058823529411765, 0.5411764705882353, 0.7647058823529411),
          (0.6509803921568628, 0.8470588235294118, 0.32941176470588235), (1.0, 0.8509803921568627, 0.1843137254901961),
          (0.8980392156862745, 0.7686274509803922, 0.5803921568627451)]
mpl.rcParams['axes.color_cycle'] = colors

font_params = {"axes.labelsize": 8,
                "axes.titlesize": 9,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
                "legend.fontsize": 8,
                'font.family': 'serif',
                'text.usetex': True,
                'font.serif': 'Computer Modern Roman',
                }
mpl.rcParams.update(font_params)

ticksize = 3.
tickwidth = .5
ax_params = {"axes.grid": False,
             "axes.facecolor": "white",
             "axes.edgecolor": "black",
             "axes.linewidth": 1,
             "xtick.direction": "in",
             "ytick.direction": "in",
             "xtick.major.width": tickwidth,
             "ytick.major.width": tickwidth,
             "xtick.minor.width": tickwidth,
             "xtick.minor.width": tickwidth,
             "xtick.major.size": ticksize,
             "xtick.minor.size": ticksize / 2,
             "ytick.major.size": ticksize,
             "ytick.minor.size": ticksize / 2}
mpl.rcParams.update(ax_params)

mpl.rcParams.update({
    "lines.linewidth": 1.1 ,
    "patch.linewidth": .1 ,
    "xtick.major.pad": 3.5 ,
    "ytick.major.pad": 3.5 ,
    })

# Set the constant defaults
mpl.rc("legend", frameon=False, numpoints=1)
mpl.rc("lines", markeredgewidth=0, solid_capstyle="round")
mpl.rc("figure", figsize=(2.95,2.95))
mpl.rc("image", cmap="jet")
mpl.rcParams['text.latex.preamble']=[r"\usepackage{siunitx}",]


#helpers
pta = r'p_{\text{T}}^{\text{asso}}'
ptt = r'p_{\text{T}}^{\text{trig}}'
pthresh = r'p_{\text{T}}^{\text{tresh}}'

f_effs = root_open('/home/christian/msc/results/finals/only_pions/output_effs.root', 'read')
f_mc = root_open('/home/christian/msc/results/finals/only_pions/output_MC.root', 'read')
f_re = root_open('/home/christian/msc/results/finals/only_pions/output_reconstructed.root', 'read')

plt.close('all')
figsize= (2.95,2.95)
eclass = 0
t = 0.0
plt.close('all')
fig = plt.figure(figsize=figsize)
h = asrootpy(f_effs.processed.eff_from_ty.get('eff_yield_class_' + str(eclass) + '_thresh_' + str(t)))
h.title = 'MC closure of total associated yield $Y$,\nHard events ('+pthresh+' 4.0GeV), $0-20\%,$'
h.Rebin2D(); h.Scale(1.0/4.0)
mega_plot(h, exclude_peak=False,)#yrange_px=(1.003, 1.015), range_2d=(0.98, 1.02), xrange_py=(0.996, 1.016))
