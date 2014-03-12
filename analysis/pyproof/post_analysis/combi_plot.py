from rootpy.io import root_open as r_open
from rootpy import asrootpy
from rootpy.plotting import root2matplotlib as rplt
from rootpy.plotting import Hist
from ROOT import Double

import numpy as np
from math import pi

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import AutoMinorLocator,  MaxNLocator

from extract_yield import extract_yield

#import rootnotes
#%config InlineBackend.figure_format = 'png'
from rootpy.plotting.style import set_style
set_style('ATLAS', mpl=True)
import matplotlib as mpl
# does not handle unicode minus properly!
mpl.rcParams['axes.unicode_minus'] = False


def project_onto_deta_with_sections(hist_2D):
    hists = []
    # near side ridge
    tmp = asrootpy(hist_2D.ProjectionY("nearside_px",
                                       hist_2D.get_xaxis().FindBin(-pi/3),
                                       hist_2D.get_xaxis().FindBin(pi/3)))
    scaling = (hist_2D.get_xaxis().FindBin(pi/3)
               - hist_2D.get_xaxis().FindBin(-pi/3) + 1)
    tmp.Scale(1.0/scaling)
    tmp.SetTitle('near side')
    hists.append(tmp)

    # far side ridge
    tmp = asrootpy(hist_2D.ProjectionY("farside_px",
                                       hist_2D.get_xaxis().FindBin(-pi/3 + pi),
                                       hist_2D.get_xaxis().FindBin(pi/3 + pi)))
    scaling = (hist_2D.get_xaxis().FindBin(pi/3 + pi)
               - hist_2D.get_xaxis().FindBin(-pi/3 + pi) + 1)
    tmp.Scale(1.0/scaling)
    tmp.SetTitle('far side')
    hists.append(tmp)

    # left overs
    tmp1 = asrootpy(hist_2D.ProjectionY("remaining1_px",
                                        0, hist_2D.get_xaxis().FindBin(-pi/3)))
    tmp2 = asrootpy(hist_2D.ProjectionY("remaining2_px",
                                        hist_2D.get_xaxis().FindBin(pi/3),
                                        hist_2D.get_xaxis().FindBin(2*pi/3)))
    tmp3 = asrootpy(hist_2D.ProjectionY("remaining3_px",
                                        hist_2D.get_xaxis().FindBin(4*pi/3),
                                        hist_2D.get_xaxis().FindBin(3*pi/2)))
    tmp1.Add(tmp2)
    tmp1.Add(tmp3)
    
    scaling = (hist_2D.get_xaxis().FindBin(-pi/3)
               - hist_2D.get_xaxis().FindBin(-2*pi/3))  # + 1 No underflow
                                                        #-> don't scale for that bin
    scaling += (hist_2D.get_xaxis().FindBin(2*pi/3)
                - hist_2D.get_xaxis().FindBin(pi/3) + 1)
    scaling += (hist_2D.get_xaxis().FindBin(3*pi/2)
                - hist_2D.get_xaxis().FindBin(4*pi/3))   # +1 No overflow
                                                        #-> don't scale for that bin
    tmp1.Scale(1.0/scaling)
    tmp1.SetTitle('remaining')
    hists.append(tmp1)
    return hists


def mega_plot(hist_2D, ax=None, range_2d=None, yrange_px=None, xrange_py=None,
              details=None, fit=False, exclude_peak=True, plot_fit=False, dphi_extra=None):
    # start with a rectangular Figure
    if ax:
        fig = plt.gcf()
    else:
        fig = plt.figure(1, figsize=(8, 10))
        #ax = plt.axes()
        ax = plt.axes([0.18, 0.18, .64, .64], axisbg='w')
        #ax_debug = plt.axes([0,0,1,1], axisbg='r')
        ax.axis('off')

    frame = 0.0
    fibo = 1.681

    heat_x = heat_y = 1.0 / fibo
    cbar_dist = 0.01
    extra_left_legend = 0.08  # may over extend to right but not left...
    xpro_x, xpro_y = heat_x, 1.0 - 1.0 / fibo
    cbar_x, cbar_y = 0.04, heat_y
    ypro_x, ypro_y = xpro_y - cbar_x - cbar_dist, heat_y
    rect_2d = [frame, frame, heat_x, heat_y]
    rect_histx = [frame, frame + heat_y, xpro_x, xpro_y]
    rect_histy = [frame + heat_x, frame, ypro_x, ypro_y]
    rect_cbar = [frame + heat_x + ypro_x + cbar_dist, frame, cbar_x, cbar_y]
    rect_legend = [frame + heat_x + extra_left_legend, frame + heat_y,
                   ypro_x + cbar_dist + cbar_x, xpro_y]
    # collect plotted lines for legend
    lns = list()

    ax_2d = _add_subplot_axes(ax, rect_2d)
    ax_2d.set_xlabel(r'$\Delta\varphi$')
    ax_2d.set_ylabel(r'$\Delta\eta$')
    if isinstance(range_2d, (list, tuple)):
        im = rplt.imshow(hist_2D, vmin=range_2d[0], vmax=range_2d[1])
    else:
        im = rplt.imshow(hist_2D)

    # x projection:
    axHistx = _add_subplot_axes(ax, rect_histx)
    if exclude_peak:
        lable = '$\Delta \phi$ projection\nexcl. peak'
    else:
        lable = '$\Delta \phi$ projection'
    _plot_x_projection(hist_2D, yrange_px=yrange_px, lns=lns,
                       lable=lable, ax=axHistx, exclude_peak=exclude_peak,
                       plot_fit=plot_fit, dphi_extra=dphi_extra)

    # y axis projection
    sub_py = project_onto_deta_with_sections(hist_2D)
    axHisty = _add_subplot_axes(ax, rect_histy)
    axHisty.set_ylim(bottom=ax_2d.get_ylim()[0], top=ax_2d.get_ylim()[1])
    # cannot use rplt since x and y are swapped
    mini, maxi = 99, -99
    for hist in sub_py:
        mini = hist.GetMinimum() if hist.GetMinimum() < mini else mini
        maxi = hist.GetMaximum() if hist.GetMaximum() > maxi else maxi
        lns.append(plt.errorbar(x=list(hist.y()), y=list(hist.x()),
                        xerr=np.array([list(hist.yerrh()), list(hist.yerrl())]),
                        label=hist.GetTitle()))
    if isinstance(xrange_py, (list, tuple)):
        axHisty.set_xlim(xrange_py)
    else:
        dist = (maxi - mini)*0.1
        axHisty.set_xlim([mini-dist, maxi+dist])
    
    axHisty.set_xlabel(r'$\frac{1}{N_{tr}}\frac{dN_{as}}{d\Delta \eta}$ [rad$^{-1}$]',
                       labelpad=0)
    
    # fix labels and ticks
    nullfmt = NullFormatter()         # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHistx.yaxis.set_major_locator(MaxNLocator(prune='lower', nbins=4))
    axHistx.xaxis.set_minor_locator(AutoMinorLocator())
    axHistx.yaxis.set_minor_locator(AutoMinorLocator())

    axHisty.yaxis.set_major_formatter(nullfmt)
    axHisty.xaxis.set_minor_locator(AutoMinorLocator())
    axHisty.yaxis.set_minor_locator(AutoMinorLocator())
    axHisty.xaxis.set_major_locator(MaxNLocator(prune='lower', nbins=4))
    # get lables for the py projection and rotate them
    axHisty.set_xticklabels(axHisty.get_xticks(),
                            rotation=35, ha='right')
    #plt.setp(axHisty.get_xmajorticklabels(), rotation=35)
    
    # colorbar axes
    cbar_ax = _add_subplot_axes(ax, rect_cbar)
    cb = fig.colorbar(im, cax=cbar_ax, )
    cb.set_label(r'$\frac{1}{N_{tr}}\frac{dN^{2}_{as}}{d\Delta \eta d\Delta \varphi}$[rad$^{-1}$]',) #  labelpad=0 rotation=270)
    
    # make the legend axes:
    ax_leg = _add_subplot_axes(ax, rect_legend)
    labels = [l.get_label() for l in lns]
    #ax_leg.xaxis.set_major_formatter(nullfmt)
    #ax_leg.yaxis.set_major_formatter(nullfmt)
    if details:
        ax_leg.text(0.05, 0.75, details)
    ax_leg.legend(lns, labels, numpoints=1, loc=10)#loc=(.03, -0.03))
    ax_leg.axis('off')

    # write title at center of plot
    #ax.text(0.5, 1.1, hist_2D.title, fontsize=12, transform=axHistx.transAxes,
    #        horizontalalignment='center', verticalalignment='center')
    ax.set_title(hist_2D.title)
    return fig


def get_subtraction(fn, title=None):
    f = r_open(fn, 'read')
    h1 = asrootpy(f.Get('totalYield'+str(0)))
    h2 = asrootpy(f.Get('totalYield'+str(3)))
    sub = h2 - h1
    if title:
        sub.title = title
    return sub


def _plot_x_projection(hist_2D, yrange_px, lable, ax, lns, exclude_peak,
                       plot_fit, dphi_extra):
    sub_px = extract_yield(hist_2D, exclude_peak)[0]
    sub_px.set_title(lable)
    if dphi_extra:
        dphi_extra.color = plt.rcParams.get('axes.color_cycle')[-1]
        lns.append(rplt.errorbar(dphi_extra, xerr=False, label=dphi_extra.title,
                                 marker='s', markersize=3))
        lower = dphi_extra.FindBin(pi/2)
        upper = dphi_extra.FindBin(3*pi/2)
        x = list(dphi_extra.x())[lower-1:upper]  # index != bin
        y1 = list(dphi_extra.y())[lower-1:upper]
        y2 = list(sub_px.y())[lower-1:upper]
        plt.fill_between(x, y1, y2, color='grey', alpha='0.5')
        diff_hist = asrootpy(sub_px - dphi_extra)
        err = Double()
        ryield = diff_hist.IntegralAndError(lower, upper, err, 'width')
        print 'excess ridge yield per rad: ', ryield, err
    lns.append(rplt.errorbar(sub_px, xerr=False, marker='d', markersize=4))
    if plot_fit:
        s_cos = sub_px.GetFunction("cos")
        d_cos = sub_px.GetFunction("double_cos")
        if False:
            print 'Fit parameters:'
            print 'single cosine: ', s_cos.GetExpFormula()
            print 'normalized chi^2: ', s_cos.GetChisquare() / s_cos.GetNDF()
            print 'Parameters: ', s_cos.GetParameter(0), s_cos.GetParameter(1)
            print 'double cosine: ', d_cos.GetExpFormula()
            print 'normalized chi^2: ', d_cos.GetChisquare() / d_cos.GetNDF()
            print 'Parameters: ', d_cos.GetParameter(0), d_cos.GetParameter(1)
        xs = np.arange(list(sub_px.x())[0], list(sub_px.x())[-1], 0.05)
        lns.append(
            plt.plot(xs, [s_cos.Eval(x) for x in xs], label='single cosing',)[0])
        lns.append(
            plt.plot(xs, [d_cos.Eval(x) for x in xs], label='double cosine')[0])
    ax.set_ylabel(
        r'$\frac{1}{N_{tr}}\frac{dN_{as}}{d\Delta \varphi}$'
        '[rad$^{-1}$]')
    ax.yaxis.labelpad = -0
    if isinstance(yrange_px, (list, tuple)):
        ax.set_ylim(yrange_px)
    else:
        ax.set_ylim(0.985*min(sub_px.y()), 1.015*max(sub_px.y()))


def _add_subplot_axes(ax, rect, axisbg='w'):
    """shamelessly stolen from http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib"""
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    # x_labelsize = subax.get_xticklabels()[0].get_size()
    # y_labelsize = subax.get_yticklabels()[0].get_size()
    # x_labelsize *= rect[2]**0.5
    # y_labelsize *= rect[3]**0.5
    # subax.xaxis.set_tick_params(labelsize=x_labelsize)
    # subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
