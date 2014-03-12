from ROOT import TH2F, TF2, Double
from math import pi
#from extract_yield import extract_yield

peak = TF2("peak", ('[0]+[1]*exp(-0.5*((x-[2])/[3])**2)*' +
                    'exp(-0.5*((y-[4])/[5])**2)'),
           -2*pi/3, 1.3, -1.5, 1.5)
dijet = TF2("dijet", ('[0]+[1]*exp(-0.5*((x-[2])/[3])**2)*' +
                      'exp(-0.5*((y-[4])/[5])**2)' +
                      '+[6]*exp(-0.5*((x-[7])/[8])**2)'),
            -2*pi/3, 4*pi/3, -1.5, 1.5)

dijet.SetParameter('p0', 0.9)
dijet.SetParameter('p1', 1)
dijet.FixParameter(2, 0)
dijet.SetParameter('p3', 0.7)
dijet.FixParameter(4, 0)
dijet.SetParameter('p5', 0.3)
dijet.SetParameter('p6', 0.5)
dijet.FixParameter(7, pi)
dijet.SetParameter('p8', 0.3)

h = TH2F("h2", "from f2", 36, -2*pi/3, 4*pi/3,
         32, -1.6, 1.6)
h.FillRandom("dijet", 10000000)
h.Scale(1/10000.0)
h_px = h.ProjectionX()
h_px.Scale(1.0/32)
zyam = h_px.GetBinContent(h_px.FindBin(1.3))
far_max = h_px.GetBinContent(h_px.FindBin(pi))
peak_max = h.GetBinContent(h.FindBin(0, 0))

peak.FixParameter(0, zyam)
peak.FixParameter(1, peak_max - zyam)
peak.FixParameter(2, 0)
peak.SetParameter('p3', 0.7)
peak.FixParameter(4, 0)
peak.SetParameter('p5', 0.3)
peak.FixParameter(6, far_max - zyam)
peak.FixParameter(7, pi)
peak.SetParameter('p8', 0.3)

h.Fit("peak", "R")
h.Draw('surf1')

