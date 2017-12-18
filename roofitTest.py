#!/usr/bin/env python

from __future__ import print_function, division

import ROOT

# w = ROOT.RooWorkspace("w")
# w.factory("Gaussian::g(x[-10,10],mean[0],sigma[3])")

x = ROOT.RooRealVar("x", "x", -10., 10.)
m = ROOT.RooRealVar("m", "m", 0.)
s = ROOT.RooRealVar("s", "s", 3.)
g = ROOT.RooGaussian("g", "gauss(x, m, s)", x, m, s)

g.Print()

outputCanvas = ROOT.TCanvas("outputCanvas", "outputCanvas")
outputCanvas.cd()
frame = x.frame()
g.plotOn(frame)
frame.Draw()
outputCanvas.SaveAs("test.png")
