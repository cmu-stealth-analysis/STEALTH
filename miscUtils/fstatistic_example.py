#!/usr/bin/env python

from __future__ import print_function, division

import os, math, subprocess
import scipy.stats
import ROOT
import stealthEnv

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

outputDirectory = stealthEnv.analysisRoot + "/ftest"
if not(os.path.isdir(outputDirectory)): subprocess.check_call("mkdir -p {oD}".format(oD=outputDirectory), shell=True, executable="/bin/bash")

nPoints = 1000
xMin = 0.05
xMax = 50.
d1 = 1
d2 = 2
f_distribution = ROOT.TGraph()
f_distribution.SetName("f_dist")
f_distribution.SetTitle("F-distribution;F;density")
f_cdf = ROOT.TGraph()
f_cdf.SetName("f_cdf")
f_cdf.SetTitle("F-distribution CDF;F;")

for xIndex in range(0, 1+nPoints):
    x = xMin*math.pow(10, (xIndex/nPoints)*math.log10(xMax/xMin))
    # x = xMin + (xIndex/nPoints)*(xMax-xMin)
    pdf_value = scipy.stats.f.pdf(x, d1, d2)
    f_distribution.SetPoint(xIndex, x, pdf_value)
    cdf_value = scipy.stats.f.cdf(x, d1, d2)
    f_cdf.SetPoint(xIndex, x, cdf_value)
    # print("At x = {x}, pdf_value: {p}, cdf_value: {c}".format(x=x, p=pdf_value, c=cdf_value))

for graphToPlot in [f_distribution, f_cdf]:
    outputCanvas = ROOT.TCanvas("c_" + graphToPlot.GetName(), "")
    graphToPlot.Draw("AC")
    ROOT.gPad.SetLogx()
    outputCanvas.SaveAs(outputDirectory + "/" + graphToPlot.GetName() + ".pdf")
