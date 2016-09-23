#!/usr/bin/python

import os
import sys
import array
import numpy as np
import ROOT
import argparse
from scipy.stats import chisquare

iJtScale = 1

xMin = 1000.
xMax = 3500.
nJtMin = 2 
nJtMax = 7 
# jet index for normalization/ratio comparison
iJtScale = 0
# jet index for determining bkg scaling
iJtBkg = 1
# bin index for determining bkg scaling
iBinBkg = 2

c1 = ROOT.TCanvas("c1","c1",600,600)
c1.SetBorderSize(0);
c1.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
c1.cd()
ROOT.gPad.SetLogy()

hST = []
hFile = ROOT.TFile("hFile.root","READ")
for j in range(nJtMin,nJtMax+1):
    hST.append( ROOT.gDirectory.Get(str(j)+"jet") )

hST[iJtBkg].SetTitle("1#gamma, "+hST[iJtBkg].GetName())
hST[iJtBkg].GetYaxis().SetRangeUser(1.1,5.e+03)
hST[iJtBkg].GetXaxis().SetTitle("S_{T} [GeV]")
hST[iJtBkg].GetXaxis().SetTitleOffset(1.1)
hST[iJtBkg].SetMarkerStyle(20)
hST[iJtBkg].SetMarkerSize(.9)
hST[iJtBkg].Draw("E")
c1.Update()

leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
leg.AddEntry(hST[iJtBkg],"Bkg","LP")

# Estimate bkgrnd analytically
StBkg = []
StBkg.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x))",xMin,xMax) )
#StBkg.append( ROOT.TF1("fSt0","TMath::Sqrt([0]*[0]/TMath::Power(x/13000.,2*[1]*TMath::Log(x)))",xMin,xMax) )
StBkg.append( ROOT.TF1("fSt1","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
StBkg.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )
xMaxFit = 2500.
for iBkg in range(0,len(StBkg)):
	xMinFit = xMin
	xMaxFit = xMax
	status = int( hST[iJtBkg].Fit("fSt"+str(iBkg),"M0N","",xMinFit,xMaxFit) )

	print status
	if status != 0 and status != 4000:
		print "WARNING: fit failed with status:",status
		continue

iPar=1
i = 0
for St in StBkg:
	print "p"+str(i+1)+"="+str(St.GetParameter(iPar))+"+-"+str(St.GetParError(iPar))
	print "chiSq / ndf =",St.GetChisquare(),"/",St.GetNDF()
	St.SetLineWidth(2)
	St.SetLineColor(i+1)
	St.Draw("SAME")
	c1.Update()
	i += 1

# Draw legend
i = 0 
label = []
label.append("1/x^{p_{1}lnS_{t}}")
label.append("1/x^{p_{2}}")
label.append("1/e^{p_{3}x}")
for St in StBkg:
    leg.AddEntry(St,label[i],"LP")
    i += 1
leg.SetBorderSize(0);
leg.Draw("SAME")
c1.Update()
