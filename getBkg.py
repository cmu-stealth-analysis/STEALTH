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
iJtBkg = 0

## MAIN ##
def main():
	'''
	c1 = ROOT.TCanvas("c1","c1",600,600)
	c1.SetBorderSize(0);
	c1.SetFrameBorderMode(0)
	ROOT.gStyle.SetTitleBorderSize(0)
	ROOT.gStyle.SetOptStat(0)
	c1.cd()
	ROOT.gPad.SetLogy()
	'''

	hST = []
	hFile = ROOT.TFile("hFile.root","READ")
	for j in range(nJtMin,nJtMax+1):
			hST.append( ROOT.gDirectory.Get("h"+str(j)+"jet") )

	'''
	hST[iJtBkg].SetTitle("1#gamma, "+hST[iJtBkg].GetName())
	#hST[iJtBkg].GetYaxis().SetRangeUser(1.1,5.e+03)
	hST[iJtBkg].GetXaxis().SetTitle("S_{T} [GeV]")
	hST[iJtBkg].GetXaxis().SetTitleOffset(1.1)
	hST[iJtBkg].SetMarkerStyle(20)
	hST[iJtBkg].SetMarkerSize(.9)
	hST[iJtBkg].Draw("E")
	c1.Update()

	leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
	leg.AddEntry(hST[iJtBkg],"Bkg","LP")
	'''
	# Estimate bkgrnd analytically
	StBkgs = []
	#StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x/13000.))",xMin,xMax) )
	StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x))",xMin,xMax) )
	StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
	StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )
	StBkgs.append( ROOT.TF1("fSt3","[0]/TMath::Exp([1]*x/13000. + [2]*pow(x,3.))",xMin,xMax) )
	#StBkgs.append( ROOT.TF1("fSt3","[0]/TMath::Exp([1]*pow(x,3.))",xMin,xMax) )
	i = 0
	for StBkg in StBkgs:
		xMinFit = xMin
		xMaxFit = xMax
		#xMaxFit = 2500.
		#if i == 3:
		#	xMaxFit = 3500.
		status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0N","",xMinFit,xMaxFit) )

		print status
		if status != 0 and status != 4000:
			print "WARNING: fit failed with status:",status
			#continue
		print "chiSq / ndf =",StBkg.GetChisquare(),"/",StBkg.GetNDF()
		print "====================="
		i += 1

	fScales = 'ScaleFactor.txt'
	with open(fScales, "r") as f:
	  scales = f.readlines()
	i = 0
	for h in hST:
		#if i > 0:
		#	continue
		print scales[i]
		plotHistovFit(h,StBkgs,float(scales[i]))
		os.rename("bkgfit.png","DATA/bkgFit_"+str(i+2)+"jet.png")
		i += 1

def plotHistovFit(hST,gFits,scale):

	c = ROOT.TCanvas("c","c",600,600)
	c.SetBorderSize(0);
	c.SetFrameBorderMode(0)
	ROOT.gStyle.SetTitleBorderSize(0)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gPad.SetLogy()
	c.cd()

	hST.SetTitle("1#gamma, "+hST.GetName())
	hST.GetYaxis().SetRangeUser(1.1e-01,4.e+03)
	hST.GetXaxis().SetTitle("S_{T} [GeV]")
	hST.GetXaxis().SetTitleOffset(1.1)
	hST.SetMarkerStyle(20)
	hST.SetMarkerSize(.9)
	hST.Draw("E")
	#c.Update()

	leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
	leg.AddEntry(hST,"Bkg","LP")

	iPar=0
	i = 0
	StBkgs = []
	for gFit in gFits:
		StBkgs.append(gFit.Clone())
		StBkgs[-1].SetParameter(iPar,float(StBkgs[-1].GetParameter(iPar))/scale)
		StBkgs[-1].SetLineWidth(2)
		StBkgs[-1].SetLineColor(i+1)
		#StBkgs[-1].Draw("SAME")
		c.Update()
		i += 1

	# Draw legend
	i = 0 
	label = []
	label.append("1/x^{p_{1}lnS_{t}}")
	#label.append("1/x^{p_{1}ln(x)}")
	label.append("1/x^{p_{2}}")
	label.append("1/e^{p_{3}x}")
	label.append("1/e^{p_{4}x + p_{5}S_{T}^{3}}")
	for StBkg in StBkgs:
			leg.AddEntry(StBkg,label[i],"LP")
			i += 1
	leg.SetBorderSize(0);
	#leg.Draw("SAME")
	c.Update()

	c.Print("bkgfit.png")

#_____ Call main() ______#
if __name__ == '__main__':
  main()