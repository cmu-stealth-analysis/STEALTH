#!/usr/bin/python

import os
import sys
import array
import numpy as np
import ROOT
import argparse
from scipy.stats import chisquare

iJtScale = 1

#xMin = 1000.
#xMax = 3500.
xMin = 1250.
xMax = 3750.
nJtMin = 2 
nJtMax = 7 
# jet index for normalization/ratio comparison
iJtScale = 0
# jet index for determining bkg scaling
iJtBkg = 0

## MAIN ##
#def main():

hST = []
hFile = ROOT.TFile("hFile.root","READ")
for j in range(nJtMin,nJtMax+1):
		hST.append( ROOT.gDirectory.Get("h"+str(j)+"jet") )
		#hST[-1].Sumw2()

# Estimate bkgrnd analytically
StBkgs = []
#StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Power(x/13000,TMath::Log(x/13000.)*[3]))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Log(x/13000.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000. + [2]*TMath::Log(x/13000.)))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Power(x/13000,[3]))",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000,[1]) + [2]/TMath::Exp([3]*x/13000.)",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x/13000.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","sqrt([0]*[0]/TMath::Power(x/13000,2*[1]) + [2]*[2]/TMath::Exp(2*[3]*x/13000.))",xMin,xMax) )
#StBkgs[-1].FixParameter(0,9.4e04)
#StBkgs[-1].FixParameter(1,40)
#StBkgs[-1].FixParameter(2,9.6e-06)
#StBkgs[-1].FixParameter(3,1.54e-07)
#StBkgs[-1].SetParLimits(0,1.e-03,1.e01)
#StBkgs[-1].SetParLimits(1,1.e-08,9)
#StBkgs[-1].SetParLimits(2,0.,9e+05)
#StBkgs[-1].SetParLimits(3,0.,9e+01)
#StBkgs[-1].SetParLimits(0,1.e04,9.e05)
#StBkgs[-1].SetParLimits(1,1.e01,9e+01)
#StBkgs[-1].SetParLimits(2,1.e-06,9.e-04)
#StBkgs[-1].SetParLimits(3,1.e-11,9.e-07)
#StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt3","[0]/TMath::Exp([1]*x/13000. + [2]*pow(x,3.))",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000. + [2]*pow(x,3.))",xMin,xMax) )
#StBkgs[-1].SetParameter(0,58762.6)
#StBkgs[-1].SetParameter(1,39.0475)
#StBkgs[-1].SetParameter(2,29.1897)
#StBkgs[-1].SetParameter(0,55080.1)
#StBkgs[-1].SetParameter(1,38.2665)
#StBkgs[-1].SetParameter(2,38.2665)

i = 0
for StBkg in StBkgs:
	xMinFit = xMin
	xMaxFit = xMax
	#xMinFit = 1250.
	#xMaxFit = 3750.
	#xMaxFit = 2500.
	if i == 0:
		pass
		#xMinFit = 1750. 
		#xMaxFit = 3200.
		#xMinFit = 1200. 
		#xMaxFit = 2000.
		xMaxFit = 2800.
	if i == 2:
		pass
		#xMinFit = 1250. 
		xMaxFit = 2800.
	#status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0N","",xMinFit,xMaxFit) )
	status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0NEI","",xMinFit,xMaxFit) )
	print status
	if status != 0 and status != 4000:
		print "WARNING: fit failed with status:",status
		#continue
	print "chiSq / ndf =",StBkg.GetChisquare(),"/",StBkg.GetNDF()
	print "====================="
	i += 1

c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetLogy()

hST[0].SetTitle("1#gamma, "+hST[0].GetName())
hST[0].GetYaxis().SetRangeUser(1.1e-01,9.e+03)
hST[0].GetXaxis().SetTitle("S_{T} [GeV]")
hST[0].GetXaxis().SetTitleOffset(1.1)
hST[0].SetMarkerStyle(20)
hST[0].SetMarkerSize(.9)
c.cd()
hST[0].Draw("E")

iPar=0
j = 0
for St in StBkgs:
	#St.SetParameter(iPar,float(St.GetParameter(iPar))/scales[i])
	St.SetLineWidth(2)
	St.SetLineColor(j+2)
	c.cd()
	St.Draw("SAME")
	#c.Update()
	j += 1

#_____ Call main() ______#
#if __name__ == '__main__':
#	  main()
