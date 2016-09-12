#!/usr/bin/python

import os
import sys
import array
import numpy as np
import ROOT
import argparse
 
parser = argparse.ArgumentParser(description='Provide file to process')
parser.add_argument('-i','--input', help='Input file name',required=True)
args = parser.parse_args()

ggIn = ROOT.TChain("ggNtuplizer/EventTree")
ggIn.Add(args.input+".root")

nEntries = ggIn.GetEntries()
print "nEntries="+str(nEntries)
nBins = 5
xMin = 750.
xMax = 3750.
nJtMin = 2
nJtMax = 7
iJtScale = 0
iJtBkg = 0

hST = []
for j in range(nJtMin,nJtMax+1):
	hST.append( ROOT.TH1F(str(j)+"jet",str(j)+"jet",nBins,xMin,xMax) )

for h in hST:
	h.Sumw2()

for jEvt in range(nEntries):

	# Initialize event
	iEvt = ggIn.LoadTree(jEvt)
	if iEvt < 0:
		break
	nb = ggIn.GetEntry(jEvt)
	if nb <= 0:
		continue
	if jEvt % 10000 == 0:
		print "Processing entry " + str(jEvt)

	#nJets = ggIn.b_nJets
	#evtSt = ggIn.b_evtSt
	nJets = ggIn.b_nJet
	evtSt = ggIn.b_evtST

	# Fill St for each jet multiplicity
	for iJet in range(0,nJtMax-nJtMin+1):
		if nJets >= nJtMax-iJet:
			hST[nJtMax-nJtMin-iJet].Fill(evtSt)
			break

for iJet in range(0,nJtMax-nJtMin+1):
	print str(iJet+2)+"-jet evts: "+str(hST[iJet].Integral())

## ==== DRAW PLOTS ==== ##

c1 = ROOT.TCanvas("c1","c1",600,600)
c1.SetBorderSize(0);
c1.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)

pUp = ROOT.TPad("upperPad", "upperPad",.005, .270, .995, .995)
pDn = ROOT.TPad("lowerPad", "lowerPad",.005, .005, .995, .270)

pUp.Draw()
pDn.Draw()

pUp.SetMargin(12.e-02,3.e-02,5.e-03,2.e-02)
pDn.SetMargin(12.e-02,3.e-02,29.e-02,4.e-02)

## Draw upper pad

pUp.cd()
pUp.SetTicky()

ROOT.gPad.SetLogy()

# Get normalization
norm = []
nStLo = []
print "==============="
for h in hST:
	norm.append(0.)
	for iBin in range (1,nBins+1):
		norm[-1] += h.GetBinContent(iBin)
		print "sT="+str( h.GetXaxis().GetBinCenter(iBin) )+" : N="+str(h.GetBinContent(iBin))
		if iBin == 1:
			nStLo.append(h.GetBinContent(iBin))
	print "==============="

# Get background estimates
iJt = 0
jtScale = []
print "==============="
for h in hST:
	jtScale.append(nStLo[iJtBkg] / nStLo[iJt])
	print "scale factor: "+str(nStLo[iJtBkg])+" / "+str(nStLo[iJt])+" = "+str(jtScale[-1])
	for iBin in range (1,nBins+1):
		print "sT="+str( h.GetXaxis().GetBinCenter(iBin) )+" : N="+str(h.GetBinContent(iBin))+" -> "+str(h.GetBinContent(iBin)*jtScale[-1])
	print "==============="
	iJt += 1

# Normalize distributions 
maxSTs = []
iJt = 0
for h in hST:
	h.Scale(1./norm[iJt])
	maxSTs.append(h.GetMaximum())
	iJt += 1

'''
print "==============="
for h in hST:
	for iBin in range (1,nBins):
		print "sT="+str( h.GetXaxis().GetBinCenter(iBin) )+" : N="+str(h.GetBinContent(iBin))
	print "==============="
'''

#hST[0].GetYaxis().SetRangeUser(2.e-04,1.4*max(maxSTs))
hST[0].GetYaxis().SetRangeUser(2.e-04,1.)
hST[0].GetXaxis().SetTitle("S_{T} [GeV]")
hST[0].GetXaxis().SetTitleOffset(1.1)
#hST[0].SetTitle("1#gamma, Ht > 700 GeV")
hST[0].SetTitle("")


hST[0].SetLineColor(1)
hST[0].Draw("E")
#hST[0].Draw("")
count = 0
for h in hST:
	if count > 0:
		h.SetLineColor(count+1)
		h.Draw("SAME E")
		c1.Update()
		#hST[i].Draw("SAME")
	count += 1

# Draw legend
count = 0
label = ''
leg = ROOT.TLegend(0.75,0.65,0.86,0.92)
for h in hST:
	if count < nJtMax - nJtMin:
		label = str(count+2)+" jets"
	else:
		label = "#geq"+str(count+2)+" jets"
	#label = str(count+2)+" jets"
	leg.AddEntry(h,label,"LP")
	count += 1
#ROOT.gStyle.SetBorderSize(0);
leg.SetBorderSize(0);
leg.Draw()
#Draw labels
tex = ROOT.TLatex()
tex.DrawLatexNDC(0.19,0.16,"N="+str(nEntries));
tex.DrawLatexNDC(0.19,0.10,"1#gamma, Ht > 700 GeV");

## Draw lower pad

pDn.cd()
pDn.SetTicky()
pDn.SetGridy()

fUnity = ROOT.TF1("fUnity","[0]",xMin,xMax)
fUnity.SetParameter( 0,1. )
fUnity.GetYaxis().SetRangeUser(0.2,2.8)

fUnity.GetXaxis().SetTitle("S_{T} [GeV]")
fUnity.GetXaxis().SetLabelSize(0.1)
fUnity.GetXaxis().SetTitleSize(0.12)
fUnity.GetXaxis().SetTickLength(0.1)
fUnity.GetXaxis().SetTitleOffset(1.14)
fUnity.GetYaxis().SetTitle("n_{j} / n_{j} = 3")
fUnity.GetYaxis().SetNdivisions(305)
fUnity.GetYaxis().SetLabelSize(0.1)
fUnity.GetYaxis().SetTitleSize(0.11)
fUnity.GetYaxis().SetTickLength(0.04)
fUnity.GetYaxis().SetTitleOffset(0.5)
fUnity.SetLineColor(2)
fUnity.SetLineWidth(1)
fUnity.SetLineStyle(4)
fUnity.SetTitle("")

fUnity.Draw()

gRatio = []
'''
iJt = 0
for h in hST:
	print "hST[" + str(iJt) + "]..." 
	for iBin in range (1,nBins+1):
		print "  >> iBin="+str(iBin)+": "+str(h.GetBinContent(iBin))
	iJt += 1
'''
# Fill ratio and errors
iJt = 0
for h in hST:
	if iJt == iJtScale:
		iJt += 1
		continue
	gRatio.append(ROOT.TGraphErrors())
	print "hST[" + str(iJt) + "]..." 
	for iBin in range (1,nBins+1):
		if not (hST[1].GetBinContent(iBin) > 0.):
			continue
		ratioX = h.GetXaxis().GetBinCenter(iBin)
		ratioY = h.GetBinContent(iBin)/hST[iJtScale].GetBinContent(iBin) 
		gRatio[-1].SetPoint( iBin, ratioX, ratioY )
		if not (h.GetBinContent(iBin) > 0.):
			continue
		errY = ratioY*np.sqrt( (h.GetBinError(iBin)/h.GetBinContent(iBin))**2 + (hST[iJtScale].GetBinError(iBin)/hST[1].GetBinContent(iBin))**2 )
		gRatio[-1].SetPointError( iBin, (xMax-xMin)/(2.*nBins), errY )
		#print " >>"+str(h.GetBinContent(iBin))+" "+str(hST[1].GetBinContent(iBin))
	gRatio[-1].SetMarkerStyle(20)
	gRatio[-1].SetMarkerColor(iJt+1)
	gRatio[-1].SetLineColor(iJt+1)
	gRatio[-1].SetLineWidth(1)
	gRatio[-1].SetMarkerSize(.9)
	gRatio[-1].Draw("P SAME")
	c1.Update()
	iJt += 1

c1.Print("sT.png")
