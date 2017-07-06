#!/usr/bin/python

import os
import sys
import ROOT
#import argparse

# x-axis range for ST plot
#xMin = 1000.
#xMin = 900.
xMin = 1200.
xMax = 3500.
xMin = 1300.
xMax = 3700.
# nJet distributions to plot
nJtMin = 2
nJtMax = 3
#nJtMax = 5
nJtMax = 6
# jet index for normalization/ratio comparison
iJtScale = 0
# jet index for determining bkg scaling
iJtBkg = iJtScale
# output directory for plots
#outDir = 'BKG/DoubleEG_ReminiAOD'
outDir = 'BKG/JetHT_ReminiAOD'
#outDir = 'BKG/DiPhoton_selA'
#outDir = "BKG/GJet_selA/"
#outDir = "BKG/SinglePhoton/"
# y-axis range
yMin = 1.1
#yMin = 1.1e-02
#yMin = 1.1e-06
#yMax = yMin*9.e+03
yMax = yMin*9.e+04
#yMax = yMin*9.e+05

## MAIN ##
def main():

    hST = []
    #hFile = ROOT.TFile("hFile.root","READ")
    hFile = ROOT.TFile("hSTs.root","READ")
    for j in range(nJtMin,nJtMax+1):
        hST.append( ROOT.gDirectory.Get("h"+str(j)+"jet") )

    # Estimate bkgrnd analytically
    StBkgs = []
    # old scheme
    #StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]) + [2]/TMath::Exp([3]*x/13000.)",xMin,xMax) )
    #StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )
    #StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000. + [2]*pow(x,3.))",xMin,xMax) )
    # new scheme
    StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
    StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x))",xMin,xMax) )
    StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )

    i = 0
    for StBkg in StBkgs:
        #xMinFit = 1300
        xMinFit = xMin
        xMaxFit = xMax
        if i == 0:
            pass
        #xMinFit = 1300. # Data B, ijt=0
        if i == 1 and iJtBkg == 0:
            pass
        #xMaxFit = 2900
        if i == 2:
            pass
        #status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0N","",xMinFit,xMaxFit) )
        status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0NEI","",xMinFit,xMaxFit) )

        print status
        if status != 0 and status != 4000:
            print "WARNING: fit failed with status:",status
          #continue
        print "chiSq / ndf =",StBkg.GetChisquare(),"/",StBkg.GetNDF()
        print "====================="
        i += 1

    with open('fit_choice_params.dat','w+') as f:
        for st in StBkgs:
            f.write('%E %E %E\n'%(st.GetParameter(0),st.GetParameter(1),st.GetParError(1)))

    fScales = 'normRatios.txt'
    with open(fScales, "r") as f:
        scales = f.readlines()
    i = 0
    for h in hST:
        print "N -> N/%f"%(float(scales[i]))
        plotHistovFit(h,StBkgs,float(scales[i]))
        os.rename("bkgfit.png","%s/bkgFit_%djet.png"%(outDir,i+2))
        i += 1

def plotHistovFit(hST,gFits,scale):

    c = ROOT.TCanvas("c","c",600,600)
    c.SetBorderSize(0);
    c.SetFrameBorderMode(0)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogy()
    c.cd()

    hST.SetTitle("2#gamma, "+hST.GetName())
    hST.GetYaxis().SetRangeUser(yMin,yMax)
    hST.GetXaxis().SetTitle("S_{T} [GeV]")
    hST.GetXaxis().SetTitleOffset(1.1)
    hST.GetXaxis().SetRangeUser(xMin,xMax)
    hST.SetMarkerStyle(20)
    hST.SetMarkerSize(.9)
    hST.Draw("E")
    #c.Update()

    leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
    leg.AddEntry(hST,"Bkg","LP")

    iPar=0
    i = 0
    StBkgs = []
    lcolor = [4,1,2]
    for gFit in gFits:
        StBkgs.append(gFit.Clone())
        StBkgs[-1].SetParameter(iPar,float(StBkgs[-1].GetParameter(iPar))/scale)
        #if i == 0:
        #	pass
        #	StBkgs[-1].SetParameter(iPar+2,float(StBkgs[-1].GetParameter(iPar+2))/scale)
        StBkgs[-1].SetLineWidth(2)
        StBkgs[-1].SetLineColor(lcolor[i])
        #StBkgs[-1].SetLineColor(i+1)
        StBkgs[-1].Draw("SAME")
        c.Update()
        i += 1

    # Draw legend
    i = 0 
    label = []
    label.append("1/x^{p_{0}}")
    label.append("1/x^{p_{1}lnS_{t}}")
    label.append("1/e^{p_{2}x}")
    # old scheme
    #label.append("1/e^{p_{1}x} #oplus 1/x^{p_{2}}")
    #label.append("1/e^{p_{3}x}")
    #label.append("1/e^{p_{4}x + p_{5}S_{T}^{3}}")
    for StBkg in StBkgs:
        leg.AddEntry(StBkg,label[i],"LP")
        i += 1
    leg.SetBorderSize(0);
    leg.Draw("SAME")
    c.Update()

    c.Print("bkgfit.png")

#_____ Call main() ______#
if __name__ == '__main__':
    main()
