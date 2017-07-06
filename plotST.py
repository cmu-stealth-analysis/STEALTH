#!/usr/bin/python

import os
import ROOT
import argparse
import re

# Register command line options
parser = argparse.ArgumentParser(description='ST processing options.')
parser.add_argument('-i','--input', nargs='+', help='Input file/s.',required=True, type=str)
parser.add_argument('-l','--stmin', default=1000., help='Min ST to plot.',type=float)
parser.add_argument('-r','--stmax', default=3500., help='Max ST to plot.',type=float)
parser.add_argument('-b','--nbins', default=5, help='Number of bins over which to plot ST.',type=int)
parser.add_argument('-H','--HTcut', default=60, help='HT cut.',type=int)
parser.add_argument('-n','--iNorm', default=1, help='HT cut.',type=int)
args = parser.parse_args()

# Set plotting parameters
# nBins for ST plot
nBins = args.nbins
# x-axis range for ST plot
xMin = args.stmin
xMax = args.stmax
# nJet distributions to plot
nJtMin = 2
#nJtMax = 3
#nJtMax = 5
#nJtMax = 7
nJtMax = 6
# jet index used as denominator for ratio plots (0:2jt, 1:3jt,...)
iJetRatio = 0 
# jet index used as control for bkg normalization (0:2jt, 1:3jt,...)
iJetBkg = iJetRatio
# histogram bin range [iBinBkgLo,iBinBkgHi+1) used as control for bkg normalization (0:underflow, 1:xMin included, ...)
#iBinBkgLo = 1
iBinBkgLo = args.iNorm
#iBinBkgLo = 7
iBinBkgHi = iBinBkgLo
# nPho and HT (only used for labels)
nPho  = 1
evtHT = args.HTcut 
if args.HTcut < 100.:
    nJtMax = 6 
    #nJtMax = 3 
    nPho = 2 
print " >> Plotting ST range: [",xMin,"->",xMax,"), in nBins:",nBins
print " >> nJets:",nJtMin,"->",nJtMax
print " >> Control sample:",iJetBkg+2,"jets, from ST bins: [",iBinBkgLo,"->",iBinBkgHi+1,")"
print " >> Denominator in ratio plots:",iJetRatio+2,"jets"

# Load input files
ggIn = ROOT.TChain("ggNtuplizer/EventTree")
for infile in args.input:
    infile = re.sub('[,\[\]]','',infile)
    infile += ".root"
    print " >> Adding input file:",infile
    ggIn.Add(infile)
nEvts = ggIn.GetEntries()
#nEvts = 500000
print " >> nEvts:",nEvts

## MAIN ##
def main():

    # Initialize nJet histograms
    hST = []
    for j in range(nJtMin,nJtMax+1):
        hST.append( ROOT.TH1F("h"+str(j)+"jet",str(j)+"jet",nBins,xMin,xMax) )
        hST[-1].Sumw2()

    #_____ BIN DATA BY ST,NJETS _____#

    # Loop over entries
    nAcc = 0
    for jEvt in range(nEvts):

        # Initialize event
        if jEvt > nEvts:
            break
        treeStatus = ggIn.LoadTree(jEvt)
        if treeStatus < 0:
            break
        evtStatus = ggIn.GetEntry(jEvt)
        if evtStatus <= 0:
            continue
        if jEvt % 1000 == 0:
            print " .. Processing entry",jEvt

        # Load nJet,ST branches
        nJets = ggIn.b_nJets
        #nJets = ggIn.b_nJet
        evtST = ggIn.b_evtST
        evtWgt = 1. 
        if not ggIn.isData:
            evtWgt = ggIn.b_evtWgt_1_pb
        if nJets < nJtMin:
            continue

        # Fill the appropriate jet distn
        nJetBins = nJtMax-nJtMin
        for iJet in range(nJetBins+1):
            # start from highest jet mutiplicity and iterate downward
            if nJets >= nJtMax-iJet: # inclue >= nJetMax
                #if nJets == nJtMax-iJet: # stop exactly nJetMax
                hST[nJetBins-iJet].Fill(evtST,evtWgt)
                #if evtST > xMin and evtST < xMax and evtST > 1000.:
                if evtST > xMin and evtST < xMax:
                    nAcc += 1
                break

    # Print out integrals for each jet distn
    #nAcc = 0
    print " >> Histogram integrals:"
    for iJet in range(nJtMax-nJtMin+1):
        #nAcc += hST[iJet].Integral()
        print " .. "+str(iJet+2)+"-jets: "+str(hST[iJet].Integral())

    #_____ OUTPUT FOR BKG ESTIMATION _____#

    # Write out histos
    hFile = ROOT.TFile("hSTs.root","RECREATE")
    for h in hST:
        h.Write()
    hFile.Close()

    # Write out normalization ratios
    bkgNorm = []
    # Get nEvts in control bin/s
    for h in hST:
        bkgNorm_ = 0
        for iBin in range(iBinBkgLo,iBinBkgHi+1):
            bkgNorm_ += h.GetBinContent(iBin)
        bkgNorm.append(bkgNorm_)
    # Write ratios to file
    iJet = 0
    ratioFile = open("normRatios.txt", "w")
    print " >> Debugging bin entries:"
    print " >> ==============="
    for h in hST:
        print "bkg norm:", bkgNorm[iJet]
        if bkgNorm[iJet] > 0:
            ratio = bkgNorm[iJetBkg]/bkgNorm[iJet]
            ratioFile.write( "%f\n" % ratio )
            print " .. ratio of norms: "+str(bkgNorm[iJetBkg])+" / "+str(bkgNorm[iJet])+" = "+str(ratio)
            for iBin in range (1,nBins+1):
                print " .. ST="+str( h.GetXaxis().GetBinLowEdge(iBin) )+": N="+str(h.GetBinContent(iBin))+" -> "+str(h.GetBinContent(iBin)*ratio)
        print " >> ==============="
        iJet += 1
    ratioFile.close()

    # Renormalize iJet histograms 
    iJet = 0
    maxSTs = []
    for h in hST:
        if bkgNorm[iJet] > 0.:
            h.Scale(1./bkgNorm[iJet])
            #h.Scale(1./h.Integral())
            maxSTs.append(h.GetMaximum())
        iJet += 1

    #______ DRAWING PLOTS _____#

    # Initialize canvas and pads
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

    ##### ST histos on upper pad #####
    pUp.cd()
    pUp.SetTicky()
    ROOT.gPad.SetLogy()

    # Draw histos
    linecolors = [1,2,8,4,49,5]
    hST[0].GetYaxis().SetRangeUser(2.e-04,9.)
    #hST[0].GetYaxis().SetRangeUser(2.e-03,9.e01)
    #hST[0].GetYaxis().SetRangeUser(2.e-03,9.e02)
    if args.HTcut < 100.:
        pass
    #hST[0].GetYaxis().SetRangeUser(2.e-04,9.e02)
    hST[0].GetXaxis().SetTitle("S_{T} [GeV]")
    hST[0].GetXaxis().SetTitleOffset(1.1)
    hST[0].SetTitle("")
    hST[0].SetLineColor(1)
    hST[0].Draw("E")
    iJet = 0
    for h in hST:
        if iJet > 0:
            h.SetLineColor(linecolors[iJet])
            h.Draw("SAME E")
            c1.Update()
        iJet += 1

    # Draw legend
    iJet = 0
    label = ''
    leg = ROOT.TLegend(0.75,0.65,0.86,0.92)
    for h in hST:
        if iJet < nJtMax - nJtMin:
            label = str(iJet+2)+" jets"
        else:
            #label = str(iJet+2)+" jets"
            label = "#geq"+str(iJet+2)+" jets"
        #label = str(iJet+2)+" jets"
        leg.AddEntry(h,label,"LP")
        iJet += 1
    leg.SetBorderSize(0);
    leg.Draw()
    #Draw labels
    tex = ROOT.TLatex()
    #tex.DrawLatexNDC(0.19,0.16,"N="+str(nEvts));
    tex.DrawLatexNDC(0.19,0.16,"N=%d"%(nAcc));
    tex.DrawLatexNDC(0.19,0.10,str(nPho)+"#gamma, HT > "+str(evtHT)+" GeV");

    ##### Ratio plots on lower pad #####
    pDn.cd()
    pDn.SetTicky()
    pDn.SetGridy()

    # Draw r = 1 and axis labels
    fUnity = ROOT.TF1("fUnity","[0]",xMin,xMax)
    fUnity.SetParameter( 0,1. )
    fUnity.GetYaxis().SetRangeUser(0.2,2.8)
    fUnity.GetXaxis().SetTitle("S_{T} [GeV]")
    fUnity.GetXaxis().SetLabelSize(0.1)
    fUnity.GetXaxis().SetTitleSize(0.12)
    fUnity.GetXaxis().SetTickLength(0.1)
    fUnity.GetXaxis().SetTitleOffset(1.14)
    fUnity.GetYaxis().SetTitle("n_{j} / n_{j} = "+str(iJetRatio+2))
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

    # Draw ratios and errors
    iJet = 0
    gRatio = []
    for h in hST:
        if iJet == iJetRatio:
            iJet += 1
            continue
        gRatio.append( make_ratio_graph("g"+str(iJet+2)+"jt_"+str(iJetRatio+2)+"jt", h, hST[iJetRatio]) )
        gRatio[-1].SetMarkerStyle(20)
        gRatio[-1].SetMarkerColor(linecolors[iJet])
        gRatio[-1].SetLineColor(linecolors[iJet])
        gRatio[-1].SetLineWidth(1)
        gRatio[-1].SetMarkerSize(.9)
        gRatio[-1].Draw("P SAME")
        c1.Update()
        iJet += 1

    c1.Print("sT.png")
    c1.Print("sT.eps")

#______ Draw Copper-Pearson (asymmetric) errors _____#
# Contributed by Marc W.
def make_ratio_graph(g_name, h_num, h_den):
    print " >> Getting ratio errors for:",g_name
    gae = ROOT.TGraphAsymmErrors()
    gae.SetName(g_name)
    h_rat = h_num.Clone()
    h_rat.Divide(h_den)
    for i in range(1, h_rat.GetNbinsX() + 1):
        n_rat = 0
        r_high = 0

        # tail = (1 - cl) / 2; for 95% CL, tail = (1 - 0.95) / 2 = 0.025
        tail = 0.16
        if h_num.GetBinError(i) == 0. or h_den.GetBinError(i) == 0.:
            continue
        n_num = pow(h_num.GetBinContent(i) / h_num.GetBinError(i), 2)
        n_den = pow(h_den.GetBinContent(i) / h_den.GetBinError(i), 2)
        q_low = ROOT.Math.fdistribution_quantile_c(1 - tail, n_num * 2,
                (n_den + 1) * 2)
        r_low = q_low * n_num / (n_den + 1)
        if n_den > 0:
            n_rat = n_num / n_den
            q_high = ROOT.Math.fdistribution_quantile_c(tail, (n_num + 1) * 2,
                    n_den * 2)
            r_high = q_high * (n_num + 1) / n_den
            gae.SetPoint(i - 1, h_rat.GetBinCenter(i), h_rat.GetBinContent(i))
            gae.SetPointError(
                    i - 1, h_rat.GetBinWidth(i) / 2, h_rat.GetBinWidth(i) / 2,
                    n_rat - r_low, r_high - n_rat
                    )
            gae.GetXaxis().SetRangeUser(h_rat.GetBinLowEdge(1),
                    h_rat.GetBinLowEdge(h_rat.GetNbinsX()) +
                    h_rat.GetBinWidth(h_rat.GetNbinsX()))
            fUnity = ROOT.TF1("fUnity","[0]",xMin,xMax)
    fUnity.SetParameter( 0,1. )
    gae.Fit("fUnity","M0","",xMin,xMax)
    print " >> chiSq / ndf =",fUnity.GetChisquare(),"/",fUnity.GetNDF()
    print " >> ====================="
    return gae

#_____ Call main() ______#
if __name__ == '__main__':
    main()
