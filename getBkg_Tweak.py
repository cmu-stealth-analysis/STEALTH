#!/usr/bin/python

import ROOT
import argparse

parser = argparse.ArgumentParser(description='Bkg plotting options.')
parser.add_argument('-j','--ijtbkg', default=0, help='0: 2jt, 1:3jt, ...', type=int)
args = parser.parse_args()

# x-axis range for ST plot
#xMin = 900.
xMin = 1200.
xMax = 3500.

#xMin = 1250.
xMin = 1300.
xMax = 3700.
xMax = 4100.
# nJet distributions to plot
nJtMin = 2 
#nJtMax = 3 
#nJtMax = 5 
nJtMax = 6 
# jet index for normalization/ratio comparison
#iJtScale = 2 
# jet index for determining bkg scaling
iJtBkg = args.ijtbkg 
# y-axis range
yMin = 1.1e-02
yMin = 1.1e-00
yMin = 1.1e-06
#yMax = yMin*9.e+03
yMax = yMin*9.e+04
#yMax = yMin*9.e+05

## MAIN ##

hST = []
hFile = ROOT.TFile("hSTs.root","READ")
for j in range(nJtMin,nJtMax+1):
    hST.append( ROOT.gDirectory.Get("h"+str(j)+"jet") )
    #hST[-1].Sumw2()

# Estimate bkgrnd analytically
StBkgs = []
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Power(x/13000,TMath::Log(x/13000.)*[3]))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Log(x/13000.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000. + [2]*TMath::Log(x/13000.)))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]*(1./TMath::Exp([1]*x/13000.) + [2]/TMath::Power(x/13000,[3]))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x/13000.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","sqrt([0]*[0]/TMath::Power(x/13000,2*[1]) + [2]*[2]/TMath::Exp(2*[3]*x/13000.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*pow(x,3.))",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000,[1]) + [2]/TMath::Exp([3]*x/13000.)",xMin,xMax) )
#StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000. + [2]*pow(x,2.))",xMin,xMax) )

StBkgs.append( ROOT.TF1("fSt0","[0]/TMath::Power(x/13000.,[1])",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt1","[0]/TMath::Power(x/13000.,[1]*TMath::Log(x))",xMin,xMax) )
StBkgs.append( ROOT.TF1("fSt2","[0]/TMath::Exp([1]*x/13000.)",xMin,xMax) )

for i,StBkg in enumerate(StBkgs):
  #xMinFit = xMin_
  #xMinFit = 1300
  xMinFit = xMin
  xMaxFit = xMax
  if i == 1 and (iJtBkg == 4):
    pass
    xMaxFit = 3500
  '''
  if i == 0 and (iJtBkg == 1 or iJtBkg == 2):
    pass
    xMaxFit = 3500
  if i == 0 and (iJtBkg == 2):
    pass
    xMaxFit = 3300
  if i == 0 and (iJtBkg == 4):
    pass
    #xMaxFit = 3400
  #if i == 1 and (iJtBkg == 0 or iJtBkg == 2 or iJtBkg >= 4):
  if i == 1 and (iJtBkg == 0 or iJtBkg == 1 or iJtBkg == 2):
    pass
    xMaxFit = 2900
  if i == 1 and (iJtBkg == 3):
    pass
    xMaxFit = 3300
  if i == 1 and (iJtBkg == 4):
    pass
    #xMaxFit = 3475
    #xMinFit = 1400 # Data B, ijt=0,2,4,5
    #xMinFit = 1200 # Data B, ijt=0,2,4,5
  if i == 2: #lower
    pass
  '''
  #status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0N","",xMinFit,xMaxFit) )
  status = int( hST[iJtBkg].Fit("fSt"+str(i),"M0NEI","",xMinFit,xMaxFit) )
  print status
  if status != 0 and status != 4000:
    print "WARNING: fit failed with status:",status
    #continue
  print "chiSq / ndf =",StBkg.GetChisquare(),"/",StBkg.GetNDF()
  print "====================="

inm = 1
with open('st_scaling_params.dat','a+') as f:
  f.write("%E %E\n"%(StBkgs[inm].GetParameter(0),StBkgs[inm].GetParameter(1)))

#hFits = ROOT.TFile("hFits.root","RECREATE")
#for st in StBkgs:
#	st.Write()
#hFits.Close()

c = ROOT.TCanvas("c","c",600,600)
c.SetBorderSize(0);
c.SetFrameBorderMode(0)
ROOT.gStyle.SetTitleBorderSize(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetLogy()

hST[iJtBkg].SetTitle("2#gamma, "+hST[iJtBkg].GetName())
hST[iJtBkg].GetYaxis().SetRangeUser(yMin,yMax)
hST[iJtBkg].GetXaxis().SetTitle("S_{T} [GeV]")
hST[iJtBkg].GetXaxis().SetTitleOffset(1.1)
hST[iJtBkg].GetXaxis().SetRangeUser(xMinFit,xMax)
hST[iJtBkg].SetMarkerStyle(20)
hST[iJtBkg].SetMarkerSize(.9)
#hST[iJtBkg].SetLineColor(0)
c.cd()
#hST[iJtBkg].Draw("HIST")
hST[iJtBkg].Draw("E")

iPar=0
j = 0
lcolors = [4,1,2]
for St in StBkgs:
  #St.SetParameter(iPar,float(St.GetParameter(iPar))/scales[i])
  St.SetLineWidth(2)
  St.SetLineColor(lcolors[j])
  c.cd()
  St.Draw("SAME")
  #c.Update()
  j += 1

# Draw legend
leg = ROOT.TLegend(0.72,0.65,0.86,0.85)
leg.AddEntry(hST[iJtBkg],"Bkg","LP")
i = 0
label = []
label.append("1/x^{p_{0}}")
label.append("1/x^{p_{1}lnS_{t}}")
label.append("1/e^{p_{2}x}")
for StBkg in StBkgs:
	leg.AddEntry(StBkg,label[i],"LP")
	i += 1
leg.SetBorderSize(0);
leg.Draw("SAME")
c.Update()

c.Print("bkgfit.png")

#_____ Call main() ______#
#if __name__ == '__main__':
#	  main()
