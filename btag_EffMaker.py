
# coding: utf-8

# In[1]:


#!/usr/bin/env python
import ROOT as ROOT
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TFile, TGraphAsymmErrors,TLatex,TLine,gStyle,TLegend,TH2D
import os, datetime
import random
import sys, optparse
import array
import math
import getpass
import socket
import json
datestr = datetime.datetime.now().strftime("%I%p%Y%m%d")


# In[2]:


inputFilename = 'BtagEff_Output.root'
inputFile = TFile(inputFilename, 'READ')
outputFilename = 'btag_eff_forweight.root'
outputFile = TFile(outputFilename, 'RECREATE')

for partonFlavor in ['b','c','light']:
    denominatorIn = inputFile.Get('h_'+partonFlavor+'eff_den')
    numeratorIn = inputFile.Get('h_'+partonFlavor+'eff_num_pass')

    xShift = denominatorIn.GetXaxis().GetBinWidth(1)/2.
    yShift = denominatorIn.GetYaxis().GetBinWidth(1)/2.

    denominatorOut = TH2D('denominator_' + partonFlavor, '', 10,-2.4,2.4,20,20.,1000.)
    numeratorOut   = TH2D('numerator_' + partonFlavor, '', 10,-2.4,2.4,20,20.,1000.)
    efficiencyOut  = TH2D('efficiency_' + partonFlavor, '', 10,-2.4,2.4,20,20.,1000.)

    for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
      for j in range(1,denominatorOut.GetYaxis().GetNbins()+1):
        binXMin = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinLowEdge(i)+xShift)
        binXMax = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinUpEdge(i)-xShift)
        binYMinPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinLowEdge(j)+yShift)
        binYMaxPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinUpEdge(j)-yShift)
        binYMinNeg = denominatorIn.GetYaxis().FindBin(-denominatorOut.GetYaxis().GetBinUpEdge(j)+yShift)
        binYMaxNeg = denominatorIn.GetYaxis().FindBin(-denominatorOut.GetYaxis().GetBinLowEdge(j)-yShift)

        denominator = denominatorIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
        denominator = denominator + denominatorIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)
        numerator = numeratorIn.Integral(binXMin,binXMax,binYMinPos,binYMaxPos)
        numerator = numerator + numeratorIn.Integral(binXMin,binXMax,binYMinNeg,binYMaxNeg)

        if(i==denominatorOut.GetXaxis().GetNbins()): # also add overflow to the last bin in jet pT
          denominator = denominator + denominatorIn.Integral(binXMax+1,denominatorIn.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
          denominator = denominator + denominatorIn.Integral(binXMax+1,denominatorIn.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)
          numerator = numerator + numeratorIn.Integral(binXMax+1,numeratorIn.GetXaxis().GetNbins()+1,binYMinPos,binYMaxPos)
          numerator = numerator + numeratorIn.Integral(binXMax+1,numeratorIn.GetXaxis().GetNbins()+1,binYMinNeg,binYMaxNeg)

        denominatorOut.SetBinContent(i,j,denominator)
        numeratorOut.SetBinContent(i,j,numerator)
        if(denominator>0.): efficiencyOut.SetBinContent(i,j,numerator/denominator)

    # check if there are any bins with 0 or 100% efficiency
    for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
      for j in range(1,denominatorOut.GetYaxis().GetNbins()+1):

        efficiency = efficiencyOut.GetBinContent(i,j)
        if(efficiency==0. or efficiency==1.):
          print ('Warning! Bin(%i,%i) for %s jets has a b-tagging efficiency of %.3f'%(i,j,partonFlavor,efficiency))

    # set efficiencies in overflow bins
    for i in range(1,denominatorOut.GetXaxis().GetNbins()+1):
      efficiencyOut.SetBinContent(i, denominatorOut.GetYaxis().GetNbins()+1, efficiencyOut.GetBinContent(i, denominatorOut.GetYaxis().GetNbins()))

    for j in range(1,denominatorOut.GetYaxis().GetNbins()+2):
      efficiencyOut.SetBinContent(denominatorOut.GetXaxis().GetNbins()+1, j, efficiencyOut.GetBinContent(denominatorOut.GetXaxis().GetNbins(), j))

    outputFile.cd()
    denominatorOut.Write()
    numeratorOut.Write()
    efficiencyOut.Write()
outputFile.Close()

print ('-------------------------------------------------------------------------------------------')
print ('successfully created and stored in %s'%outputFilename)
print ('')
