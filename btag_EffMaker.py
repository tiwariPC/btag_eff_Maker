
#!/usr/bin/env python
import ROOT as ROOT
from ROOT import TCanvas, TColor, TGaxis, TH1F, TPad, TFile, TGraphAsymmErrors,TLatex,TLine,gStyle,TLegend,TH2D
import os, datetime
import random
import sys, optparse
from array import array
import math
import getpass
import socket
import json
datestr = datetime.datetime.now().strftime("%I%p%Y%m%d")

## ----- command line argument
usage = "analyzer for bb+DM (debugging) "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-inDir", "--inputDir",  dest="inputDir", default=".")
parser.add_argument("-D", "--outputdir", dest="outputdir")
parser.add_argument("-T", "--testing", action="store_true",  dest="testing")

args = parser.parse_args()

if args.inputDir:
    dirName = args.inputDir
else:
    print('Please provide an Input Directory')
    sys.exit()

if args.outputdir:
    outputdir = str(args.outputdir)
else:
    print('Please provide an Output Directory')
    sys.exit()

if args.testing == None:
    istest = False
else:
    istest = args.testing

def runbbdm(txtfile):
  inputFile = TFile(txtfile, 'READ')
  outputFilename = 'Output_'txtfile
  outputFile = TFile(outputFilename, 'RECREATE')
  outputFile.cd()
  for partonFlavor in ['btag', 'ctag', 'lighttag']:
      denominatorIn = inputFile.Get('h_'+partonFlavor+'_den')
      for wp in ['lwp', 'mwp']:
          numeratorIn = inputFile.Get('h_'+partonFlavor+'_num_pass_'+wp)
          print('wp', wp)

          xShift = denominatorIn.GetXaxis().GetBinWidth(1)/2.
          yShift = denominatorIn.GetYaxis().GetBinWidth(1)/2.
          bins_pT = [20.0, 50.0, 80.0, 120.0, 200.0, 300.0, 400.0, 500.0, 700.0, 1000.0]
          bins_eta = [-2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5]

          denominatorOut = TH2D('denominator_' + partonFlavor+'_'+wp, '', 6, array('d', bins_eta), 9, array('d', bins_pT))
          numeratorOut = TH2D('numerator_' + partonFlavor+'_'+wp, '', 6, array('d', bins_eta), 9, array('d', bins_pT))
          efficiencyOut = TH2D('efficiency_' + partonFlavor+'_'+wp, '', 6, array('d', bins_eta), 9, array('d', bins_pT))

          for i in range(1, denominatorOut.GetXaxis().GetNbins()+1):
            for j in range(1, denominatorOut.GetYaxis().GetNbins()+1):
              binXMin = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinLowEdge(i)+xShift)
              binXMax = denominatorIn.GetXaxis().FindBin(denominatorOut.GetXaxis().GetBinUpEdge(i)-xShift)

              binYMinPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinLowEdge(j)+yShift)
              binYMaxPos = denominatorIn.GetYaxis().FindBin(denominatorOut.GetYaxis().GetBinUpEdge(j)-yShift)

              denominator = denominatorIn.Integral(binXMin, binXMax, binYMinPos, binYMaxPos)
              numerator = numeratorIn.Integral(binXMin, binXMax, binYMinPos, binYMaxPos)

              # also add overflow to the last bin in jet pT
              if(i == denominatorOut.GetXaxis().GetNbins()):
                denominator = denominator + denominatorIn.Integral(binXMax+1, denominatorIn.GetXaxis().GetNbins()+1, binYMinPos, binYMaxPos)
                numerator = numerator + numeratorIn.Integral(binXMax+1, numeratorIn.GetXaxis().GetNbins()+1, binYMinPos, binYMaxPos)

              denominatorOut.SetBinContent(i, j, denominator)
              numeratorOut.SetBinContent(i, j, numerator)
              if(denominator > 0.):
                efficiencyOut.SetBinContent(i, j, numerator/denominator)

          # check if there are any bins with 0 or 100% efficiency
          for i in range(1, denominatorOut.GetXaxis().GetNbins()+1):
            for j in range(1, denominatorOut.GetYaxis().GetNbins()+1):
              efficiency = efficiencyOut.GetBinContent(i, j)
              if(efficiency == 0. or efficiency == 1.):
                print('Warning! Bin(%i,%i) for %s jets has a b-tagging efficiency of %.3f' %(i, j, partonFlavor, efficiency))
          # set efficiencies in overflow bins
          for i in range(1, denominatorOut.GetXaxis().GetNbins()+1):
            efficiencyOut.SetBinContent(i, denominatorOut.GetYaxis().GetNbins()+1, efficiencyOut.GetBinContent(i, denominatorOut.GetYaxis().GetNbins()))

          for j in range(1, denominatorOut.GetYaxis().GetNbins()+2):
            efficiencyOut.SetBinContent(denominatorOut.GetXaxis().GetNbins()+1, j, efficiencyOut.GetBinContent(denominatorOut.GetXaxis().GetNbins(), j))

          print('writing eff : ', 'h_'+partonFlavor+'_num_pass_'+wp)
          efficiencyOut.Write()
  outputFile.Close()
  print('-------------------------------------------------------------------------------------------')
  print('successfully created and stored in %s' % outputFilename)
  print('')

if __name__ == '__main__':
  files = [f for f in os.listdir(infile) if f.endswith(".root")]
  n = mp.cpu_count()  
  final = [files[i * n:(i + 1) * n] for i in range((len(files) + n - 1) // n)]
  if istest:
    runbbdm, final[0]
  else:
    for i in range(len(final)):
      try:
        pool = mp.Pool(mp.cpu_count())
        pool.map(runbbdm, final[i])
        pool.close()
        pool.join()
      except Exception as e:
        print(e)
        print("Corrupt file inside input txt file is detected! Skipping this txt file:  ", final[i])
        continue
