#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain,TGraphAsymmErrors, TMath, TH2D, TLorentzVector, TF1, AddressOf
import ROOT as ROOT
import os
import random
import sys, optparse
from array import array
import math

ROOT.gROOT.SetBatch(True)

ROOT.gROOT.LoadMacro("Loader.h+")

usage = "usage: %prog [options] arg1 arg2"
parser = optparse.OptionParser(usage)


parser.add_option("-i", "--inputfile",  dest="inputfile")
parser.add_option("-o", "--outputfile", dest="outputfile")
parser.add_option("-D", "--outputdir", dest="outputdir")
parser.add_option("-F", "--farmout", action="store_true",  dest="farmout")

(options, args) = parser.parse_args()


if options.farmout==None:
    isfarmout = False
else:
    isfarmout = options.farmout


inputfilename = options.inputfile
outputdir = options.outputdir


pathlist = inputfilename.split("/")
sizeoflist = len(pathlist)
#print ('sizeoflist = ',sizeoflist)
rootfile='tmphist'
rootfile = pathlist[sizeoflist-1]
textfile = rootfile+".txt"

if outputdir!='.': os.system('mkdir -p '+outputdir)

if options.outputfile is None or options.outputfile==rootfile:
    if not isfarmout:
        outputfilename = "/Output_"+rootfile
    else:
        outputfilename = "/Output_"+rootfile.split('.')[0]+".root"
else:
    outputfilename = "/"+options.outputfile



outfilename = outputdir + outputfilename
#else:
#    outfilename = options.outputfile

print "Input:",options.inputfile, "; Output:", outfilename


#outfilename= 'SkimmedTree.root'
skimmedTree = TChain("tree/treeMaker")

if isfarmout:
    infile = open(inputfilename)
    failcount=0
    for ifile in infile:
        try:
            f_tmp = TFile.Open(ifile.rstrip(),'READ')
            if f_tmp.IsZombie():            # or fileIsCorr(ifile.rstrip()):
                failcount += 1
                continue
            skimmedTree.Add(ifile.rstrip())
        except:
            failcount += 1
    if failcount>0: print "Could not read %d files. Skipping them." %failcount

if not isfarmout:
    skimmedTree.Add(inputfilename)

def AnalyzeDataSet():

    outfile = TFile(outfilename,'RECREATE')

    NEntries = skimmedTree.GetEntries()
    NEntries = 1000000
    print 'NEntries = '+str(NEntries)
    npass = 0


    h_beff_num=TH2D("h_beff_num","h_beff_num_",10,-2.4,2.4,20,20.,1000.)
    h_beff_den=TH2D("h_beff_den","h_beff_den",10,-2.4,2.4,20,20.,1000.)
    h_lighteff_num=TH2D("h_lighteff_num","h_lighteff_num",10,-2.4,2.4,20,20.,1000.)
    h_lighteff_den=TH2D("h_lighteff_den","h_lighteff_den",10,-2.4,2.4,20,20.,1000.)

    CSVMWP = 0.8484
    for ievent in range(NEntries):
        if ievent%100==0: print "Processed "+str(ievent)+" of "+str(NEntries)+" events."
        skimmedTree.GetEntry(ievent)
        ## Get all relevant branches
        try:

            nTHINJets                  = skimmedTree.__getattr__('THINnJet')
            thinjetP4                  = skimmedTree.__getattr__('THINjetP4')
            thinJetCSV                 = skimmedTree.__getattr__('THINjetCISVV2')
            passThinJetTightID         = skimmedTree.__getattr__('THINjetPassIDTight')
            THINjetHadronFlavor        = skimmedTree.__getattr__('THINjetHadronFlavor')

        except Exception as e:
            print e
            print "Corrupt file detected! Skipping 1 event."
            continue



        for nb in range(nTHINJets):
            if THINjetHadronFlavor[nb]==5:
                h_beff_den.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())
                h_lighteff_den.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())
                if thinJetCSV[nb] > CSVMWP:
                    h_beff_num.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())
                else:
                    h_lighteff_num.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())

            #if THINjetHadronFlavor[nb]==4 or THINjetHadronFlavor[nb]==0:
                #h_lighteff_den.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())
                #if thinJetCSV[nb] > CSVMWP:
                    #h_lighteff_num.Fill(thinjetP4[nb].Eta(),thinjetP4[nb].Pt())

    h_beff_num.Write()
    h_beff_den.Write()
    h_lighteff_num.Write()
    h_lighteff_den.Write()
    outfile.Write()


    print "ROOT file written to", outfilename

    print "Completed."


if __name__ == "__main__":
    AnalyzeDataSet()
