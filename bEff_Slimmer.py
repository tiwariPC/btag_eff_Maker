#!/usr/bin/env python
from ROOT import TFile, TTree, TH1F, TH1D, TH1, TCanvas, TChain, TGraphAsymmErrors, TMath, TH2D, TLorentzVector, AddressOf, gROOT, TNamed
import ROOT as ROOT
import os,traceback
import sys, optparse,argparse
from array import array
import math
import numpy as numpy
import pandas
from root_pandas import read_root
from pandas import  DataFrame, concat
from pandas import Series
import time
import glob

## for parallel threads in interactive running
from multiprocessing import Process
import multiprocessing as mp


isCondor = False

## user packages
## in local dir
sys.path.append('skim_configs')
import  triggers as trig
import variables as branches
import filters as filters
#import genPtProducer as GenPtProd

## from commonutils
if isCondor:sys.path.append('ExoPieUtils/commonutils/')
else:sys.path.append('../ExoPieUtils/commonutils/')
import MathUtils as mathutil
from MathUtils import *
import BooleanUtils as boolutil


## from analysisutils
if isCondor:sys.path.append('ExoPieUtils/analysisutils/')
else:sys.path.append('../ExoPieUtils/analysisutils/')
import analysis_utils as anautil
import genPtProducer as GenPtProd
######################################################################################################
## All import are done before this
######################################################################################################


runInteractive = False
runOn2016 = False
runOn2017 = False
runOn2018 = False

## ----- start if clock

start = time.clock()


## ----- command line argument
usage = "analyzer for bb+DM (debugging) "
parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--inputfile",  dest="inputfile",
                    default="myfiles.txt")
parser.add_argument("-inDir", "--inputDir",  dest="inputDir", default=".")
parser.add_argument("-runOnTXT", "--runOnTXT",
                    action="store_true", dest="runOnTXT")
parser.add_argument("-o", "--outputfile",
                    dest="outputfile", default="out.root")
parser.add_argument("-D", "--outputdir", dest="outputdir")
parser.add_argument("-F", "--farmout", action="store_true",  dest="farmout")
parser.add_argument("-y", "--year", dest="year", default="Year")

#parser.add_argument("runOnTXT", "--runOnTXT",dest="runOnTXT")
#parser.set_defaults(runOnTXT=False)
## add argument for debug, default to be false

args = parser.parse_args()

if args.farmout == None:
    isfarmout = False
else:
    isfarmout = args.farmout

if args.inputDir and isfarmout:
    dirName = args.inputDir

if args.year == '2016':
    runOn2016 = True
elif args.year == '2017':
    runOn2017 = True
elif args.year == '2018':
    runOn2018 = True
else:
    print('Please provide on which year you want to run?')
    sys.exit()

runOnTxt = False
if args.runOnTXT:
    runOnTxt = True

if isfarmout:
    infile = args.inputfile

else:
    print "No file is provided for farmout"


outputdir = '.'
if args.outputdir:
    outputdir = str(args.outputdir)

infilename = "ExoPieElementTuples.root"

debug_ = False

if runOn2016 or runOn2017 or runOn2018:
    from TheaCorrection import TheaCorrection_2016 as TheaCorrection
#please provide the latest recommended working points from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation
if runOn2016:
    filter_list = filters.filters2016
    LWP = 0.2217
    MWP = 0.6321
    TWP = 0.8953
elif runOn2017:
    filter_list = filters.filters2017
    LWP = 0.1522
    MWP = 0.4941
    TWP = 0.8001
elif runOn2018:
    filter_list = filters.filters2018
    LWP = 0.1241
    MWP = 0.4184
    TWP = 0.7527
def whichsample(filename):
    print "SkimTree:whichsample:-- file name is ", filename
    sample = -999
    if "TTT" in filename:
        sample = 6
    elif "WJetsToLNu_HT" in filename:
        sample = 24
    elif "ZJetsToNuNu_HT" in filename or "DYJetsToLL" in filename:
        sample = 23
    return sample

def TextToList(textfile):
    return([iline.rstrip() for iline in open(textfile)])

## the input file list and key is caught in one variable as a python list,
#### first element is the list of rootfiles
#### second element is the key, user to name output.root

def runbbdm(txtfile):
    infile_ = []
    outfilename = ""
    prefix = "Skimmed_"
    ikey_ = ""

    if not runInteractive:
        print "running for ", txtfile
        infile_ = TextToList(txtfile)
        outfile = txtfile.split('/')[-1].replace('.txt', '.root')
        #key_=txtfile[1]


        outfilename = outfile  # prefix+key_+".root"
        print "outfilename", outfilename

    if runInteractive:
        infile_ = TextToList(txtfile)
        print "infile_ = ", infile_
        #print "running code for ",infile_
        prefix_ = ''  # '/eos/cms/store/group/phys_exotica/bbMET/2017_skimmedFiles/locallygenerated/'
        if outputdir != '.':
            prefix_ = outputdir+'/'
        #print "prefix_", prefix_
        outfilename = prefix_+txtfile.split('/')[-1].replace('.txt', '.root')
        print 'outfilename',  outfilename

    samplename = whichsample(outfilename)
    print ("samplename = ", samplename)

    if runOn2016:
        triglist = trig.trigger2016
        eletrig = trig.Electrontrigger2016
        muontrig = trig.Muontrigger2016
        mettrig = trig.METtrigger2016
        photontrig = trig.Photontrigger2016
    elif runOn2017:
        triglist = trig.trigger2017
        eletrig = trig.Electrontrigger2017
        muontrig = trig.Muontrigger2017
        mettrig = trig.METtrigger2017
        photontrig = trig.Photontrigger2017
    elif runOn2018:
        triglist = trig.trigger2018
        eletrig = trig.Electrontrigger2018
        muontrig = trig.Muontrigger2018
        mettrig = trig.METtrigger2018
        photontrig = trig.Photontrigger2018

    bins_pT     = [20.0,50.0,80.0,120.0,200.0,300.0,400.0,500.0,700.0,1000.0]
    bins_eta    = [-2.5,-1.5,-0.5,0.0,0.5,1.5,2.5]

    h_btag_num_pass_mwp    =  TH2D("h_btag_num_pass_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_btag_num_fail_mwp    =  TH2D("h_btag_num_fail_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_btag_den             =  TH2D("h_btag_den","",6,array('d',bins_eta),9,array('d',bins_pT))

    h_ctag_num_pass_mwp    =  TH2D("h_ctag_num_pass_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_ctag_num_fail_mwp    =  TH2D("h_ctag_num_fail_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_ctag_den             =  TH2D("h_ctag_den","",6,array('d',bins_eta),9,array('d',bins_pT))

    h_lighttag_num_pass_mwp=  TH2D("h_lighttag_num_pass_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_lighttag_num_fail_mwp=  TH2D("h_lighttag_num_fail_mwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_lighttag_den         =  TH2D("h_lighttag_den","",6,array('d',bins_eta),9,array('d',bins_pT))

    h_btag_num_pass_lwp    = TH2D("h_btag_num_pass_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_btag_num_fail_lwp    = TH2D("h_btag_num_fail_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))

    h_ctag_num_pass_lwp    = TH2D("h_ctag_num_pass_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_ctag_num_fail_lwp    = TH2D("h_ctag_num_fail_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))

    h_lighttag_num_pass_lwp= TH2D("h_lighttag_num_pass_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))
    h_lighttag_num_fail_lwp= TH2D("h_lighttag_num_fail_lwp","",6,array('d',bins_eta),9,array('d',bins_pT))

    passfilename = open("skim_configs/outfilename.txt", "w")

    passfilename.write(outfilename)
    passfilename.close()


    ## following can be moved to outputtree.py if we manage to change the name of output root file.
    outfilenameis = outfilename
    outfile = TFile(outfilenameis,'RECREATE')

    ## following can be moved to outputtree.py if we manage to change the name of output root file.

    if runOn2016:
        jetvariables = branches.allvars2016
    elif runOn2017:
        jetvariables = branches.allvars2017
    elif runOn2018:
        jetvariables = branches.allvars2018

    filename = infile_

    ieve = 0;icount = 0
    #print "running on", filename
    for df in read_root(filename, 'tree/treeMaker', columns=jetvariables, chunksize=125000):
        if runOn2016:
            var_zip = zip(df.runId, df.lumiSection, df.eventId, df.isData, df.mcWeight,
                          df.prefiringweight, df.prefiringweightup, df.prefiringweightdown,
                          df.pu_nTrueInt, df.pu_nPUVert,
                          df.hlt_trigName, df.hlt_trigResult, df.hlt_filterName, df.hlt_filterResult,
                          df.pfpatMet_smear, df.pfMetCorrPt, df.pfMetCorrPhi, df.pfMetCorrUnc,
                          df.pfMetCorrSig, df.pfpatCaloMETPt, df.pfpatCaloMETPhi, df.pfTRKMETPt_, df.pfTRKMETPhi_,
                          df.nEle, df.elePx, df.elePy, df.elePz, df.eleEnergy, df.eleIsPassVeto, df.eleIsPassLoose, df.eleIsPassTight, df.eleD0, df.eleDz,
                          df.eleCharge, df.nPho, df.phoPx, df.phoPy, df.phoPz, df.phoEnergy, df.phoIsPassLoose, df.phoIsPassTight,
                          df.nMu, df.muPx, df.muPy, df.muPz, df.muEnergy, df.isLooseMuon, df.isTightMuon, df.PFIsoLoose, df.PFIsoMedium, df.PFIsoTight, df.PFIsoVeryTight, df.muCharge,
                          df.HPSTau_n, df.HPSTau_Px, df.HPSTau_Py, df.HPSTau_Pz, df.HPSTau_Energy, df.disc_decayModeFinding, df.disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byTightIsolationMVArun2017v2DBoldDMwLT2017,
                          df.disc_againstMuonLoose3, df.disc_againstMuonTight3, df.disc_againstElectronLooseMVA6, df.disc_againstElectronMediumMVA6, df.disc_againstElectronTightMVA6,
                          df.nGenPar, df.genParId, df.genMomParId, df.genParSt, df.genParPx, df.genParPy, df.genParPz, df.genParE,
                          df.THINnJet, df.THINjetPx, df.THINjetPy, df.THINjetPz, df.THINjetEnergy, df.THINbRegNNResolution, df.THINbRegNNCorr,df.THINisPUJetIDLoose,df.THINisPUJetIDMedium,df.THINisPUJetIDTight,
                          df.THINjetPassIDLoose, df.THINjetDeepCSV_b, df.THINjetHadronFlavor, df.THINjetCEmEF, df.THINjetCHadEF, df.THINjetNEmEF, df.THINjetNHadEF, df.THINjetCMulti, df.THINjetNMultiplicity, df.THINjetCorrUncUp, df.THINjetNPV,
                          df.FATnJet, df.FATjetPx, df.FATjetPy, df.FATjetPz, df.FATjetEnergy, df.FATgenjetpx, df.FATgenjetpy, df.FATgenjetpz, df.FATgenjetE, df.FATjetPassIDLoose,
                          df.FATjet_DoubleSV, df.FATjet_probQCDb, df.FATjet_probHbb, df.FATjet_probQCDc, df.FATjet_probHcc, df.FATjet_probHbbc,
                          df.FATjet_prob_bbvsLight, df.FATjet_prob_ccvsLight, df.FATjet_prob_TvsQCD, df.FATjet_prob_WvsQCD, df.FATjet_prob_ZHbbvsQCD,
                          df.FATjetSDmass, df.FATN2_Beta1_, df.FATN2_Beta2_, df.FATjetCHSPRmassL2L3Corr, df.FATjetCHSSDmassL2L3Corr, df.FATjetTau1, df.FATjetTau2)
        elif runOn2017:
            var_zip = zip(df.runId, df.lumiSection, df.eventId, df.isData, df.mcWeight,
                          df.prefiringweight, df.prefiringweightup, df.prefiringweightdown,
                          df.pu_nTrueInt, df.pu_nPUVert,
                          df.hlt_trigName, df.hlt_trigResult, df.hlt_filterName, df.hlt_filterResult,
                          df.pfpatmodifiedMet_smear, df.pfmodifiedMetCorrPt, df.pfmodifiedMetCorrPhi, df.pfmodifiedMetCorrUnc,
                          df.pfmodifiedMetCorrSig, df.pfpatCaloMETPt, df.pfpatCaloMETPhi, df.pfTRKMETPt_, df.pfTRKMETPhi_,
                          df.nEle, df.elePx, df.elePy, df.elePz, df.eleEnergy, df.eleIsPassVeto, df.eleIsPassLoose, df.eleIsPassTight, df.eleD0, df.eleDz,
                          df.eleCharge, df.nPho, df.phoPx, df.phoPy, df.phoPz, df.phoEnergy, df.phoIsPassLoose, df.phoIsPassTight,
                          df.nMu, df.muPx, df.muPy, df.muPz, df.muEnergy, df.isLooseMuon, df.isTightMuon, df.PFIsoLoose, df.PFIsoMedium, df.PFIsoTight, df.PFIsoVeryTight, df.muCharge,
                          df.HPSTau_n, df.HPSTau_Px, df.HPSTau_Py, df.HPSTau_Pz, df.HPSTau_Energy, df.disc_decayModeFinding, df.disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byTightIsolationMVArun2017v2DBoldDMwLT2017,
                          df.disc_againstMuonLoose3, df.disc_againstMuonTight3, df.disc_againstElectronLooseMVA6, df.disc_againstElectronMediumMVA6, df.disc_againstElectronTightMVA6,
                          df.nGenPar, df.genParId, df.genMomParId, df.genParSt, df.genParPx, df.genParPy, df.genParPz, df.genParE,
                          df.THINnJet, df.THINjetPx, df.THINjetPy, df.THINjetPz, df.THINjetEnergy, df.THINbRegNNResolution, df.THINbRegNNCorr,df.THINisPUJetIDLoose,df.THINisPUJetIDMedium,df.THINisPUJetIDTight,
                          df.THINjetPassIDTight, df.THINjetDeepCSV_b, df.THINjetHadronFlavor, df.THINjetCEmEF, df.THINjetCHadEF, df.THINjetNEmEF, df.THINjetNHadEF, df.THINjetCMulti, df.THINjetNMultiplicity, df.THINjetCorrUncUp, df.THINjetNPV,
                          df.FATnJet, df.FATjetPx, df.FATjetPy, df.FATjetPz, df.FATjetEnergy, df.FATgenjetpx, df.FATgenjetpy, df.FATgenjetpz, df.FATgenjetE, df.FATjetPassIDTight,
                          df.FATjet_DoubleSV, df.FATjet_probQCDb, df.FATjet_probHbb, df.FATjet_probQCDc, df.FATjet_probHcc, df.FATjet_probHbbc,
                          df.FATjet_prob_bbvsLight, df.FATjet_prob_ccvsLight, df.FATjet_prob_TvsQCD, df.FATjet_prob_WvsQCD, df.FATjet_prob_ZHbbvsQCD,
                          df.FATjetSDmass, df.FATN2_Beta1_, df.FATN2_Beta2_, df.FATjetCHSPRmassL2L3Corr, df.FATjetCHSSDmassL2L3Corr, df.FATjetTau1, df.FATjetTau2)
        elif runOn2018:
            df['prefiringweight'] = 1.0
            df['prefiringweightup'] = 1.0
            df['prefiringweightdown'] = 1.0
            var_zip = zip(df.runId, df.lumiSection, df.eventId, df.isData, df.mcWeight,
                          df.prefiringweight, df.prefiringweightup, df.prefiringweightdown,
                          df.pu_nTrueInt, df.pu_nPUVert,
                          df.hlt_trigName, df.hlt_trigResult, df.hlt_filterName, df.hlt_filterResult,
                          df.pfpatMet_smear, df.pfMetCorrPt, df.pfMetCorrPhi, df.pfMetCorrUnc,
                          df.pfMetCorrSig, df.pfpatCaloMETPt, df.pfpatCaloMETPhi, df.pfTRKMETPt_, df.pfTRKMETPhi_,
                          df.nEle, df.elePx, df.elePy, df.elePz, df.eleEnergy, df.eleIsPassVeto, df.eleIsPassLoose, df.eleIsPassTight, df.eleD0, df.eleDz,
                          df.eleCharge, df.nPho, df.phoPx, df.phoPy, df.phoPz, df.phoEnergy, df.phoIsPassLoose, df.phoIsPassTight,
                          df.nMu, df.muPx, df.muPy, df.muPz, df.muEnergy, df.isLooseMuon, df.isTightMuon, df.PFIsoLoose, df.PFIsoMedium, df.PFIsoTight, df.PFIsoVeryTight, df.muCharge,
                          df.HPSTau_n, df.HPSTau_Px, df.HPSTau_Py, df.HPSTau_Pz, df.HPSTau_Energy, df.disc_decayModeFinding, df.disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017, df.disc_byTightIsolationMVArun2017v2DBoldDMwLT2017,
                          df.disc_againstMuonLoose3, df.disc_againstMuonTight3, df.disc_againstElectronLooseMVA6, df.disc_againstElectronMediumMVA6, df.disc_againstElectronTightMVA6,
                          df.nGenPar, df.genParId, df.genMomParId, df.genParSt, df.genParPx, df.genParPy, df.genParPz, df.genParE,
                          df.THINnJet, df.THINjetPx, df.THINjetPy, df.THINjetPz, df.THINjetEnergy, df.THINbRegNNResolution, df.THINbRegNNCorr,df.THINisPUJetIDLoose,df.THINisPUJetIDMedium,df.THINisPUJetIDTight,
                          df.THINjetPassIDTight, df.THINjetDeepCSV_b, df.THINjetHadronFlavor, df.THINjetCEmEF, df.THINjetCHadEF, df.THINjetNEmEF, df.THINjetNHadEF, df.THINjetCMulti, df.THINjetNMultiplicity, df.THINjetCorrUncUp, df.THINjetNPV,
                          df.FATnJet, df.FATjetPx, df.FATjetPy, df.FATjetPz, df.FATjetEnergy, df.FATgenjetpx, df.FATgenjetpy, df.FATgenjetpz, df.FATgenjetE, df.FATjetPassIDTight,
                          df.FATjet_DoubleSV, df.FATjet_probQCDb, df.FATjet_probHbb, df.FATjet_probQCDc, df.FATjet_probHcc, df.FATjet_probHbbc,
                          df.FATjet_prob_bbvsLight, df.FATjet_prob_ccvsLight, df.FATjet_prob_TvsQCD, df.FATjet_prob_WvsQCD, df.FATjet_prob_ZHbbvsQCD,
                          df.FATjetSDmass, df.FATN2_Beta1_, df.FATN2_Beta2_, df.FATjetCHSPRmassL2L3Corr, df.FATjetCHSSDmassL2L3Corr, df.FATjetTau1, df.FATjetTau2)
        for run, lumi, event, isData, mcWeight_,\
                prefiringweight_, prefiringweightup_, prefiringweightdown_,\
                pu_nTrueInt_, pu_nPUVert_,\
                trigName_, trigResult_, filterName, filterResult,\
                met_smear, met_, metphi_, metUnc_,\
                metCorrSig, patCaloMETPt, patCaloMETPhi, TRKMETPt_, TRKMETPhi_,\
                nele_, elepx_, elepy_, elepz_, elee_, elevetoid_, elelooseid_, eletightid_, eleD0_, eleDz_,\
                eleCharge_, npho_, phopx_, phopy_, phopz_, phoe_, pholooseid_, photightID_,\
                nmu_, mupx_, mupy_, mupz_, mue_, mulooseid_, mutightid_, muisoloose, muisomedium, muisotight, muisovtight, muCharge_,\
                nTau_, tau_px_, tau_py_, tau_pz_, tau_e_, tau_dm_, tau_isLoose_, tau_isoMedium_, tau_isoTight_,\
                Taudisc_againstLooseMuon, Taudisc_againstTightMuon, Taudisc_againstLooseElectron, Taudisc_againstMediumElectron, Taudisc_againstTightElectron,\
                nGenPar_, genParId_, genMomParId_, genParSt_, genpx_, genpy_, genpz_, gene_,\
                nak4jet_, ak4px_, ak4py_, ak4pz_, ak4e_, ak4bRegNNResolution, ak4bRegNNCorr,ak4PUJetIDLoose,ak4PUJetIDMedium,ak4PUJetIDTight,\
                ak4PassID_, ak4deepcsv_, ak4flavor_, ak4CEmEF_, ak4CHadEF_, ak4NEmEF_, ak4NHadEF_, ak4CMulti_, ak4NMultiplicity_, ak4JEC_, ak4NPV_,\
                fatnJet, fatjetPx, fatjetPy, fatjetPz, fatjetEnergy, fatgenjetPx, fatgenjetPy, fatgenjetPz, fatgenjetEnergy, fatjetPassID,\
                fatjet_DoubleSV, fatjet_probQCDb, fatjet_probHbb, fatjet_probQCDc, fatjet_probHcc, fatjet_probHbbc,\
                fatjet_prob_bbvsLight, fatjet_prob_ccvsLight, fatjet_prob_TvsQCD, fatjet_prob_WvsQCD, fatjet_prob_ZHbbvsQCD,\
                fatjetSDmass, fatN2_Beta1_, fatN2_Beta2_, fatjetCHSPRmassL2L3Corr, fatjetCHSSDmassL2L3Corr, fatjetTau1, fatjetTau2\
                in var_zip:


            if ieve % 1000 == 0:
                print ("Processed", ieve, "Events")
            ieve = ieve + 1
            # -------------------------------------------------
            # MC Weights
            # -------------------------------------------------
            # mcweight[0] = 0.0
            # if isData == 1:
            #     mcweight[0] = 1.0
            # if not isData:
            #     if mcWeight_ < 0:
            #         mcweight[0] = -1.0
            #     if mcWeight_ > 0:
            #         mcweight[0] = 1.0
            # h_total.Fill(1.)
            # h_total_mcweight.Fill(1., mcweight[0])

            # -------------------------------------------------
            ## Trigger selection
            # -------------------------------------------------

            eletrigdecision = False
            mudecision = False
            metdecision = False
            phodecision = False

            eletrigstatus = [(anautil.CheckFilter(
                trigName_, trigResult_, eletrig[itrig])) for itrig in range(len(eletrig))]
            mutrigstatus = [(anautil.CheckFilter(
                trigName_, trigResult_, muontrig[itrig])) for itrig in range(len(muontrig))]
            mettrigstatus = [(anautil.CheckFilter(
                trigName_, trigResult_, mettrig[itrig])) for itrig in range(len(mettrig))]
            photrigstatus = [(anautil.CheckFilter(
                trigName_, trigResult_, photontrig[itrig])) for itrig in range(len(photontrig))]

            eletrigdecision = boolutil.logical_OR(eletrigstatus)
            mutrigdecision = boolutil.logical_OR(mutrigstatus)
            mettrigdecision = boolutil.logical_OR(mettrigstatus)
            photrigdecision = boolutil.logical_OR(photrigstatus)

            # ------------------------------------------------------
            ## Filter selection
            # ------------------------------------------------------
            filterdecision = False
            filterstatus = [False for ifilter in range(len(filter_list))]
            filterstatus = [anautil.CheckFilter(
                filterName, filterResult, filter_list[ifilter]) for ifilter in range(len(filter_list))]

            filterdecision = boolutil.logical_AND(filterstatus)
            if filterdecision == False and isData:
                continue

            # ------------------------------------------------------
            ## PFMET Selection
            # --------------------------------------------------------
            pfmetstatus = (met_ > 200.0)


            '''
            ****   *      ****
            *      *      *
            ***    *      ***
            *      *      *
            ****   ****   ****
            '''
            elept = getPt(elepx_, elepy_)
            eleeta = getEta(elepx_, elepy_, elepz_)
            elephi = getPhi(elepx_, elepy_)

            ele_pt10_eta2p5_vetoID = boolutil.logical_and3((elept > 10.0), (elevetoid_),  numpy.logical_and(
                numpy.logical_or(numpy.abs(eleeta) > 1.566, numpy.abs(eleeta) < 1.4442), (numpy.abs(eleeta) < 2.5)))

            ele_pt10_eta2p5_looseID = boolutil.logical_and3((elept > 10.0), (elelooseid_),  numpy.logical_and(
                numpy.logical_or(numpy.abs(eleeta) > 1.566, numpy.abs(eleeta) < 1.4442), (numpy.abs(eleeta) < 2.5)))

            ele_pt30_eta2p5_tightID = boolutil.logical_and3((elept > 30.0), (eletightid_),  numpy.logical_and(numpy.logical_or(boolutil.logical_and3(numpy.abs(eleeta) > 1.566, numpy.abs(
                eleD0_) < 0.10, numpy.abs(eleDz_) < 0.20), boolutil.logical_and3(numpy.abs(eleeta) < 1.4442, numpy.abs(eleD0_) < 0.05, numpy.abs(eleDz_) < 0.10)), (numpy.abs(eleeta) < 2.5)))

            pass_ele_veto_index = boolutil.WhereIsTrue(ele_pt10_eta2p5_vetoID)
            pass_ele_loose_index = boolutil.WhereIsTrue(
                ele_pt10_eta2p5_looseID)

            '''
            **     *  *     *
            * *  * *  *     *
            *  *   *  *     *
            *      *  *     *
            *      *   *****
            '''
            mupt = getPt(mupx_, mupy_)
            mueta = getEta(mupx_, mupy_, mupz_)
            muphi = getPhi(mupx_, mupy_)
            mu_pt10_eta2p4_looseID_looseISO = boolutil.logical_and4(
                mupt > 10.0, numpy.abs(mueta) < 2.4,  mulooseid_, muisoloose)
            mu_pt30_eta2p4_tightID_tightISO = boolutil.logical_and4(
                (mupt > 30.0), (numpy.abs(mueta) < 2.4), (mutightid_), (muisotight))

            pass_mu_index = boolutil.WhereIsTrue(
                mu_pt10_eta2p4_looseID_looseISO)

            '''
            *******   *      *   ******
            *     *   *      *  *      *
            *******   ********  *      *
            *         *      *  *      *
            *         *      *   ******
            '''

            phopt = getPt(phopx_, phopy_)
            phoeta = getEta(phopx_, phopy_, phopz_)
            phophi = getPhi(phopx_, phopy_)

            pho_pt15_eta2p5_looseID = boolutil.logical_and3(
                (phopt > 15.0),   (numpy.abs(phoeta) < 2.5),  (pholooseid_))
            pass_pho_index = boolutil.WhereIsTrue(pho_pt15_eta2p5_looseID)

            cleanedPho_ag_ele = []; cleanedPho_ag_mu = [];pass_pho_index_cleaned=[]
            if npho_ > 0: #and ep_nEle > 0:
                cleanedPho_ag_ele = anautil.jetcleaning(pho_pt15_eta2p5_looseID, ele_pt10_eta2p5_looseID, phoeta, eleeta, phophi, elephi, 0.4)
                cleanedPho_ag_mu  = anautil.jetcleaning(pho_pt15_eta2p5_looseID, mu_pt10_eta2p4_looseID_looseISO, phoeta, mueta, phophi, muphi, 0.4)
                cleanedPhoton     = boolutil.logical_AND_List2(cleanedPho_ag_ele,cleanedPho_ag_mu)
                pass_pho_index_cleaned = boolutil.WhereIsTrue(cleanedPhoton)

            ## Fill variables for the CRs which require lepton.
            WenuRecoil = -1
            WenuRecoilSmearPt = -1
            Wenumass = -1
            WenuPhi = -10

            WmunuRecoil = -1
            WmunuRecoilSmearPt = -1
            Wmunumass = -1
            WmunuPhi = -10

            ZeeRecoil = -1
            ZeeRecoilSmear = -1
            Zeemass = -1
            ZeePhi = -10

            ZmumuRecoil = -1
            ZmumuRecoilSmear = -1
            Zmumumass = -1
            ZmumuPhi = -10

            GammaRecoil = -1
            GammaRecoilSmearPt = -1
            GammaPhi = -10
            if debug_:
                print ('Reached Fill variables')

            # ------------------
            # Z CR
            # ------------------
            ## for dielectron
            if len(pass_ele_loose_index) == 2:
                iele1 = pass_ele_loose_index[0]
                iele2 = pass_ele_loose_index[1]
                if eleCharge_[iele1]*eleCharge_[iele2] < 0:
                    ee_mass = InvMass(elepx_[iele1], elepy_[iele1], elepz_[iele1], elee_[
                                      iele1], elepx_[iele2], elepy_[iele2], elepz_[iele2], elee_[iele2])
                    zeeRecoilPx = -(met_*math.cos(metphi_) +
                                    elepx_[iele1] + elepx_[iele2])
                    zeeRecoilPy = -(met_*math.sin(metphi_) +
                                    elepy_[iele1] + elepy_[iele2])
                    ZeeRecoilPt = math.sqrt(zeeRecoilPx**2 + zeeRecoilPy**2)
                    if ee_mass > 60.0 and ee_mass < 120.0 and ZeeRecoilPt > 200.0:
                        ZeeRecoil = ZeeRecoilPt
                        Zeemass = ee_mass
                        ZeePhi = mathutil.ep_arctan(
                            zeeRecoilPx, zeeRecoilPy)
                    zeeRecoilSmearPx = - \
                        (met_*math.cos(metphi_) +
                         elepx_[iele1] + elepx_[iele2])
                    zeeRecoilSmearPy = - \
                        (met_*math.sin(metphi_) +
                         elepy_[iele1] + elepy_[iele2])
                    ZeeRecoilSmearPt = math.sqrt(
                        zeeRecoilSmearPx**2 + zeeRecoilSmearPy**2)
                    if ee_mass > 60.0 and ee_mass < 120.0 and ZeeRecoilSmearPt > 200.0:
                        ZeeRecoilSmear = ZeeRecoilSmearPt
            ## for dimu
            if len(pass_mu_index) == 2:
                imu1 = pass_mu_index[0]
                imu2 = pass_mu_index[1]
                if muCharge_[imu1]*muCharge_[imu2] < 0:
                    mumu_mass = InvMass(mupx_[imu1], mupy_[imu1], mupz_[imu1], mue_[
                                        imu1], mupx_[imu2], mupy_[imu2], mupz_[imu2], mue_[imu2])
                    zmumuRecoilPx = -(met_*math.cos(metphi_) +
                                      mupx_[imu1] + mupx_[imu2])
                    zmumuRecoilPy = -(met_*math.sin(metphi_) +
                                      mupy_[imu1] + mupy_[imu2])
                    ZmumuRecoilPt = math.sqrt(
                        zmumuRecoilPx**2 + zmumuRecoilPy**2)
                    if mumu_mass > 60.0 and mumu_mass < 120.0 and ZmumuRecoilPt > 200.0:
                        ZmumuRecoil = ZmumuRecoilPt
                        Zmumumass = mumu_mass
                        ZmumuPhi = mathutil.ep_arctan(
                            zmumuRecoilPx, zmumuRecoilPy)
                    zmumuRecoilSmearPx = - \
                        (met_*math.cos(metphi_) +
                         mupx_[imu1] + mupx_[imu2])
                    zmumuRecoilSmearPy = - \
                        (met_*math.sin(metphi_) +
                         mupy_[imu1] + mupy_[imu2])
                    ZmumuRecoilSmearPt = math.sqrt(
                        zmumuRecoilSmearPx**2 + zmumuRecoilSmearPy**2)
                    if mumu_mass > 60.0 and mumu_mass < 120.0 and ZmumuRecoilSmearPt > 200.0:
                        ZmumuRecoilSmear = ZmumuRecoilSmearPt
            if len(pass_ele_loose_index) == 2:
                ZRecoilstatus = (ZeeRecoil > 200.0) or (
                    ZeeRecoilSmear > 200.0)
            elif len(pass_mu_index) == 2:
                ZRecoilstatus = (ZmumuRecoil > 200.0) or (
                    ZmumuRecoilSmear > 200.0)
            else:
                ZRecoilstatus = False
            if debug_:
                print 'Reached Z CR'

            # ------------------
            # W CR
            # ------------------
            ## for Single electron
            if len(pass_ele_loose_index) == 1:
                ele1 = pass_ele_loose_index[0]
                # transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}
                e_mass = MT(elept[ele1], met_, DeltaPhi(elephi[ele1], metphi_))
                WenuRecoilPx = -(met_*math.cos(metphi_) + elepx_[ele1])
                WenuRecoilPy = -(met_*math.sin(metphi_) + elepy_[ele1])
                WenuRecoilPt = math.sqrt(WenuRecoilPx**2 + WenuRecoilPy**2)
                if WenuRecoilPt > 200.0:
                   WenuRecoil = WenuRecoilPt
                   Wenumass = e_mass
                   WenuPhi = mathutil.ep_arctan(WenuRecoilPx, WenuRecoilPy)
                WenuRecoilSmearPx = - \
                    (met_*math.cos(metphi_) + elepx_[ele1])
                WenuRecoilSmearPy = - \
                    (met_*math.sin(metphi_) + elepy_[ele1])
                WenuRecoilSmearPt = math.sqrt(
                    WenuRecoilSmearPx**2 + WenuRecoilSmearPy**2)

            ## for Single muon
            if len(pass_mu_index) == 1:
                mu1 = pass_mu_index[0]
                # transverse mass defined as sqrt{2pT*MET*(1-cos(dphi)}
                mu_mass = MT(mupt[mu1], met_, DeltaPhi(muphi[mu1], metphi_))
                WmunuRecoilPx = -(met_*math.cos(metphi_) + mupx_[mu1])
                WmunuRecoilPy = -(met_*math.sin(metphi_) + mupy_[mu1])
                WmunuRecoilPt = math.sqrt(WmunuRecoilPx**2 + WmunuRecoilPy**2)
                if WmunuRecoilPt > 200.0:
                   WmunuRecoil = WmunuRecoilPt
                   Wmunumass = mu_mass
                   WmunuPhi = mathutil.ep_arctan(
                       WmunuRecoilPx, WmunuRecoilPy)
                WmunuRecoilSmearPx = - \
                    (met_*math.cos(metphi_) + mupx_[mu1])
                WmunuRecoilSmearPy = - \
                    (met_*math.sin(metphi_) + mupy_[mu1])
                WmunuRecoilSmearPt = math.sqrt(
                    WmunuRecoilSmearPx**2 + WmunuRecoilSmearPy**2)

            if len(pass_ele_loose_index) == 1:
                WRecoilstatus = (WenuRecoil > 200.0) or (
                    WenuRecoilSmearPt > 200.0)
            elif len(pass_mu_index) == 1:
                WRecoilstatus = (WmunuRecoil > 200.0) or (
                    WmunuRecoilSmearPt > 200.0)
            else:
                WRecoilstatus = False
            if debug_:
                print 'Reached W CR'

            # ------------------
            # Gamma CR
            # ------------------
            ## for Single photon
            if len(pass_pho_index) >= 1:
                pho1 = pass_pho_index[0]
                GammaRecoilPx = -(met_*math.cos(metphi_) + phopx_[pho1])
                GammaRecoilPy = -(met_*math.sin(metphi_) + phopy_[pho1])
                GammaRecoilPt = math.sqrt(GammaRecoilPx**2 + GammaRecoilPy**2)
                if GammaRecoilPt > 200.0:
                    GammaRecoil = GammaRecoilPt
                    GammaPhi = mathutil.ep_arctan(
                        GammaRecoilPx, GammaRecoilPy)
                GammaRecoilSmearPx = - \
                    (met_*math.cos(metphi_) + phopx_[pho1])
                GammaRecoilSmearPy = - \
                    (met_*math.sin(metphi_) + phopy_[pho1])
                GammaRecoilSmearPt = math.sqrt(
                    GammaRecoilSmearPx**2 + GammaRecoilSmearPy**2)

            GammaRecoilStatus = (GammaRecoil > 200.0) or (
                GammaRecoilSmearPt > 200.0)
            if debug_:
                print 'Reached Gamma CR'

            if pfmetstatus == False and ZRecoilstatus == False and WRecoilstatus == False and GammaRecoilStatus == False:
                continue

            '''
            *******   *****   *******
               *      *          *
               *      ****       *
               *      *          *
            ***       *****      *
            '''
            '''
            ak4pt = [getPt(ak4px_[ij], ak4py_[ij]) for ij in range(nak4jet_)]
            ak4eta = [getEta(ak4px_[ij], ak4py_[ij], ak4pz_[ij]) for ij in range(nak4jet_)]
            ak4phi = [getPhi(ak4px_[ij], ak4py_[ij]) for ij in range(nak4jet_)]
            '''
            ak4pt = getPt(ak4px_, ak4py_)
            ak4eta = getEta(ak4px_, ak4py_, ak4pz_)
            ak4phi = getPhi(ak4px_, ak4py_)
            if runOn2016:
                #ak4PassID_Calc = [jetID_(ak4CEmEF_[ij],ak4CHadEF_[ij],ak4NEmEF_[ij],ak4NHadEF_[ij],ak4CMulti_[ij],ak4NMultiplicity_[ij],ak4eta[ij])[0] for ij in range(nak4jet_)]
                ak4PassID_Calc = ak4PassID_
            if runOn2017:
                #ak4PassID_Calc = [jetID_(ak4CEmEF_[ij],ak4CHadEF_[ij],ak4NEmEF_[ij],ak4NHadEF_[ij],ak4CMulti_[ij],ak4NMultiplicity_[ij],ak4eta[ij])[1] for ij in range(nak4jet_)]
                ak4PassID_Calc = ak4PassID_
            if runOn2018:
                #ak4PassID_Calc = [jetID_(ak4CEmEF_[ij],ak4CHadEF_[ij],ak4NEmEF_[ij],ak4NHadEF_[ij],ak4CMulti_[ij],ak4NMultiplicity_[ij],ak4eta[ij])[1] for ij in range(nak4jet_)]
                ak4PassID_Calc = ak4PassID_

            #ak4_pt30_eta4p5_IDT  =  ( (ak4pt[ij] > 30.0) and (abs(ak4eta[ij]) < 4.5) and (ak4PassID_Calc[ij] ) ) for ij in range(nak4jet_)]
            ak4_pt30_eta4p5_IDT = boolutil.logical_and3(
                (ak4pt > 30.0), (numpy.abs(ak4eta) < 4.5), (ak4PassID_Calc))

            ##--- jet cleaning
            jetCleanAgainstEle = []
            jetCleanAgainstMu = []
            pass_jet_index_cleaned = []


            if len(ak4_pt30_eta4p5_IDT) > 0:
                DRCut = 0.4
                jetCleanAgainstEle = anautil.jetcleaning(
                    ak4_pt30_eta4p5_IDT, ele_pt10_eta2p5_vetoID, ak4eta, eleeta, ak4phi, elephi, DRCut)
                jetCleanAgainstMu = anautil.jetcleaning(
                    ak4_pt30_eta4p5_IDT, mu_pt10_eta2p4_looseID_looseISO, ak4eta, mueta, ak4phi, muphi, DRCut)
                jetCleaned = boolutil.logical_AND_List3(
                    ak4_pt30_eta4p5_IDT, jetCleanAgainstEle, jetCleanAgainstMu)
                pass_jet_index_cleaned = boolutil.WhereIsTrue(jetCleaned)
                if debug_:
                    print "pass_jet_index_cleaned = ", pass_jet_index_cleaned, "nJets= ", len(ak4px_)

            '''
            ********    *        *       *
               *      *    *     *       *
               *     *      *    *       *
               *     ********    *       *
               *     *      *    *       *
               *     *      *     *******
            '''
            taupt = getPt(tau_px_, tau_py_)
            taueta = getEta(tau_px_, tau_py_, tau_pz_)
            tauphi = getPhi(tau_px_, tau_py_)

            tau_eta2p3_iDLdm_pt18 = boolutil.logical_AND_List4(
                (taupt > 18.0), (numpy.abs(taueta) < 2.3), (tau_isLoose_), (tau_dm_))

            if debug_:
                print "tau_eta2p3_iDLdm_pt18 = ", tau_eta2p3_iDLdm_pt18

            tau_eta2p3_iDLdm_pt18_looseEleVeto_looseMuVeto = boolutil.logical_and6((taupt > 18.0), (numpy.abs(
                taueta) < 2.3), (tau_isLoose_), (tau_dm_), (Taudisc_againstLooseElectron), (Taudisc_againstLooseMuon))
            tau_eta2p3_iDLdm_pt18_looseEleVeto_tightMuVeto = boolutil.logical_and6((taupt > 18.0), (numpy.abs(
                taueta) < 2.3), (tau_isLoose_), (tau_dm_), (Taudisc_againstLooseElectron), (Taudisc_againstTightMuon))
            tau_eta2p3_iDLdm_pt18_mediumEleVeto_looseMuVeto = boolutil.logical_and6((taupt > 18.0), (numpy.abs(
                taueta) < 2.3), (tau_isLoose_), (tau_dm_), (Taudisc_againstMediumElectron), (Taudisc_againstLooseMuon))
            tau_eta2p3_iDLdm_pt18_tightEleVeto_looseMuVeto = boolutil.logical_and6((taupt > 18.0), (numpy.abs(
                taueta) < 2.3), (tau_isLoose_), (tau_dm_), (Taudisc_againstTightElectron), (Taudisc_againstLooseMuon))
            tau_eta2p3_iDLdm_pt18_tightEleVeto_tightMuVeto = boolutil.logical_and6((taupt > 18.0), (numpy.abs(
                taueta) < 2.3), (tau_isLoose_), (tau_dm_), (Taudisc_againstTightElectron), (Taudisc_againstTightMuon))

            tau_eta2p3_iDLdm_pt18_looseEleVeto_looseMuVeto_index = boolutil.WhereIsTrue(
                tau_eta2p3_iDLdm_pt18_looseEleVeto_looseMuVeto)
            tau_eta2p3_iDLdm_pt18_looseEleVeto_tightMuVeto_index = boolutil.WhereIsTrue(
                tau_eta2p3_iDLdm_pt18_looseEleVeto_tightMuVeto)
            tau_eta2p3_iDLdm_pt18_mediumEleVeto_looseMuVeto_index = boolutil.WhereIsTrue(
                tau_eta2p3_iDLdm_pt18_mediumEleVeto_looseMuVeto)
            tau_eta2p3_iDLdm_pt18_tightEleVeto_looseMuVeto_index = boolutil.WhereIsTrue(
                tau_eta2p3_iDLdm_pt18_tightEleVeto_looseMuVeto)
            tau_eta2p3_iDLdm_pt18_tightEleVeto_tightMuVeto_index = boolutil.WhereIsTrue(
                tau_eta2p3_iDLdm_pt18_tightEleVeto_tightMuVeto)

            tauCleanAgainstEle = []
            tauCleanAgainstMu = []
            pass_tau_index_cleaned_DRBased = []
            if len(tau_eta2p3_iDLdm_pt18) > 0:
                DRCut = 0.4
                tauCleanAgainstEle = anautil.jetcleaning(
                    tau_eta2p3_iDLdm_pt18, ele_pt10_eta2p5_looseID,         taueta, eleeta, tauphi, elephi, DRCut)
                tauCleanAgainstMu = anautil.jetcleaning(
                    tau_eta2p3_iDLdm_pt18, mu_pt10_eta2p4_looseID_looseISO, taueta, mueta,  tauphi, muphi,  DRCut)
                tauCleaned = boolutil.logical_AND_List3(
                    tau_eta2p3_iDLdm_pt18, tauCleanAgainstEle, tauCleanAgainstMu)
                pass_tau_index_cleaned_DRBased = boolutil.WhereIsTrue(
                    tauCleaned)
                if debug_:
                    print "pass_tau_index_cleaned_DRBased", pass_tau_index_cleaned_DRBased

            for ithinjet in pass_jet_index_cleaned:
                #==================== for flavor 5 ==================================
                if ak4flavor_[ithinjet]==5:
                    h_btag_den.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > LWP:
                        h_btag_num_pass_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_btag_num_fail_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > MWP:
                        h_btag_num_pass_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_btag_num_fail_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])

                #==================== for flavor 4 ==================================
                if ak4flavor_[ithinjet]==4:
                    h_ctag_den.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > LWP:
                        h_ctag_num_pass_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_ctag_num_fail_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > MWP:
                        h_ctag_num_pass_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_ctag_num_fail_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])

                #==================== for light flavor  ==================================

                if ak4flavor_[ithinjet]!=4 and ak4flavor_[ithinjet]!=5:
                    h_lighttag_den.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > LWP:
                        h_lighttag_num_pass_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_lighttag_num_fail_lwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    if ak4deepcsv_[ithinjet] > MWP:
                        h_lighttag_num_pass_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
                    else:
                        h_lighttag_num_fail_mwp.Fill(ak4eta[ithinjet],ak4pt[ithinjet])
    outfile.cd()
    h_btag_num_pass_mwp.Write()
    h_btag_num_fail_mwp.Write()
    h_btag_num_pass_lwp.Write()
    h_btag_num_fail_lwp.Write()
    h_btag_den.Write()
    h_ctag_num_pass_mwp.Write()
    h_ctag_num_fail_mwp.Write()
    h_ctag_num_pass_lwp.Write()
    h_ctag_num_fail_lwp.Write()
    h_ctag_den.Write()
    h_lighttag_num_pass_mwp.Write()
    h_lighttag_num_fail_mwp.Write()
    h_lighttag_num_pass_lwp.Write()
    h_lighttag_num_fail_lwp.Write()
    h_lighttag_den.Write()
    outfile.Write()
    print "output written to ", outfilename
    end = time.clock()
    print "%.4gs" % (end-start)

if __name__ == '__main__':
    if not runInteractive:
        txtFile = infile

        runbbdm(txtFile)

    if runInteractive and runOnTxt:
	filesPath = dirName+'/*txt'
	files = glob.glob(filesPath)
        n = 1  # submit n txt files at a time, make equal to cores
        final = [files[i * n:(i + 1) * n]
                 for i in range((len(files) + n - 1) // n)]
        print 'final', final
        for i in range(len(final)):
            print 'first set', final[i]

            try:
                pool = mp.Pool(1)
                pool.map(runbbdm, final[i])
                pool.close()
                pool.join()
	    except Exception as e:
		print e
		print "Corrupt file inside input txt file is detected! Skipping this txt file:  ", final[i]
		continue

    if runInteractive and not runOnTxt:
        ''' following part is for interactive running. This is still under testing because output file name can't be changed at this moment '''
        inputpath = "/eos/cms/store/group/phys_exotica/bbMET/ExoPieElementTuples/MC_2017miniaodV2_V1/"

        os.system('rm dirlist.txt')
        os.system("ls -1 "+inputpath+" > dirlist.txt")

        allkeys = [idir.rstrip() for idir in open('dirlist.txt')]
        alldirs = [inputpath+"/"+idir.rstrip() for idir in open('dirlist.txt')]

        pool = mp.Pool(6)
        allsample = []
        for ikey in allkeys:
            dirpath = inputpath+"/"+ikey
            txtfile = ikey+".txt"
            os.system("find "+dirpath +
                      "  -name \"*.root\" | grep -v \"failed\"  > "+txtfile)
            fileList = TextToList(txtfile)
            ## this is the list, first element is txt file with all the files and second element is the ikey (kind of sample name identifier)
            sample_ = [txtfile, ikey]
            ## push information about one sample into global list.
            allsample.append(sample_)
        print allsample
        pool.map(runbbdm, allsample)
        ## this works fine but the output file name get same value becuase it is done via a text file at the moment, need to find a better way,
