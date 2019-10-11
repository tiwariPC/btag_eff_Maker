import numpy

def getPt(px_, py_):
    return  numpy.sqrt( px_**2 + py_**2)


def GenPtProducer(sample,nGenPar, genParId, genMomParId, genParSt,st_genParPx,st_genParPy):
    pt_list=[]
    #################
    # WJets
    #################
    if sample==24:
        goodLepID = []
        for ig in range(nGenPar):
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]
            if ( (abs(PID) != 11) & (abs(PID) != 12) &  (abs(PID) != 13) & (abs(PID) != 14) &  (abs(PID) != 15) &  (abs(PID) != 16) ): continue
            if ( ( (status != 1) & (abs(PID) != 15)) | ( (status != 2) & (abs(PID) == 15)) ): continue
            if ( (abs(momPID) != 24) & (momPID != PID) ): continue
            goodLepID.append(ig)

        if len(goodLepID) == 2 :
            pt_list.append(getPt((st_genParPx[goodLepID[0]]+st_genParPx[goodLepID[1]]),(st_genParPy[goodLepID[0]]+st_genParPy[goodLepID[1]])))
    #################
    #ZJets
    #################
    if sample == 23:
        goodLepID = []
        for ig in range(nGenPar):
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]


            if ( (abs(PID) != 12) &  (abs(PID) != 14) &  (abs(PID) != 16) ) : continue
            if ( status != 1 ) : continue
            if ( (momPID != 23) & (momPID != PID) ) : continue
            goodLepID.append(ig)

        if len(goodLepID) == 2 :
            pt_list.append(getPt((st_genParPx[goodLepID[0]]+st_genParPx[goodLepID[1]]),(st_genParPy[goodLepID[0]]+st_genParPy[goodLepID[1]])))
    #################
    #TTBar
    #################
    if (sample==6):
        goodLepID = []
        for ig in range(nGenPar):
            PID    = genParId[ig]
            momPID = genMomParId[ig]
            status = genParSt[ig]
            if ( abs(PID) == 6) :
                goodLepID.append(ig)
        if(len(goodLepID)==2):
            pt_list.append(getPt(st_genParPx[goodLepID[0]],st_genParPy[goodLepID[0]]))
            pt_list.append(getPt(st_genParPx[goodLepID[1]],st_genParPy[goodLepID[1]]))

    return pt_list
