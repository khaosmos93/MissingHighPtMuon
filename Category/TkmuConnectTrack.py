import numpy
import sys
import ROOT
import math
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Events, Handle
from ROOT import TVector3,TLorentzVector,TFile,TH1D,TRandom3
DeltaR = ROOT.Math.VectorUtil.DeltaR
DeltaPhi = ROOT.Math.VectorUtil.DeltaPhi
DeltaR2 = lambda a, b: DeltaR(a.p4(), b.p4())       # for reco::Candidates
DeltaPhi2 = lambda a, b: DeltaPhi(a.p4(), b.p4())   # for reco::Candidates

dR = lambda a, b : math.sqrt( (a.eta() - b.eta())**2 + (a.phi()-b.phi())**2 )


Muhandle = Handle ("std::vector<reco::Muon>")
Mulabel = ("muons","","RECO")
Tkhandle = Handle ("std::vector<reco::Track>")
Tklabel = ("generalTracks","","RECO")
METhandle = Handle ("std::vector<reco::PFMET>")
METlabel = ("pfMet","","RECO")

vertices, vertexLabel = Handle("std::vector<reco::Vertex>"), "offlinePrimaryVertices"


files = [
          ('RunBv2','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunBv2.root'),
          ('RunC','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunC.root'),
          ('RunD','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunD.root'),
          ('RunE','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunE.root'),
          ('RunF','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunF.root'),
          ('RunG','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunG.root'),
          ('RunHv2','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunHv2.root'),
          ('RunHv3','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_RunHv3.root'),

          #('DY400_800','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_DY400_800.root'),
          #('DY800_1400','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_DY800_1400.root'),
          #('DY1400_2300','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_DY1400_2300.root'),
          #('DY2300_3500','file:/cms/home/msoh/MissingHighPtMuon/Files/TagMuonSelectionEvents_DY2300_3500.root'),

          #'file:/cms/home/msoh/MissingHighPtMuon/PE/Files/Mu600/SelectedMuons_RunHv3.root',
          #'file:/cms/home/msoh/MissingHighPtMuon/PE/Files/Mu600/SelectedMuons_RunHv2.root',
          #'file:/cms/home/msoh/MissingHighPtMuon/PE/Files/Mu600/SelectedMuons_DY1400_2300.root',
        ]

#vervos = False
vervos = True
#printKinematics = False
printKinematics = True
OnlyBarrel = False
MaxEvPerFile = 10000000

pi = 3.141592
dPhim = pi - 0.3*pi
dPhiM = pi + 0.3*pi

print "Start!"

for name, fileName in files:
  print "\n\n\n",name," : ", fileName
  if not printKinematics: fn='TrackerMuonTrackConnectEvents_'+name+'.txt'
  elif printKinematics: fn='TrackerMuonTrackConnectKinematics_'+name+'.txt'
  fout = open(fn, 'wt')
  eventCount=0
  events = Events(fileName)
  for event in events:
    if eventCount > MaxEvPerFile: break
    eventCount = eventCount + 1
    if (vervos): print "\n eventCount : ", eventCount


    try: event.getByLabel (Mulabel, Muhandle)
    except RuntimeError: print "No muon info"
    muons=Muhandle.product()

    try: event.getByLabel (Tklabel, Tkhandle)
    except RuntimeError: print "No Track info"
    tracks=Tkhandle.product()

    try: event.getByLabel (METlabel, METhandle)
    except RuntimeError: print "No MET info"
    MET=METhandle.product()

    try: event.getByLabel(vertexLabel, vertices)
    except RuntimeError: print "No Vertex info"
    vtxs=vertices.product()

    ### Vertices ###
    numGoodVTXs=0
    goodVTXs=[]
    for vtx in vtxs:
      if ( vtx.ndof() > 4 and
           abs(vtx.z() <= 24) and
           abs(vtx.position().rho()) <=2 ):
        numGoodVTXs+=1
        goodVTXs.append(vtx)
    if (vervos): print "number of good VTX in event: ",numGoodVTXs
    if numGoodVTXs==0: continue
    else: PV=goodVTXs[0]

    #if len(vertices.product()) == 0 or vertices.product()[0].ndof() < 4: continue
    #else: PV = vertices.product()[0]

    ### Tag Muons ###
    numMuPassed=0
    selectedMuons=[]
    for mu in muons:
      if mu.pt() < 53 : continue
      if not (mu.isGlobalMuon() and mu.isTrackerMuon()): continue
      #if abs(mu.eta())>2.4: continue
      #if abs(mu.dB())>0.2: continue
      if not (mu.muonBestTrack().isNonnull()): continue
      if not (mu.globalTrack().isNonnull()): continue
      if not (mu.innerTrack().isNonnull()): continue
      if not(abs(mu.muonBestTrack().dxy(PV.position())) < 0.2): continue
      if not((mu.muonBestTrack().ptError()/mu.muonBestTrack().pt()) < 0.3): continue
      if not((mu.isolationR03().sumPt / mu.innerTrack().pt()) < 0.10): continue
      if not(mu.globalTrack().hitPattern().trackerLayersWithMeasurement() > 5): continue
      if not(mu.globalTrack().hitPattern().numberOfValidPixelHits() > 0): continue
      if not(mu.globalTrack().hitPattern().numberOfValidMuonHits() > 0): continue
      if not( mu.numberOfMatchedStations() > 1 or
             (mu.numberOfMatchedStations() == 1 and not(mu.stationMask() == 1 or mu.stationMask() == 16)) or
             (mu.numberOfMatchedStations() == 1 and (mu.stationMask() == 1 or mu.stationMask() == 16) and mu.numberOfMatchedRPCLayers()>2)
             ): continue
      numMuPassed+=1
      selectedMuons.append(mu)
    if (vervos): print "number of muons passing Z' selection in event: ",numMuPassed

    if numMuPassed != 1:
      if (vervos): print "not one high pT Z' muon in event"
      continue

    isBarrel = True
    if selectedMuons[0].pt() < 600:
      if (vervos): print "Tag pT < 600!"
      continue
    if abs(selectedMuons[0].eta()) > 0.9:
      #if (vervos): print "Tag muon is not in barrel!"
      isBarrel = False
      #if OnlyBarrel: continue
    #if (OnlyBarrel and not isBarrel): continue


    ### No another muon ###
    numAnotherMuon=0
    #AnotherMuons=[]
    for muon in muons:
      if (muon.pt() == selectedMuons[0].pt() and
          muon.eta() == selectedMuons[0].eta() and
          muon.phi() == selectedMuons[0].phi()
         ): continue
      if not muon.isGlobalMuon(): continue
      if not muon.pt() > 53: continue
      if not( muon.numberOfMatchedStations() > 1 or
             (muon.numberOfMatchedStations() == 1 and not(muon.stationMask() == 1 or muon.stationMask() == 16)) or
             (muon.numberOfMatchedStations() == 1 and (muon.stationMask() == 1 or muon.stationMask() == 16) and muon.numberOfMatchedRPCLayers()>2)
             ): continue
      numAnotherMuon+=1
      #AnotherMuons.append(muon)
    if (vervos): print "number of another muon in event: ",numAnotherMuon
    if numAnotherMuon > 0: continue


    ### Probe Tracks ###
    numTkPassed=0
    selectedTracks=[]
    for tk in tracks:
      if tk.pt() < 10: continue  # small pT cut on probe track
      if abs(tk.eta()) > 0.9: continue  # track in Barrel
      if not ( abs( tk.phi()-selectedMuons[0].phi() ) < dPhiM and abs( tk.phi()-selectedMuons[0].phi() ) > dPhim): continue  # b2b in phi
      if not ( tk.hitPattern().trackerLayersWithMeasurement() > 5 ): continue
      if not ( tk.hitPattern().numberOfValidPixelHits() > 0 ): continue
      numTkPassed+=1
      selectedTracks.append(tk)

    if (vervos): print "number of tracks passing pre-selection in event: ",numTkPassed
    if numTkPassed < 1: continue


    ### Category : isTracker muon connected to the general track ###
    numTkCategory=0
    categoryTracks=[]
    for selTk in selectedTracks:
      numTkMuMatch = 0
      for muu in muons:
        if not muu.isTrackerMuon(): continue
        if ( muu.innerTrack().pt() == selTk.pt() and
             muu.innerTrack().eta() == selTk.eta() and
             muu.innerTrack().phi() == selTk.phi()
            ): numTkMuMatch+=1
      if numTkMuMatch>0:
        numTkCategory+=1
        categoryTracks.append(selTk)

    if (vervos): print "number of tracks passing category selection in event: ",numTkCategory
    if numTkCategory < 1: continue



    if (printKinematics):
      for selMu in selectedMuons:
        muPt = selMu.pt()
        muEta = selMu.eta()
        muPhi = selMu.phi()
        muDxy = selMu.muonBestTrack().dxy(PV.position())
        muDz = selMu.muonBestTrack().dz(PV.position())
        MuonInfo = '''
Muon :    pT= %(muPt)s    eta= %(muEta)s    phi= %(muPhi)s    isBarrel= %(isBarrel)s    dxy= %(muDxy)s    dz= %(muDz)s'''
        fout.write(MuonInfo % locals())
        #print "muon : pT=", selMu.pt(), "\t eta=", selMu.eta(), "\t phi=", selMu.phi(), "\t isBarrel : ", isBarrel, "\tdxy=", selMu.muonBestTrack().dxy(PV.position()), "\tdz=", selMu.muonBestTrack().dz(PV.position())

      for catTk in categoryTracks:
        selectedTkSumPtOverPt=-999
        TkSumPt = 0
        for tki in tracks:
          if (tki.pt() == catTk.pt() and
              tki.eta() == catTk.eta() and
              tki.phi() == catTk.phi()
              ): continue
          if abs(tki.dz(PV.position())) > 0.2: continue
          if dR(tki,catTk) < 0.3:
            TkSumPt += tki.pt()
        selectedTkSumPtOverPt=TkSumPt/catTk.pt()

        tkPt = catTk.pt()
        tkEta = catTk.eta()
        tkPhi = catTk.phi()
        tkDxy = catTk.dxy(PV.position())
        tkDz = catTk.dz(PV.position())
        tkDPhi = abs( catTk.phi()-selectedMuons[0].phi() )
        TrackInfo = '''
Track :   pT= %(tkPt)s    eta= %(tkEta)s    phi= %(tkPhi)s
          dxy= %(tkDxy)s,    dz= %(tkDz)s,    dPhi= %(tkDPhi)s,    SumPt/Pt= %(selectedTkSumPtOverPt)s'''
        fout.write(TrackInfo % locals())
        #print "catTk :", "\tpT=", catTk.pt(), "\teta=", catTk.eta(), "\tphi=", catTk.phi(), "\n\t", "\tdxy=", catTk.dxy(PV.position()), "\tdz=", catTk.dz(PV.position()), "\tdPhi = ", abs( catTk.phi()-selectedMuons[0].phi() ), "\t SumPt/Pt : ", selectedTkSumPtOverPt

      for met in MET:
        metPt = met.pt()
        metEta = met.eta()
        metPhi = met.phi()
        METInfo = '''
MET :     pT= %(metPt)s,    eta= %(metEta)s,    phi= %(metPhi)s'''
        fout.write(METInfo % locals() + '\n' )
        #print "MET : pt=", met.pt(), "\t eta=", met.eta(), "\t phi=", met.phi()

    run = event.eventAuxiliary().run()
    lumi = event.eventAuxiliary().luminosityBlock()
    ev = event.eventAuxiliary().event()
    EventID = '%(run)s:%(lumi)s:%(ev)s'  ###   [%(eventCount)s]'
    fout.write(EventID % locals() + '\n')
    #print "Run:Lumi:Event = ", event.eventAuxiliary().run(), ":", event.eventAuxiliary().luminosityBlock(), ":", event.eventAuxiliary().event()  , "\t [", eventCount, "]"
    if (printKinematics): fout.write('\n')

  fout.close()

print "finished!"



