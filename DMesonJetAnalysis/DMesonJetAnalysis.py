#!/usr/bin/env python
#python program to perform a D meson jet analysis

import ROOT
import math
import DMesonJetProjectors

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, projector):
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        
    def Start(self):
        self.fProjector.GetInvMassHisograms(self.fTrigger, self.fDMeson, self.fJetDefinitions, 
                                            self.fBinSet, self.fNMassBins, self.fMinMass, self.fMaxMass)
        self.PlotInvMassPlots()
        
    def PlotInvMassPlots(self):
        for name,bins in self.fBinSet.fBins.iteritems():
            self.PlotInvMassPlotsBinSet(name,bins)
    
    def GenerateInvMassCanvas(self, name, n):
        rows = int(math.floor(math.sqrt(n)))
        cols = int(math.ceil(float(n) / rows))
        c = ROOT.TCanvas(name, name, cols*400, rows*400)
        c.Divide(cols, rows)
        globalList.append(c)
        return c
    
    def PlotInvMassPlotsBinSet(self, name, bins):
        c = self.GenerateInvMassCanvas(name, len(bins))
        for i,bin in enumerate(bins):
            if bin.fInvMassHisto:
                c.cd(i+1)
                bin.fInvMassHisto.Draw()

class DMesonJetAnalysis:
    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, config):
        binSet = DMesonJetProjectors.BinSet()
        for binLists in config["bins"]:
            binSet.AddBins(name = binLists["name"], jetPtLimits = binLists["jet_pt"], DPtLimits = binLists["d_pt"], ZLimits = binLists["d_z"])
         
        eng = DMesonJetAnalysisEngine(config["trigger"], config["d_meson"], 
                                     binSet, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                     config["jets"], self.fProjector)
        self.fAnalysisEngine.append(eng)
        eng.Start()
        