#!/usr/bin/env python
#python script to do EMCal jet trigger turn on curves

import ROOT
import argparse
import IPython
import copy
import math

globalList = []
globalCanvasList = []

def SaveCanvases(path, label):
    for canvas in globalCanvasList:
        fname = "{0}/{1}{2}.pdf".format(path, label, canvas.GetName())
        canvas.SaveAs(fname)

class JetDefinition:
    def __init__(self, _type, _radius):
        self.fType = _type
        self.fRadius = _radius
    
    def GetBranchName(self):
        if self.fType == "Charged":
            return "Jet_AKTCharged{0}_tracks_pT0150_pt_scheme".format(self.fRadius)
        elif self.fType == "Full":
            return "Jet_AKTFull{0}_tracks_pT0150_caloClusters_E0300_pt_scheme".format(self.fRadius)
    
    def GetName(self):
        return "{0}_{1}".format(self.fType, self.fRadius)
    
    def GetTitle(self):
        return "{0} {1}".format(self.fType, self.fRadius)
        
    def AddObservable(self, observable):
        if not hasattr(self, "fObservableList"):
            self.fObservableList = []
        self.fObservableList.append(observable)
        
    def GenerateCanvasList(self, label="", logY=False):
        canvasList = []
        
        for obs in self.fObservableList:
            if label:
                cname = "{0}_{1}_{2}".format(label, self.GetName(), obs.fName)
                ctitle = "{0} {1} {2}".format(label, self.GetTitle(), obs.fTitle)
            else:
                cname = "{0}_{1}".format(self.GetName(), obs.fName)
                ctitle = "{0} {1}".format(self.GetTitle(), obs.fTitle)
            print("Generating canvas {0} for observable {1}".format(cname, obs.fName))
            canvas = ROOT.TCanvas(cname, ctitle)
            if logY:
                canvas.SetLogy(obs.fLogy)
            canvasList.append(canvas)
            globalCanvasList.append(canvas)
            
        return canvasList
    
    def DrawLegends(self, legend, ratioLegend):
        for canvas in self.fCanvasList:
            canvas.cd()
            legend.Draw()
        for canvas in self.fRatioCanvasList:
            canvas.cd()
            ratioLegend.Draw()

    def DrawHisto(self, canvasList, histoName, opt, marker, color):
        for obs,canvas in zip(self.fObservableList, canvasList):
            canvas.cd()
            
            if not hasattr(obs, histoName):
                continue
            
            histo = getattr(obs, histoName)
            
            resetMax = True
            for histInList in canvas.GetListOfPrimitives():
                if not isinstance(histInList, ROOT.TH1):
                    continue
                if histInList.GetMaximum() > histo.GetMaximum():
                    resetMax = False
                    break
            if resetMax:
                firstHist = canvas.GetListOfPrimitives().At(0)
                if hasattr(firstHist, "SetMaximum"):
                    if obs.fLogy:
                        firstHist.SetMaximum(histo.GetMaximum()*5)
                    else:
                        firstHist.SetMaximum(histo.GetMaximum()*1.2)
            if not obs.fLogy:
                histo.SetMinimum(0)
            histo.SetMarkerSize(1)
            histo.SetMarkerStyle(marker)
            histo.SetMarkerColor(color)
            histo.SetLineColor(color)
            histo.Draw(opt)

class Observable:
    def __init__(self, name, title, bins, min, max, logY=False, norm=False):
        self.fName = name
        self.fTitle = title
        self.fBins = bins
        self.fMin = min
        self.fMax = max
        self.fHistogram = None
        self.fLogy = logY
        self.fNorm = norm
        
    def CreateHistogram(self, hname, htitle, xtitle="", ytitle=""):
        hname = "{0}_{1}".format(hname, self.fName)
        if not xtitle:
            xtitle = self.fTitle
        if not ytitle:
            ytitle = "entries"
        htitle = "{0} {1};{2};{3}".format(htitle, self.fTitle, xtitle, ytitle)
        return ROOT.TH1F(hname, htitle, self.fBins, self.fMin, self.fMax)
    
    def AddCut(self, obs):
        if not hasattr(self, "fListOfCuts"):
            self.fListOfCuts = []
        self.fListOfCuts.append(copy.deepcopy(obs))
        
class TreeProjector:    
    def __init__(self, name, title, triggerScale, path, subdirs, filename, trigger, color, marker):
        self.CreateChain(path, subdirs, filename, trigger)
        self.fName = name
        self.fTitle = title
        self.fColor = color
        self.fMarker = marker
        self.fTriggerScale = triggerScale
        self.fProjectionOk = False
        
    def AddJetDefinition(self, jetDef):
        if not hasattr(self, "fJetDefinitionList"):
            self.fJetDefinitionList = []
        self.fJetDefinitionList.append(jetDef)
    
    def CreateChain(self, path, subdirs, filename, trigger):
        treeName = "AliAnalysisTaskEmcalJetTree_{0}_jets".format(trigger)
        self.fChain = ROOT.TChain(treeName)
    
        for dir in subdirs:
            fullFileName = "{0}/{1}/{2}".format(path, dir, filename)
            print("Adding {0}".format(fullFileName))
            self.fChain.Add(fullFileName)
            
    def ProjectTree(self, maxEntries=0):
        print("Projecting {0}, number of events: {1}".format(self.fName, self.fChain.GetEntries()))
        
        for jetDef in self.fJetDefinitionList:
            for observable in jetDef.fObservableList:
                observable.fHistogram = observable.CreateHistogram("{0}_{1}".format(self.fName, jetDef.GetName()), 
                                                                   "{0} {1}".format(self.fTitle, jetDef.GetTitle()))
                observable.fHistogram.Sumw2()

        for i, event in enumerate(self.fChain):
            if maxEntries>0 and i >= maxEntries:
                break
            if i % 100000 == 0:
                print("Entry {0}".format(i))
            for jetDef in self.fJetDefinitionList:
                jets = getattr(event, jetDef.GetBranchName())
                for jet in jets:
                    for observable in jetDef.fObservableList:
                        skip = False
                        if hasattr(observable, "fListOfCuts"):
                            for cut in observable.fListOfCuts:
                                val = getattr(jet, cut.fName)
                                if val > cut.fMax or val < cut.fMin:
                                    skip = True
                                    break
                        if skip:
                            continue
                        val = getattr(jet, observable.fName)
                        #print dir(jet)
                        #print("Value {0} for observable {1}".format(val, observable.fName))
                        observable.fHistogram.Fill(val)
                        
        for jetDef in self.fJetDefinitionList:
            for observable in jetDef.fObservableList:               
                if observable.fNorm:
                    integral = observable.fHistogram.Integral()
                    if integral > 1e-10:
                        observable.fHistogram.Scale(1. / integral)
                        
        self.fProjectionOk = True
            
    def GenerateTriggerTurnOn(self, maxEntries, treeMB):
        if not self.fProjectionOk:
            self.ProjectTree(maxEntries)
        
        if not treeMB:
            return
        
        if not treeMB.fProjectionOk:
           treeMB.ProjectTree(maxEntries) 
        
        for jetDefRare, jetDefMB in zip(self.fJetDefinitionList, treeMB.fJetDefinitionList):
            for obsRare, obsMB in zip(jetDefRare.fObservableList, jetDefMB.fObservableList):
                obsRare.fRatio = obsMB.CreateHistogram("TurnOn_{0}_{1}".format(jetDefRare.GetName(), self.fName), 
                                                       "TurnOn {0} {1}".format(jetDefRare.GetTitle(), self.fTitle),
                                                       "", "ratio")
                obsRare.fRatio.Divide(obsRare.fHistogram, obsMB.fHistogram)
            
    def Draw(self, doRatios=False, treeProj=None, legend=None, ratioLegend=None):
        if treeProj:
            for jetDef1, jetDef2 in zip(self.fJetDefinitionList, treeProj.fJetDefinitionList):
                if hasattr(jetDef2, "fCanvasList") and not hasattr(jetDef1, "fCanvasList"):
                    jetDef1.fCanvasList = jetDef2.fCanvasList
                if hasattr(jetDef2, "fRatioCanvasList") and not hasattr(jetDef1, "fRatioCanvasList"):
                    jetDef1.fRatioCanvasList = jetDef2.fRatioCanvasList

        for jetDef in self.fJetDefinitionList:
            print("Drawing {0} {1}...".format(self.fTitle, jetDef.GetTitle()))
            if hasattr(jetDef, "fCanvasList") and jetDef.fCanvasList:
                print("Using pre-existing canvas list...")
                opt = "same"
            else:
                print("Generating canvas list...")
                jetDef.fCanvasList = jetDef.GenerateCanvasList("", True)
                opt = ""          
        
            jetDef.DrawHisto(jetDef.fCanvasList, "fHistogram", opt, self.fMarker, self.fColor)
        
            if doRatios:
                if hasattr(jetDef, "fRatioCanvasList") and jetDef.fRatioCanvasList:
                    opt = "same"
                else:
                    jetDef.fRatioCanvasList = jetDef.GenerateCanvasList("Ratio", False)
                    opt = ""
            
                jetDef.DrawHisto(jetDef.fRatioCanvasList, "fRatio", opt, self.fMarker, self.fColor)
            
            jetDef.DrawLegends(legend, ratioLegend)
    
    def SaveRoot(self, file):
        file.cd()
        for jetDef in self.fJetDefinitionList:
            for observable in jetDef.fObservableList:
                if hasattr(observable, "fHistogram"):
                    observable.fHistogram.Write()
                if hasattr(observable, "fRatio"):
                    observable.fRatio.Write()

class TriggerTurnOn:
    def __init__(self, triggerList, triggerScaleList, path, subdirs, filename, entries):
        self.fListOfTreeProjectors = dict()
        colorList = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+1, ROOT.kOrange+1]
        markerList = [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenDiamond, ROOT.kOpenCross]
        for trigger,triggerScale,color,marker in zip(triggerList,triggerScaleList,colorList,markerList):
            if triggerScale:
                print("Adding trigger {0} which will be scaled using {1}".format(trigger,triggerScale))
                triggerScaleObj = self.fListOfTreeProjectors[triggerScale]
            else:
                print("Adding trigger {0} which will not be scaled".format(trigger))
                triggerScaleObj = None
            self.fListOfTreeProjectors[trigger] = TreeProjector(trigger, trigger, triggerScaleObj, path, subdirs, filename, trigger, color, marker)
            print("Number of triggers is now {0}".format(len(self.fListOfTreeProjectors)))
        self.fEntries = entries
        
    def AddJetDefinition(self, jetDef):
        for proj in self.fListOfTreeProjectors.itervalues():
            proj.AddJetDefinition(copy.deepcopy(jetDef))

    def GenerateTriggerTurnOn(self):
        for proj in self.fListOfTreeProjectors.itervalues():
            proj.GenerateTriggerTurnOn(self.fEntries, proj.fTriggerScale)
    
    def GenerateLegends(self, projList):
        legend = ROOT.TLegend(0.5, 0.9, 0.9, 0.7)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(43)
        legend.SetTextSize(12)
        for proj in projList:
            entry = legend.AddEntry(None, proj.fTitle, "pe")
            entry.SetMarkerSize(1)
            entry.SetMarkerColor(proj.fColor)
            entry.SetMarkerStyle(proj.fMarker)
            entry.SetLineColor(proj.fColor)
        return legend
    
    def GenerateRatioLegends(self, projList):
        legend = ROOT.TLegend(0.5, 0.9, 0.9, 0.7)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.SetTextFont(43)
        legend.SetTextSize(12)
        for proj in projList:
            if not proj.fTriggerScale: 
                continue
            entry = legend.AddEntry(None, "{0}/{1}".format(proj.fTitle, proj.fTriggerScale.fTitle), "pe")
            entry.SetMarkerSize(1)
            entry.SetMarkerColor(proj.fColor)
            entry.SetMarkerStyle(proj.fMarker)
            entry.SetLineColor(proj.fColor)
        return legend
            
    def Draw(self):
        self.fLegend = self.GenerateLegends(self.fListOfTreeProjectors.values())
        self.fRatioLegend = self.GenerateRatioLegends(self.fListOfTreeProjectors.values())
        
        prev = None
        for proj in self.fListOfTreeProjectors.itervalues():
            proj.Draw(True, prev, self.fLegend, self.fRatioLegend)
            prev = proj
        
    def SaveAll(self, path, label=""):
        file = ROOT.TFile("{0}/{1}TriggerTurnOnAnalysis.root".format(path, label), "recreate")
        for proj in self.fListOfTreeProjectors.itervalues():
            proj.SaveRoot(file)
            
        file.Close()

def AddCuts(obs, ptcut):
    if ptcut > 0:
        obs.AddCut(Observable("fPt","#it{p}_{T} (GeV/#it{c})", 10, ptcut, 1000, True, False))

def main(train, period, trainNumbers, inputPath, filename, triggers, entries, ptcut):
    ROOT.gSystem.Load("libCGAL")
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)
    
    TrainNumberList = trainNumbers.split(",")
    TrainNumbers = []
    for TrainNumberRange in TrainNumberList:
        Range = TrainNumberRange.split(":")
        TrainNumbers.extend(range(int(Range[0]), int(Range[len(Range)-1])+1))
    
    for TrainNumber in TrainNumbers:
        train += "_{0}".format(TrainNumber)
    
    periodList = []
    
    if period == "LHC12x":
        periodList = ["LHC12a", "LHC12b", "LHC12c", "LHC12d", "LHC12e", 
                      "LHC12f", "LHC12g", "LHC12h", "LHC12i"]
    else:
        periodList.append(period)
        
    print periodList
    
    path = "{0}/{1}".format(inputPath, train)
    
    triggerList = []
    triggerScaleList = []
    for trigger in triggers.split(","):
        pair = trigger.split(":")
        triggerList.append(pair[0])
        if len(pair) == 1:
            triggerScaleList.append("")
        else:
            triggerScaleList.append(pair[1])
    
    triggerTurnOn = TriggerTurnOn(triggerList, triggerScaleList, path, periodList, filename, entries)
    chargedJetR040 = JetDefinition("Charged", "R040")
    chargedJetR040.AddObservable(Observable("fPt","#it{p}_{T} (GeV/#it{c})", 17, 5, 90, True, False))
    etaObs = Observable("fEta","#it{#eta}", 15, -0.9, 0.9, False, True)
    AddCuts(etaObs, ptcut)
    chargedJetR040.AddObservable(etaObs)
    phiObs = Observable("fPhi","#it{#phi}", 50, 0, 2*math.pi, False, True)
    AddCuts(phiObs, ptcut)
    chargedJetR040.AddObservable(phiObs)
    #chargedJetR040.AddObservable(Observable("fNConstituents","jet constituents", 50, -0.5, 49.5, True))
    zObs = Observable("fZLeading","#it{z}^{leading}_{||}", 12, 0, 1.2, False, True)
    AddCuts(zObs, ptcut)
    chargedJetR040.AddObservable(zObs)
    
    chargedJetR020 = copy.deepcopy(chargedJetR040)
    chargedJetR020.fRadius = "R020"

    fullJetR040 = copy.deepcopy(chargedJetR040)
    fullJetR040.fType = "Full"
    nefObs = Observable("fNEF","NEF", 12, 0, 1.2, False, True)
    AddCuts(nefObs, ptcut)
    fullJetR040.AddObservable(nefObs)
    
    fullJetR020 = copy.deepcopy(fullJetR040)
    fullJetR020.fRadius = "R020"
    
    triggerTurnOn.AddJetDefinition(chargedJetR040)
    triggerTurnOn.AddJetDefinition(chargedJetR020)
    triggerTurnOn.AddJetDefinition(fullJetR040)
    triggerTurnOn.AddJetDefinition(fullJetR020)
    triggerTurnOn.GenerateTriggerTurnOn()
    triggerTurnOn.Draw()
    label = ""
    if ptcut > 0:
        label = "Cut{0}GeV".format(ptcut)
    triggerTurnOn.SaveAll(path, label)
    SaveCanvases(path, label)
    
    globalList.append(triggerTurnOn)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal jet trigger turn on curves.')
    parser.add_argument('trainNumbers', metavar='trainNumbers',
                        help='Comma separated train numbers; use colon to define ranges.')
    parser.add_argument('--trigger', metavar='trigger',
                        default='EMCEJE',
                        help='Rare trigger name (e.g. EMCEJE)')
    parser.add_argument('--train', metavar='train',
                        default='Jets_EMC_pp',
                        help='Train name to be analyzed (e.g. Jets_EMC_pp)')
    parser.add_argument('--period', metavar='period',
                        default='LHC12x',
                        help='Period to be analyzed (e.g. LHC12x)')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Volumes/DATA/ALICE/JetResults",
                        help='Input path')
    parser.add_argument('--file-name', metavar='file-name',
                        default="AnalysisResults.root",
                        help='File name (e.g. AnalysisResults.root')
    parser.add_argument('--entries', metavar='entries',
                        default=0, type=int,
                        help='Input path')
    parser.add_argument('--pt-cut', metavar='ptcut',
                        default=0, type=int,
                        help='Jet pT cut')
    args = parser.parse_args()
    
    
    main(args.train, args.period, args.trainNumbers, args.input_path, args.file_name, args.trigger, args.entries, args.pt_cut)
    
    IPython.embed()