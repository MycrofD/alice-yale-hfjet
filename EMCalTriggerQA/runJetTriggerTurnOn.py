#!/usr/bin/env python
#python script to do EMCal jet trigger turn on curves

import ROOT
import argparse
import IPython
import copy

globalList = []

class JetDefinition:
    def __init__(self, _type, _radius):
        self.fType = _type
        self.fRadius = _radius
    
    def GetBranchName(self):
        if self.fType == "Charged":
            return "Jet_AKTCharged{0}_tracks_pT0150_pt_scheme".format(self.fRadius)
        elif self.fType == "Full":
            return "Jet_AKTFull{0}_tracks_pT0150_caloClusters_E0300_pt_scheme".format(self.fRadius)

class Observable:
    def __init__(self, name, title, bins, min, max, logY=False):
        self.fName = name
        self.fTitle = title
        self.fBins = bins
        self.fMin = min
        self.fMax = max
        self.fHistogram = None
        self.fLogy = logY
        
    def CreateHistogram(self, hname, htitle, xtitle="", ytitle=""):
        hname = "{0}_{1}".format(hname, self.fName)
        if not xtitle:
            xtitle = self.fTitle
        if not ytitle:
            ytitle = "entries"
        htitle = "{0} {1};{2};{3}".format(htitle, self.fTitle, xtitle, ytitle)
        return ROOT.TH1F(hname, htitle, self.fBins, self.fMin, self.fMax)
        
class TreeProjector:    
    def __init__(self, name, title, path, subdirs, filename, trigger, color, marker):
        self.CreateChain(path, subdirs, filename, trigger)
        self.fName = name
        self.fTitle = title
        self.fColor = color
        self.fMarker = marker
        
    def SetJetDefinition(self, jetDef):
        self.fJetDefinition = jetDef
    
    def CreateChain(self, path, subdirs, filename, trigger):
        treeName = "AliAnalysisTaskEmcalJetTree_{0}_jets".format(trigger)
        self.fChain = ROOT.TChain(treeName)
    
        for dir in subdirs:
            fullFileName = "{0}/{1}/{2}".format(path, dir, filename)
            print("Adding {0}".format(fullFileName))
            self.fChain.Add(fullFileName)
    
    def AddObservable(self, observable):
        if not hasattr(self, "fObservableList"):
            self.fObservableList = []
        self.fObservableList.append(observable)
            
    def ProjectTree(self, maxEntries=0):
        print("Projecting {0}".format(self.fName))
        
        for observable in self.fObservableList:
            observable.fHistogram = observable.CreateHistogram(self.fName, self.fTitle)
            observable.fHistogram.Sumw2()

        for i, event in enumerate(self.fChain):
            if maxEntries>0 and i >= maxEntries:
                break
            if i % 100000 == 0:
                print("Entry {0}".format(i))
            for observable in self.fObservableList:
                jets = getattr(event, self.fJetDefinition.GetBranchName())
                for jet in jets:
                    observable.fHistogram.Fill(getattr(jet, observable.fName))
            
    def GenerateCanvasList(self, label="", noLogY=False):
        canvasList = []
        for obs in self.fObservableList:
            if label:
                cname = "{0}_{1}".format(label, obs.fName)
                ctitle = "{0} {1}".format(label, obs.fTitle)
            else:
                cname = obs.fName
                ctitle = obs.fTitle
            print("Generating canvas {0} for observable {1}".format(cname, obs.fName))
            canvas = ROOT.TCanvas(cname, ctitle)
            if not noLogY:
                canvas.SetLogy(obs.fLogy)
            canvasList.append(canvas)
            globalList.append(canvas)
            
        return canvasList
                
    def Draw(self, canvasList=None, ratioCanvasList=None, doRatios=False):
        print("Drawing {0}...".format(self.fTitle))
        if canvasList:
            print("Using pre-existing canvas list...")
            self.fListOfCanvases = canvasList;
            opt = "same"
        else:
            print("Generating canvas list...")
            self.fListOfCanvases = self.GenerateCanvasList()
            opt = ""          
        
        self.DrawHisto(self.fListOfCanvases, "fHistogram", opt)
        
        if doRatios:
            if len(ratioCanvasList) > 0:
                self.ListOfRatioCanvases = ratioCanvasList
                opt = "same"
            else:
                self.ListOfRatioCanvases = self.GenerateCanvasList("Ratio", True)
                opt = ""
            
            self.DrawHisto(self.ListOfRatioCanvases, "fRatio", opt)
    
    def DrawHisto(self, canvasList, histoName, opt):
        for obs,canvas in zip(self.fObservableList, canvasList):
            canvas.cd()
            histo = getattr(obs, histoName)
            histo.SetMarkerSize(1)
            histo.SetMarkerStyle(self.fMarker)
            histo.SetMarkerColor(self.fColor)
            histo.SetLineColor(self.fColor)
            histo.Draw(opt)

class TriggerTurnOn:
    def __init__(self, triggerList, MBtrigger, path, subdirs, filename, entries):
        self.fListOfTreeProjectors = []
        colorList = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+1, ROOT.kOrange+1]
        markerList = [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenDiamond, ROOT.kOpenCross]
        for trigger,color,marker in zip(triggerList,colorList,markerList):
            self.fListOfTreeProjectors.append(TreeProjector(trigger, trigger, path, subdirs, filename, trigger, color, marker))
        self.fMBTreeProjector = TreeProjector(MBtrigger, MBtrigger, path, subdirs, filename, MBtrigger, ROOT.kBlack, ROOT.kFullCircle)
        self.fEntries = entries
        
    def SetJetDefinition(self, jetDef):
        for proj in self.fListOfTreeProjectors:
            proj.SetJetDefinition(jetDef)
        self.fMBTreeProjector.SetJetDefinition(jetDef)
    
    def AddObservable(self, observable):
        for proj in self.fListOfTreeProjectors:
            proj.AddObservable(copy.copy(observable))
        self.fMBTreeProjector.AddObservable(copy.copy(observable))
    
    def GenerateTriggerTurnOn(self):
        self.fMBTreeProjector.ProjectTree(self.fEntries)
        for proj in self.fListOfTreeProjectors:
            proj.ProjectTree(self.fEntries)
            self.GenerateTriggerTurnOnForProj(proj)
            
    def GenerateTriggerTurnOnForProj(self, proj):
        for obsRare, obsMB in zip(proj.fObservableList, self.fMBTreeProjector.fObservableList):
            obsRare.fRatio = obsMB.CreateHistogram("TurnOn_{0}".format(proj.fName), 
                                                   "TurnOn {0}".format(proj.fTitle),
                                                   "", "ratio")
            obsRare.fRatio.Divide(obsRare.fHistogram, obsMB.fHistogram)
            
    def Draw(self):
        self.fMBTreeProjector.Draw()
        self.fListOfCanvases = self.fMBTreeProjector.fListOfCanvases
        self.fListOfRatioCanvases = []
        for proj in self.fListOfTreeProjectors:
            proj.Draw(self.fListOfCanvases, self.fListOfRatioCanvases, True)

def main(train, period, trainNumbers, inputPath, filename, triggers, MBtrigger, entries):
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
    
    triggerTurnOn = TriggerTurnOn(triggers.split(","), MBtrigger, "{0}/{1}".format(inputPath, train), periodList, filename, entries)
    triggerTurnOn.SetJetDefinition(JetDefinition("Full", "R040"))
    triggerTurnOn.AddObservable(Observable("fPt","#it{p}_{T} (GeV/#it{c})", 19, 5, 100, True))
    triggerTurnOn.AddObservable(Observable("fNEF","#it{p}_{T} (GeV/#it{c})", 10, 0, 1, True))
    triggerTurnOn.GenerateTriggerTurnOn()
    triggerTurnOn.Draw()
    
    globalList.append(triggerTurnOn)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal jet trigger turn on curves.')
    parser.add_argument('trainNumbers', metavar='trainNumbers',
                        help='Comma separated train numbers; use colon to define ranges.')
    parser.add_argument('--trigger', metavar='trigger',
                        default='EMCEJE',
                        help='Rare trigger name (e.g. EMCEJE)')
    parser.add_argument('--mb-trigger', metavar='MBtrigger',
                        default='AnyINT',
                        help='MB trigger name (e.g. AnyINT)')
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
    args = parser.parse_args()
    
    
    main(args.train, args.period, args.trainNumbers, args.input_path, args.file_name, args.trigger, args.mb_trigger, args.entries)
    
    IPython.embed()