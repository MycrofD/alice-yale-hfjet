#!/usr/bin/env python
#python script to do EMCal jet spectra

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
    def __init__(self, name, title, tree, color, marker):
        self.fName = name
        self.fTitle = title
        self.fColor = color
        self.fMarker = marker
        self.fProjectionOk = False
        self.fTree = tree
        
    def AddJetDefinition(self, jetDef):
        if not hasattr(self, "fJetDefinitionList"):
            self.fJetDefinitionList = []
        self.fJetDefinitionList.append(jetDef)
            
    def ProjectTree(self, maxEntries=0):
        print("Projecting {0}, number of events: {1}".format(self.fName, self.fTree.GetEntries()))
        
        for jetDef in self.fJetDefinitionList:
            for observable in jetDef.fObservableList:
                observable.fHistogram = observable.CreateHistogram("{0}_{1}".format(self.fName, jetDef.GetName()), 
                                                                   "{0} {1}".format(self.fTitle, jetDef.GetTitle()))
                observable.fHistogram.Sumw2()

        for i, event in enumerate(self.fTree):
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
                        print("Value {0} for observable {1}".format(val, observable.fName))
                        observable.fHistogram.Fill(val)
                        
        for jetDef in self.fJetDefinitionList:
            for observable in jetDef.fObservableList:               
                if observable.fNorm:
                    integral = observable.fHistogram.Integral()
                    if integral > 1e-10:
                        observable.fHistogram.Scale(1. / integral)
                        
        self.fProjectionOk = True
            
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

def AddCuts(obs, ptcut):
    if ptcut > 0:
        obs.AddCut(Observable("fPt","#it{p}_{T} (GeV/#it{c})", 10, ptcut, 1000, True, False))

def main(train, radius, jetType, entries=0, filename="AnalysisResults.root", inputPath="/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/anaDev"):
    ROOT.gSystem.Load("libCGAL")
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.TH1.AddDirectory(False)
    
    filePath = "{0}/{1}/{2}".format(inputPath, train, filename)
    
    treeName = "AliAnalysisTaskEmcalJetTree_jets"
    
    file = ROOT.TFile.Open(filePath)
    
    tree = file.Get(treeName)
    
    treeProj = TreeProjector("JetSpectra", "Jet Spectra", tree, ROOT.kBlue+2, ROOT.kFullCircle)

    jetDef = JetDefinition(jetType, radius)
    jetDef.AddObservable(Observable("fPt","#it{p}_{T} (GeV/#it{c})", 17, 5, 90, True, False))
    jetDef.AddObservable(Observable("fEta","#it{#eta}", 15, -0.9, 0.9, False, True))
    jetDef.AddObservable(Observable("fPhi","#it{#phi}", 50, 0, 2*math.pi, False, True))
    jetDef.AddObservable(Observable("fNConstituents","jet constituents", 50, -0.5, 49.5, True))
    jetDef.AddObservable(Observable("fZLeading","#it{z}^{leading}_{||}", 12, 0, 1.2, False, True))
    
    treeProj.AddJetDefinition(jetDef)
    treeProj.ProjectTree(entries)
    treeProj.Draw()
    #SaveCanvases(path)
    
    globalList.append(treeProj)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal jet trigger turn on curves.')
    parser.add_argument('--train', metavar='train',
                        default='',
                        help='Train name to be analyzed (e.g. Jets_EMC_pp)')
    parser.add_argument('--file-name', metavar='AnalysisResults.root',
                        default='AnalysisResults.root',
                        help='File name')
    parser.add_argument('--input-path', metavar='path/to/file',
                        default="/Users/sa639/Documents/Work/ALICE/alice-yale-hfjet/anaDev",
                        help='Input path')
    parser.add_argument('--entries', metavar='entries',
                        default=0, type=int,
                        help='Input path')
    parser.add_argument('--radius', metavar='R',
                        default='R040', 
                        help='Jet radius')
    parser.add_argument('--jet-type', metavar='Charged|Full',
                        default='Charged', 
                        help='Jet type')
    args = parser.parse_args()
    
    
    main(args.train, args.radius, args.jet_type, args.entries, args.file_name, args.input_path)
    
    IPython.embed()