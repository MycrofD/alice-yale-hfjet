#!/usr/bin/env python
#python script to do EMCal cluster trigger turn on curves

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math
from collections import OrderedDict
import yaml
import IPython
from cProfile import label

globalList = []

class TriggerTurnOnAnalysis:
    def __init__(self, detector, file, normalize, outputPath, triggerList, baseTrigger, clusters, cells):
        self.fNormalize = normalize
        self.fDetector = detector
        self.fName = "{0}_Cluster".format(detector)
        self.fTitle = "{0} Cluster".format(detector)
        self.fOutputPath = outputPath
        self.fColorList = dict()
        self.fTriggerList = triggerList
        self.fBaseTrigger = baseTrigger
        self.fColorList[self.fBaseTrigger] = ROOT.kBlack
        colors = [ROOT.kRed+2, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kOrange+2]
        for i,trigger in enumerate(self.fTriggerList):
            if trigger == self.fBaseTrigger:
                continue
            self.fColorList[trigger] = colors[i]
         
        self.fHistogramName = "fHistClusPhiEtaEnergy_0"
        self.fListName = "AliAnalysisTaskEmcalJetQA_{clusName}_{cellName}_{triggerName}_histos/histosAliAnalysisTaskEmcalJetQA_{clusName}_{cellName}_{triggerName}/{clusName}".format(clusName=clusters, cellName=cells, triggerName="{0}")
        self.fFile = file
        
    def GetHistogram(self, trigger):
        listName = self.fListName.format(trigger)
        
        chain = listName.split("/")
        list = None
        for name in chain:
            if list:
                list = list.FindObject(name)
            else:
                list = self.fFile.Get(name)
                if list:
                    hEventsName = "fHistEventCount"
                    hevents = list.FindObject(hEventsName)
                    if not hevents:
                        print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hEventsName, name))
                        list.Print()
                        exit(1) 
    
                    self.fEvents[trigger] = hevents.GetBinContent(1)
                    print("Total number of events for {0} is {1}".format(trigger, self.fEvents[trigger]))
                
        if not list:
            self.fFile.ls()
            print("Could not get list '{0}' from file '{1}'! Aborting...".format(listName, self.fFile.GetName()))
            exit(1)
                
        hist = list.FindObject(self.fHistogramName)
        if not hist:
            list.Print()
            print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(self.fHistogramName, listName))
            exit(1)
                           
        self.fHistograms[trigger] = hist.Clone("Clusters_{0}".format(trigger))
        self.fHistograms[trigger].Sumw2()
    
    def GetHistograms(self):
        self.fHistograms = dict()
        self.fEvents = dict()
        for trigger in self.fTriggerList:
            self.GetHistogram(trigger)
            
    def ProjectHistograms(self):
        self.fProjections = dict()
        for trigger,hist in self.fHistograms.iteritems():
            if self.fDetector == "EMCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, 0, hist.GetYaxis().FindBin(4), "e")
            elif self.fDetector == "DCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, hist.GetYaxis().FindBin(4), -1, "e")
            elif self.fDetector == "EMCal+DCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, 0, -1, "e")
            if self.fNormalize:
                self.fProjections[trigger].Scale(1. / self.fEvents[trigger], "width")
                self.fProjections[trigger].GetYaxis().SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{E}_{cluster}}")
            else:
                self.fProjections[trigger].Scale(1, "width")
                self.fProjections[trigger].GetYaxis().SetTitle("#frac{d#it{N}}{d#it{E}_{cluster}}")
            
            
    def Ratios(self):
        self.fRatios = dict()
        for trigger,hist in self.fProjections.iteritems():
            if trigger == self.fBaseTrigger:
                continue
            ratio = hist.Clone("Ratio_{0}".format(trigger))
            ratio.Divide(self.fProjections[self.fBaseTrigger])
            if self.fNormalize:
                ratio.GetYaxis().SetTitle("#frac{#it{N}_{evt, {" + self.fBaseTrigger + "}}{#it{N}_{evt, L1} #frac{d#it{N}_{L1}}{d#it{E}_{cluster}} / #frac{d#it{N}_{" + self.fBaseTrigger + "}}{d#it{E}_{cluster}}")
            else:
                ratio.GetYaxis().SetTitle("#frac{d#it{N}_{L1}}{d#it{E}_{cluster}} / #frac{d#it{N}_{" + self.fBaseTrigger + "}}{d#it{E}_{cluster}}")
            self.fRatios[trigger] = ratio
            globalList.append(ratio)
    
    def Plot(self):
        min = 1e-6
        max = 50
        maxR = 100
        if not self.fNormalize:
            min *= self.fEvents[self.fBaseTrigger]
            max *= self.fEvents[self.fBaseTrigger]
            maxR = 1.2
        self.DoPlot(self.fProjections, "{0}_Spectra".format(self.fName), "{0} Spectra".format(self.fTitle), True, min, max)
        self.DoPlot(self.fRatios, "{0}_Ratios".format(self.fName), "{0} Ratios".format(self.fTitle),  False, 0, maxR)
        
    def DoPlot(self, list, name, title, logY, min, max):
        canvas = ROOT.TCanvas(name, title)
        canvas.cd()
        if logY:
            canvas.SetLogy()
        leg = ROOT.TLegend(0.6, 0.7, 0.85, 0.85)
        leg.SetTextFont(43)
        leg.SetTextSize(12)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        i = 0
        for label,hist in list.iteritems():
            hist.SetLineColor(self.fColorList[label])
            hist.SetMarkerColor(self.fColorList[label])
            hist.SetMarkerStyle(ROOT.kOpenCircle)
            hist.SetMarkerSize(1)
            hist.GetXaxis().SetRangeUser(0,30)
            leg.AddEntry(hist, label, "pe")
            if i == 0:
                hist.Draw()
                hist.GetYaxis().SetRangeUser(min,max)
            else:
                hist.Draw("same")
            i += 1
            
        leg.Draw()
        globalList.append(leg)
        globalList.append(canvas)
        canvas.SaveAs("{0}/{1}.pdf".format(self.fOutputPath, name))
        
def main(yamlConfig):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = "{0}/{1}/AnalysisResults.root".format(yamlConfig["input_path"], yamlConfig["train"])
    outputPath = "{0}/{1}".format(yamlConfig["input_path"], yamlConfig["train"])
    
    file = ROOT.TFile.Open(fileName);
    if not file or file.IsZombie():
        print("Could not open file '{0}'! Aborting...".format(fileName))
        exit(1)

    analysis = []
    
    for config in yamlConfig["configs"]:
        if not config["active"]:
            continue
        analysis.append(TriggerTurnOnAnalysis(config["detector"], file, False, outputPath, config["triggers"], config["base_trigger"], yamlConfig["clusters"], yamlConfig["cells"]))

    for ana in analysis:
        ana.GetHistograms()
        ana.ProjectHistograms()
        ana.Ratios()
        ana.Plot()
    
    globalList.append(ana)
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config)
    
    IPython.embed()
    