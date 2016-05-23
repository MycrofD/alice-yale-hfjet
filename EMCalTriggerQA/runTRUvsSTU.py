#!/usr/bin/env python
#python script to do TRU vs STU correlation studies

import ROOT
import argparse
import math
import yaml
import IPython

globalList = []

class TRUvsSTUanalysis:
    def __init__(self, runlist, train, triggerlist, path):
        self.fRunList = runlist
        self.fTrain = train
        self.fInputPath = path
        self.fTriggerList = triggerlist
        
        self.fColors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+1, ROOT.kOrange+1, ROOT.kMagenta+1, ROOT.kTeal+1, ROOT.kYellow+1]
        
    def LoadHistograms(self):
        self.fHistograms = dict()
        
        for run in self.fRunList:
            self.LoadHistogramsForRunNumber(run)
            
    def LoadHistogramsForRunNumber(self, run):
        self.fHistograms[run] = dict()
        
        fileName = "{0}/{1}/{2}/AnalysisResults.root".format(self.fInputPath, self.fTrain, run)
        file = ROOT.TFile.Open(fileName)
        if not file:
            print "Could not open file '" + fileName + "'! Aborting..."
            exit(1)
            
        if file.IsZombie():
            print "Could not open file '" + fileName + "'! Aborting..."
            exit(1)
            
        for trigger in self.fTriggerList:
            listName = "AliEmcalTriggerQATask_" + trigger + "_histos"
            hlistName = "histosAliEmcalTriggerQATask_" + trigger
        
            list = file.Get(listName)
            if not list:
                file.ls()
                print "Could not get list '" + listName + "' from file '" + fileName + "'! Aborting..."
                exit(1)
            
            hlist = list.FindObject(hlistName)
            if not hlist:
                list.Print()
                print "Could not get hash list '" + hlistName + "' in list '" + listName + "' from file '" + fileName + "'! Aborting..."
                exit(1)
                
            trulist = hlist.FindObject("ByTRU")
            if not trulist:
                hlist.Print()
                print("Could not get hash list 'ByTRU' in list '{0}'! Aborting...".format(hlist.GetName()))
                exit(1)
                
            outhname = "TRUvsSTU_{0}_{1}_{2}".format(run, trigger, "{0}")
            self.fHistograms[run][trigger] = self.LoadHistogramsFromList(trulist, outhname)
        
        file.Close()
            
    def LoadHistogramsFromList(self, list, outhname):
        inhname = "EMCTRQA_histFastORAmpSTUVsTRU{0}"
        res = dict()
        for i in range(0, 52):
            h = list.FindObject(inhname.format(i))
            if h:
                res[i] = h.Clone(outhname.format(i))
                
        return res

    def RunAnalysis(self):
        for trigger in self.fTriggerList:
            if "EMC" in trigger:
                trulist = range(0, 32)
            else:
                trulist = [32, 33, 36, 37, 38, 39, 42, 43, 44, 45, 48, 49, 50, 51]
        
            for tru in trulist:
                self.PlotTRUvsSTU(trigger, tru)
            
    def PlotTRUvsSTU(self, trigger, tru):
        cname = "{0}_TRUvsSTU_TRU{1}".format(trigger, tru)
        ctitle = "{0} TRUvsSTU TRU={1}".format(trigger, tru)
        canvas = ROOT.TCanvas(cname, ctitle)
        canvas.SetLogz()
        legend = ROOT.TLegend(0.6, 0.7, 0.95, 0.95)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        
        first = True;
        for i,run in enumerate(self.fRunList):
            h = self.fHistograms[run][trigger][tru]
            h.GetXaxis().SetRangeUser(0,200)
            h.GetYaxis().SetRangeUser(0,200)
            h.SetMarkerColor(self.fColors[i])
            h.SetFillColor(self.fColors[i])
            h.SetLineColor(self.fColors[i])
            if first:
                h.Draw()
                first = False
            else:
                h.Draw("same")
                
            legend.AddEntry(h, "Run {0}".format(run), "f")
    
        legend.Draw()
        
        globalList.append(canvas)

def main(config):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    ana = TRUvsSTUanalysis(config["runlist"], config["train"], config["triggerlist"], config["path"])
    ana.LoadHistograms()
    ana.RunAnalysis()
    globalList.append(ana)
    
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='TRU vs STU analysis.')
    parser.add_argument('configYAML', metavar='configYAML',
                        help='YAML format coniguration file')
    args = parser.parse_args()
    
    f = open(args.configYAML, 'r')
    config = yaml.load(f)
    f.close()
    
    main(config)
    
    IPython.embed()