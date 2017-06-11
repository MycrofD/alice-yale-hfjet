#!/usr/local/bin/python
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
        self.fSlopes = dict()
        self.fSlopeErrors = dict()
        
        for run in self.fRunList:
            self.LoadHistogramsForRunNumber(run)
            
    def LoadHistogramsForRunNumber(self, run):
        self.fHistograms[run] = dict()
        self.fSlopes[run] = dict()
        self.fSlopeErrors[run] = dict()
        
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
                
            hevents = list.FindObject("fHistEventCount")
            if hevents:
                print("Run {0}: number of events ({1}) = {2}".format(run, trigger, hevents.GetBinContent(1)))
            
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
            self.fSlopes[run][trigger] = dict()
            self.fSlopeErrors[run][trigger] = dict()
        
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
        self.fTRUlist = dict()
        
        for trigger in self.fTriggerList:
            if "EMC" in trigger:
                self.fTRUlist[trigger] = range(0, 32)
            else:
                self.fTRUlist[trigger] = [32, 33, 36, 37, 38, 39, 42, 43, 44, 45, 48, 49, 50, 51]
        
            for tru in self.fTRUlist[trigger]:
                self.PlotTRUvsSTU(trigger, tru)
                
        self.PlotSlopes()
        self.SaveRootFile()
            
    def PlotTRUvsSTU(self, trigger, tru):
        cname = "{0}_TRUvsSTU_TRU{1}".format(trigger, tru)
        ctitle = "{0} TRUvsSTU TRU={1}".format(trigger, tru)
        canvas = ROOT.TCanvas(cname, ctitle, 600, 600)
        canvas.SetLogz()
        legend = ROOT.TLegend(0.1, 0.65, 0.5, 0.90)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        
        first = True;
        for i,run in enumerate(self.fRunList):
            funct = ROOT.TF1("{0}_{1}_TRU{2}_fit".format(trigger, run, tru), "pol1", 20, 200)
            funct.SetLineColor(self.fColors[i])
            globalList.append(funct)
            h = self.fHistograms[run][trigger][tru]
            h.SetMarkerColor(self.fColors[i])
            h.SetFillColor(self.fColors[i])
            h.SetLineColor(self.fColors[i])
            
            h.GetXaxis().SetRangeUser(20,200)
            h.GetYaxis().SetRangeUser(20,200)
            h.Fit(funct, "0")
            
            h.GetXaxis().SetRangeUser(0,200)
            h.GetYaxis().SetRangeUser(0,200)
            if first:
                h.Draw()
                first = False
            else:
                h.Draw("same")
                
            funct.Draw("same")
            self.fSlopes[run][trigger][tru] = funct.GetParameter(1)
            self.fSlopeErrors[run][trigger][tru] = funct.GetParError(1)
            legend.AddEntry(h, "Run {0}, slope {1:.2f} #pm {2:.2f}".format(run, funct.GetParameter(1), funct.GetParError(1)), "f")
    
        legend.Draw()
        globalList.append(legend)
        globalList.append(canvas)
        canvas.SaveAs("{0}/{1}/{2}.pdf".format(self.fInputPath, self.fTrain, canvas.GetName()))
        
    def PlotSlopes(self):
        nTRUs = 0;
        for trigger in self.fTriggerList:
            nTRUs += len(self.fTRUlist[trigger])
            
        cname = "TRUvsSTU_Slopes"
        ctitle = "TRUvsSTU slopes"
        canvas = ROOT.TCanvas(cname, ctitle, 1200, 400)
        canvas.SetLeftMargin(0.06)
        canvas.SetRightMargin(0.04)
        legend = ROOT.TLegend(0.85, 0.25, 0.98, 0.90)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
            
        self.fTRUvsSTUgraphs = dict()
        
        first = True
        for irun,run in enumerate(self.fRunList):
            g = ROOT.TGraphErrors(nTRUs)
            g.SetName("SlopeTRUvsSTU_Run{0}".format(run))
            g.SetMarkerColor(self.fColors[irun])
            g.SetLineColor(self.fColors[irun])
            g.SetMarkerStyle(ROOT.kOpenCircle)
            g.SetMarkerSize(1)
            globalList.append(g)
            self.fTRUvsSTUgraphs[run] = g
            legend.AddEntry(g, "Run {0}".format(run), "pe")
            totpoints = 0
            for trigger in self.fTriggerList:
                for ipoint,tru in enumerate(self.fTRUlist[trigger]):
                    g.SetPoint(ipoint+totpoints, tru, self.fSlopes[run][trigger][tru])
                    g.SetPointError(ipoint+totpoints, 0, self.fSlopeErrors[run][trigger][tru])
                totpoints = len(self.fTRUlist[trigger])
            
            if first:
                g.Draw("ap")
                first = False
            else:
                g.Draw("p")
            
            g.GetXaxis().SetTitle("TRU #")
            g.GetYaxis().SetTitle("TRU vs STU slope")
            g.GetYaxis().SetTitleOffset(0.6)
                
        legend.Draw()
        globalList.append(legend)
        globalList.append(canvas)
        
        canvas.SaveAs("{0}/{1}/{2}.pdf".format(self.fInputPath, self.fTrain, canvas.GetName()))
        
    def SaveRootFile(self):
        out = ROOT.TFile.Open("{0}/{1}/TRUvsSTU.root".format(self.fInputPath, self.fTrain), "recreate")
        out.cd()
        
        for run in self.fRunList:
            for trigger in self.fTriggerList:
                for tru in self.fTRUlist[trigger]:
                    self.fHistograms[run][trigger][tru].Write()
                    fname = "{0}_{1}_TRU{2}_fit".format(trigger, run, tru)
                    f = self.fHistograms[run][trigger][tru].GetFunction(fname)
                    if f:
                        f.Write()
            self.fTRUvsSTUgraphs[run].Write()
            
        out.Close()

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
