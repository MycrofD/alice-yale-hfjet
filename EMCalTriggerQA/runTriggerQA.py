#!/usr/bin/env python
#python script to do EMCal trigger QA plots

import ROOT
import argparse
from ROOT import gROOT

import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4  # 0.075882352

def Plot1D(hlist, hname, thresholds, nevents, color):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '", hname, "' in list '", hlist.GetName(), "'. Skipping..."
        return
    
    canvas = ROOT.TCanvas(hname, hname)
    canvas.SetLogy()
    canvas.cd()
    
    hist.Sumw2()
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(1)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.Draw()
    globalList.append(canvas)
    
    print hname
    totEntries = hist.GetEntries()
    print "Total number of entries is", totEntries
    
    for th in thresholds:
        thADC = th / ADCtoGeV
        entries = hist.Integral(hist.GetXaxis().FindBin(thADC), -1)
        suppression = entries / nevents
        #print "Entries above " + str(th) + " GeV (" + str(thADC) + " ADC): " + str(entries) + ". Suppression = " + str(suppression)
        print str(th) + ";" + str(thADC) + ";" + str(entries) + ";" + str(suppression)
    

def main(train, trigger="EMC7", inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = inputPath + "/" + train + "/AnalysisResults.root"
    listName = "AliEmcalTriggerQATaskPP_" + trigger + "_histos"
    #hlistName = "histos" + "AliEmcalTriggerQATaskPP_" + trigger + "_AliEmcalTriggerQAAP_Cent0"
    hlistName = "histos" + "AliEmcalTriggerQATaskPP_" + trigger + "_Cent0"
    
    file = ROOT.TFile.Open(fileName);
    if not file:
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
        
    if file.IsZombie():
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
    
    list = file.Get(listName)
    if not list:
        print "Could not get list '" + listName + "from file '" + fileName + "'! Aborting..."
        exit(1)
        
    hlist = list.FindObject(hlistName)
    if not hlist:
        print "Could not get hash list '" + hlistName + "' in list '" + listName + "from file '" + fileName + "'! Aborting..."
        exit(1)
    
    thresholds_GA = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    thresholds_JE = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
    
    hEventsName = "fHistEventCount"
    hevents = list.FindObject(hEventsName)
    if not hevents:
        print "Could not find histogram '" + hEventsName + "' in list '" + listName + "'. Aborting..."
        exit(1) 
    
    nevents = hevents.GetBinContent(1)
    
    print "Total number of events is", nevents
    
    Plot1D(hlist, "EMCTRQA_histEMCalMaxPatchAmpEMCGAHRecalc", thresholds_GA, nevents, ROOT.kRed+1)
    Plot1D(hlist, "EMCTRQA_histEMCalMaxPatchAmpEMCJEHRecalc", thresholds_JE, nevents, ROOT.kBlue+1)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--trigger', metavar='trigger',
                        default="EMC7",
                        help='Trigger')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.trigger, args.input_path)
    
    IPython.embed()