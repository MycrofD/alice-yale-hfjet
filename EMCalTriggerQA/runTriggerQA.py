#!/usr/bin/env python
#python script to do EMCal trigger QA plots

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum

import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4  # 0.075882352

class Axis(Enum):
    GeV = 1
    ADC = 2
    
def Plot2D(hlist, hname, thresholds1, thresholds2, nevents, axis):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '", hname, "' in list '", hlist.GetName(), "'. Skipping..."
        return
    
    canvas = ROOT.TCanvas(hname, hname)
    canvas.SetLogz()
    canvas.cd()
    
    max = 1000
    
    if axis is Axis.GeV:
        max *= ADCtoGeV
        hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
        hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
        hist.GetYaxis().SetLimits(hist.GetYaxis().GetXmin() * ADCtoGeV, hist.GetYaxis().GetXmax() * ADCtoGeV)
        hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() + " (GeV)")
    
    hist.Sumw2()
    hist.GetXaxis().SetRangeUser(0,max)
    hist.GetYaxis().SetRangeUser(0,max)
    hist.Draw("colz")
    globalList.append(canvas)
    
    print hname
    totEntries = hist.GetEntries()
    print "Total number of entries is", totEntries
    
    print "Partial number of events is", nevents
    
    for th1,th2 in zip(thresholds1, thresholds2):
        if axis is Axis.ADC:
            th1 /= ADCtoGeV
            th2 /= ADCtoGeV
            
        entries = hist.Integral(hist.GetXaxis().FindBin(th1), -1, hist.GetXaxis().FindBin(th2), -1)
        suppression = entries / nevents
        
        print str(th1) + ";" + str(th2) + ";" + str(entries) + ";" + str(suppression)

def Plot1D(hlist, hname, thresholds, zeroth, nevents, color, axis):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '", hname, "' in list '", hlist.GetName(), "'. Skipping..."
        return
    
    canvas = ROOT.TCanvas(hname, hname)
    canvas.SetLogy()
    canvas.cd()
    
    max = 1000
    
    if axis is Axis.GeV:
        max *= ADCtoGeV
        hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
        hist.GetXaxis().SetTitle("Energy (GeV)")
    
    hist.Sumw2()
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(1)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.GetXaxis().SetRangeUser(0,max)
    hist.Draw()
    globalList.append(canvas)
    
    print hname
    totEntries = hist.GetEntries()
    print "Total number of entries is", totEntries
    
    if zeroth >= 0:
        if axis is Axis.ADC:
            zeroth /= ADCtoGeV
        print 'Calculating partial number of events above ', zeroth
        nevents = hist.Integral(hist.GetXaxis().FindBin(zeroth), -1)
    
    print "Partial number of events is", nevents
    
    for th in thresholds:
        if axis is Axis.ADC:
            th /= ADCtoGeV
        entries = hist.Integral(hist.GetXaxis().FindBin(th), -1)
        suppression = entries / nevents
        print str(th) + ";" + str(entries) + ";" + str(suppression)
        
    return nevents
    

def main(train, trigger="EMC7", type="Recalc", axis="ADC", inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = inputPath + "/" + train + "/AnalysisResults.root"
    listName = "AliEmcalTriggerQATaskPP_" + trigger + "_histos"
    hlistName = "histos" + "AliEmcalTriggerQATaskPP_" + trigger
    #hlistName = "histos" + "AliEmcalTriggerQATaskPP_" + trigger
    
    file = ROOT.TFile.Open(fileName);
    if not file:
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
        
    if file.IsZombie():
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
    
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
    
    thresholds_GA = [3, 4,  5,  6,  7,  8,  9, 10]
    thresholds_JE = [6, 8, 10, 12, 14, 16, 18, 20]
    
    hEventsName = "fHistEventCount"
    hevents = list.FindObject(hEventsName)
    if not hevents:
        print "Could not find histogram '" + hEventsName + "' in list '" + listName + "'. Aborting..."
        exit(1) 
    
    nevents = hevents.GetBinContent(1)
    
    print "Total number of events is", nevents
    
    eaxis = Axis.ADC
    if axis == "GeV":
        eaxis = Axis.GeV
    
    nevents = Plot1D(hlist, "EMCTRQA_histEMCalMaxPatchAmpEMCGAH" + type, thresholds_GA, 2.5, -1, ROOT.kRed+1, eaxis)
    Plot1D(hlist, "EMCTRQA_histEMCalMaxPatchAmpEMCJEH" + type, thresholds_JE, -1, nevents, ROOT.kBlue+1, eaxis)
    
    Plot2D(hlist, "EMCTRQA_histEMCalEMCGAHMaxVsEMCJEHMax" + type, thresholds_JE, thresholds_GA, nevents, eaxis)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--trigger', metavar='trigger',
                        default="EMC7",
                        help='Trigger')
    parser.add_argument('--type', metavar='type',
                        default="Recalc",
                        help='Trigger')
    parser.add_argument('--axis', metavar='axis',
                        default="ADC",
                        help='Trigger')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.trigger, args.type, args.axis, args.input_path)
    
    IPython.embed()