#!/usr/bin/env python
#python script to do EMCal trigger QA plots

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math
from collections import OrderedDict
import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4 

class Axis(Enum):
    GeV = 1
    ADC = 2
    
class TriggerConfiguration:
    def __init__(self, _lab, _th, _hlist):
        self.fLabel = _lab
        self.fThreshold = _th
        self.fHistogramList = _hlist
        self.fHistogramNameAxis = dict()
        self.fHistograms = dict()
        self.fSharedEntries = OrderedDict()
        self.fSharedEntriesErr = OrderedDict()
        self.fShares = OrderedDict()
        self.fSharesErr = OrderedDict()
        
    def AddHistogram(self, _hname, _axis):
        self.fHistogramNameAxis[_hname] = _axis
        
    def RetrieveHistograms(self):
        for hname in self.fHistogramNameAxis.iterkeys():
            hist = self.fHistogramList.FindObject(hname)
            if not hist:
                print("Could not get histogram '{0}' in list '{1}'. Skipping...".format(hname, self.fHistogramList.GetName()))
                self.fHistogramList.Print()
                break
            
            print("Total number of entries in {0} is {1}.".format(hname, hist.GetEntries()))
            if not "GeV" in hist.GetXaxis().GetTitle():
                hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
                hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
                hist.GetYaxis().SetLimits(hist.GetYaxis().GetXmin() * ADCtoGeV, hist.GetYaxis().GetXmax() * ADCtoGeV)
                hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() + " (GeV)")
                
            self.fHistograms[hname] = hist
            
    def CalculateSuppression(self, nevents):
        hist = self.fHistograms.values()[0]
        
        self.fEvents = nevents
        if self.fHistogramNameAxis[hist.GetName()] == "x":
            self.fEntries = hist.Integral(hist.GetXaxis().FindBin(self.fThreshold), -1, 0, -1)
        elif self.fHistogramNameAxis[hist.GetName()] == "y":
            self.fEntries = hist.Integral(0, -1, hist.GetYaxis().FindBin(self.fThreshold), -1)
        else:
            print("Error no axis defined!!!")

        self.fEntriesErr = math.sqrt(self.fEntries)
        self.fEventsErr = math.sqrt(self.fEvents)
        
        self.fSuppression = self.fEntries / self.fEvents
        self.fErrSuppression= self.fSuppression * (1. / self.fEntriesErr + 1. / self.fEventsErr) 
    
    def FindCommonHistogram(self, other):
        for key in self.fHistogramNameAxis.iterkeys():
            if other.fHistogramNameAxis.has_key(key):
                return key
        return
    
    def CalculateShare(self, other):
        hname = self.FindCommonHistogram(other)
        
        if hname:
            print("Calculating share between {0} and {1} ({2})".format(self.fLabel, other.fLabel, hname))
        else:
            print("Skipping {0} with {1}".format(self.fLabel, other.fLabel))
            return 
        
        if other.fHistogramNameAxis[hname] == self.fHistogramNameAxis[hname]:
            #self trigger is a subset of the other
            if other.fThreshold >= self.fThreshold:
                #the other trigger has a higher threshold than self so other is fully contained in self
                other.fSharedEntries[self.fLabel] = other.fEntries
            else:
                #the other trigger has a lower threshold than self so applying both triggers is equivalent to applying only self
                other.fSharedEntries[self.fLabel] = self.fEntries
        else:
            #the two triggers are applied to different variables
            if other.fHistogramNameAxis[hname] == "x":
                thx = other.fThreshold
            elif other.fHistogramNameAxis[hname] == "y":
                thy = other.fThreshold
                
            if self.fHistogramNameAxis[hname] == "x":
                thx = self.fThreshold
            elif self.fHistogramNameAxis[hname] == "y":
                thy = self.fThreshold
        
            print("{0}: applying {1} to the x axis ({2}) and {3} to the y axis ({4})".format(self.fLabel, thx, self.fHistograms[hname].GetXaxis().GetTitle(), thy, self.fHistograms[hname].GetYaxis().GetTitle()))
        
            other.fSharedEntries[self.fLabel] = self.fHistograms[hname].Integral(self.fHistograms[hname].GetXaxis().FindBin(thx), -1, self.fHistograms[hname].GetYaxis().FindBin(thy), -1)
        
        other.fSharedEntriesErr[self.fLabel] = math.sqrt(other.fSharedEntries[self.fLabel])
        
        if other.fEntries > 0:
            other.fShares[self.fLabel] = other.fSharedEntries[self.fLabel] / other.fEntries
        
            if other.fSharedEntriesErr[self.fLabel] > 0 and other.fEntriesErr > 0:
                other.fSharesErr[self.fLabel] = other.fShares[self.fLabel] * (1. / other.fSharedEntriesErr[self.fLabel] + 1. / other.fEntriesErr)
            else:
                other.fSharesErr[self.fLabel] = 1e9
        else:
            other.fShares[self.fLabel] = 0
        
        print("The shared entries are {0}. Total entries of {1} are {2}, of {3} are {4}".format(other.fSharedEntries[self.fLabel], self.fLabel, self.fEntries, other.fLabel, other.fEntries))
        
    def Print(self, keys):
        res = "|  *{0: ^9}*  ".format(self.fLabel)
        for key in keys:
            if self.fShares.has_key(key):
                share = self.fShares[key]
            else:
                share = -1
            res += "|  {0: ^11.3f}  ".format(share)
        res += "|"
        print(res)
        
class TriggerCombination:
    def __init__(self, label, rate):
        self.fLabel = label
        self.fRate = rate
        
    def Print(self):
        print("|  *{0: ^20}*  |  {1:.4f}  |".format(self.fLabel, self.fRate))
        
class TriggerAnalysis:
    fTriggerList = OrderedDict()
    
    def AddTrigger(self, trigger):
        self.fTriggerList[trigger.fLabel] = trigger
        return trigger
    
    def DoAnalysis(self, nevents):
        for trigger in self.fTriggerList.itervalues():
            trigger.RetrieveHistograms()
            trigger.CalculateSuppression(nevents)
        
        for trigger1 in self.fTriggerList.itervalues():
            for trigger2 in self.fTriggerList.itervalues():
                trigger1.CalculateShare(trigger2)
        
        res = "|  *Triggers *  "
        for labels in self.fTriggerList.iterkeys():
            res += "|  *{0: ^9}*  ".format(labels)
        res += "|"
        print(res)
            
        for trigger in self.fTriggerList.itervalues():
            trigger.Print(self.fTriggerList.keys())
        
    def PrintTriggerSuppression(self, list, baseTrigger, subtract):
        label = ""    
        rate = 0
        
        for i,trigger in enumerate(list):
            if label:
                label += "+"
            label += trigger
            #print("Rate for trigger {0} is {1:.3f}".format(trigger, self.fTriggerList[baseTrigger].fShares[trigger]))
            rate += self.fTriggerList[baseTrigger].fShares[trigger] 
            for j in range(i+1,len(list)):
                #print("Subtracting share with {0}: {1:.3f} * {2:.3f} = {3:.3f}".format(list[j], self.fTriggerList[baseTrigger].fShares[trigger], self.fTriggerList[trigger].fShares[list[j]], self.fTriggerList[baseTrigger].fShares[trigger] * self.fTriggerList[trigger].fShares[list[j]]))
                rate -= self.fTriggerList[baseTrigger].fShares[trigger] * self.fTriggerList[trigger].fShares[list[j]]
        
        rate -= subtract
        #print("Subtracting {0:.3f}".format(subtract))
        
        #print("{0} = {1:.3f}".format(label, rate))
        
        return TriggerCombination(label, rate)
        
def CalculateTriggerSuppression(hlist, type, nevents, thresholds):
    triggerAna = TriggerAnalysis()
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("EG1", thresholds["G1"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal", "EMCal", type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCGAHMax{2}".format("EMCal","DCal",type), "y")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("EG2", thresholds["G2"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal", "EMCal", type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCGAHMax{2}".format("EMCal","DCal",type), "y")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("EJ1", thresholds["J1"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal", "EMCal", type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCJEHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal","EMCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("EJ2", thresholds["J2"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal", "EMCal", type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCJEHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal","EMCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("DG1", thresholds["G1"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal", "DCal", type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal","EMCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCGAHMax{2}".format("EMCal","DCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("DG2", thresholds["G2"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal", "DCal", type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal","EMCal",type), "y")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCGAHMax{2}".format("EMCal","DCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("DJ1", thresholds["J1"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal", "DCal", type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCJEHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("DJ2", thresholds["J2"], hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal", "DCal", type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCJEHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal","DCal",type), "x")
    
    trgConf = triggerAna.AddTrigger(TriggerConfiguration("EMC7+DMC7", 0, hlist))
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("DCal", "DCal", type), "x")
    trgConf.AddHistogram("EMCTRQA_hist{0}EMCGAHMaxVs{1}EMCJEHMax{2}".format("EMCal", "EMCal", type), "x")
    
    triggerAna.DoAnalysis(nevents)
    
    highTriggers = ["EG1", "EJ1", "DG1", "DJ1"]
    lowTriggers = ["EG2", "EJ2", "DG2", "DJ2"]
    
    highTriggersGammaOnly = ["EG1", "DG1"]
    lowTriggersGammaOnly = ["EG2", "DG2"]
    
    highTriggersJetOnly = ["EJ1", "DJ1"]
    lowTriggersJetOnly = ["EJ2", "DJ2"]
    
    EMCalHighTriggers = ["EG1", "EJ1"]
    EMCalLowTriggers = ["EG2", "EJ2"]
    
    DCalHighTriggers = ["DG1", "DJ1"]
    DCalLowTriggers = ["DG2", "DJ2"]
    
    summary = OrderedDict()
    
    summary["highTriggers"] = triggerAna.PrintTriggerSuppression(highTriggers, "EMC7+DMC7", 0)
    summary["lowTriggers"] = triggerAna.PrintTriggerSuppression(lowTriggers, "EMC7+DMC7", summary["highTriggers"].fRate)
    
    summary["highTriggersGammaOnly"] = triggerAna.PrintTriggerSuppression(highTriggersGammaOnly, "EMC7+DMC7", 0)
    summary["lowTriggersGammaOnly"] = triggerAna.PrintTriggerSuppression(lowTriggersGammaOnly, "EMC7+DMC7", summary["highTriggersGammaOnly"].fRate)

    summary["highTriggersJetOnly"] = triggerAna.PrintTriggerSuppression(highTriggersJetOnly, "EMC7+DMC7", 0)
    summary["lowTriggersJetOnly"] = triggerAna.PrintTriggerSuppression(lowTriggersJetOnly, "EMC7+DMC7", summary["highTriggersJetOnly"].fRate)

    summary["EMCalHighTriggers"] = triggerAna.PrintTriggerSuppression(EMCalHighTriggers, "EMC7+DMC7", 0)
    summary["EMCalLowTriggers"] = triggerAna.PrintTriggerSuppression(EMCalLowTriggers, "EMC7+DMC7", summary["EMCalHighTriggers"].fRate)
    
    summary["DCalHighTriggers"] = triggerAna.PrintTriggerSuppression(DCalHighTriggers, "EMC7+DMC7", 0)
    summary["DCalLowTriggers"] = triggerAna.PrintTriggerSuppression(DCalLowTriggers, "EMC7+DMC7", summary["DCalHighTriggers"].fRate)
    
    print("|  *{0: ^20}*  |  *{1}*  |".format("Trigger combination", "Suppression over L0"))
    for trigger in summary.itervalues():
        trigger.Print()
    
def Plot2D(hlist, hname, trigLab1, trigLab2,
           thresholds1, thresholds2, nevents, axis, fname, suffix, maxX, maxY, rebin):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return
    
    canvas = ROOT.TCanvas(hname, hname)
    canvas.SetLogz()
    canvas.cd()
    
    if axis is Axis.GeV:
        hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
        hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
        hist.GetYaxis().SetLimits(hist.GetYaxis().GetXmin() * ADCtoGeV, hist.GetYaxis().GetXmax() * ADCtoGeV)
        hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() + " (GeV)")
    else:
        max /= ADCtoGeV
    
    hist.Sumw2()
    hist.GetXaxis().SetRangeUser(0,maxX)
    hist.GetYaxis().SetRangeUser(0,maxY)
    hist.Draw("colz")
    globalList.append(canvas)
    
    print hname
    totEntries = hist.GetEntries()
    print "Total number of entries is", totEntries
    
    colors = (ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1, ROOT.kOrange+1, ROOT.kGreen+1, ROOT.kMagenta+1, ROOT.kCyan+1, ROOT.kPink+1, ROOT.kTeal+1, ROOT.kYellow+1, ROOT.kSpring+1)
    projList = []
    ratioList = []
    canvasProj = ROOT.TCanvas(hname + "_proj", hname + "_proj")
    canvasProj.SetLogy()
    canvasRatio = ROOT.TCanvas(hname + "_ratio", hname + "_ratio")
    globalList.append(projList)
    globalList.append(ratioList)
    globalList.append(canvasProj)
    globalList.append(canvasRatio)
    
    legend = ROOT.TLegend(0.7, 0.67, 0.88, 0.88, "", "NBNDC")
    legend.SetBorderSize(0)
    globalList.append(legend)
    
    ith = 0
    for th1 in thresholds1:
        if axis is Axis.ADC:
            th1 /= ADCtoGeV
        proj = hist.ProjectionX(hname + "_proj_" + str(th1), hist.GetYaxis().FindBin(th1), hist.GetNbinsY()+1)
        if rebin > 1:
            proj.Rebin(rebin)
        proj.SetTitle("{0} {1} GeV".format(trigLab1, th1))
        projList.append(proj)
        canvasProj.cd()
        proj.SetLineColor(colors[ith])
        proj.SetLineWidth(2)
        proj.SetLineStyle(ith+1)
        
        legend.AddEntry(proj, proj.GetTitle(), "l")

        if ith == 0: proj.Draw("lhist")
        else: proj.Draw("lhistsame")
        
        ratio = proj.Clone(proj.GetName() + "_ratio")
        ratio.Divide(projList[0])
        ratioList.append(ratio)
        canvasRatio.cd()
        
        if ith == 0: ratio.Draw("lhist")
        else: ratio.Draw("lhistsame")
        
        for th2 in thresholds2:
            if axis is Axis.ADC:
                th2 /= ADCtoGeV
            #print("Applying th {0} to {1} and th {2} to {3}".format(th2, hist.GetXaxis().GetTitle(), th1, hist.GetYaxis().GetTitle()))
            entriesBoth = hist.Integral(hist.GetXaxis().FindBin(th2), -1, hist.GetYaxis().FindBin(th1), -1)
            entries1 = hist.Integral(0, -1, hist.GetYaxis().FindBin(th1), -1)
            entries2 = hist.Integral(hist.GetXaxis().FindBin(th2), -1, 0, -1)
            errentriesBoth = math.sqrt(entriesBoth)
            errentries1 = math.sqrt(entries1)
            errentries2 = math.sqrt(entries2)
            share1 = entriesBoth / entries1
            share2 = entriesBoth / entries2
            suppression1 = entries1 / nevents
            suppression2 = entries2 / nevents
            errsuppression1 = suppression1 * (1. / math.sqrt(entries1) + 1. / math.sqrt(nevents))
            errsuppression2 = suppression2 * (1. / math.sqrt(entries2) + 1. / math.sqrt(nevents))
            errshare1 = share1 * (1. / math.sqrt(entriesBoth) + 1. / math.sqrt(entries1))
            errshare2 = share2 * (1. / math.sqrt(entriesBoth) + 1. / math.sqrt(entries2))
            print("{0}:{1}, {2}:{3}, "
                  "ev. {0} = {4:.1f} ({5:.1f}), "
                  "ev. {2} = {8:.1f} ({9:.1f}), "
                  "suppr. {0} = {6:.4f} ({7:.4f}), "
                  "suppr. {2} = {10:.4f} ({11:.4f}), "
                  "share {0} = {12:.4f} ({13:.4f}), "
                  "share {2} = {14:.4f} ({15:.4f})".format(trigLab1, th1, trigLab2, th2, 
                                                           entries1, errentries1,
                                                           suppression1, errsuppression1,
                                                           entries2, errentries2,
                                                           suppression2, errsuppression2,
                                                           share1, errshare1,
                                                           share2, errshare2))
            
        ith += 1
        
    canvasRatio.cd()
    legend.Draw()
    canvasRatio.SaveAs("{0}_Ratio{1}".format(fname, suffix))
    
    canvasProj.cd()
    legend.Draw()
    canvasProj.SaveAs("{0}_Proj{1}".format(fname, suffix))
    
    canvas.SaveAs("{0}{1}".format(fname, suffix))
    
    return canvas

def Plot1D(hlist, hname, thresholds, nevents, color, marker, axis, canvas, opt, max):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return
    
    canvas.cd()
    
    if axis is Axis.GeV:
        hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
        hist.GetXaxis().SetTitle("Energy (GeV)")
    else:
        max /= ADCtoGeV
    
    hist.Sumw2()
    hist.SetMarkerStyle(marker)
    hist.SetMarkerSize(1)
    hist.SetMarkerColor(color)
    hist.SetLineColor(color)
    hist.GetXaxis().SetRangeUser(0,max)
    hist.Draw(opt)
    
    print hname
    totEntries = hist.GetEntries()
    print "Total number of entries is", totEntries
    
    for th in thresholds:
        if axis is Axis.ADC:
            th /= ADCtoGeV
        entries = hist.Integral(hist.GetXaxis().FindBin(th), -1)
        suppression = entries / nevents
        print str(th) + ";" + str(entries) + ";" + str(suppression)
    
    return hist
        
def PlotPatchAmp(hname, variable, offline, recalc, patchtype, patchsize, hlist, nevents, thresholds, eaxis, max):
    canvas = ROOT.TCanvas(variable + patchtype + patchsize, variable + patchtype + patchsize)
    canvas.SetLogy()
    
    opt = ""
    
    legend = ROOT.TLegend(0.7, 0.67, 0.88, 0.88, "", "NBNDC")
    legend.SetBorderSize(0)
    
    if offline >= 0:
        hist = Plot1D(hlist, hname + variable + "EMC" + patchtype + "H" + "Offline", thresholds, nevents, offline, ROOT.kFullCircle, eaxis, canvas, opt, max)
        legend.AddEntry(hist, "Offline", "pe")
        opt = "same"
    if recalc >= 0:
        hist = Plot1D(hlist, hname + variable + "EMC" + patchtype + "H" + "Recalc", thresholds, nevents, recalc, ROOT.kOpenCircle, eaxis, canvas, opt, max)
        legend.AddEntry(hist, "Recalc", "pe")
        
    canvas.cd()
    legend.Draw()
    globalList.append(canvas)
    globalList.append(legend)
        
    return canvas

class PlotNHitsRes:
    def __init__(self, canvas, histogram, badchannels):
        self.fCanvas = canvas
        self.fHistogram = histogram
        self.fBadChannels = badchannels
        
def PlotNHits(hlist, hname, nevents, nhitsTh, label, bins, minBin, maxBin):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return
    
    canvas = ROOT.TCanvas("{0}_{1}".format(label,hname), "{0} {1}".format(label, hist.GetYaxis().GetTitle()))
    canvas.SetLogy()
    
    histNHits = ROOT.TH1D("FastORNHits{0}".format(label), "FastORNHits{0};number of hits per event;probability".format(label), bins, minBin, maxBin)
    
    f = open("badchannels_{0}.txt".format(hname), 'w')

    badchannels = []

    for i in range(1, hist.GetNbinsX()):
        nhits = hist.GetBinContent(i) / nevents
        histNHits.Fill(nhits)
        if nhits > nhitsTh:
            badchannels.append(i-1)
            f.write(str(i-1) + '\n')
        
    f.close()
    canvas.cd()
    histNHits.GetXaxis().SetRangeUser(0,1.01)
    histNHits.Draw()
    globalList.append(canvas)
    globalList.append(histNHits)
    
    return PlotNHitsRes(canvas, histNHits, badchannels)

def GeneratePedestal(hlist, hname, nevents, nhitsTh):
    canvas = ROOT.TCanvas("PedestalCanvas", "PedestalCanvas")
    canvas2 = ROOT.TCanvas("BadChannels", "BadChannels")
    canvas2.SetLogy()
    
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return

    histPedestal = ROOT.TH1D("Pedestal", "Pedestal;FastOR abs. ID;ADC counts", 5000, 0, 5000)
    
    f = open('pedestal.txt', 'w')
    
    iplots = 0;
    
    hist.Sumw2()
    proj = hist.ProjectionY(hist.GetName() + "_projAll", 0, -1)
    canvas2.cd()
    proj.SetMarkerColor(ROOT.kBlack)
    proj.SetLineColor(ROOT.kBlack)
    proj.SetMarkerStyle(ROOT.kFullCircle)
    proj.SetMarkerSize(0.9)
    proj.Scale(1. / proj.Integral())
    proj.Draw()
    proj.GetXaxis().SetRangeUser(0, 150)
    globalList.append(proj)
    
    for i in range(1, hist.GetNbinsX()):
        ped = 0
        for j in range(1, hist.GetNbinsY()):
            nhits = hist.GetBinContent(i, j) / nevents
            if nhits > nhitsTh:
                ped = hist.GetYaxis().GetBinLowEdge(j)
        histPedestal.SetBinContent(i, ped)
        f.write(str(ped) + '\n')

    f.close()
    canvas.cd()
    histPedestal.GetXaxis().SetRangeUser(0, 4500)
    histPedestal.Draw()
    globalList.append(canvas)
    globalList.append(canvas2)
    globalList.append(histPedestal)
    
    return canvas

def ParseThresholds(strTh):
    listTh = []
    listStrTh = strTh.split(",")
    for th in listStrTh:
        listTh.append(int(th))
    return listTh

class TriggerSuppressionCurve:
    def __init__(self, name, title, det, triggerPatch, patchType, threholds):
        self.fName = name
        self.fTitle = title
        self.fDetector = det
        self.fTriggerPatch = triggerPatch
        self.fTriggerPatchType = patchType
        self.fThresholds = threholds
    
    def LoadHistogram(self, file):
        if self.fDetector == "EMCal":
            self.fBaseTrigger = "CEMC7"
            self.fMarkerStyle = ROOT.kFullCircle
        elif self.fDetector == "DCal":
            self.fBaseTrigger = "CDMC7"
            self.fMarkerStyle = ROOT.kOpenCircle
            
        if self.fTriggerPatch == "EMCGAH":
            self.fColor = ROOT.kRed+1
        elif self.fTriggerPatch == "EMCJEH":
            self.fColor = ROOT.kBlue+1
        
        listName = "AliEmcalTriggerQATask_{0}_histos".format(self.fBaseTrigger)
        hlistName = "histosAliEmcalTriggerQATask_{0}".format(self.fBaseTrigger)
        hname = "EMCTRQA_hist{0}MaxPatchAmp{1}{2}".format(self.fDetector, self.fTriggerPatch, self.fTriggerPatchType)
        
        list = file.Get(listName)
        hlist = list.FindObject(hlistName)
        self.fHistogram = hlist.FindObject(hname)
        
    def CalculateSuppression(self):
        self.fGraph = ROOT.TGraphErrors(len(self.fThresholds))
        self.fGraph.SetName(self.fName)
        self.fGraph.SetTitle(self.fTitle)
        self.fGraph.SetMarkerStyle(self.fMarkerStyle)
        self.fGraph.SetMarkerSize(1)
        self.fGraph.SetMarkerColor(self.fColor)
        self.fGraph.SetLineColor(self.fColor)
        
        self.fTotalEvents = self.fHistogram.GetEntries()
        print("Total events are {0}".format(self.fTotalEvents))
        
        for i,th in enumerate(self.fThresholds):
            events = self.fHistogram.Integral(self.fHistogram.GetXaxis().FindBin(th/ADCtoGeV), -1)
            sup = events / self.fTotalEvents
            err = (1. / math.sqrt(events) + 1. / math.sqrt(self.fTotalEvents)) * sup
            print("Suppression for th {0} is {1} ({2})".format(th, sup, err))
            self.fGraph.SetPoint(i, th, sup)
            self.fGraph.SetPointError(i, 0, err)
            
def PlotTriggerSuppression(file, triggerList):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    leg = ROOT.TLegend(0.53,0.58,0.88,0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    canvas = ROOT.TCanvas("TriggerSuppression", "Trigger Suppression")
    blankHist = ROOT.TH1C("blankHist", "blankHist;Threshold (GeV); Suppression over EMCal/DCal L0", 25, 0, 25)
    blankHist.SetMinimum(2e-3)
    blankHist.SetMaximum(6e-1)
    canvas.SetLogy()
    blankHist.Draw()
    globalList.append(blankHist)
    
    for trigger in triggerList:
        trigger.LoadHistogram(file)
        trigger.CalculateSuppression()
        trigger.fGraph.Draw("same p")
        leg.AddEntry(trigger.fGraph, trigger.fTitle, "pe")
        globalList.append(trigger.fGraph)
        
    leg.Draw()
    globalList.append(leg)
    globalList.append(canvas)
        

def main(train, trigger="EMC7", det="EMCal", offline=True, recalc=True, 
         GApatch=True, JEpatch =True, L0vsJEpatch=True, GAvsJEpatch=True, L0vsGApatch=True, 
         pedestal=True, badchannels=True, level0=True, level1=False, triggerSuppression=True,
         plotTriggerSuppression=False, gammaTh="4,5,6,7,8,9,10", jetTh="14,15,16,17,18,19,20",
         axis="ADC", run="", jetsize="16x16", inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    suffix = ""
    if run:
        suffix += "_{0}".format(run)
    if trigger:
        suffix += "_{0}".format(trigger)
    suffix += ".pdf"
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = "{0}/{1}/{2}/AnalysisResults.root".format(inputPath, train, run)
    
    if trigger:
        listName = "AliEmcalTriggerQATask_" + trigger + "_histos"
        hlistName = "histosAliEmcalTriggerQATask_" + trigger
    else:
        listName = "AliEmcalTriggerQATask_histos"
        hlistName = "histosAliEmcalTriggerQATask"
    
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
    
    thresholds_L0 = [ 2.5 ]
    thresholds_GA = ParseThresholds(gammaTh)
    thresholds_JE = ParseThresholds(jetTh)
    
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

    if GApatch:
        if offline:
            offlineColor = ROOT.kRed+1
        else:
            offlineColor = -1
        
        if recalc:
            recalcColor = ROOT.kOrange+1
        else:
            recalcColor = -1
        
        canvas = PlotPatchAmp("EMCTRQA_hist{0}".format(det), "MaxPatchAmp", offlineColor, recalcColor, "GA", "2x2", hlist, nevents, thresholds_GA, eaxis, 50)
        canvas.SaveAs("{0}{1}".format(canvas.GetName(), suffix))
        
        canvas = PlotPatchAmp("EMCTRQA_hist{0}".format(det), "PatchAmp", offlineColor, recalcColor, "GA", "2x2", hlist, nevents, [], eaxis, 50)
        canvas.SaveAs("{0}{1}".format(canvas.GetName(), suffix))
        
    if JEpatch:
        if offline:
            offlineColor = ROOT.kBlue+1
        else:
            offlineColor = -1
        
        if recalc:
            recalcColor = ROOT.kGreen+1
        else:
            recalcColor = -1
            
        canvas = PlotPatchAmp("EMCTRQA_hist{0}".format(det), "MaxPatchAmp", offlineColor, recalcColor, "JE", jetsize, hlist, nevents, thresholds_JE, eaxis, 100)
        canvas.SaveAs("{0}{1}".format(canvas.GetName(), suffix))
        
        canvas = PlotPatchAmp("EMCTRQA_hist{0}".format(det), "PatchAmp", offlineColor, recalcColor, "JE", jetsize, hlist, nevents, [], eaxis, 100)
        canvas.SaveAs("{0}{1}".format(canvas.GetName(), suffix))
    
    if triggerSuppression:
        thresholdsStr = triggerSuppression.split(",")
        thresholds = dict()
        for thStr in thresholdsStr:
            thresholds[thStr.split(":")[0]] = int(thStr.split(":")[1])
        print(thresholds)
        if offline:
            CalculateTriggerSuppression(hlist, "Offline", nevents, thresholds)
            
        if recalc:
            CalculateTriggerSuppression(hlist, "Recalc", nevents, thresholds)
            
    if plotTriggerSuppression:
        if recalc:
            triggerList = []
            triggerList.append(TriggerSuppressionCurve("EMCal_GA", "EMCal GA (4x4)", "EMCal", "EMCGAH", "Recalc", thresholds_GA))
            triggerList.append(TriggerSuppressionCurve("DCal_GA", "DCal GA (4x4)", "DCal", "EMCGAH", "Recalc", thresholds_GA))
            triggerList.append(TriggerSuppressionCurve("EMCal_JE", "EMCal JE (32x32)", "EMCal", "EMCJEH", "Recalc", thresholds_JE))
            triggerList.append(TriggerSuppressionCurve("DCal_JE", "DCal JE (32x32)", "DCal", "EMCJEH", "Recalc", thresholds_JE))
        
            PlotTriggerSuppression(file, triggerList)
    
    if GAvsJEpatch:
        if offline:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCGAHMaxVsEMCJEHMaxOffline".format(det), "GA", "JE", thresholds_GA, thresholds_JE, nevents, eaxis, "MaxGA2x2vsJE{0}".format(jetsize), "_Offline"+suffix, 100, 100, 1)
            
        if recalc:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCGAHMaxVsEMCJEHMaxRecalc".format(det), "GA", "JE", thresholds_GA, thresholds_JE, nevents, eaxis, "MaxGA2x2vsJE{0}".format(jetsize), "_Recalc"+suffix, 100, 100, 1)
            
    if L0vsJEpatch:
        if offline:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCL0MaxVsEMCJEHMaxOffline".format(det), "L0", "JE", thresholds_L0, thresholds_JE, nevents, eaxis, "MaxL02x2vsJE{0}".format(jetsize), "_Offline"+suffix, 100, 100, 1)
            
        if recalc:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCL0MaxVsEMCJEHMaxRecalc".format(det), "L0", "JE", thresholds_L0, thresholds_JE, nevents, eaxis, "MaxL02x2vsJE{0}".format(jetsize), "_Recalc"+suffix, 100, 100, 1)
            
    if L0vsGApatch:
        if offline:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCL0MaxVsEMCGAHMaxOffline".format(det), "L0", "GA", thresholds_L0, thresholds_GA, nevents, eaxis, "MaxL02x2vsGA2x2", "_Offline"+suffix, 20, 20, 1)
            
        if recalc:
            canvas = Plot2D(hlist, "EMCTRQA_hist{0}EMCL0MaxVsEMCGAHMaxRecalc".format(det), "L0", "GA", thresholds_L0, thresholds_GA, nevents, eaxis, "MaxL02x2vsGA2x2", "_Recalc"+suffix, 20, 20, 1)
            
    if badchannels:
        if level0:
            L0badchannelsRes = PlotNHits(hlist, "EMCTRQA_histFastORL0", nevents, 3e-3, "Level0", 1500, 0, 1.5)
            L0badchannelsRes.fCanvas.SaveAs("{0}{1}".format(L0badchannelsRes.fCanvas.GetName(), suffix))
            
            L0badchannelsLargeRes = PlotNHits(hlist, "EMCTRQA_histLargeAmpFastORL0", nevents, 7e-7, "Level0", 1000, 0, 0.002)
            L0badchannelsLargeRes.fCanvas.SaveAs("{0}{1}".format(L0badchannelsLargeRes.fCanvas.GetName(), suffix))
            
            L0badchannels = set(L0badchannelsRes.fBadChannels)
            L0badchannels.update(L0badchannelsLargeRes.fBadChannels)
            
            f = open("L0badchannels.txt", 'w')
            
            for i in L0badchannels:
                 f.write(str(i) + '\n')
                 
        if level1:
            L1badchannelsRes = PlotNHits(hlist, "EMCTRQA_histFastORL1", nevents, 3e-3, "Level1", 1500, 0, 1.5)
            L1badchannelsRes.fCanvas.SaveAs("{0}{1}".format(L1badchannelsRes.fCanvas.GetName(), suffix))
            
            L1badchannelsLargeRes = PlotNHits(hlist, "EMCTRQA_histLargeAmpFastORL1", nevents, 7e-7, "Level1", 1000, 0, 0.002)
            L1badchannelsLargeRes.fCanvas.SaveAs("{0}{1}".format(L1badchannelsLargeRes.fCanvas.GetName(), suffix))
            
            L1badchannels = set(L1badchannelsRes.fBadChannels)
            L1badchannels.update(L1badchannelsLargeRes.fBadChannels)
            
            f = open("L1badchannels.txt", 'w')
            
            for i in L1badchannels:
                 f.write(str(i) + '\n')
        
            
    if pedestal:
        canvas = GeneratePedestal(hlist, "EMCTRQA_histFastORL0Amp", nevents, 0.01)
        canvas.SaveAs("{0}{1}".format(canvas.GetName(), suffix))
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--trigger', metavar='trigger',
                        default="EMC7",
                        help='Trigger')
    parser.add_argument('--detector', metavar='trigger',
                        default="EMCal",
                        help='Detector (EMCal or DCal)')
    parser.add_argument('--offline', action='store_const',
                        default=False, const=True,
                        help='Offline')
    parser.add_argument('--recalc', action='store_const',
                        default=False, const=True,
                        help='Offline')
    parser.add_argument('--GApatch', action='store_const',
                        default=False, const=True,
                        help='GA patch')
    parser.add_argument('--JEpatch', action='store_const',
                        default=False, const=True,
                        help='JE patch')
    parser.add_argument('--GAvsJEpatch', action='store_const',
                        default=False, const=True,
                        help='GA vs JE patch')
    parser.add_argument('--L0vsJEpatch', action='store_const',
                        default=False, const=True,
                        help='L0 vs JE patch')
    parser.add_argument('--L0vsGApatch', action='store_const',
                        default=False, const=True,
                        help='L0 vs GA patch')
    parser.add_argument('--trigger-suppression', metavar='triggerSuppression',
                        default="",
                        help='Trigger suppression analysis: thresholds G1,J1,G2,J2')
    parser.add_argument('--plot-trigger-suppression', action='store_const',
                        default=False, const=True,
                        help='Plot Trigger suppression as a funtion of threshold')
    parser.add_argument('--gamma', metavar='gamma',
                        default="",
                        help='Threholds for gamma patches')
    parser.add_argument('--jet', metavar='jet',
                        default="",
                        help='Threholds for jet patches')
    parser.add_argument('--pedestal', action='store_const',
                        default=False, const=True,
                        help='Run pedestal analysis')
    parser.add_argument('--badchannels', action='store_const',
                        default=False, const=True,
                        help='Run bad channel analysis')
    parser.add_argument('--level0', action='store_const',
                        default=False, const=True,
                        help='Run bad channel analysis for L0')
    parser.add_argument('--level1', action='store_const',
                        default=False, const=True,
                        help='Run bad channel analysis for L1')
    parser.add_argument('--axis', metavar='axis',
                        default="ADC",
                        help='Axis type (GeV or ADC)')
    parser.add_argument('--run', metavar='run',
                        default="",
                        help='Run number')
    parser.add_argument('--size', metavar='size',
                        default="16x16",
                        help='Jet patch size')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.trigger, args.detector, args.offline, args.recalc, 
         args.GApatch, args.JEpatch, args.L0vsJEpatch, args.GAvsJEpatch, args.L0vsGApatch, 
         args.pedestal, args.badchannels, args.level0, args.level1, args.trigger_suppression,
         args.plot_trigger_suppression, args.gamma, args.jet,
         args.axis, args.run, args.size, args.input_path)
    
    IPython.embed()