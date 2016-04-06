#!/usr/bin/env python
#python script to do HLT EMCal trigger QA plots

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math

import IPython

globalList = []

def GetHistogramListOfflineAnalysis(fileName, listTaskName="AliEmcalTriggerQATask_histos", listName="AliEmcalTriggerQATask"):
    file = ROOT.TFile(fileName)
    if not file or file.IsZombie():
        print("Could not find file '{0}".format(filename))
        return
    listTask = file.Get(listTaskName)
    if not listTask:
        print("Could not find list '{0}'".format(listTaskName))
        return
    list = listTask.FindObject(listName)
    if not list:
        print("Could not find list '{0}'".format(listName))
        return
    
    histList = dict()
    next = ROOT.TIter(list)
    obj = next()
    while obj:
        histList[obj.GetName()] = obj
        obj = next()
    
    return histList

def PlotbySM(histList, histNamePattern, minSM, maxSM, events, det, plotTitle, plotName, period, logz):
    n = maxSM - minSM + 1 #number of canvases
    ncol = 2
    nrow = int(math.ceil(n / ncol))
    canvas = ROOT.TCanvas("{0}_{1}_{2}".format(period, plotName, det), "{0} - {1} - {2}".format(period, plotTitle, det), 780, 180*nrow)

    canvas.Divide(ncol, nrow)
    
    for i in range(minSM, maxSM+1):
        padcol = i % 2
        padrow = (maxSM - i) / 2
        padid = padcol + padrow * 2 + 1
        pad = canvas.cd(padid)
        pad.SetLogz(logz)
        pad.SetLeftMargin(0.09)
        pad.SetRightMargin(0.18)
        pad.SetBottomMargin(0.15)
        pad.SetTopMargin(0.1)
        histName = "{0}_SM{1}".format(histNamePattern, i)
        hist = histList[histName]
        if not hist:
            print("Could not find histogram '{0}'".format(histName))
            continue
        hist.GetXaxis().SetTitleFont(43)
        hist.GetXaxis().SetTitleSize(11)
        hist.GetXaxis().SetLabelFont(43)
        hist.GetXaxis().SetLabelSize(11)
        hist.GetXaxis().SetTitleOffset(3.6)
        hist.GetYaxis().SetTitleFont(43)
        hist.GetYaxis().SetTitleSize(11)
        hist.GetYaxis().SetLabelFont(43)
        hist.GetYaxis().SetLabelSize(11)
        hist.GetYaxis().SetTitleOffset(2.0)
        hist.GetZaxis().SetTitleFont(43)
        hist.GetZaxis().SetTitleSize(11)
        hist.GetZaxis().SetLabelFont(43)
        hist.GetZaxis().SetLabelSize(11)
        hist.GetZaxis().SetTitleOffset(3.4)
        hist.Scale(1. / events)
        hist.Draw("colz")
        
        t = ROOT.TPaveText(0.05, 0.9, 0.95, 1.0, "NB NDC")
        t.AddText("SM {0}".format(i))
        t.SetTextAlign(22)
        t.SetTextFont(43)
        t.SetTextSize(16)
        t.SetFillStyle(0)
        t.SetBorderSize(0)
        t.Draw()
        globalList.append(t)
    
    globalList.append(canvas)
    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))
    return canvas

def GetEvents(histList, histName):
    hist = histList[histName]
    if not hist:
        print("Could not find histogram '{0}'".format(histName))
        return
    return hist.GetBinContent(1)

def Plot1D(histList, histName, events, th, plotTitle, plotName, period):
    hist = histList[histName]
    if not hist:
        print("Could not find histogram '{0}'".format(histName))
        return
    
    hist.Sumw2()
    hist.Scale(1. / events)
    
    canvas = ROOT.TCanvas("{0}_{1}".format(period, plotName), plotTitle, 800, 600)
    canvas.cd()
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.SetMarkerSize(0.8)
    hist.SetMarkerColor(ROOT.kBlue+1)
    hist.SetLineColor(ROOT.kBlue+1)
    hist.Draw()
    
    absIdList = []
    
    for ibin in range(1, hist.GetXaxis().GetNbins()+1):
        if hist.GetBinContent(ibin) > th:
            absIdList.append(ibin-1)
            
    pave = ROOT.TPaveText(0.1, 0.7, 0.9, 0.2, "NB NDC")
    pave.SetTextAlign(13)
    pave.SetTextFont(43)
    pave.SetTextSize(12)
    pave.SetFillStyle(0)
    pave.SetTextColor(ROOT.kRed)
    pave.SetBorderSize(0)
    
    absIdText = ""
    
    for absId in absIdList:
        if absIdText:
            absIdText = "{0}, {1}".format(absIdText, absId)
        else:
            absIdText = "{0}".format(absId)
        if len(absIdText) > 110:
            pave.AddText(absIdText)
            print absIdText
            absIdText = ""
        
    pave.Draw()
    
    globalList.append(pave)
    globalList.append(canvas)
    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))
    return canvas
        
def main(fileName, offline, period):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    if offline:
        histList = GetHistogramListOfflineAnalysis(fileName)
    else:
        histList = GetHistogramListOnlineAnalysis(fileName)
        
    events = GetEvents(histList, "EMCTRQA_histEvents")
    
    PlotbySM(histList, "EMCTRQA_histFastORL0", 0, 5, events, "EMCAL1", "Hit frequency", "HitFrequency", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0", 6, 11, events, "EMCAL2", "Hit frequency", "HitFrequency", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0", 12, 19, events, "DCAL", "Hit frequency", "HitFrequency", period, False)
    
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 0, 5, events, "EMCAL1", "Hit frequency (amp > 400)", "HitFrequencyLargeAmp", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 6, 11, events, "EMCAL2", "Hit frequency (amp > 400)", "HitFrequencyLargeAmp", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 12, 19, events, "DCAL", "Hit frequency (amp > 400)", "HitFrequencyLargeAmp", period, False)
    
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 0, 5, events, "EMCAL1", "Average amplitude", "AverageAmplitude", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 6, 11, events, "EMCAL2", "Average amplitude", "AverageAmplitude", period, False)
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 12, 19, events, "DCAL", "Average amplitude", "AverageAmplitude", period, False)
    
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 0, 5, events, "EMCAL1", "FEE vs TRU", "FEEvsTRU", period, True)
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 6, 11, events, "EMCAL2", "FEE vs TRU", "FEEvsTRU", period, True)
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 12, 19, events, "DCAL", "FEE vs TRU", "FEEvsTRU", period, True)
    
    Plot1D(histList, "EMCTRQA_histFastORL0", events, 3e-3, "Hit frequency", "HitFrequency", period)
    Plot1D(histList, "EMCTRQA_histFastORL0LargeAmp", events, 7e-7, "Hit frequency (amp > 400)", "HitFrequencyLargeAmp", period)
    Plot1D(histList, "EMCTRQA_histFastORL0Amp", events, 10000, "Average amplitude", "AverageAmplitude", period)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='HLT EMCal trigger QA.')
    parser.add_argument('fileName', metavar='fileName',
                        help='Name of the output file')
    parser.add_argument('--offline', action='store_const',
                        default=False, const=True,
                        help='Offline')
    parser.add_argument('-p', '--period',
                        help='Run period (e.g. LHC16b)')
    args = parser.parse_args()
    
    main(args.fileName, args.offline, args.period)
    
    IPython.embed()