#!/usr/local/bin/python
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

def GetHistogramListOnlineAnalysis(fileName):
    file = ROOT.TFile(fileName)
    if not file or file.IsZombie():
        print("Could not find file '{0}".format(filename))
        return
    
    histList = dict()
    for key in file.GetListOfKeys():
        histList[key.GetName()] = file.Get(key.GetName())
    
    return histList

def PlotbySM(histList, histNamePattern, minSM, maxSM, events, det, plotTitle, plotName, period, logz,
             minX = 0, maxX = 0, minY = 0, maxY = 0):
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
        
        if minX < maxX:
            hist.GetXaxis().SetRangeUser(minX, maxX)
        if minY < maxY:
            hist.GetYaxis().SetRangeUser(minY, maxY)
        
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
    
    canvas = ROOT.TCanvas("{0}_{1}".format(period, plotName), "{0} - {1}".format(period, plotTitle), 800, 600)
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

def PlotMaxPatch(histList, patchType, triggerType, events, plotTitle, plotName, period):
    detectors = ["EMCal", "DCal"]
    colors = [ROOT.kRed+1, ROOT.kBlue+1]
    markers = [ROOT.kFullCircle, ROOT.kOpenCircle]
    options = ["", "same"]
    canvas = ROOT.TCanvas("{0}_{1}".format(period, plotName), "{0} - {1}".format(period, plotTitle))
    canvas.SetLogy(True)
    canvas.cd()
    
    legend = ROOT.TLegend(0.6, 0.9, 0.9, 0.7)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    for det,col,marker,opt in zip(detectors,colors,markers,options):
        histName = "EMCTRQA_hist{0}MaxPatchAmp{1}{2}".format(det, triggerType, patchType)
        hist = histList[histName]
        hist.Sumw2()
        hist.SetMarkerSize(0.8)
        hist.SetMarkerStyle(marker)
        hist.SetLineColor(col)
        hist.SetMarkerColor(col)
        hist.Scale(1./events)
        hist.GetYaxis().SetTitle("entries / events")
        legend.AddEntry(hist, det, "pe")
        hist.Draw(opt)
        
    legend.Draw()
    globalList.append(legend)
    globalList.append(canvas)
    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))   
    return canvas

def Plot2D(histList, histNamePrefix, patchType, triggerType, events, plotTitle, plotName, period):
    canvas = ROOT.TCanvas("{0}_{1}".format(period, plotName), "{0} - {1}".format(period, plotTitle))
    canvas.cd()
    
    histName = "{0}{1}{2}".format(histNamePrefix, triggerType, patchType)
    hist = histList[histName]
    hist.Scale(1./events)
    hist.GetZaxis().SetTitle("entries / events")
    hist.Draw("colz")
    globalList.append(canvas)
    canvas.SaveAs("{0}.pdf".format(canvas.GetName()))   
    return canvas
    
def SaveRootFile(histList, period):
    fname = "AnalysisResults_{0}.root".format(period)
    file = ROOT.TFile(fname, "recreate")
    file.cd()
    for name, hist in histList.iteritems():
        hist.Write()
    file.Close()
    
def main(fileName, offline, period, saveRoot):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.TH1.AddDirectory(False)
    
    if offline:
        histList = GetHistogramListOfflineAnalysis(fileName)
    else:
        histList = GetHistogramListOnlineAnalysis(fileName)
    
    if saveRoot:
        SaveRootFile(histList, period)
        
    events = GetEvents(histList, "EMCTRQA_histEvents")
    
    for i in range(0,2):
        PlotbySM(histList, "EMCTRQA_histFastORL{0}".format(i), 0, 5, events, "EMCAL1", "L{0} hit frequency by SM".format(i), "L{0}HitFrequency".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}".format(i), 6, 11, events, "EMCAL2", "L{0} hit frequency by SM".format(i), "L{0}HitFrequency".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}".format(i), 12, 19, events, "DCAL", "L{0} hit frequency by SM".format(i), "L{0}HitFrequency".format(i), period, False)
        
        PlotbySM(histList, "EMCTRQA_histFastORL{0}LargeAmp".format(i), 0, 5, events, "EMCAL1", "L{0} hit frequency (amp > 400) by SM".format(i), "L{0}HitFrequencyLargeAmp".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}LargeAmp".format(i), 6, 11, events, "EMCAL2", "L{0} hit frequency (amp > 400) by SM".format(i), "L{0}HitFrequencyLargeAmp".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}LargeAmp".format(i), 12, 19, events, "DCAL", "L{0} hit frequency (amp > 400) by SM".format(i), "L{0}HitFrequencyLargeAmp".format(i), period, False)
        
        PlotbySM(histList, "EMCTRQA_histFastORL{0}Amp".format(i), 0, 5, events, "EMCAL1", "L{0} Average amplitude by SM".format(i), "L{0}AverageAmplitude".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}Amp".format(i), 6, 11, events, "EMCAL2", "L{0} Average amplitude by SM".format(i), "L{0}AverageAmplitude".format(i), period, False)
        PlotbySM(histList, "EMCTRQA_histFastORL{0}Amp".format(i), 12, 19, events, "DCAL", "L{0} Average amplitude by SM".format(i), "L{0}AverageAmplitude".format(i), period, False)
        
        Plot1D(histList, "EMCTRQA_histFastORL{0}".format(i), events, 3e-3, "L{0} Hit frequency".format(i), "L{0}HitFrequency".format(i), period)
        Plot1D(histList, "EMCTRQA_histFastORL{0}LargeAmp".format(i), events, 7e-7, "L{0} Hit frequency (amp > 400)".format(i), "L{0}HitFrequencyLargeAmp".format(i), period)
        Plot1D(histList, "EMCTRQA_histFastORL{0}Amp".format(i), events, 10000, "L{0} Average amplitude".format(i), "L{0}AverageAmplitude".format(i), period)
 
    
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 0, 5, events, "EMCAL1", "FEE vs TRU", "FEEvsTRU", period, True, 0, 250, 0, 20)
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 6, 11, events, "EMCAL2", "FEE vs TRU", "FEEvsTRU", period, True, 0, 250, 0, 20)
    PlotbySM(histList, "EMCTRQA_histFEEvsTRU", 12, 19, events, "DCAL", "FEE vs TRU", "FEEvsTRU", period, True, 0, 250, 0, 20)

    PlotbySM(histList, "EMCTRQA_histFEEvsSTU", 0, 5, events, "EMCAL1", "FEE vs STU", "FEEvsSTU", period, True, 0, 250, 0, 20)
    PlotbySM(histList, "EMCTRQA_histFEEvsSTU", 6, 11, events, "EMCAL2", "FEE vs STU", "FEEvsSTU", period, True, 0, 250, 0, 20)
    PlotbySM(histList, "EMCTRQA_histFEEvsSTU", 12, 19, events, "DCAL", "FEE vs STU", "FEEvsSTU", period, True, 0, 250, 0, 20)
   
    PlotMaxPatch(histList, "Recalc", "EMCL0", events, "L0 trigger max patch (recalc)", "L0MaxPatchRecalc", period)
    PlotMaxPatch(histList, "Recalc", "EMCGAH", events, "L1 GA trigger max patch (recalc)", "L1GAHMaxPatchRecalc", period)
    PlotMaxPatch(histList, "Recalc", "EMCJEH", events, "L1 JE trigger max patch (recalc)", "L1JEHMaxPatchRecalc", period)

    Plot2D(histList, "EMCTRQA_histMaxEdgePos", "Recalc", "EMCL0", events, "L0 trigger max patch position (recalc)", "L0MaxPatchRecalcPosition", period)
    Plot2D(histList, "EMCTRQA_histMaxEdgePos", "Recalc", "EMCGAH", events, "L1 GA trigger max patch position (recalc)", "L1GAHMaxPatchRecalcPosition", period)
    Plot2D(histList, "EMCTRQA_histMaxEdgePos", "Recalc", "EMCJEH", events, "L1 JE trigger max patch position (recalc)", "L1JEHMaxPatchRecalcPosition", period)

    Plot2D(histList, "EMCTRQA_histAmpEdgePos", "Recalc", "EMCL0", events, "L0 trigger average amplitude (recalc)", "L0AverageAmpRecalcPosition", period)
    Plot2D(histList, "EMCTRQA_histAmpEdgePos", "Recalc", "EMCGAH", events, "L1 GA trigger average amplitude (recalc)", "L1GAHAverageAmpRecalcPosition", period)
    Plot2D(histList, "EMCTRQA_histAmpEdgePos", "Recalc", "EMCJEH", events, "L1 JE trigger average amplitude (recalc)", "L1JEHAverageAmpRecalcPosition", period)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='HLT EMCal trigger QA.')
    parser.add_argument('fileName', metavar='fileName',
                        help='Name of the output file')
    parser.add_argument('--offline', action='store_const',
                        default=False, const=True,
                        help='Offline')
    parser.add_argument('-p', '--period',
                        default='LHC15j',
                        help='Run period (e.g. LHC16b)')
    parser.add_argument('--save-root', action='store_const',
                        default=False, const=True,
                        help='Save a ROOT file with the histograms.')
    args = parser.parse_args()
    
    main(args.fileName, args.offline, args.period, args.save_root)
    
    IPython.embed()
