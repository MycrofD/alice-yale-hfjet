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

def PlotbySM(histList, histNamePattern, minSM, maxSM, det, plotTitle):
    n = maxSM - minSM + 1 #number of canvases
    ncol = 2
    nrow = int(math.ceil(n / ncol))
    canvas = ROOT.TCanvas("{0}_{1}".format(histNamePattern, det), "{0} - {1}".format(plotTitle, det), 780, 180*nrow)

    canvas.Divide(ncol, nrow)
    
    for i in range(minSM, maxSM+1):
        padcol = i % 2
        padrow = (maxSM - i) / 2
        padid = padcol + padrow * 2 + 1
        pad = canvas.cd(padid)
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
        
def main(fileName, offline):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    if offline:
        histList = GetHistogramListOfflineAnalysis(fileName)
    else:
        histList = GetHistogramListOnlineAnalysis(fileName)
    
    PlotbySM(histList, "EMCTRQA_histFastORL0", 0, 5, "EMCAL1", "Hit frequency")
    PlotbySM(histList, "EMCTRQA_histFastORL0", 6, 11, "EMCAL2", "Hit frequency")
    PlotbySM(histList, "EMCTRQA_histFastORL0", 12, 19, "DCAL", "Hit frequency")
    
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 0, 5, "EMCAL1", "Hit frequency (amp > 400)")
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 6, 11, "EMCAL2", "Hit frequency (amp > 400)")
    PlotbySM(histList, "EMCTRQA_histFastORL0LargeAmp", 12, 19, "DCAL", "Hit frequency (amp > 400)")
    
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 0, 5, "EMCAL1", "Integrated amplitude")
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 6, 11, "EMCAL2", "Integrated amplitude")
    PlotbySM(histList, "EMCTRQA_histFastORL0Amp", 12, 19, "DCAL", "Integrated amplitude")
    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='HLT EMCal trigger QA.')
    parser.add_argument('fileName', metavar='fileName',
                        help='Name of the output file')
    parser.add_argument('--offline', action='store_const',
                        default=False, const=True,
                        help='Offline')
    args = parser.parse_args()
    
    main(args.fileName, args.offline)
    
    IPython.embed()