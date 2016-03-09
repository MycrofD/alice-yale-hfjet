#!/usr/bin/env python
#python script to do EMCal trigger QA plots

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array

import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4 

class Axis(Enum):
    GeV = 1
    ADC = 2
    
def Plot2D(hlist, hname, thresholds1, thresholds2, nevents, axis):
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return
    
    canvas = ROOT.TCanvas(hname, hname)
    canvas.SetLogz()
    canvas.cd()
    
    max = 100
    
    if axis is Axis.GeV:
        hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
        hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
        hist.GetYaxis().SetLimits(hist.GetYaxis().GetXmin() * ADCtoGeV, hist.GetYaxis().GetXmax() * ADCtoGeV)
        hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() + " (GeV)")
    else:
        max /= ADCtoGeV
    
    hist.Sumw2()
    hist.GetXaxis().SetRangeUser(0,max)
    hist.GetYaxis().SetRangeUser(0,max)
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
    
    ith = 0
    for th1 in thresholds1:
        if axis is Axis.ADC:
            th1 /= ADCtoGeV
        
        proj = hist.ProjectionX(hname + "_proj_" + str(th1), hist.GetYaxis().FindBin(th1), hist.GetNbinsY()+1)
        proj.Rebin(11)
        proj.SetTitle("L0 " + str(th1) + " GeV")
        projList.append(proj)
        canvasProj.cd()
        proj.SetLineColor(colors[ith])
        proj.SetLineWidth(2)
        proj.SetLineStyle(ith+1)

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
                
            entries = hist.Integral(hist.GetXaxis().FindBin(th2), -1, hist.GetYaxis().FindBin(th1), -1)
            suppression = entries / nevents
            print str(th1) + ";" + str(th2) + ";" + str(entries) + ";" + str(suppression)
            
        ith += 1
    
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

def PlotNHits(hlist, hname, nevents, nhitsTh):
    canvas = ROOT.TCanvas("FastORNHits", "FastORNHits")
    canvas.SetLogy()
    hist = hlist.FindObject(hname)
    if not hist:
        print "Could not get histogram '" + hname + "' in list '" + hlist.GetName() + "'. Skipping..."
        hlist.Print()
        return
    
    histNHits = ROOT.TH1D("FastORNHits", "FastORNHits;number of hits per event;probability", 1500, 0, 1.5)
    
    f = open('badchannels.txt', 'w')

    for i in range(1, hist.GetNbinsX()):
        nhits = hist.GetBinContent(i) / nevents
        histNHits.Fill(nhits)
        if nhits > nhitsTh:
            f.write(str(i-1) + '\n')
        
    f.close()
    canvas.cd()
    histNHits.GetXaxis().SetRangeUser(0,1.01)
    histNHits.Draw()
    globalList.append(canvas)
    globalList.append(histNHits)
    
    return canvas

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

def main(train, trigger="EMC7", offline=True, recalc=True, GApatch=True, JEpatch =True, L0vsJEpatch=True, GAvsJEpatch=True, pedestal=True,
         axis="ADC", run="237673", jetsize="16x16", inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = inputPath + "/" + train + "/AnalysisResults"
    if run:
        fileName += "_" + run 
    if jetsize:
        fileName += "_" + jetsize    
    fileName += ".root"
    
    if trigger:
        listName = "AliEmcalTriggerQATaskPP_" + trigger + "_histos"
        hlistName = "histosAliEmcalTriggerQATaskPP_" + trigger
    else:
        listName = "AliEmcalTriggerQATaskPP_histos"
        hlistName = "histosAliEmcalTriggerQATaskPP"
    
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
    
    thresholds_GA = [0, 3, 4, 6, 8]
    thresholds_JE = [5, 10, 15, 20]
    
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
        
        canvas = PlotPatchAmp("EMCTRQA_histEMCal", "MaxPatchAmp", offlineColor, recalcColor, "GA", "2x2", hlist, nevents, thresholds_GA, eaxis, 50)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
        
        canvas = PlotPatchAmp("EMCTRQA_histEMCal", "PatchAmp", offlineColor, recalcColor, "GA", "2x2", hlist, nevents, [], eaxis, 50)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
        
    if JEpatch:
        if offline:
            offlineColor = ROOT.kBlue+1
        else:
            offlineColor = -1
        
        if recalc:
            recalcColor = ROOT.kGreen+1
        else:
            recalcColor = -1
            
        canvas = PlotPatchAmp("EMCTRQA_histEMCal", "MaxPatchAmp", offlineColor, recalcColor, "JE", jetsize, hlist, nevents, thresholds_JE, eaxis, 100)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
        
        canvas = PlotPatchAmp("EMCTRQA_histEMCal", "PatchAmp", offlineColor, recalcColor, "JE", jetsize, hlist, nevents, [], eaxis, 100)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
    
    if GAvsJEpatch:
        if offline:
            canvas = Plot2D(hlist, "EMCTRQA_histEMCalEMCGAHMaxVsEMCJEHMaxOffline", thresholds_GA, thresholds_JE, nevents, eaxis)
            canvas.SaveAs("MaxGA2x2vsJE" + jetsize + "_"+run+"_"+trigger+"_Offline.pdf")
            
        if recalc:
            canvas = Plot2D(hlist, "EMCTRQA_histEMCalEMCGAHMaxVsEMCJEHMaxRecalc", thresholds_GA, thresholds_JE, nevents, eaxis)
            canvas.SaveAs("MaxGA2x2vsJE" + jetsize + "_"+run+"_"+trigger+"_Recalc.pdf")
            
    if L0vsJEpatch:
        if offline:
            canvas = Plot2D(hlist, "EMCTRQA_histEMCalEMCL0MaxVsEMCJEHMaxOffline", thresholds_JE, thresholds_GA, nevents, eaxis)
            canvas.SaveAs("MaxGA2x2vsJE" + jetsize + "_"+run+"_"+trigger+"_Offline.pdf")
            
        if recalc:
            canvas = Plot2D(hlist, "EMCTRQA_histEMCalEMCL0MaxVsEMCJEHMaxRecalc", thresholds_JE, thresholds_GA, nevents, eaxis)
            canvas.SaveAs("MaxGA2x2vsJE" + jetsize + "_"+run+"_"+trigger+"_Recalc.pdf")
            
    if pedestal:
        canvas = PlotNHits(hlist, "EMCTRQA_histFastORL0", nevents, 0.01)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
        
        canvas = GeneratePedestal(hlist, "EMCTRQA_histFastORL0Amp", nevents, 0.01)
        canvas.SaveAs(canvas.GetName() + "_" + run + "_" + trigger + ".pdf")
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--trigger', metavar='trigger',
                        default="EMC7",
                        help='Trigger')
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
    parser.add_argument('--pedestal', action='store_const',
                        default=False, const=True,
                        help='Run pedestal analysis')
    parser.add_argument('--axis', metavar='axis',
                        default="ADC",
                        help='Axis type (GeV or ADC)')
    parser.add_argument('--run', metavar='run',
                        default="",
                        help='Run number')
    parser.add_argument('--size', metavar='size',
                        default="",
                        help='Jet patch size')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.trigger, args.offline, args.recalc, args.GApatch, args.JEpatch, args.L0vsJEpatch, args.GAvsJEpatch, args.pedestal, args.axis, args.run, args.size, args.input_path)
    
    IPython.embed()