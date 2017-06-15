#!/usr/bin/env python
#python script to do study effect of L0 threshold on the L1 triggers

import ROOT
import argparse
import math
import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4 

def MakeL1vsL0(title, hist, L0thresholds):
    baseline = hist.ProjectionX("{0}_nocut".format(title))
    baseline.SetTitle("L0 patch > 0 GeV")
    baseline.GetYaxis().SetTitle("counts")
    baseline.GetXaxis().SetRangeUser(0,30)
    #baseline.Rebin(3)
    baseline.Sumw2()
    histograms = []
    for th in L0thresholds:
        h = hist.ProjectionX("{0}_{1}".format(title, int(th*10)), hist.GetYaxis().FindBin(th), hist.GetYaxis().GetNbins())
        #h.Rebin(3)
        h.Sumw2()
        h.GetXaxis().SetRangeUser(0,30)
        h.SetTitle("L0 patch > {0:.2f} GeV".format(th))
        h.GetYaxis().SetTitle("counts")
        histograms.append(h)
    result = CompareSpectra(baseline, histograms, title, "", "P", "L0 threshold X / no L0 threshold")
    for obj in result:
        globalList.append(obj)
        if isinstance(obj, ROOT.TCanvas):
            obj.SaveAs("~/{0}.pdf".format(obj.GetName()))

def main(filename):
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)

    taskName = "AliEmcalTriggerQATask"
    triggerClassName = "INT7"
    histogramNames = {"FEE_L0vsJE": "EMCTRQA_histEMCalEMCL0MaxVsEMCalEMCJEHMaxOffline", 
                      "FOR_GAvsJE": "EMCTRQA_histEMCalEMCGAHMaxVsEMCalEMCJEHMaxRecalc"}
    L0thresholds = [2.5, 3.5, 4.0, 4.5, 5.0]

    file = ROOT.TFile(filename)
    if not file or file.IsZombie():
        print("Could not open file {0}!".format(filename))
        exit(1)
    else:
        print("File {0} successfully open!".format(filename))

    listName = "_".join([taskName, triggerClassName, "histos"])
    datalist = file.Get(listName)
    if not datalist:
        print("Could not get list {0} from file {1}!".format(listName, filename))
        file.ls()
        exit(1)
    else:
        print("List {0} successfully retrieved from file {1}!".format(listName, filename))

    hlistName = "histos{0}_{1}".format(taskName, triggerClassName)
    hdatalist = datalist.FindObject(hlistName)
    if not hdatalist:
        print("Could not get list {0}/{1} from file {2}!".format(listName, hlistName, filename))
        datalist.ls()
        exit(1)
    else:
        print("List {0} successfully retrieved from list {1}!".format(hlistName, listName))

    for title, histName in histogramNames.iteritems():
        hist = hdatalist.FindObject(histName)
        if not hist:
            print("Could not get histogram {0} from list {1}!".format(histName, hlistName))
            datalist.ls()
            exit(1)
        else:
            print("Histogram {0} successfully retrieved from list {1}!".format(histName, hlistName))
        if not "GeV" in hist.GetXaxis().GetTitle():
            hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
            hist.GetXaxis().SetTitle(hist.GetXaxis().GetTitle() + " (GeV)")
            hist.GetYaxis().SetLimits(hist.GetYaxis().GetXmin() * ADCtoGeV, hist.GetYaxis().GetXmax() * ADCtoGeV)
            hist.GetYaxis().SetTitle(hist.GetYaxis().GetTitle() + " (GeV)")
        MakeL1vsL0(title, hist, L0thresholds)

def CompareSpectra(baseline, spectra, comparisonName, opt="", optRatio="", yaxisRatio="ratio", doSpectra="logy", doRatio="lineary", c=None, cRatio=None, leg=None, legRatio=None, styles=None):
    results = []
    baselineRatio = None
    mainRatioHist = None
    mainHist = None
    
    print("CompareSpectra: {0}".format(comparisonName))
    print(baseline.GetName())
    for s in spectra:
        print(s.GetName())

    if styles:
        colors = styles["colors"]
        markers = styles["markers"]
        lines = styles["lines"]        
    else:
        colors = [ROOT.kBlack, ROOT.kBlue+2, ROOT.kRed+2, ROOT.kGreen+2, ROOT.kOrange+2, ROOT.kAzure+2, ROOT.kMagenta+2, ROOT.kCyan+2]
        markers = [ROOT.kOpenCircle, ROOT.kFullCircle, ROOT.kFullSquare, ROOT.kFullTriangleUp, ROOT.kFullTriangleDown, ROOT.kFullDiamond, ROOT.kFullStar]
        lines = [1, 2, 9, 5, 7, 10, 4]

    if doSpectra:
        if not c:
            c = ROOT.TCanvas(comparisonName, comparisonName)

        results.append(c)
        c.cd()
        if doSpectra == "logy":
            c.SetLogy()

        if leg:
            leg.SetY1(leg.GetY1()-0.04*(len(spectra)+1))
        else:
            leg = ROOT.TLegend(0.55, 0.87-0.04*(len(spectra)+1), 0.85, 0.87)
            leg.SetName("{0}_legend".format(c.GetName()))
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextFont(43)
            leg.SetTextSize(16)

        results.append(leg)

        if "hist" in opt:
            baseline.SetLineColor(colors[0])
            baseline.SetLineWidth(2)
            baseline.SetLineStyle(lines[0])
            leg.AddEntry(baseline, baseline.GetTitle(), "l")
        else: 
            baseline.SetMarkerColor(colors[0])
            baseline.SetLineColor(colors[0])
            baseline.SetMarkerStyle(markers[0])
            baseline.SetMarkerSize(1.2)
            leg.AddEntry(baseline, baseline.GetTitle(), "pe")

        baseline.Draw(opt)
        if "frac" in baseline.GetYaxis().GetTitle():
            c.SetLeftMargin(0.12)
            baseline.GetYaxis().SetTitleOffset(1.4)
        if not "same" in opt:
            opt += "same"

        for obj in c.GetListOfPrimitives():
            if isinstance(obj, ROOT.TH1):
                mainHist = obj
                mainHist.SetMinimum(-1111)
                mainHist.SetMaximum(-1111)
                print("Main histogram is: {0}".format(mainHist.GetName()))
                break
        minY = min(mainHist.GetMinimum(0), baseline.GetMinimum(0))
        maxY = max(mainHist.GetMaximum(), baseline.GetMaximum())

    if doRatio:
        cname = "{0}_Ratio".format(comparisonName)
        if not cRatio:
            cRatio = ROOT.TCanvas(cname, cname)
            cRatio.SetGridx()
            cRatio.SetGridy()
        results.append(cRatio)
        cRatio.cd()
        if doRatio == "logy":
            cRatio.SetLogy()

        if legRatio:
            legRatio.SetY1(legRatio.GetY1()-0.04*(len(spectra)))
        else:
            legRatio = ROOT.TLegend(0.55, 0.87-0.04*(len(spectra)+1), 0.85, 0.87)
            legRatio.SetName("{0}_legend".format(cRatio.GetName()))
            legRatio.SetFillStyle(0)
            legRatio.SetBorderSize(0)
            legRatio.SetTextFont(43)
            legRatio.SetTextSize(16)

        results.append(legRatio)
    first = True
    for color, marker, line, h in zip(colors[1:], markers[1:], lines[1:], spectra):
        if doSpectra:
            c.cd()
            if minY > h.GetMinimum(0):
                minY = h.GetMinimum(0)
            if maxY < h.GetMaximum():
                maxY = h.GetMaximum()
            h.Draw(opt)
            if "hist" in opt:
                h.SetLineColor(color)
                h.SetLineWidth(3)
                h.SetLineStyle(line)
                leg.AddEntry(h, h.GetTitle(), "l")
            else:
                h.SetMarkerColor(color)
                h.SetLineColor(color)
                h.SetMarkerStyle(marker)
                h.SetMarkerSize(1.2)
                leg.AddEntry(h, h.GetTitle(), "pe")

        if doRatio:
            cRatio.cd()
            #hRatio = h.Clone("{0}_Ratio".format(h.GetName()))
            #hRatio.Divide(baseline)
            hRatio = ROOT.TGraphAsymmErrors(h, baseline, "cl=0.683 b(1,1) mode")
            hRatio.SetName("{0}_Ratio".format(h.GetName()))
            hRatio.SetTitle("{0} Ratio".format(h.GetTitle()))
            if first: hRatio.Draw(optRatio+"A")
            else: hRatio.Draw(optRatio)
            first = False
            print(hRatio.GetName())
            for ibin in range(0,hRatio.GetN(),5):
                if hRatio.GetX()[ibin] < 10:
                    continue
                if hRatio.GetX()[ibin] > 30:
                    break
                print("\\item JE patch = ${0:.2f}$ GeV $\\rightarrow$ ratio = ${1:.3f} + {2:.4f} - {3:.4f}$".format(hRatio.GetX()[ibin], hRatio.GetY()[ibin], hRatio.GetErrorYhigh(ibin), hRatio.GetErrorYlow(ibin)))
            if not baselineRatio:
                baselineRatio = hRatio
            if "hist" in optRatio:
                hRatio.SetLineColor(color)
                hRatio.SetLineWidth(3)
                hRatio.SetLineStyle(line)
                legRatio.AddEntry(hRatio, h.GetTitle(), "l")
            else:
                hRatio.SetMarkerColor(color)
                hRatio.SetLineColor(color)
                hRatio.SetMarkerStyle(marker)
                hRatio.SetMarkerSize(1.2)
                legRatio.AddEntry(hRatio, h.GetTitle(), "pe")
            results.append(hRatio)
            hRatio.GetYaxis().SetTitle(yaxisRatio)
            if not mainRatioHist:
                mainRatioHist = hRatio

            if not "same" in optRatio:
                optRatio += "same"

    if doRatio:
        mainRatioHist.GetXaxis().SetRangeUser(0,30)
        mainRatioHist.GetYaxis().SetRangeUser(0,1.8)
        cRatio.cd()
        legRatio.Draw()

    if doSpectra:
        if doSpectra == "logy":
            maxY *= 10
            minY /= 5
        else:
            maxY *= 1.5
            if minY < 0.2:
                minY = 0
            else:
                minY *= 0.6
        mainHist.SetMinimum(minY)
        mainHist.SetMaximum(maxY)
        c.cd()
        leg.Draw()

    return results

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Study of the effect of L0 threshold on the L1 triggers.')
    parser.add_argument('filename', metavar='FILENAME',
                        help='ROOT file name')
    args = parser.parse_args()
    
    main(args.filename)
    
    IPython.embed()
