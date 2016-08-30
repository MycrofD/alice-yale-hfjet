#!/usr/bin/env python
#python program to perform a D meson jet analysis

import ROOT
import math
import DMesonJetProjectors
from DMesonJetBase import *
import array
import copy
import DMesonJetUtils
import collections

globalList = []

class DMesonJetAnalysisEngine:
    def __init__(self, figTitle, collision, trigger, dmeson, binSet, nMassBins, minMass, maxMass, jets, spectra, projector):
        self.fFigureTitle = figTitle
        self.fCollision = collision
        self.fTrigger = trigger
        self.fDMeson = dmeson
        self.fJetDefinitions = jets
        self.fProjector = projector
        self.fBinMultiSet = binSet
        self.fNMassBins = nMassBins
        self.fMinMass = minMass
        self.fMaxMass = maxMass
        self.fSpectra = collections.OrderedDict()
        self.fCanvases = []
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            name = "{0}_{1}".format(self.fDMeson, s["name"])
            self.fSpectra[s["name"]] = Spectrum(s, name, self.fBinMultiSet.fBinSets)
            
    def CompareSpectra(self):
        spectraToCompare = []
        axisBaseline = None
        for name,s in self.fSpectra.iteritems():
            if not s.fNormHistogram:
                continue
            if len(s.fAxis) != 1:
                continue
            if axisBaseline:
                if axisBaseline.fName != s.fAxis[0].fName:
                    continue
            else:
                axisBaseline = s.fAxis[0]
            h = s.fNormHistogram.Clone("{0}_copy".format(s.fNormHistogram.GetName()))
            if s.fTitle:
                h.SetTitle(s.fTitle)
            globalList.append(h)
            spectraToCompare.append(h)
        if len(spectraToCompare) < 2:
            return
        results = DMesonJetUtils.CompareSpectra(spectraToCompare[0], spectraToCompare[1:], "{0}_SpectraComparison".format(self.fDMeson), "", "hist")
        for obj in results:
            if obj and isinstance(obj, ROOT.TCanvas):
                self.fCanvases.append(obj)
            globalList.append(obj)        

    def SaveRootFile(self, file):
        file.cd()
        for rlist in self.fBinMultiSet.GenerateInvMassRootLists():
            rlist.Write("{0}_{1}".format(self.fDMeson, rlist.GetName()), ROOT.TObject.kSingleKey)
        for s in self.fSpectra.itervalues():
            if s.fNormHistogram:
                s.fNormHistogram.Write()
            slist = s.GenerateRootList()
            slist.Write(slist.GetName(), ROOT.TObject.kSingleKey)

    def SavePlots(self, path, format):
        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(path, c.GetName(), format))

    def CreateMassFitter(self, name):
        if "D0" in self.fDMeson:
            minMass = self.fMinMass
            maxMass = self.fMaxMass
            minFitRange = self.fMinMass
            maxFitRange = self.fMaxMass
            startingSigma = 0.01
            startingSigmaBkg = 0.01
            massFitTypeSig = ROOT.MassFitter.kGaus
            massFitTypeBkg = ROOT.MassFitter.kExpo
        elif self.fDMeson == "DStar":
            print("Not implemented for DStar!")

        fitter = ROOT.MassFitter(name, massFitTypeSig, massFitTypeBkg, minMass, maxMass)
        fitter.GetFitFunction().SetParameter(1, startingSigmaBkg)
        fitter.GetFitFunction().SetParameter(4, startingSigma)
        fitter.SetFitRange(minFitRange, maxFitRange)
        if "SignalOnly" in self.fDMeson:
            fitter.DisableBkg()
        elif "BackgroundOnly" in self.fDMeson:
            fitter.DisableSig()
        globalList.append(fitter)

        return fitter

    def DoProjections(self):
        self.fProjector.GetInvMassHisograms(self.fTrigger, self.fDMeson, self.fJetDefinitions, 
                                            self.fBinMultiSet, self.fNMassBins, self.fMinMass, self.fMaxMass)

        self.fEvents = self.fProjector.fTotalEvents

    def Start(self, engines):
        self.fEngines = engines
        if not "MCTruth" in self.fDMeson:
            self.FitInvMassPlots()
            self.PlotInvMassPlots()
        
        if not "BackgroundOnly" in self.fDMeson:
            self.GenerateSpectra()
            self.PlotSpectra()

    def PlotSpectra(self):
        for s in self.fSpectra.itervalues():
            if len(s.fAxis) == 1:
                self.PlotSpectrum1D(s)
            elif len(s.fAxis) == 2:
                self.PlotSpectrum2D(s)
            elif len(s.fAxis) == 3:
                self.PlotSpectrum3D(s)
            else:
                print("Not able to plot spectra with dim > 3!")

    def PlotSpectrum1D(self, s):
        # Spectrum
        c = ROOT.TCanvas("{0}_canvas".format(s.fNormHistogram.GetName()), s.fNormHistogram.GetName())
        c.SetLogy()
        self.fCanvases.append(c)
        c.cd()
        h = s.fNormHistogram.DrawCopy()

        h.SetMarkerColor(ROOT.kBlue+2)
        h.SetMarkerStyle(ROOT.kFullCircle)
        h.SetMarkerSize(0.9)
        h.SetLineColor(ROOT.kBlue+2)

        pave = ROOT.TPaveText(0.60, 0.88, 0.9, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fFigureTitle)
        pave.AddText(self.fCollision)
        pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
        pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
        pave.AddText(s.fTitle)
        pave.Draw()

        globalList.append(c)
        globalList.append(h)
        globalList.append(pave)
        
        # Uncertainty
        c = ROOT.TCanvas("{0}_canvas".format(s.fUncertainty.GetName()), s.fUncertainty.GetName())
        self.fCanvases.append(c)
        c.cd()
        h = s.fUncertainty.DrawCopy("hist")

        h.SetLineColor(ROOT.kBlue+2)
        h.SetLineWidth(2)
        h.SetFillColorAlpha(ROOT.kBlue+2, 0.25)
        h.SetFillStyle(1001)

        pave = ROOT.TPaveText(0.10, 0.88, 0.4, 0.55, "NB NDC")
        pave.SetFillStyle(0)
        pave.SetBorderSize(0)
        pave.SetTextFont(43)
        pave.SetTextSize(15)
        pave.AddText(self.fFigureTitle)
        pave.AddText(self.fCollision)
        pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4")
        pave.AddText("with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
        pave.AddText(s.fTitle)
        pave.Draw()

        globalList.append(c)
        globalList.append(h)
        globalList.append(pave)
        
        # Mass
        if s.fMass:
            c = ROOT.TCanvas("{0}_canvas".format(s.fMass.GetName()), s.fMass.GetName())
            self.fCanvases.append(c)
            c.cd()
            h = s.fMass.DrawCopy()
    
            h.SetMarkerColor(ROOT.kBlue+2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue+2)
            h.GetYaxis().SetRangeUser(1.84, 1.91)
    
            pave = ROOT.TPaveText(0.10, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText("{0} {1}".format(self.fFigureTitle, self.fCollision))
            pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
            pave.AddText(s.fTitle)
            pave.Draw()
            
            line = ROOT.TLine(h.GetXaxis().GetBinLowEdge(1), 1.86484, h.GetXaxis().GetBinUpEdge(h.GetXaxis().GetNbins()), 1.86484)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
    
            globalList.append(c)
            globalList.append(h)
            globalList.append(pave)
            globalList.append(line)
        
        # Mass width
        if s.fMassWidth:
            c = ROOT.TCanvas("{0}_canvas".format(s.fMassWidth.GetName()), s.fMassWidth.GetName())
            self.fCanvases.append(c)
            c.cd()
            h = s.fMassWidth.DrawCopy()
    
            h.SetMarkerColor(ROOT.kBlue+2)
            h.SetMarkerStyle(ROOT.kFullCircle)
            h.SetMarkerSize(0.9)
            h.SetLineColor(ROOT.kBlue+2)
            #h.GetYaxis().SetRangeUser(1.84, 1.91)
    
            pave = ROOT.TPaveText(0.10, 0.88, 0.8, 0.68, "NB NDC")
            pave.SetFillStyle(0)
            pave.SetBorderSize(0)
            pave.SetTextFont(43)
            pave.SetTextSize(15)
            pave.AddText("{0} {1}".format(self.fFigureTitle, self.fCollision))
            pave.AddText("Charged Jets, Anti-#it{k}_{T}, #it{R} = 0.4 with D^{0} #rightarrow K^{-}#pi^{+} and c.c.")
            pave.AddText(s.fTitle)
            pave.Draw()
    
            globalList.append(c)
            globalList.append(h)
            globalList.append(pave)
            globalList.append(line)

    def PlotSpectrum2D(self, s):
        c = ROOT.TCanvas("{0}_canvas".format(s.fNormHistogram.GetName()), s.fNormHistogram.GetName())
        c.SetRightMargin(0.18)
        c.SetLogz()
        self.fCanvases.append(c)
        c.cd()
        h = s.fNormHistogram.DrawCopy("colz")
        h.GetZaxis().SetTitleOffset(1.4)

        globalList.append(c)
        globalList.append(h)       

    def PlotSpectrum3D(self, s):
        print("PlotSpectrum3D not implemented!")

    def GenerateSpectra(self):
        for s in self.fSpectra.itervalues():
            if len(s.fAxis) == 1:
                self.GenerateSpectrum1D(s)
            elif len(s.fAxis) == 2:
                self.GenerateSpectrum2D(s)
            elif len(s.fAxis) == 3:
                self.GenerateSpectrum3D(s)
            else:
                print("Not able to generate spectra with dim > 3!")
            s.GenerateNormalizedSpectrum(self.fEvents)

    def BuildSpectrum1D(self, s, name, yaxis):
        hist = ROOT.TH1D(name, name, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins))
        hist.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(yaxis)
        hist.Sumw2()
        return hist

    def GenerateSpectrum1DInvMassFit(self, s):
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fMass = self.BuildSpectrum1D(s, "{0}_Mass".format(s.fName), "D^{0} mass (GeV/#it{c}^{2})")
        s.fMassWidth = self.BuildSpectrum1D(s, "{0}_MassWidth".format(s.fName), "D^{0} mass width (GeV/#it{c}^{2})")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{m} - <#it{m}>| < 2#sigma")
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                if "SignalOnly" in self.fDMeson:
                    signal = bin.fInvMassHisto.Integral()
                    signal_unc = math.sqrt(bin.fInvMassHisto.Integral())
                else:
                    signal = bin.fMassFitter.GetSignal()
                    signal_unc = bin.fMassFitter.GetSignalError()

                s.fHistogram.SetBinContent(xbin, signal)
                s.fHistogram.SetBinError(xbin, signal_unc)
                s.fBackground.SetBinContent(xbin, bin.fMassFitter.GetBackground(2))
                s.fBackground.SetBinError(xbin, bin.fMassFitter.GetBackgroundError(2))
                s.fUncertainty.SetBinContent(xbin, signal_unc/signal) 
                s.fMass.SetBinContent(xbin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, bin.fMassFitter.GetSignalMeanError())
                s.fMassWidth.SetBinContent(xbin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum1DLikeSignMethod(self, s):
        if "SignalOnly" in self.fDMeson:
            return self.GenerateSpectrum1DInvMassFit(s)
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fMass = self.BuildSpectrum1D(s, "{0}_Mass".format(s.fName), "D^{0} mass (GeV/#it{c}^{2})")
        s.fMassWidth = self.BuildSpectrum1D(s, "{0}_MassWidth".format(s.fName), "D^{0} mass width (GeV/#it{c}^{2})")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{m} - <#it{m}>| < 2#sigma")

        s.fSignalPlusBackground = self.BuildSpectrum1D(s, "{0}_SignalPlusBackground".format(s.fName), "counts")
        s.fLikeSignBackground = self.BuildSpectrum1D(s, "{0}_LikeSignBackground".format(s.fName), "counts")
        s.fSignalMinusLikeSign = self.BuildSpectrum1D(s, "{0}_SignalMinusLikeSign".format(s.fName), "counts")
        eng = None
        for engTest in self.fEngines:
            if engTest.fTrigger == self.fTrigger and engTest.fDMeson == s.fLikeSignTree:
                eng = engTest
                break
        if not eng:
            print("Could not find engine with trigger '{0}' and tree name '{1}'".format(self.fTrigger, s.fLikeSignTree))
            return
        binSetName = s.fBins[0]
        s.fLikeSignSubtractedBinSet = copy.deepcopy(self.fBinMultiSet.fBinSets[binSetName])
        s.fLikeSignSubtractedBinSet.fName = "{0}_LikeSignSubtracted".format(s.fLikeSignSubtractedBinSet.fName)
        for i,bin in enumerate(s.fLikeSignSubtractedBinSet.fBins):
            binSig_1 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
            binSig_2 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())

            sigBkgErr = ROOT.Double(0)
            sigBkg = bin.fInvMassHisto.IntegralAndError(binSig_1, binSig_2, sigBkgErr)
            s.fSignalPlusBackground.SetBinContent(i+1, sigBkg)
            s.fSignalPlusBackground.SetBinError(i+1, sigBkgErr)

            lsErr = ROOT.Double(0)
            ls = eng.fBinMultiSet.fBinSets[binSetName].fBins[i].fInvMassHisto.IntegralAndError(binSig_1, binSig_2, lsErr)
            s.fLikeSignBackground.SetBinContent(i+1, ls)
            s.fLikeSignBackground.SetBinError(i+1, lsErr)

            sig = sigBkg - ls
            sigErr = math.sqrt(sigBkgErr**2 + lsErr**2)
            s.fSignalMinusLikeSign.SetBinContent(i+1, sig)
            s.fSignalMinusLikeSign.SetBinError(i+1, sigErr)

            bin.fInvMassHisto.Add(eng.fBinMultiSet.fBinSets[binSetName].fBins[i].fInvMassHisto, -1)
            bin.fMassFitter = None

        self.FitInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins, s.fLikeSignSubtractedBinSet.fFitOptions)
        self.PlotInvMassPlotsBinSet(s.fLikeSignSubtractedBinSet.fName, s.fLikeSignSubtractedBinSet.fBins)

        for i,bin in enumerate(s.fLikeSignSubtractedBinSet.fBins):
            if bin.fMassFitter is None:
                print("The bin printed below does not have a mass fitter!")
                bin.Print()
                continue
            xbin = i+1
            signal = bin.fMassFitter.GetSignal()
            signal_unc = bin.fMassFitter.GetSignalError()

            s.fHistogram.SetBinContent(xbin, signal)
            s.fHistogram.SetBinError(xbin, signal_unc)
            s.fBackground.SetBinContent(xbin, bin.fMassFitter.GetBackground(2))
            s.fBackground.SetBinError(xbin, bin.fMassFitter.GetBackgroundError(2))
            s.fUncertainty.SetBinContent(xbin, signal_unc/signal) 
            s.fMass.SetBinContent(xbin, bin.fMassFitter.GetSignalMean())
            s.fMass.SetBinError(xbin, bin.fMassFitter.GetSignalMeanError())
            s.fMassWidth.SetBinContent(xbin, bin.fMassFitter.GetSignalWidth())
            s.fMassWidth.SetBinError(xbin, bin.fMassFitter.GetSignalWidthError())

    def GenerateSpectrum1DSideBandMethod(self, s):
        s.fSideBandHistograms = []
        s.fSignalHistograms = []
        s.fMass = None
        s.fMassWidth = None
        s.fHistogram = self.BuildSpectrum1D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum1D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fBackground = self.BuildSpectrum1D(s, "{0}_Bkg".format(s.fName), "background |#it{{m}} - <#it{{m}}>| < {0}#sigma".format(int(s.fBinCountSignalSigmas)))
        s.fSideBandWindowTotalHistogram = self.BuildSpectrum1D(s, "{0}_SideBandWindowTotal".format(s.fName), "counts")
        s.fSignalWindowTotalHistogram = self.BuildSpectrum1D(s, "{0}_SignalWindowTotal".format(s.fName), "counts")
        print("Generating spectrum {0}".format(s.fName))
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if "SignalOnly" in self.fDMeson:
                    sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), 0, -1, "e")
                    sig.SetTitle(bin.GetTitle())
                    s.fSignalHistograms.append(sig)
                    s.fSignalWindowTotalHistogram.Add(sig)
                else:
                    binSBL_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBL_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
                    binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
                    binSig_1 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
                    binSig_2 = bin.fBinCountAnalysisHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + s.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())

                    if binSBL_1 < 1:
                        binSBL_1 = 1
                    if binSBR_2 > bin.fBinCountAnalysisHisto.GetXaxis().GetNbins():
                        binSBR_2 = bin.fBinCountAnalysisHisto.GetXaxis().GetNbins()
                    sbL = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowL_{1}".format(s.fName, bin.GetName()), binSBL_1, binSBL_2, "e")
                    sbR = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SideBandWindowR_{1}".format(s.fName, bin.GetName()), binSBR_1, binSBR_2, "e")
                    sig = bin.fBinCountAnalysisHisto.ProjectionY("{0}_SignalWindow_{1}".format(s.fName, bin.GetName()), binSig_1, binSig_2, "e")
                    sbTotal = sbL.Clone("{0}_SideBandWindow_{1}".format(s.fName, bin.GetName()))
                    sbTotal.Add(sbR)

                    bkgErrSigWindow = ROOT.Double(0.)
                    bkgSigWindow = bin.fMassFitter.GetBackgroundAndError(bkgErrSigWindow, s.fBinCountSignalSigmas)
                    print("Bin: {0}".format(bin.GetTitle()))
                    print("The background in side bands is: {0} + {1} = {2}".format(sbL.Integral(0,-1), sbR.Integral(0,-1), sbTotal.Integral(0,-1)))
                    print("The estimated background in the signal window is {0} +/- {1}".format(bkgSigWindow, bkgErrSigWindow))
                    print("The total signal+background is {0}, which is the same from the invariant mass plot {1} or summing signal and background {2}".format(sig.Integral(0,-1), bin.fInvMassHisto.Integral(binSig_1, binSig_2), bin.fMassFitter.GetSignal()+bin.fMassFitter.GetBackground(s.fBinCountSignalSigmas)))

                    sbTotalIntegral = sbTotal.Integral(0, -1)
                    if sbTotalIntegral > 0:
                        sbTotal.Scale(1. / sbTotalIntegral)
    
                        for xbin in range(0, sbTotal.GetNbinsX()+2):
                            if sbTotal.GetBinContent(xbin) == 0:
                                continue
                            error = math.sqrt(sbTotal.GetBinError(xbin)**2/sbTotal.GetBinContent(xbin)**2+bkgErrSigWindow**2/bkgSigWindow**2)*sbTotal.GetBinContent(xbin)*bkgSigWindow
                            cont = sbTotal.GetBinContent(xbin)*bkgSigWindow
                            sbTotal.SetBinError(xbin, error)
                            sbTotal.SetBinContent(xbin, cont)

                    integralError = ROOT.Double(0.)
                    integral = sbTotal.IntegralAndError(0, -1, integralError)
                    print("The total normalized side-band background is {0} +/- {1}".format(integral, integralError))

                    sbTotal.SetTitle(bin.GetTitle())
                    s.fSideBandHistograms.append(sbTotal)
                    s.fSideBandWindowTotalHistogram.Add(sbTotal)

                    sig.SetTitle(bin.GetTitle())
                    s.fSignalHistograms.append(sig)
                    s.fSignalWindowTotalHistogram.Add(sig)

        s.fHistogram.Add(s.fSignalWindowTotalHistogram)
        s.fHistogram.Add(s.fSideBandWindowTotalHistogram, -1)
        s.fBackground.Add(s.fSideBandWindowTotalHistogram)

        for xbin in range(0, s.fHistogram.GetNbinsX()+2):
            if s.fHistogram.GetBinContent(xbin) > 0:
                s.fUncertainty.SetBinContent(xbin, s.fHistogram.GetBinError(xbin) / s.fHistogram.GetBinContent(xbin))
            else:
                s.fUncertainty.SetBinContent(xbin, 0)

    def GenerateSpectrum1D(self, s):
        if s.fAnalysisType == AnalysisType.SideBand:
            self.GenerateSpectrum1DSideBandMethod(s)
        elif s.fAnalysisType == AnalysisType.LikeSign:
            self.GenerateSpectrum1DLikeSignMethod(s)
        else:
            self.GenerateSpectrum1DInvMassFit(s)

    def BuildSpectrum2D(self, s, name, zaxis):
        hist = ROOT.TH2D(name, name, len(s.fAxis[0].fBins)-1, array.array('d',s.fAxis[0].fBins), len(s.fAxis[1].fBins)-1, array.array('d',s.fAxis[1].fBins))
        hist.GetXaxis().SetTitle(s.fAxis[0].GetTitle())
        hist.GetYaxis().SetTitle(s.fAxis[1].GetTitle())
        hist.GetZaxis().SetTitle(zaxis)
        hist.Sumw2()
        return hist

    def GenerateSpectrum2D(self, s):
        s.fHistogram = self.BuildSpectrum2D(s, s.fName, "counts")
        s.fUncertainty = self.BuildSpectrum2D(s, "{0}_Unc".format(s.fName), "relative statistical uncertainty")
        s.fMass = self.BuildSpectrum2D(s, "{0}_Mass".format(s.fName), "mass")
        s.fBackground = self.BuildSpectrum2D(s, "{0}_Bkg".format(s.fName), "background |#it{m} - <#it{m}>| < 3#sigma")
        for binSetName in s.fBins:
            for bin in self.fBinMultiSet.fBinSets[binSetName].fBins:
                if bin.fMassFitter is None:
                    print("The bin printed below does not have a mass fitter!")
                    bin.Print()
                    continue
                xbin = s.fHistogram.GetXaxis().FindBin(bin.GetBinCenter(s.fAxis[0].fName))
                ybin = s.fHistogram.GetYaxis().FindBin(bin.GetBinCenter(s.fAxis[1].fName))
                if "SignalOnly" in self.fDMeson:
                    signal = bin.fInvMassHisto.Integral()
                    signal_unc = math.sqrt(bin.fInvMassHisto.Integral())
                else:
                    signal = bin.fMassFitter.GetSignal()
                    signal_unc = bin.fMassFitter.GetSignalError()

                s.fHistogram.SetBinContent(xbin, ybin, signal)
                s.fHistogram.SetBinError(xbin, ybin, signal_unc)
                s.fBackground.SetBinContent(xbin, ybin, bin.fMassFitter.GetBackground())
                s.fBackground.SetBinError(xbin, ybin, bin.fMassFitter.GetBackgroundError())
                s.fUncertainty.SetBinContent(xbin, ybin, signal_unc/signal)
                s.fMass.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalMean())
                s.fMass.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalMeanError())
                s.fMassWidth.SetBinContent(xbin, ybin, bin.fMassFitter.GetSignalWidth())
                s.fMassWidth.SetBinError(xbin, ybin, bin.fMassFitter.GetSignalWidthError())
                
    def GenerateSpectrum3D(self, s):
        print("GenerateSpectrum3D not implemented!")
                
    def DrawFitResults(self,bin):
        if bin.fMassFitter is None:
            return
        
        bin.fMassFitter.Draw("same");
        
        fitStatus = int(bin.fMassFitter.GetFitStatus())
        if fitStatus == 0:
            chi2Text = bin.fMassFitter.GetChisquareString().Data()
        else:
            chi2Text = "Fit failed"

        paveSig = ROOT.TPaveText(0.165, 0.795, 0.490, 0.92, "NB NDC")
        globalList.append(paveSig)
        paveSig.SetBorderSize(0)
        paveSig.SetFillStyle(0)
        paveSig.SetTextFont(43)
        paveSig.SetTextSize(14)
        paveSig.SetTextAlign(13)
        paveSig.AddText("{0}, {1}".format(bin.fMassFitter.GetSignalString().Data(), 
                                          bin.fMassFitter.GetBackgroundString().Data()))
        paveSig.AddText("{0}, {1}, {2}".format(bin.fMassFitter.GetSignalOverBackgroundString().Data(), 
                                          bin.fMassFitter.GetSignalOverSqrtSignalBackgroundString().Data(),
                                          chi2Text))
        paveSig.Draw()

        paveFit = ROOT.TPaveText(0.48, 0.51, 0.97, 0.77, "NB NDC")
        globalList.append(paveFit)
        paveFit.SetBorderSize(0)
        paveFit.SetFillStyle(0)
        paveFit.SetTextFont(43)
        paveFit.SetTextSize(14)
        paveFit.SetTextAlign(23)

        paveFit.AddText(bin.fMassFitter.GetSignalMeanString().Data())
        paveFit.AddText(bin.fMassFitter.GetSignalWidthString().Data())
        paveFit.AddText(bin.fMassFitter.GetBkgPar1String().Data())
        paveFit.AddText(bin.fMassFitter.GetTotalEntriesString().Data())
        paveFit.Draw()

    def FitInvMassPlots(self):
        for binSet in self.fBinMultiSet.fBinSets.itervalues():
            self.FitInvMassPlotsBinSet(binSet.fName,binSet.fBins,binSet.fFitOptions)
            
    def FitInvMassPlotsBinSet(self,name,bins,fitOptions):
        pdgMass = ROOT.TDatabasePDG.Instance().GetParticle(421).Mass()
        
        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            fitter = self.CreateMassFitter("{0}_{1}_{2}_fitter".format(self.fDMeson, name, bin.GetName()))
            bin.SetMassFitter(fitter)
            integral = bin.fInvMassHisto.Integral(1, bin.fInvMassHisto.GetXaxis().GetNbins())
            fitter.GetFitFunction().FixParameter(0, integral) # total integral is fixed
            fitter.GetFitFunction().SetParameter(2, integral / 100) # signal integral (start with very small signal)
            fitter.GetFitFunction().SetParLimits(2, 0, integral) # signal integral has to be contained in the total integral
            fitter.GetFitFunction().SetParameter(3, pdgMass) # start fitting using PDG mass
            print("Fitting bin {0}".format(bin.GetTitle()))
            
            fitter.Fit(bin.fInvMassHisto, fitOptions);
            
    def PlotInvMassPlots(self):
        for binSet in self.fBinMultiSet.fBinSets.itervalues():
            self.PlotInvMassPlotsBinSet(binSet.fName,binSet.fBins)

    def PlotInvMassLikeSign(self, name, bin):
        spectrum = None
        for s in self.fSpectra.itervalues():
            for binSetName in s.fBins:
                if binSetName == name and s.fAnalysisType == AnalysisType.LikeSign:
                    spectrum = s
                    break
        if not spectrum:
            return None

        eng = None
        for engTest in self.fEngines:
            if engTest.fTrigger == self.fTrigger and engTest.fDMeson == s.fLikeSignTree:
                eng = engTest
                break
        if not eng:
            return None

        binLS = None
        for binTest in eng.fBinMultiSet.fBinSets[name].fBins:
            if binTest.IsSameOf(bin):
                binLS = binTest
                break

        hls = binLS.fInvMassHisto.DrawCopy("hist")
        hls.SetLineColorAlpha(ROOT.kCyan+2, 0.4)
        hls.SetFillColorAlpha(ROOT.kCyan+2, 0.4)
        hls.SetFillStyle(1001)
        globalList.append(hls)

        return hls

    def PlotInvMassSideBands(self, name, bin):
        spectrum = None
        for s in self.fSpectra.itervalues():
            for binSetName in s.fBins:
                if binSetName == name and s.fAnalysisType == AnalysisType.SideBand:
                    spectrum = s
                    break
        if not spectrum:
            return None

        binSBL_1 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - spectrum.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
        binSBL_2 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - spectrum.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
        binSBR_1 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + spectrum.fSideBandMinSigmas*bin.fMassFitter.GetSignalWidth())
        binSBR_2 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + spectrum.fSideBandMaxSigmas*bin.fMassFitter.GetSignalWidth())
        binSig_1 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() - spectrum.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())
        binSig_2 = bin.fInvMassHisto.GetXaxis().FindBin(bin.fMassFitter.GetSignalMean() + spectrum.fBinCountSignalSigmas*bin.fMassFitter.GetSignalWidth())

        bin.fSideBandWindowHisto = bin.fInvMassHisto.Clone(bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSBWindow"))
        bin.fSignalWindowHisto = bin.fInvMassHisto.Clone(bin.fInvMassHisto.GetName().replace("InvMass", "InvMassSigWindow"))

        for xbin in range(0,bin.fInvMassHisto.GetNbinsX()+2):
            if not ((xbin >= binSBL_1 and xbin < binSBL_2) or (xbin >= binSBR_1 and xbin < binSBR_2)):
                bin.fSideBandWindowHisto.SetBinContent(xbin, 0)
                bin.fSideBandWindowHisto.SetBinError(xbin, 0)
            if not (xbin >= binSig_1 and xbin < binSig_2):
                bin.fSignalWindowHisto.SetBinContent(xbin, 0)
                bin.fSignalWindowHisto.SetBinError(xbin, 0)
        hsb = bin.fSideBandWindowHisto.DrawCopy("hist")
        hsb.SetFillColorAlpha(ROOT.kGreen+2, 0.4)
        hsb.SetFillStyle(1001)
        hsb.SetLineColorAlpha(ROOT.kGreen+2, 0.4)
        hsig = bin.fSignalWindowHisto.DrawCopy("hist same")
        hsig.SetFillColorAlpha(ROOT.kRed+2, 0.4)
        hsig.SetFillStyle(1001)
        hsig.SetLineColorAlpha(ROOT.kRed+2, 0.4)
        globalList.append(hsb)
        globalList.append(hsig)

        return hsb

    def PlotInvMassPlotsBinSet(self, name, bins):
        cname = "{0}_{1}".format(self.fDMeson, name)
        c = DMesonJetUtils.GenerateMultiCanvas(cname, len(bins))
        self.fCanvases.append(c)
        globalList.append(c)
        for i,bin in enumerate(bins):
            if not bin.fInvMassHisto:
                continue
            pad = c.cd(i+1)
            pad.SetLeftMargin(0.12)
            pad.SetRightMargin(0.05)
            pad.SetTopMargin(0.08)
            pad.SetBottomMargin(0.13)
            if not "SignalOnly" in self.fDMeson:
                SB = self.PlotInvMassSideBands(name, bin)
                LS = self.PlotInvMassLikeSign(name, bin)
            else:
                SB = None
                LS = None
            if SB:
                h = bin.fInvMassHisto.DrawCopy("same")
                globalList.append(h)
                SB.SetMaximum(h.GetMaximum()*1.8)
                h = SB
            elif LS:
                h = bin.fInvMassHisto.DrawCopy("same")
                globalList.append(h)
                LS.SetMaximum(h.GetMaximum()*1.8)
                h = LS
            else:
                h = bin.fInvMassHisto.DrawCopy()
                globalList.append(h)
                h.SetMaximum(h.GetMaximum()*1.8)
            h.GetXaxis().SetTitleFont(43)
            h.GetXaxis().SetTitleOffset(2.3)
            h.GetXaxis().SetTitleSize(19)
            h.GetXaxis().SetLabelFont(43)
            h.GetXaxis().SetLabelOffset(0.009)
            h.GetXaxis().SetLabelSize(18)
            h.GetYaxis().SetTitleFont(43)
            h.GetYaxis().SetTitleOffset(2.3)
            h.GetYaxis().SetTitleSize(19)
            h.GetYaxis().SetLabelFont(43)
            h.GetYaxis().SetLabelOffset(0.009)
            h.GetYaxis().SetLabelSize(18)
            htitle = ROOT.TPaveText(0.12, 0.99, 0.95, 0.93, "NB NDC")
            htitle.SetBorderSize(0)
            htitle.SetFillStyle(0)
            htitle.SetTextFont(43)
            htitle.SetTextSize(18)
            htitle.AddText(bin.GetTitle())
            htitle.Draw()
            globalList.append(htitle)
            self.DrawFitResults(bin)

class DMesonJetAnalysis:
    def __init__(self, name):
        self.fName = name
        self.fAnalysisEngine = []
        self.fCanvases = []

    def SetProjector(self, projector):
        self.fProjector = projector
        
    def StartAnalysis(self, figTitle, collision, config):
        binMultiSet = DMesonJetProjectors.BinMultiSet()
        for binLists in config["binLists"]:
            if "active" in binLists and not binLists["active"]:
                continue
            limitSetList = []
            for name, binList in binLists["bins"].iteritems():
                limitSetList.append((name, binList))
            if "cuts" in binLists:
                cuts = binLists["cuts"]
            else:
                cuts = []
            if "bin_count_analysis" in binLists:
                bin_count_analysis = binLists["bin_count_analysis"]
            else:
                bin_count_analysis = None
            if "efficiency" in binLists and binLists["efficiency"]:
                effWeight = DMesonJetProjectors.EfficiencyWeightCalculator("{0}/{1}".format(self.fProjector.fInputPath, binLists["efficiency"]["file_name"]), binLists["efficiency"]["list_name"], binLists["efficiency"]["object_name"])
                fitOptions = "0 WL S"
            else:
                effWeight = DMesonJetProjectors.SimpleWeight()
                fitOptions = "0 L S"
            binMultiSet.AddBinSet(binLists["name"], limitSetList, cuts, bin_count_analysis, effWeight, fitOptions)

        for trigger in config["trigger"]:
            for d_meson in config["d_meson"]:
                binset_copy = copy.deepcopy(binMultiSet)
                eng = DMesonJetAnalysisEngine(figTitle, collision, trigger, d_meson, 
                                              binset_copy, config["n_mass_bins"], config["min_mass"], config["max_mass"],
                                              config["jets"], config["spectra"], self.fProjector)
                self.fAnalysisEngine.append(eng)
                eng.DoProjections()

        for eng in self.fAnalysisEngine:
            if not "LikeSign" in eng.fDMeson:
                eng.Start(self.fAnalysisEngine)
                eng.CompareSpectra()

        self.CompareSpectra(config["spectra"])

    def CompareSpectra(self, spectra):
        if len(self.fAnalysisEngine) <= 1:
            return
        for s in spectra:
            if "active" in s and not s["active"]:
                continue
            if "axis" in s and len(s["axis"]) != 1:
                continue
            baseline = self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(self.fAnalysisEngine[0].fSpectra[s["name"]].fHistogram.GetName()))
            baseline.SetTitle(self.fAnalysisEngine[0].fDMeson)
            globalList.append(baseline)
            spectraToCompare = []
            for eng in self.fAnalysisEngine[1:]:
                if not eng.fSpectra[s["name"]].fHistogram:
                    continue
                h = eng.fSpectra[s["name"]].fHistogram.Clone("{0}_copy".format(eng.fSpectra[s["name"]].fHistogram.GetName()))
                h.SetTitle(eng.fDMeson)
                globalList.append(h)
                spectraToCompare.append(h)

            results = DMesonJetUtils.CompareSpectra(baseline, spectraToCompare, "{0}_SpectraComparison".format(s["name"]))
            for obj in results:
                if isinstance(obj, ROOT.TCanvas):
                    self.fCanvases.append(obj)
                globalList.append(obj)

    def SaveRootFile(self, path):
        file = ROOT.TFile.Open("{0}/{1}.root".format(path, self.fName), "recreate")
        file.cd()
        for eng in self.fAnalysisEngine:
            eng.SaveRootFile(file)
        file.Close()
        
    def SavePlots(self, path, format):
        fullPath = "{0}/{1}/{2}".format(path, self.fName, format)
        if not os.path.isdir(fullPath):
            os.makedirs(fullPath)
        for eng in self.fAnalysisEngine:
            eng.SavePlots(fullPath, format)

        for c in self.fCanvases:
            c.SaveAs("{0}/{1}.{2}".format(fullPath, c.GetName(), format))
