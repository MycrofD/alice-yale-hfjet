#!/usr/local/bin/python
#python script to do EMCal cluster trigger turn on curves

import ROOT
import argparse
import array
import math
import yaml
import IPython

globalList = []

ADCtoGeV = 0.018970588 * 4 

class TriggerTurnOnAnalysis:
    def __init__(self, name, detector, file, outputPath, triggerList, baseTrigger):
        self.fBaseName = name
        self.fFile = file
        self.fMaxObservableValue = 40
        self.fNormalize = False
        self.fDetector = detector
        self.fOutputPath = outputPath
        self.fColorList = dict()
        self.fMarkerList = dict()
        self.fTriggerList = triggerList
        self.fBaseTrigger = baseTrigger
        self.fColorList[self.fBaseTrigger] = ROOT.kBlack
        self.fMarkerList[self.fBaseTrigger] = ROOT.kFullCircle
        colors = [ROOT.kRed+2, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kOrange+2]
        markers = [ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenTriangleDown, ROOT.kOpenDiamond, ROOT.kOpenCross, ROOT.kOpenStar]
        for color,marker,trigger in zip(colors, markers, self.fTriggerList):
            if trigger == self.fBaseTrigger:
                continue
            self.fColorList[trigger] = color
            self.fMarkerList[trigger] = marker

    def GetHistograms(self):
        self.fHistograms = dict()
        self.fEvents = dict()
        for trigger in self.fTriggerList:
            self.GetHistogram(trigger)
            
    def Ratios(self):
        self.fRatios = dict()
        for trigger,hist in self.fProjections.iteritems():
            if trigger == self.fBaseTrigger:
                continue
            ratio = hist.Clone("Ratio_{0}".format(trigger))
            ratio.Divide(self.fProjections[self.fBaseTrigger])
            if self.fNormalize:
                ratio.GetYaxis().SetTitle("#frac{#it{N}_{evt, {" + self.fBaseTrigger + "}}{#it{N}_{evt, rare} #frac{d#it{N}_{rare}}{d" + self.fObservableName + "} / #frac{d#it{N}_{" + self.fBaseTrigger + "}}{d" + self.fObservableName + "}")
            else:
                ratio.GetYaxis().SetTitle("#it{N}_{trigger} / #it{N}_{" + self.fBaseTrigger + "}")
            self.fRatios[trigger] = ratio
            globalList.append(ratio)
    
    def Plot(self):
        min = 1e-6
        max = 1e3
        maxR = 100
        if not self.fNormalize:
            min *= self.fEvents[self.fBaseTrigger]
            max *= self.fEvents[self.fBaseTrigger]
            maxR = 1.2
        self.DoPlot(self.fProjections, "{0}_Spectra".format(self.fName), "{0} Spectra".format(self.fTitle), True, min, max)
        self.DoPlot(self.fRatios, "{0}_Ratios".format(self.fName), "{0} Ratios".format(self.fTitle),  False, 0, maxR)
        
    def DoPlot(self, list, name, title, logY, min, max):
        canvas = ROOT.TCanvas(name, title)
        canvas.cd()
        if logY:
            canvas.SetLogy()
        leg = ROOT.TLegend(0.6, 0.7, 0.85, 0.85)
        leg.SetTextFont(43)
        leg.SetTextSize(12)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        i = 0
        for label,hist in list.iteritems():
            hist.SetLineColor(self.fColorList[label])
            hist.SetMarkerColor(self.fColorList[label])
            hist.SetMarkerStyle(self.fMarkerList[label])
            hist.SetMarkerSize(1)
            hist.GetXaxis().SetRangeUser(0,self.fMaxObservableValue)
            leg.AddEntry(hist, label, "pe")
            if i == 0:
                hist.Draw()
                hist.GetYaxis().SetRangeUser(min,max)
            else:
                hist.Draw("same")
            i += 1
            
        leg.Draw()
        globalList.append(leg)
        globalList.append(canvas)
        canvas.SaveAs("{0}/{1}.pdf".format(self.fOutputPath, name))

class ClusterTriggerTurnOnAnalysis(TriggerTurnOnAnalysis):
    def __init__(self, name, detector, file, outputPath, triggerList, baseTrigger, clusters, cells):
        TriggerTurnOnAnalysis.__init__(self, name, detector, file, outputPath, triggerList, baseTrigger)

        self.fName = "{0}_Cluster".format(self.fBaseName)
        self.fTitle = "{0} Cluster".format(self.fBaseName)
        self.fObservableName = "#it{E}_{cluster}"
        self.fListName = "AliAnalysisTaskEmcalJetQA_{clusName}_{cellName}_{triggerName}_histos/histosAliAnalysisTaskEmcalJetQA/{clusName}".format(clusName=clusters, cellName=cells, triggerName="{0}")
        if "SM" in detector:
            self.fListName += "/BySM"
            self.fHistogramName = "fHistClusEnergy_{0}_0".format(detector)
        else:
            self.fHistogramName = "fHistClusPhiEtaEnergy_0" 
        
    def GetHistogram(self, trigger):
        listName = self.fListName.format(trigger)
        
        chain = listName.split("/")
        list = None
        scaleFactor = 1
        for name in chain:
            if list:
                list = list.FindObject(name)
            else:
                list = self.fFile.Get(name)
                if list:
                    hEventsName = "fHistEventCount"
                    hevents = list.FindObject(hEventsName)
                    if not hevents:
                        print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hEventsName, name))
                        list.Print()
                        exit(1) 
    
                    self.fEvents[trigger] = hevents.GetBinContent(1)
                    print("Total number of events for {0} is {1}".format(trigger, self.fEvents[trigger]))
                
        if not list:
            self.fFile.ls()
            print("Could not get list '{0}' from file '{1}'! Aborting...".format(listName, self.fFile.GetName()))
            exit(1)
                
        hist = list.FindObject(self.fHistogramName)
        if not hist:
            list.Print()
            print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(self.fHistogramName, listName))
            exit(1)
                           
        self.fHistograms[trigger] = hist.Clone("Clusters_{0}".format(trigger))
        self.fHistograms[trigger].Sumw2()
            
    def ProjectHistograms(self):
        self.fProjections = dict()
        for trigger,hist in self.fHistograms.iteritems():
            if self.fDetector == "EMCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, 0, hist.GetYaxis().FindBin(4), "e")
            elif self.fDetector == "DCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, hist.GetYaxis().FindBin(4), -1, "e")
            elif self.fDetector == "EMCal+DCal":
                self.fProjections[trigger] = hist.ProjectionZ("ClusterEnergy_{0}".format(trigger), 0, -1, 0, -1, "e")
            else: #by SM
                self.fProjections[trigger] = hist.Clone("ClusterEnergy_{0}".format(trigger))
            self.fProjections[trigger].GetXaxis().SetTitle(self.fObservableName + " (GeV)")
            if self.fNormalize:
                self.fProjections[trigger].Scale(1. / self.fEvents[trigger], "width")
                self.fProjections[trigger].GetYaxis().SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d" + self.fObservableName + "} (GeV)^{-1}")
            else:
                self.fProjections[trigger].GetYaxis().SetTitle("counts")
                
class PatchTriggerTurnOnAnalysis(TriggerTurnOnAnalysis):
    def __init__(self, name, detector, file, outputPath, triggerList, baseTrigger, triggerPatchName, patchType):
        TriggerTurnOnAnalysis.__init__(self, name, detector, file, outputPath, triggerList, baseTrigger)

        self.fName = "{0}_{1}_{2}_Patch".format(self.fBaseName, triggerPatchName, patchType)
        self.fTitle = "{0} {1} {2} Patch".format(self.fBaseName, triggerPatchName, patchType)
        self.fObservableName = "#it{E}_{patch}"
        self.fTriggerPatchName = triggerPatchName
        self.fTriggerPatchType = patchType
        self.fListName = "AliEmcalTriggerQATask_{0}_histos/histosAliEmcalTriggerQATask_{0}"
        self.fHistogramName = "EMCTRQA_hist{detectorName}PatchAmp{triggerPatch}{type}".format(detectorName="{detector}", triggerPatch=triggerPatchName, type=patchType) 
        
    def GetHistogram(self, trigger):
        listName = self.fListName.format(trigger)
        
        chain = listName.split("/")
        list = None
        for name in chain:
            if list:
                list = list.FindObject(name)
            else:
                list = self.fFile.Get(name)
                if list:
                    hEventsName = "fHistEventCount"
                    hevents = list.FindObject(hEventsName)
                    if not hevents:
                        print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hEventsName, name))
                        list.Print()
                        exit(1) 
    
                    self.fEvents[trigger] = hevents.GetBinContent(1)
                    print("Total number of events for {0} is {1}".format(trigger, self.fEvents[trigger]))
                
        if not list:
            self.fFile.ls()
            print("Could not get list '{0}' from file '{1}'! Aborting...".format(listName, self.fFile.GetName()))
            exit(1)
            
        if self.fDetector == "EMCal" or self.fDetector == "DCal":
            hname = self.fHistogramName.format(detector=self.fDetector)
            hist = list.FindObject(hname)
            if not hist:
                list.Print()
                print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hname, listName))
                exit(1)
                        
        elif self.fDetector == "EMCal+DCal":
            hname = self.fHistogramName.format(detector="EMCal")
            hist1 = list.FindObject(hname)
            if not hist1:
                list.Print()
                print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hname, listName))
                exit(1)
                
            hname = self.fHistogramName.format(detector="DCal")
            hist2 = list.FindObject(hname)
            if not hist2:
                list.Print()
                print("Could not get histogram '{0}' from list '{1}'! Aborting...".format(hname, listName))
                exit(1)
                
            hist = hist1
            hist.Add(hist2)
        
        else:
            print("Detector '{0}' unknown! Aborting...".format(self.fDetector))
            exit(1)
        
        self.fHistograms[trigger] = hist.Clone("{0}_{1}_Patches_{2}".format(self.fTriggerPatchName, self.fTriggerPatchType, trigger))
        self.fHistograms[trigger].Sumw2()
            
    def ProjectHistograms(self):
        self.fProjections = dict()
        for trigger,hist in self.fHistograms.iteritems():
            hist.GetXaxis().SetLimits(hist.GetXaxis().GetXmin() * ADCtoGeV, hist.GetXaxis().GetXmax() * ADCtoGeV)
            hist.GetXaxis().SetTitle(self.fObservableName + " (GeV)")
            self.fProjections[trigger] = hist
            globalList.append(self.fProjections[trigger])
            if self.fNormalize:
                self.fProjections[trigger].Scale(1. / self.fEvents[trigger], "width")
                self.fProjections[trigger].GetYaxis().SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d" + self.fObservableName + "} (GeV)^{-1}")
            else:
                self.fProjections[trigger].GetYaxis().SetTitle("counts")


def CompareSpectra(anaList, search, baseLine):
    canvases = dict()
    ratioCanvases = dict()
    baseLineHistos = dict()
    
    for ana in anaList:
        if ana.fName == baseLine:
            for trigger in ana.fTriggerList:
                canvases[trigger] = ROOT.TCanvas("{0}_Comparison".format(trigger), "{0} Comparison".format(trigger))
                canvases[trigger].SetLogy()
                baseLineHistos[trigger] = ana.fProjections[trigger].Clone()
                baseLineHistos[trigger].SetLineColor(ROOT.kBlack)
                baseLineHistos[trigger].SetMarkerColor(ROOT.kBlack)
                #baseLineHistos[trigger].Rebin(3)
                baseLineHistos[trigger].Scale(1. / baseLineHistos[trigger].Integral())
                baseLineHistos[trigger].GetXaxis().SetRangeUser(0,20)
                globalList.append(baseLineHistos[trigger])
                
                canvases[trigger].cd()
                baseLineHistos[trigger].Draw()
                
                ratioCanvases[trigger] = ROOT.TCanvas("{0}_Ratio".format(trigger), "{0} Ratio".format(trigger))
                
                ratio = baseLineHistos[trigger].Clone("{0}_ratio".format(baseLineHistos[trigger].GetName()))
                ratio.Divide(baseLineHistos[trigger])
                ratioCanvases[trigger].cd()
                ratio.Draw("hist")
                ratio.GetYaxis().SetRangeUser(0,3)
                globalList.append(ratio)
                
            break

    colors = [ROOT.kRed+2, ROOT.kBlue+2, ROOT.kGreen+2, ROOT.kMagenta+2, ROOT.kCyan+2, ROOT.kOrange+2, ROOT.kTeal+2, ROOT.kYellow+2, ROOT.kPink+1, ROOT.kSpring, ROOT.kViolet+7, ROOT.kAzure+1]
    iSM = 0
    for ana in anaList:
        if not search in ana.fName or ana.fName == baseLine:
            continue
        for trigger in ana.fTriggerList:
            canvases[trigger].cd()
            hist = ana.fProjections[trigger].Clone()
            hist.SetLineColor(colors[iSM])
            hist.SetMarkerColor(colors[iSM])
            #hist.Rebin(3)
            hist.Scale(1. / hist.Integral())
            hist.Draw("same")
            
            ratioCanvases[trigger].cd()
            ratio = hist.Clone("{0}_ratio".format(hist.GetName()))
            ratio.Divide(baseLineHistos[trigger])
            ratio.Draw("same")
            
            globalList.append(hist)
            globalList.append(ratio)
            
        first = False
        iSM += 1
        
    globalList.append(canvases)
    globalList.append(ratioCanvases)

def main(yamlConfig, sm):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = "{0}/{1}/AnalysisResults.root".format(yamlConfig["input_path"], yamlConfig["train"])
    outputPath = "{0}/{1}".format(yamlConfig["input_path"], yamlConfig["train"])
    
    file = ROOT.TFile.Open(fileName);
    if not file or file.IsZombie():
        print("Could not open file '{0}'! Aborting...".format(fileName))
        exit(1)

    analysis = []
    
    for config in yamlConfig["configs"]:
        if not config["active"]:
            continue
        if "name" in config:
            name = config["name"]
        else:
            name = config["detector"]

        ana = ClusterTriggerTurnOnAnalysis(name, config["detector"], file, outputPath, config["triggers"], config["base_trigger"], yamlConfig["clusters"], yamlConfig["cells"])
        if "max_cluster_e" in config: ana.fMaxObservableValue = config["max_cluster_e"]
        analysis.append(ana)
        if "patches" in config:
            for patch_config in config["patches"]:
                ana = PatchTriggerTurnOnAnalysis(name, config["detector"], file, outputPath, config["triggers"], config["base_trigger"], patch_config["trigger"], patch_config["type"])
                if "max_patch_e" in patch_config: ana.fMaxObservableValue = patch_config["max_patch_e"]
                analysis.append(ana)

        if sm:
            if config["detector"] == "EMCal":
                SMlist = range(0,12)
            elif config["detector"] == "DCal":
                SMlist = range(12,18)
            else:
                SMlist = []
            for iSM in SMlist:
                analysis.append(ClusterTriggerTurnOnAnalysis("SM{0}".format(iSM), file, outputPath, config["triggers"], config["base_trigger"], yamlConfig["clusters"], yamlConfig["cells"]))

    for ana in analysis:
        ana.GetHistograms()
        ana.ProjectHistograms()
        ana.Ratios()
        ana.Plot()

    if sm: CompareSpectra(analysis, "SM", "SM1_Cluster")
    
    globalList.append(analysis)
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='EMCal trigger QA.')
    parser.add_argument('yaml', metavar='config.yaml',
                        help='YAML configuration file')
    parser.add_argument("--sm", action='store_const',
                        default=False, const = True,
                        help='Plot trigger turn-on curves SM by SM.')
    args = parser.parse_args()
    
    f = open(args.yaml, 'r')
    config = yaml.load(f)
    f.close()

    main(config, args.sm)
    
    IPython.embed()
