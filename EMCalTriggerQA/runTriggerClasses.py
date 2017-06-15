#!/usr/bin/env python
#python script to do trigger class analysis

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math
import yaml
from collections import OrderedDict

import IPython

globalList = []

class TriggerCombination:
    def __init__(self, label, rate):
        self.fLabel = label
        self.fRate = rate
        
    def Print(self):
        print("|  *{0: ^20}*  |  {1:.4f}  |".format(self.fLabel, self.fRate))
        
        
def PrintTriggerSuppression(shares, baseTrigger, subtract):
    label = ""    
    rate = 0
    
    for i,trigger in enumerate(list):
        if label:
            label += "+"
        label += trigger
        #print("Rate for trigger {0} is {1:.3f}".format(trigger, self.fTriggerList[baseTrigger].fShares[trigger]))
        rate += shares[baseTrigger].fShares[trigger] 
        for j in range(i+1,len(list)):
            #print("Subtracting share with {0}: {1:.3f} * {2:.3f} = {3:.3f}".format(list[j], self.fTriggerList[baseTrigger].fShares[trigger], self.fTriggerList[list[j]].fShares[trigger], self.fTriggerList[baseTrigger].fShares[trigger] * self.fTriggerList[list[j]].fShares[trigger]))
            rate -= self.fTriggerList[baseTrigger].fShares[trigger] * self.fTriggerList[trigger].fShares[list[j]]
    
    rate -= subtract
    
    #print("{0} = {1:.3f}".format(label, rate))
    
    return TriggerCombination(label, rate)

def AnalyzeTriggerClass(trigger, triggerList, file):
    
    listName = "AliEmcalTriggerQATask_" + trigger["label"] + "_histos"
    hlistName = "histosAliEmcalTriggerQATask_" + trigger["label"]
    
    list = file.Get(listName)
    if not list:
        file.ls()
        print("Could not get list '{0}' from file '{1}'! Aborting...".format(listName, file.GetName()))
        exit(1)
        
    hlist = list.FindObject(hlistName)
    if not hlist:
        list.Print()
        print("Could not get hash list '{0}' in list '{1}' from file '{2}'! Aborting...".format(hlistName, listName, fileName))
        exit(1)

    htrigger = list.FindObject("fHistTriggerClasses")
    if not htrigger:
        print("Could not find histogram fHistTriggerClasses in list {0}".format(list.GetName()))
        return
    
    hevents = list.FindObject("fHistEventCount")
    if not hevents:
        print("Could not find histogram fHistEventCount in list {0}".format(list.GetName()))
        return
    
    events = dict()
    shares = dict()
    
    events[trigger["label"]] = hevents.GetBinContent(1)
    
    #print("Analyzing {0}".format(trigger["class"]))
    
    for triggerComp in triggerList:
        if len(triggerComp["class"]) != 1:
            continue
        if triggerComp["label"] == trigger["label"]:
            continue
        triggerCompClass = triggerComp["class"][0]
        events[triggerComp["label"]] = htrigger.GetBinContent(htrigger.GetXaxis().FindBin(triggerCompClass))
        #print("Events in {0} is {1:.2f} ({2:.2f})".format(triggerComp["label"], events[triggerComp["label"]], math.sqrt(events[triggerComp["label"]])))
    
    res = "|  *{0: ^10}*  ".format(trigger["label"])   
    
    for triggerComp in triggerList:
        if not events.has_key(triggerComp["label"]):
            continue
        if events[trigger["label"]] > 0:
            shares[triggerComp["label"]] = events[triggerComp["label"]] / events[trigger["label"]]
        else:
            shares[triggerComp["label"]] = 0
            
        err = shares[triggerComp["label"]] * (1. / math.sqrt(events[triggerComp["label"]]) + 1. / math.sqrt(events[trigger["label"]]))
        #print("Fraction of events in {0} is {1:.4f} ({2:.4f})".format(triggerComp["label"], shares[triggerComp["label"]], err))    
        res += "|  {0: ^12.4f}  ".format(shares[triggerComp["label"]])
    res += "|"
    print(res)
    return shares

class TriggerCombination:
    def __init__(self, label, rate):
        self.fLabel = label
        self.fRate = rate
        
    def Print(self):
        print("|  *{0: ^20}*  |  {1:.4f}  |".format(self.fLabel, self.fRate))

def PrintTriggerSuppression(shares, list, baseTrigger, subtract):
    label = ""    
    rate = 0
    
    for i,trigger in enumerate(list):
        if label:
            label += "+"
        label += trigger
        #print("Rate for trigger {0} is {1:.3f}".format(trigger, self.fTriggerList[baseTrigger].fShares[trigger]))
        rate += shares[baseTrigger][trigger] 
        for j in range(i+1,len(list)):
            #print("Subtracting share with {0}: {1:.3f} * {2:.3f} = {3:.3f}".format(list[j], self.fTriggerList[baseTrigger].fShares[trigger], self.fTriggerList[list[j]].fShares[trigger], self.fTriggerList[baseTrigger].fShares[trigger] * self.fTriggerList[list[j]].fShares[trigger]))
            rate -= shares[baseTrigger][trigger] * shares[trigger][list[j]]
    
    rate -= subtract
    
    #print("{0} = {1:.3f}".format(label, rate))
    
    return TriggerCombination(label, rate)

def main(train, confFile, inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    f = open(confFile, 'r')
    config = yaml.load(f)
    f.close()
    
    fileName = "{0}/{1}/AnalysisResults.root".format(inputPath, train)
    
    file = ROOT.TFile.Open(fileName);
    if not file:
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
        
    if file.IsZombie():
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
    
    triggerList = config["triggers"]
    res = "|  * Triggers *  "
    for trigger in triggerList:
        res += "|  *{0: ^10}*  ".format(trigger["label"])
    res += "|"
    print(res)
    
    shares = dict()
    
    for trigger in triggerList:
        shares[trigger["label"]] = AnalyzeTriggerClass(trigger, triggerList, file)
        
    highTriggers = ["CEMC7EG1", "CEMC7EJ1", "CDMC7DG1", "CDMC7DJ1"]
    lowTriggers = ["CEMC7EG2", "CEMC7EJ2", "CDMC7DG2", "CDMC7DJ2"]
    
    highTriggersGammaOnly = ["CEMC7EG1", "CDMC7DG1"]
    lowTriggersGammaOnly = ["CEMC7EG2", "CDMC7DG2"]
    
    highTriggersJetOnly = ["CEMC7EJ1", "CDMC7DJ1"]
    lowTriggersJetOnly = ["CEMC7EJ2", "CDMC7DJ2"]
    
    EMCalHighTriggers = ["CEMC7EG1", "CEMC7EJ1"]
    EMCalLowTriggers = ["CEMC7EG2", "CEMC7EJ2"]
    
    DCalHighTriggers = ["CDMC7DG1", "CDMC7DJ1"]
    DCalLowTriggers = ["CDMC7DG2", "CDMC7DJ2"]
    
    summary = OrderedDict()
    
    summary["highTriggers"] = PrintTriggerSuppression(shares, highTriggers, "CEMC7CDMC7", 0)
    summary["lowTriggers"] = PrintTriggerSuppression(shares, lowTriggers, "CEMC7CDMC7", summary["highTriggers"].fRate)
    
    summary["highTriggersGammaOnly"] = PrintTriggerSuppression(shares, highTriggersGammaOnly, "CEMC7CDMC7", 0)
    summary["lowTriggersGammaOnly"] = PrintTriggerSuppression(shares, lowTriggersGammaOnly, "CEMC7CDMC7", summary["highTriggersGammaOnly"].fRate)

    summary["highTriggersJetOnly"] = PrintTriggerSuppression(shares, highTriggersJetOnly, "CEMC7CDMC7", 0)
    summary["lowTriggersJetOnly"] = PrintTriggerSuppression(shares, lowTriggersJetOnly, "CEMC7CDMC7", summary["highTriggersJetOnly"].fRate)

    summary["EMCalHighTriggers"] = PrintTriggerSuppression(shares, EMCalHighTriggers, "CEMC7CDMC7", 0)
    summary["EMCalLowTriggers"] = PrintTriggerSuppression(shares, EMCalLowTriggers, "CEMC7CDMC7", summary["EMCalHighTriggers"].fRate)
    
    summary["DCalHighTriggers"] = PrintTriggerSuppression(shares, DCalHighTriggers, "CEMC7CDMC7", 0)
    summary["DCalLowTriggers"] = PrintTriggerSuppression(shares, DCalLowTriggers, "CEMC7CDMC7", summary["DCalHighTriggers"].fRate)
    
    print("|  *{0: ^20}*  |  *{1}*  |".format("Trigger combination", "Suppression over L0"))
    for trigger in summary.itervalues():
        trigger.Print()

        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Trigger classes analysis.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--conf', metavar='conf',
                        help='YAML configuration file')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.conf, args.input_path)
    
    IPython.embed()
