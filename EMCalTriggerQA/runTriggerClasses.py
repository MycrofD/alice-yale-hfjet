#!/usr/bin/env python
#python script to do trigger class analysis

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math
import yaml

import IPython

globalList = []

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
    
    events = dict()
    shares = dict()
    
    #print("Analyzing {0}".format(trigger["class"]))
    
    for triggerComp in triggerList:
        events[triggerComp["label"]] = 0
        for triggerCompClass in triggerComp["class"]:
            events[triggerComp["label"]] += htrigger.GetBinContent(htrigger.GetXaxis().FindBin(triggerCompClass))
        #print("Events in {0} is {1:.2f} ({2:.2f})".format(triggerComp["label"], events[triggerComp["label"]], math.sqrt(events[triggerComp["label"]])))
    
    res = "|  *{0: ^10}*  ".format(trigger["label"])   
    
    for triggerComp in triggerList:
        if events[trigger["label"]] > 0:
            shares[triggerComp["label"]] = events[triggerComp["label"]] / events[trigger["label"]]
        else:
            shares[triggerComp["label"]] = 0
            
        err = shares[triggerComp["label"]] * (1. / math.sqrt(events[triggerComp["label"]]) + 1. / math.sqrt(events[trigger["label"]]))
        #print("Fraction of events in {0} is {1:.4f} ({2:.4f})".format(triggerComp["label"], shares[triggerComp["label"]], err))    
        res += "|  {0: ^12.3f}  ".format(shares[triggerComp["label"]])
    res += "|"
    print(res)

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
    
    for trigger in triggerList:
        AnalyzeTriggerClass(trigger, triggerList, file)
        
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