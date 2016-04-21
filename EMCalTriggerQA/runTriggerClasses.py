#!/usr/bin/env python
#python script to do trigger class analysis

import ROOT
import argparse
from ROOT import gROOT
from enum import Enum
import array
import math

import IPython

globalList = []

def AnalyzeTriggerClass(trigger, triggerList, file):
    
    if trigger:
        triggerSimple = trigger.split("-")[0]
        listName = "AliEmcalTriggerQATask_" + triggerSimple + "_histos"
        hlistName = "histosAliEmcalTriggerQATask_" + triggerSimple
    else:
        triggerSimple = ""
        listName = "AliEmcalTriggerQATask_histos"
        hlistName = "histosAliEmcalTriggerQATask"
    
    list = file.Get(listName)
    if not list:
        file.ls()
        print "Could not get list '" + listName + "' from file '" + file.GetName() + "'! Aborting..."
        exit(1)
        
    hlist = list.FindObject(hlistName)
    if not hlist:
        list.Print()
        print "Could not get hash list '" + hlistName + "' in list '" + listName + "' from file '" + fileName + "'! Aborting..."
        exit(1)
    
    
    htrigger = list.FindObject("fHistTriggerClasses")
    if not htrigger:
        print("Could not find histogram fHistTriggerClasses in list {0}".format(list.GetName()))
        return
    
    events = dict()
    shares = dict()
    
    print("Analyzing {0}".format(trigger))
    
    for triggerComp in triggerList:
        events[triggerComp] = htrigger.GetBinContent(htrigger.GetXaxis().FindBin(triggerComp))
        print("Events in {0} is {1:.2f} ({2:.2f})".format(triggerComp, events[triggerComp], math.sqrt(events[triggerComp])))
    
    for triggerComp in triggerList:
        if events[trigger] > 0:
            shares[triggerComp] = events[triggerComp] / events[trigger] * 100
        else:
            shares[triggerComp] = 0
            
        err = shares[triggerComp] * (1. / math.sqrt(events[triggerComp]) + 1. / math.sqrt(events[trigger]))
        print("Fraction of events in {0} is {1:.2f} ({2:.2f})".format(triggerComp, shares[triggerComp], err))    
    

def main(train, triggers, inputPath="/Users/sa639/Documents/Work/ALICE/TriggerQA"):
    
    ROOT.TH1.AddDirectory(False)
    ROOT.gStyle.SetOptTitle(False)
    ROOT.gStyle.SetOptStat(0)
    
    fileName = "{0}/{1}/AnalysisResults.root".format(inputPath, train)
    
    file = ROOT.TFile.Open(fileName);
    if not file:
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
        
    if file.IsZombie():
        print "Could not open file '" + fileName + "'! Aborting..."
        exit(1)
    
    triggerList = triggers.split(",")
    for trigger in triggerList:
        AnalyzeTriggerClass(trigger, triggerList, file)
        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Trigger classes analysis.')
    parser.add_argument('train', metavar='train',
                        help='Train to be analyzed')
    parser.add_argument('--trigger', metavar='trigger',
                        default="EMC7",
                        help='Trigger')
    parser.add_argument('--input-path', metavar='input-path',
                        default="/Users/sa639/Documents/Work/ALICE/TriggerQA",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.train, args.trigger, args.input_path)
    
    IPython.embed()