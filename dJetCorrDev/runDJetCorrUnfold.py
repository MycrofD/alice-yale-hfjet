#!/usr/bin/env python
#python script to use the class DJetCorrUnfold

import ROOT
import argparse
from ROOT import gROOT
from ROOT import TGaxis

import IPython

import commonFunctions
import runDJetCorrAnalysis
import runDJetCorrResponse

def main(trainData, trainResp, loadLibs = True, noSave = False, inputPath = "$JETRESULTS"):
	
	ROOT.TH1.AddDirectory(False)
	
	if loadLibs:
		commonFunctions.LoadDJetCorrClasses()

	from ROOT import DJetCorrAnalysis
	from ROOT import DJetCorrResponse
	from ROOT import DJetCorrUnfold
  
  	ana = runDJetCorrAnalysis.main(trainData, False, False, False, False, True, False, False, inputPath)
   	resp = runDJetCorrResponse.main(trainResp, False, False, False, inputPath)
  
  	unfold = DJetCorrUnfold(ana, resp)
   	unfold.SetDataParamIndex(0)
   	unfold.SetRespParamIndex(0)
   	unfold.SetSavePlots(not noSave)
   	unfold.SetSaveOutputFile(not noSave);
   	unfold.SetRegParam(2, 8, 2)
    #unfold.SetUseEfficiency(False)
    #unfold.SetUseKinEfficiency(False)
	unfold.Start()

	unfold.UpdateAllCanvases()

	return unfold

if __name__ == '__main__':
    # runDJetCorrUnfold.py executed as script
    
    parser = argparse.ArgumentParser(description='D meson jet correlation response matrix.')
    parser.add_argument('trainData', metavar='trainData',
                        help='Train with data to be analyzed')
    parser.add_argument('trainResp', metavar='trainResp',
						help='MC train used for the response matrix')
    parser.add_argument('--no-Libs', action='store_const',
                        default=False, const=True,
                        help='Load the DJetCorr libraries')
    parser.add_argument('--no-save', action='store_const',
						default=False, const=True,
						help='Do not save the output file')
    parser.add_argument('--inputPath', metavar='inputPath',
                        default="$JETRESULTS",
                        help='Input path')
    args = parser.parse_args()
    
    main(args.trainData, args.trainResp, not args.no_Libs, args.no_save, args.inputPath)
    
    IPython.embed()