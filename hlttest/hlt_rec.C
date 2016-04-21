#include <AliHLTSystem.h>
#include <AliDAQ.h>
#include <AliHLTConfiguration.h>
#include <AliReconstruction.h>
#include <AliHLTPluginBase.h>
#include <AliCDBManager.h>

#include <iostream>

void hlt_rec(const char* input="./")
{
  // For real data:
  AliCDBManager::Instance()->SetDefaultStorage("raw://");
  
  if (!gSystem->AccessPathName("galice.root")){
    std::cerr << "please delete the galice.root or run at different place." <<std:: endl;
    return;
  }
  
  if (!input) {
    std::cerr << "please specify input or run without arguments" << std::endl;
    return;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  
  AliHLTSystem* gHLT = AliHLTPluginBase::GetInstance();
  
  /*
  gHLT->LoadComponentLibraries("libAliHLTUtil.so");
  gHLT->LoadComponentLibraries("libAliHLTRCU.so");
  gHLT->LoadComponentLibraries("libAliHLTCalo.so");
  gHLT->LoadComponentLibraries("libAliHLTEMCAL.so");
  gHLT->LoadComponentLibraries("libAliHLTPHOS.so");
  gHLT->LoadComponentLibraries("libAliHLTGlobal.so");
  */
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  // New handling of the specification: Specification will be the DDL ID
  // There will be 1 publisher per ddl, but the raw analyzers receive input from ALL DDLs
  
  Int_t ddlOffset = 4608; // The DDL offset for EMCAL (for PHOS the number is 1792)

  TString arg;
  TString rps;
  
  for (int module = 0; module <= AliDAQ::NumberOfDdls("EMCAL"); module++) {
    TString publisher;
    // raw data publisher components
    publisher.Form("EMCAL-RP_%02d", module);
    arg.Form("-verbose -minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec %d ", ddlOffset + module, module);
                
    if (rps.Length()) rps += " ";
    rps += publisher;
    AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
  }
                
  // Raw analyzer
  arg = "";
  AliHLTConfiguration rawConf("EMCAL-RA", "EmcalRawCrude", rps.Data(), arg.Data());
                
  // Raw analyzer for TRU data
  arg = "";
  AliHLTConfiguration truConf("EMCAL-TRU", "EmcalTruAnalyzer", rps.Data(), arg.Data());

  // STU raw analyser
  arg = "";
  AliHLTConfiguration stuConf("EMCAL-STU", "EmcalStuAnalyzer", rps.Data(), arg.Data());

  // digit maker components
  arg = TString::Format("-sethighgainfactor 0.0153 -setlowgainfactor 0.2448 -setdigitthresholds 0.005 0.002");
  TString dmInput = "EMCAL-RA";
  AliHLTConfiguration dmConf("EMCAL-DM", "EmcalDigitMaker", dmInput.Data(), arg.Data());
            
  //arg = TString::Format("-digitthreshold 0.005 -recpointthreshold 0.1 -modulemode");
  //TString clInput = "EMCAL-DM";
  //AliHLTConfiguration clConf("EMCAL-CF", "EmcalClusterizer", clInput.Data(), arg.Data());
            
  // Tigger data merger
  arg = "";
  TString tdInput = "EMCAL-TRU EMCAL-STU";
  AliHLTConfiguration trgdata("EMCAL-TRG", "EmcalTriggerDataMaker", tdInput.Data(), arg.Data());

  arg = "";
  TString tmInput = "EMCAL-DM EMCAL-TRG";
  AliHLTConfiguration trgmaker("EMCAL-TM", "EmcalTriggerMaker", tmInput.Data(), arg.Data());

  // The call for the trigger QA
  arg = "-newTriggerBitConfig -noHistoReset -debugLevel0";
  TString tqaInput = "EMCAL-DM EMCAL-TRG EMCAL-TM";
  AliHLTConfiguration htConf("EMCAL-TQA", "EmcalTriggerQA", tqaInput, arg.Data());
  
  // ESD convertert
  //arg = "";
  //TString ecInput = "";
  //AliHLTConfiguration esdcconf("ESD-CONVERTER", "GlobalEsdConverter", ecInput.Data(), arg.Data());
  
  // EMC histo maker 
  //arg = TString::Format("-pushfraction 5 -beverbose 1");
  //TString fwInput = "";
  //AliHLTConfiguration hfConf("emcalHisto", "EmcalRawHistoMaker", fwInput.Data(), arg.Data());
  
  // Write the root file 
  arg = "-datafile roothisto.root -concatenate-events -overwrite ";
  TString rwInput = "EMCAL-TQA";
  AliHLTConfiguration rwConf("rootFileHisto","ROOTFileWriter", rwInput.Data(), arg.Data());
  
  //////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off 
  //
  /////////////////////////////////////////////////////////////////////
  
  AliReconstruction rec;
  
  // uncomment for simulation
  //rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  //rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  //rec.SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2013/OCDB");
  rec.SetDefaultStorage("raw://");
  
  rec.SetRunReconstruction(":");
  rec.SetEventRange(0,-1);
  rec.SetInput(input);
  
  //rec.SetInput("/Users/sa639/Documents/Work/ALICE/Data/2013/LHC13d/000195873/raw/13000195873000.10.root");
  //rec.SetInput("/Users/sa639/Documents/Work/ALICE/Data/2015/LHC15l/000239913/raw/15000239913030.1000.root");
  //rec.SetInput("/Users/sa639/Documents/Work/ALICE/Data/2015/LHC15i/236141/raw/15000236141030.101.root");
  //rec.SetInput("/Volumes/DATA/ALICE/Data/2015/LHC15o/000244918/raw/15000244918019.100.root");
  //rec.SetInput("/Volumes/DATA/ALICE/Data/2015/LHC15j/000237673/raw/15000237673030.101.root");
  rec.SetInput("/Volumes/DATA/ALICE/Data/2016/LHC16c/000251269/raw/16000251269030.200.root");
  //rec.SetInput("files_LHC16c_251269.txt");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunMultFinder(kFALSE);
  rec.SetRunVertexFinderTracks(kFALSE);
  rec.SetRunV0Finder(kFALSE);
  rec.SetRunCascadeFinder(kFALSE);
  
  rec.SetRunReconstruction("HLT");
  //rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  
  //TString option = "libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTEMCAL.so libAliHLTGlobal.so chains=";
  //option+="ESD-CONVERTER";
  //option += "rootFileHisto loglevel=0x5f";
  TString option = "chains=rootFileHisto loglevel=0x5f";
  rec.SetOption("HLT", option);
  
  rec.Run();
}
