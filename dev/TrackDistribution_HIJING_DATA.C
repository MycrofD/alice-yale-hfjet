void LoadLibs();
TH1* GetTrackHIJING();

const char* localTrainPath  = "$JETRESULTS";
const char* train           = "LEGO_TRAIN_MC_289";  // AOD149 LHC12a17d_fix 0-10%
const char* fname           = "AnalysisResults-GOOD.root";

const char* trainData       = "LEGO_TRAIN_374";     // AOD149 LHC11h
const char* myMacrosPath    = "/Users/saiola/Documents/Work/ALICE/aliemcaldev/saiola/MyMacros";

const Int_t nbins = 17;
Double_t bins[nbins+1] = {0,1,2,3,4,5,6,8,10,12,14,16,18,20,25,30,35,40};

void TrackDistribution_HIJING_DATA()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  LoadLibs();

  TH1::AddDirectory(kFALSE);

  TH1* trackHijing = GetTrackHIJING();
  TH1* trackData = GetTrackData();  

  TCanvas *canvas = new TCanvas("c","c");
  canvas->SetLogy();
  canvas->cd();

  TH1* myBlankHisto = new TH1F("myBlankHisto","myBlankHisto",150,0,150);
  myBlankHisto->GetYaxis()->SetRangeUser(1e-2,1e4);
  myBlankHisto->GetYaxis()->SetTitle("#frac{1}{#it{N}_{evt}} #frac{d#it{N}}{d#it{p}_{T}} (GeV/#it{c})^{-1}");
  myBlankHisto->GetXaxis()->SetRangeUser(0,10);
  myBlankHisto->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  myBlankHisto->Draw();

  trackData->SetMarkerStyle(kFullCircle);
  trackData->SetMarkerSize(0.8);
  trackData->SetMarkerColor(kBlack);
  trackData->SetLineColor(kBlack);
  trackData->Draw("same");

  Printf("Data bin 1 = %.2f, bin 2 = %.2f", trackData->GetBinContent(1), trackData->GetBinContent(2));

  Double_t integralData = 0;
  Double_t errIntegralData = 0;
  integralData = trackData->IntegralAndError(1,trackData->GetNbinsX(),errIntegralData,"width");

  Printf("Integral = %.2f +/- %.2f", integralData, errIntegralData);  

  trackHijing->SetMarkerStyle(kOpenCircle);
  trackHijing->SetMarkerSize(0.8);
  trackHijing->SetMarkerColor(kRed);
  trackHijing->SetLineColor(kRed);
  trackHijing->Draw("same");

  Printf("HIJING bin 1 = %.2f, bin 2 = %.2f", trackHijing->GetBinContent(1), trackHijing->GetBinContent(2));

  Double_t integralHijing = 0;
  Double_t errIntegralHijing = 0;
  integralHijing = trackHijing->IntegralAndError(1,trackHijing->GetNbinsX(),errIntegralHijing,"width");

  Printf("Integral = %.2f +/- %.2f", integralHijing, errIntegralHijing);  

  TLegend *leg = new TLegend(0.40, 0.75, 0.85, 0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(trackData, "Data 0-10% (LHC11h)", "pe");
  //leg->AddEntry(trackHijing, "HIJING 0-10% (LHC12a17d) - Detector level w/ emb. signal", "pe");
  leg->AddEntry(trackHijing, "HIJING 0-10% (LHC12a17d) - Particle level", "pe");
  leg->Draw();
}


TH1* GetTrackData()
{
  TString fileName(Form("%s/%s.root", myMacrosPath, trainData));
  TFile *currentFile = TFile::Open(fileName);

  if (!currentFile || currentFile->IsZombie()) {
    Printf("Could not open file %s!", fileName.Data());
    return 0;
  }
  
  TH3 *hist = dynamic_cast<TH3*>(currentFile->Get("TrackPt_0"));
  TH1 *cent = dynamic_cast<TH3*>(currentFile->Get("Centrality"));

  Double_t events = cent->Integral(cent->GetXaxis()->FindBin(0.0001), cent->GetXaxis()->FindBin(9.9999));
  
  TH1 *pt = hist->ProjectionZ("pt");

  TH1 *ptRebin = pt->Rebin(nbins,"ptrebin",bins);

  delete pt;
  pt = 0;

  ptRebin->Sumw2();
  ptRebin->Scale(1./events,"width");
  
  return ptRebin;
}


TH1* GetTrackHIJING()
{
  TString fileName(Form("%s/%s/%s",localTrainPath, train, fname));
  TFile *currentFile = TFile::Open(fileName);

  if (!currentFile || currentFile->IsZombie()) {
    Printf("Could not open file %s!", fileName.Data());
    return 0;
  }

  const UInt_t kStepReconstructed = 0;   //Reconstructed particles vs pTrec
  const UInt_t kStepSecondaries = 1;     //Reconstructed secondaries vs pTmc
  const UInt_t kStepReconstructedMC = 2; //Reconstructed particles vs pTmc
  const UInt_t kStepMCAcceptance = 3;    //All generated primaries in acceptance

  //Get MC containers
  TString dirName(Form("PWG4_HighPtSpectraCent%dTrackType%dCuts%d%s", 10, 0, 0, "kMB"));

  TString contPosName(dirName);
  contPosName += "/containerPos";
  AliCFContainer *dataPos = dynamic_cast<AliCFContainer*>(currentFile->Get(contPosName));
  
  TString contNegName(dirName);
  contNegName += "/containerNeg";
  AliCFContainer *dataNeg = dynamic_cast<AliCFContainer*>(currentFile->Get(contNegName));

  TString contListName(dirName);
  contListName += "/";
  contListName += Form("chist0HighPtSpectraCent%dTrackType%dCuts%d%s", 10, 0, 0, "kMB");
  TList *list = dynamic_cast<TList*>(currentFile->Get(contListName));

  TH1* hEvents = dynamic_cast<TH1*>(list->FindObject("fNEventSel"));
  Double_t nevents = hEvents->GetBinContent(1);

  AliCFGridSparse* gridPos = dataPos->GetGrid(kStepReconstructed);
  AliCFGridSparse* gridNeg = dataNeg->GetGrid(kStepReconstructed);

  THnSparse *thnPos = gridPos->GetGrid();
  THnSparse *thnNeg = gridNeg->GetGrid();

  TH1* ptPos = thnPos->Projection(0);
  TH1* ptNeg = thnNeg->Projection(0);

  TH1* pt = static_cast<TH1*>(ptPos->Clone("pt"));
  pt->Add(ptNeg);

  TH1 *ptRebin = pt->Rebin(nbins,"ptrebin",bins);

  ptRebin->Scale(1./nevents, "width");

  return ptRebin;
}

void LoadLibs()
{
  printf("Load libraries\n");

  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
}
