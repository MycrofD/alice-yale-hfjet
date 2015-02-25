#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <Riostream.h>
#include <TF1.h>
#include <TLegend.h>
#include <TList.h>
#include <TPaveText.h>
#include <TLegendEntry.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TSeqCollection.h>
#include <TObjArray.h>

#define EPSILON 1e-3

Int_t nBins = 0;
Double_t *bins = 0;

Double_t minX = -60;
Double_t maxX = 160;

class DeltaPtContainer {
public:

  enum EDeltaPtMethod {
    kRC=0,
    kSPE=1
  };

  DeltaPtContainer();
  DeltaPtContainer(Int_t cent, EDeltaPtMethod meth, Int_t ptcut, Double_t r);

  TString GetHistoName();
  TString GetFileName();
  TString GetPtCutLabel();
  TString GetRadiusLabel();
  TString GetEtaLabel();
  TString GetMethodLabel();
  TString GetCentLabel();
  void Prepare();
  void RecursiveFit();
  void AddLegendEntry(TLegend* leg);
  void Plot(Style_t marker, Color_t color);

  Int_t fCent;
  enum EDeltaPtMethod fMethod;
  Int_t fConstCut;
  Double_t fRadius;
  TH1* fDeltaPtHisto;
  TH1* fDeltaPtHistoSmallBins;
  TF1* fFitFunction;
  Double_t fMean;
  Double_t fSigma;
  TFile* fInputFile;
  Bool_t fShowMethod;
  Bool_t fShowPtCut;
  Bool_t fShowRadius;
  Bool_t fShowCent;

private:
  Bool_t fInitialized;
};

void GenerateBins();
void myOptions(Int_t lStat=0);
void myLegendSetUp(TLegend *currentLegend, Float_t currentTextSize);

void MakeDeltaPtPlot(Int_t rWrite = 2)
{
  // rWrite     = 0 / 1 / 2 -> unsaved / saved as .eps / saved as .pdf ;

  TH1::AddDirectory(kFALSE);

  GenerateBins();
  myOptions();

  const Int_t n = 3;
  
  /*
  DeltaPtContainer containers[2] = {DeltaPtContainer(0, DeltaPtContainer::kRC, 0, 0.2), DeltaPtContainer(0, DeltaPtContainer::kRC, 1, 0.2)};
  Style_t markers[2] = {kFullCircle, kFullSquare};
  Color_t colors[2] = {kBlue+2,kRed+2};
  */

  /*
  DeltaPtContainer containers[2] = {DeltaPtContainer(0, DeltaPtContainer::kRC, 1, 0.2), DeltaPtContainer(0, DeltaPtContainer::kRC, 1, 0.3)};
  Style_t markers[2] = {kFullCircle, kFullTriangleUp};
  Color_t colors[2] = {kBlue+2,kGreen+2};
  */

  DeltaPtContainer containers[n] = {DeltaPtContainer(0, DeltaPtContainer::kRC, 0, 0.1), DeltaPtContainer(0, DeltaPtContainer::kRC, 0, 0.2), DeltaPtContainer(0, DeltaPtContainer::kRC, 0, 0.3)};
  Style_t markers[n] = {kFullSquare, kFullCircle, kFullTriangleUp};
  Color_t colors[n] = {kMagenta+2, kBlue+2, kGreen+2};

  TCanvas *myCan = new TCanvas("c","c",800,600);
  myCan->cd();
  myCan->SetLogy();

  myCan->SetLeftMargin(0.18);
  myCan->SetBottomMargin(0.15);
  myCan->SetRightMargin(0.04);
  myCan->SetTopMargin(0.04);

  TH1* myBlankHisto = new TH1F("myBlankHisto", "myBlankHisto", 500, -250, 250);
  myBlankHisto->GetXaxis()->SetTitleOffset(1.30);
  myBlankHisto->GetYaxis()->SetTitleOffset(1.60);
  myBlankHisto->GetXaxis()->SetTitleSize(0.045);
  myBlankHisto->GetYaxis()->SetTitleSize(0.045);
  myBlankHisto->GetXaxis()->SetLabelSize(0.045);
  myBlankHisto->GetYaxis()->SetLabelSize(0.045);
  myBlankHisto->GetXaxis()->SetRangeUser(minX,maxX);
  myBlankHisto->GetYaxis()->SetRangeUser(1e-7,5e-1);
  myBlankHisto->GetXaxis()->SetTitle("#delta#it{p}_{T} (GeV/#it{c})");
  myBlankHisto->GetYaxis()->SetTitle("Probability density");

  myBlankHisto->Draw();

  TLegend *leg = new TLegend(0.56,0.32,0.90,0.69);
  myLegendSetUp(leg,0.035);

  TString radiusLabel;
  TString etaLabel;
  TString ptCutLabel;
  TString methodLabel;
  TString centLabel;

  for (Int_t i = 0; i < n; i++) {
    containers[i].fShowRadius = kTRUE;
    containers[i].fShowMethod = kFALSE;
    containers[i].fShowCent = kFALSE;
    containers[i].fShowPtCut = kFALSE;

    containers[i].Prepare();
    containers[i].Plot(markers[i],colors[i]);
    containers[i].AddLegendEntry(leg);

    radiusLabel = containers[i].GetRadiusLabel();
    etaLabel = containers[i].GetEtaLabel();
    ptCutLabel = containers[i].GetPtCutLabel();
    methodLabel = containers[i].GetMethodLabel();
    centLabel = containers[i].GetCentLabel();
  }

  TLatex *alice = new TLatex(0.25,0.90,"ALICE   Pb-Pb   #sqrt{#it{s}_{NN}} = 2.76 TeV");
  alice->SetNDC();
  alice->SetTextSize(0.045);
  alice->Draw();

  TLegend *myLegend = new TLegend(0.60,0.71,0.7,0.87);
  myLegendSetUp(myLegend,0.035);
  myLegend->AddEntry((TObject*)0,centLabel,"");
  //myLegend->AddEntry((TObject*)0,radiusLabel,"");    
  myLegend->AddEntry((TObject*)0,etaLabel,"");
  myLegend->AddEntry((TObject*)0,methodLabel,"");
  myLegend->AddEntry((TObject*)0,ptCutLabel,"");
  myLegend->Draw();

  leg->Draw();

  myCan->Update();

  if (rWrite == 1)  myCan->SaveAs("fig_DeltaPt.eps");
  else if (rWrite == 2)  myCan->SaveAs("fig_DeltaPt.pdf");
}

void GenerateBins()
{
  Double_t minSmallBins = -50;
  Double_t maxSmallBins = 50;
  Int_t nSmallBins = 40;

  Double_t maxLargeBins = 150;
  Int_t nLargeBins = 20;

  nBins = nSmallBins + nLargeBins;
  bins = new Double_t[nBins+1];

  Double_t smallBinWidth = (maxSmallBins - minSmallBins) / nSmallBins;
  Double_t largeBinWidth = (maxLargeBins - maxSmallBins) / nLargeBins;  
  
  for (Int_t i = 0; i <= nSmallBins; i++) {
    bins[i] = minSmallBins + smallBinWidth*i;
    //Printf("bins[%d] = %.1f", i, bins[i]);
  }
  for (Int_t i = nSmallBins+1; i <= nBins; i++) {
    bins[i] = maxSmallBins + largeBinWidth*(i-nSmallBins);
    Printf("bins[%d] = %.1f", i, bins[i]);
  }
}

void myOptions(Int_t lStat)
{
  // Set gStyle
  Int_t font = 42;
  // From plain
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(10);
  gStyle->SetCanvasColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(10);
  gStyle->SetStatBorderSize(1);
  gStyle->SetLegendBorderSize(1);
  //
  gStyle->SetDrawBorder(0);
  gStyle->SetTextFont(font);
  gStyle->SetStatFont(font);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.97);
  gStyle->SetStatY(0.98);
  gStyle->SetStatH(0.03);
  gStyle->SetStatW(0.3);
  gStyle->SetTickLength(0.02,"y");
  gStyle->SetEndErrorSize(3);
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelFont(font,"xyz"); 
  gStyle->SetLabelOffset(0.01,"xyz");
  gStyle->SetTitleFont(font,"xyz");  
  gStyle->SetTitleOffset(1.2,"xyz");  
  gStyle->SetTitleSize(0.045,"xyz");  
  gStyle->SetMarkerSize(1); 
  gStyle->SetPalette(1,0);
  gStyle->SetErrorX(0);
  gStyle->SetLegendBorderSize(0);

  if (lStat) {
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
  }
  else {
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
  }
}

void myLegendSetUp(TLegend *currentLegend, Float_t currentTextSize)
{
  currentLegend->SetTextFont(42);
  currentLegend->SetBorderSize(0);
  currentLegend->SetFillStyle(0);
  currentLegend->SetFillColor(0);
  currentLegend->SetMargin(0.25);
  currentLegend->SetTextSize(currentTextSize);
  currentLegend->SetEntrySeparation(0.5);
}

DeltaPtContainer::DeltaPtContainer() : 
  fCent(0), fMethod(kRC), fConstCut(0), fRadius(0.2), 
  fDeltaPtHisto(0), fDeltaPtHistoSmallBins(0), fFitFunction(0), 
  fMean(0), fSigma(0), fInputFile(0),
  fShowMethod(kTRUE), fShowPtCut(kFALSE), fShowRadius(kFALSE), fShowCent(kFALSE), fInitialized(kFALSE) {;}

DeltaPtContainer::DeltaPtContainer(Int_t cent, EDeltaPtMethod meth, Int_t ptcut, Double_t r) : 
  fCent(cent), fMethod(meth), fConstCut(ptcut), fRadius(r), 
  fDeltaPtHisto(0), fDeltaPtHistoSmallBins(0), fFitFunction(0), 
  fMean(0), fSigma(0), fInputFile(0),
  fShowMethod(kTRUE), fShowPtCut(kFALSE), fShowRadius(kFALSE), fShowCent(kFALSE), fInitialized(kFALSE) {;}

void DeltaPtContainer::Prepare()
{
  if (!fInputFile) {
    TString fname(GetFileName());

    if (fname.IsNull()) return;

    TSeqCollection *flist = gROOT->GetListOfFiles();
    TIter next(flist);
    while(TFile *f = dynamic_cast<TFile*>(next())) {
      TString n(f->GetName());
      if (n.Contains(fname)) {
        fInputFile = f;
        break;
      }
    }
    if (!fInputFile) {
      fInputFile = TFile::Open(Form("$JETRESULTS/%s",fname.Data()));
      if (!fInputFile || fInputFile->IsZombie()) {
        Printf("Unable to open file %s!", fname.Data());
        fInputFile = 0;
        return;
      }
    }
  }

  TString hname(GetHistoName());

  if (hname.IsNull()) return;

  TH2* temp = static_cast<TH2*>(fInputFile->Get(hname));
  if (!temp) {
    Printf("Unable to get histogram %s!", hname.Data());
    return;
  }

  fDeltaPtHistoSmallBins = temp->ProjectionY(Form("%s_smallbins", hname.Data()));
  delete temp;
  temp = 0;

  fDeltaPtHisto = fDeltaPtHistoSmallBins->Rebin(nBins, hname, bins);
  
  fDeltaPtHistoSmallBins->Sumw2();
  fDeltaPtHisto->Sumw2();

  fDeltaPtHisto->Scale(1./fDeltaPtHisto->Integral(), "width");

  fInitialized = kTRUE;

  RecursiveFit();
}

void DeltaPtContainer::RecursiveFit()
{
  if (!fInitialized) return;

  Printf("Now starting the recursive fit for histogram %s.", fDeltaPtHistoSmallBins->GetName());
  
  Double_t oldMean = 0;
  Double_t oldSigma = 0;
  fMean = fDeltaPtHistoSmallBins->GetMean();
  fSigma = fDeltaPtHistoSmallBins->GetRMS();
  Int_t count = 0;
  do {
    fDeltaPtHistoSmallBins->Fit("gaus","0+","", fMean - 3 * fSigma, fMean + 0.5 * fSigma);
    fFitFunction = static_cast<TF1*>(fDeltaPtHistoSmallBins->GetListOfFunctions()->At(0));
    oldMean = fMean;
    oldSigma = fSigma;
    fMean =  fFitFunction->GetParameter("Mean");
    fSigma = fFitFunction->GetParameter("Sigma");
    count++;
  } while ((TMath::Abs(oldMean-fMean) > EPSILON || TMath::Abs(oldSigma-fSigma) > EPSILON) && count < 100);
  if (count == 100) {				
    Printf("WARNING!!! The recursive fitting procedure did not converge after 100 attempts for histogram %s.", fDeltaPtHistoSmallBins->GetName());
  }
  else {
    Printf("The recursive fitting procedure converged after %d attempts for histogram %s.", count, fDeltaPtHistoSmallBins->GetName());
  }
}

void DeltaPtContainer::AddLegendEntry(TLegend* leg)
{
  if (!fInitialized) return;

  TString t;

  if (fShowMethod) {
    t = GetMethodLabel();
    if (fShowPtCut || fShowRadius || fShowCent) t += ", ";
  }

  if (fShowCent) {
    t += GetCentLabel();
    if (fShowPtCut || fShowRadius) t += ", ";
  }

  if (fShowRadius) {
    t += GetRadiusLabel();
    if (fShowPtCut) t += ", ";
  }

  if (fShowPtCut) t += GetPtCutLabel();

  leg->AddEntry(fDeltaPtHisto, t, "pe");
  TLegendEntry *entry0 = leg->AddEntry((TObject*)0, Form("RMS = %4.1f GeV/  #it{c}", fDeltaPtHisto->GetRMS()), "");
  TLegendEntry *entry1 = leg->AddEntry((TObject*)0, Form("#mu^{LHS} = %4.1f GeV/ #it{c}", fMean), "");
  TLegendEntry *entry2 = leg->AddEntry((TObject*)0, Form("#sigma^{LHS} = %4.1f GeV/ #it{c}", fSigma), "");
  entry0->SetTextSize(0.032);
  entry1->SetTextSize(0.032);
  entry2->SetTextSize(0.032);
}

void DeltaPtContainer::Plot(Style_t marker, Color_t color)
{
  if (!fInitialized) return;

  fDeltaPtHisto->SetMarkerStyle(marker);
  fDeltaPtHisto->SetMarkerSize(1.35);
  fDeltaPtHisto->SetMarkerColor(color);
  fDeltaPtHisto->SetLineColor(color);

  fDeltaPtHisto->Draw("same");
}


TString DeltaPtContainer::GetHistoName()
{
  TString hname;

  if (fMethod == kRC) hname = Form("DeltaPtRC_Cent%d",fCent);
  else if (fMethod == kSPE) hname = Form("DeltaPtSPE_Cent%d",fCent);

  return hname;
}

TString DeltaPtContainer::GetFileName()
{
  TString fname;
  if (fRadius == 0.1) {
    if (fConstCut == 0) fname = "LEGO_TRAIN_494_JetType0.root";
  }
  else if (fRadius == 0.2) {
    if (fConstCut == 0) fname = "LEGO_TRAIN_374_JetType0.root";
    else if (fConstCut == 1) fname = "LEGO_TRAIN_489_JetType0.root";
  }
  else if (fRadius == 0.3) {
    if (fConstCut == 0) fname = "LEGO_TRAIN_454_JetType0.root";
    else if (fConstCut == 1) fname = "LEGO_TRAIN_458_JetType0.root";
  }
  return fname;
}

TString DeltaPtContainer::GetPtCutLabel()
{
  TString label;
  if (fConstCut == 0) label = "#it{p}_{T,track(clus)} > 0.15 (0.30) GeV/#it{c}";
  else if (fConstCut == 1) label = "#it{p}_{T,track(clus)} > 1 GeV/#it{c}";

  return label;
}

TString DeltaPtContainer::GetRadiusLabel()
{
  TString label;
  if (fRadius == 0.1) label = "#it{R} = 0.1";
  else if (fRadius == 0.2) label = "#it{R} = 0.2";
  else if (fRadius == 0.3) label = "#it{R} = 0.3";

  return label;
}

TString DeltaPtContainer::GetEtaLabel()
{
  TString label;
  if (fRadius == 0.1) label = "| #eta_{jet}| < 0.6";
  else if (fRadius == 0.2) label = "| #eta_{jet}| < 0.5";
  else if (fRadius == 0.3) label = "| #eta_{jet}| < 0.4";

  return label;
}

TString DeltaPtContainer::GetMethodLabel()
{
  TString label;

  if (fMethod == kRC) label = "Random Cones";
  else if (fMethod == kSPE) label = "Track Embedding";

  return label;
}

TString DeltaPtContainer::GetCentLabel()
{
  TString label;

  if (fCent == 0) label = "0-10% Centrality";
  else if (fCent == 1) label = "10-30% Centrality";

  return label;
}
