#include <TFile.h>
#include <TList.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <Riostream.h>
#include <TSystem.h>

const char *rootPath = "$JETRESULTS";
const char *fileName = "LEGO_TRAIN_MC_370/MergedResults-169838.root";
const char *listName = "AliAnalysisTaskSAJF_Jet_AKTFullR020_MCParticlesSelected_pT0000_pt_scheme_EMCAL_histos";
const char *histName = "fHistTracksJetPt_0";

void MultiplyPt(TH1* hist);
void SetHistoStyle(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t color);
void Plot(const char* yaxisTitle, Double_t maxY, TH1* lo, TH1* mi, TH1* hi);

void PlotJetConstituents()
{
  TH1::AddDirectory(kFALSE);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TString fname(Form("%s/%s", rootPath, fileName));
  gSystem->ExpandPathName(fname);

  TFile *file = TFile::Open(fname);
  if (!file || file->IsZombie()) {
    Printf("Could not open file %s!", fname.Data());
  }

  TList *list = dynamic_cast<TList*>(file->Get(listName));
  if (!list) {
    Printf("Could not find list %s!", listName);
    return;
  }

  TH2 *fHistTracksJetPt = dynamic_cast<TH2*>(list->FindObject(histName));
  if (!fHistTracksJetPt) {
    Printf("Could not find histogram %s!", histName);
    return;
  }

  TH2 *fHistTracksJetPt_clone = dynamic_cast<TH2*>(fHistTracksJetPt->Clone("fHistTracksJetPt_clone"));

  delete list;
  list = 0;
  fHistTracksJetPt = 0;
  
  file->Close();
  delete file;
  file = 0;
  
  TH1* contituentsLoPt = fHistTracksJetPt_clone->ProjectionX("contituentsLoPt", 
                                                             fHistTracksJetPt_clone->GetYaxis()->FindBin(30.0001), fHistTracksJetPt_clone->GetYaxis()->FindBin(39.9999));
  TH1* contituentsMiPt = fHistTracksJetPt_clone->ProjectionX("contituentsMiPt", 
                                                             fHistTracksJetPt_clone->GetYaxis()->FindBin(60.0001), fHistTracksJetPt_clone->GetYaxis()->FindBin(69.9999));
  TH1* contituentsHiPt = fHistTracksJetPt_clone->ProjectionX("contituentsHiPt", 
                                                             fHistTracksJetPt_clone->GetYaxis()->FindBin(100.0001), fHistTracksJetPt_clone->GetYaxis()->FindBin(109.9999));

  SetHistoStyle(contituentsLoPt, kFullCircle, 1, kRed+2);
  SetHistoStyle(contituentsMiPt, kFullSquare, 1, kBlue+2);
  SetHistoStyle(contituentsHiPt, kFullTriangleUp, 1, kGreen+2);

  TH1* contituentsLoPtTimesPt = static_cast<TH1*>(contituentsLoPt->Clone("contituentsLoPtTimesPt"));
  TH1* contituentsMiPtTimesPt = static_cast<TH1*>(contituentsMiPt->Clone("contituentsMiPtTimesPt"));
  TH1* contituentsHiPtTimesPt = static_cast<TH1*>(contituentsHiPt->Clone("contituentsHiPtTimesPt"));

  contituentsLoPtTimesPt->Scale(1./contituentsLoPt->Integral(), "");
  contituentsMiPtTimesPt->Scale(1./contituentsMiPt->Integral(), "");
  contituentsHiPtTimesPt->Scale(1./contituentsHiPt->Integral(), "");

  MultiplyPt(contituentsLoPtTimesPt);
  MultiplyPt(contituentsMiPtTimesPt);
  MultiplyPt(contituentsHiPtTimesPt);

  contituentsLoPt->Rebin(2);
  contituentsMiPt->Rebin(2);
  contituentsHiPt->Rebin(2);

  contituentsLoPtTimesPt->Rebin(2);
  contituentsMiPtTimesPt->Rebin(2);
  contituentsHiPtTimesPt->Rebin(2);

  contituentsLoPtTimesPt->Scale(1., "width");
  contituentsMiPtTimesPt->Scale(1., "width");
  contituentsHiPtTimesPt->Scale(1., "width");

  contituentsLoPt->Scale(1./contituentsLoPt->Integral(), "width");
  contituentsMiPt->Scale(1./contituentsMiPt->Integral(), "width");
  contituentsHiPt->Scale(1./contituentsHiPt->Integral(), "width");

  Plot("Probability density", 1., contituentsLoPt,contituentsMiPt,contituentsHiPt);
  Plot("(Probability density)   #times  #it{p}_{T}^{const} (GeV/#it{c})", 1., contituentsLoPtTimesPt,contituentsMiPtTimesPt,contituentsHiPtTimesPt);
}

void MultiplyPt(TH1* hist)
{
  if (!hist) return;
  for (Int_t i = 1; i <= hist->GetNbinsX(); i++) {
    hist->SetBinError(i, hist->GetBinError(i)*hist->GetXaxis()->GetBinCenter(i));    
    hist->SetBinContent(i, hist->GetBinContent(i)*hist->GetXaxis()->GetBinCenter(i));
  }
}

void SetHistoStyle(TH1* hist, Style_t markerStyle, Size_t markerSize, Color_t color) 
{
  if (!hist) return;
  hist->SetMarkerStyle(markerStyle);
  hist->SetMarkerSize(markerSize);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
}

void Plot(const char* yaxisTitle, Double_t maxY, TH1* lo, TH1* mi, TH1* hi)
{
  TCanvas *canvas = new TCanvas();
  canvas->SetLogy();
  canvas->cd();
  
  TH1* myBlankHist = new TH1F("myBlankHist","myBlankHist", 1500, 0, 150);
  myBlankHist->GetXaxis()->SetRangeUser(0., 100);
  myBlankHist->GetXaxis()->SetTitle("#it{p}_{T}^{const} (GeV/#it{c})");  
  myBlankHist->GetYaxis()->SetRangeUser(1e-4, maxY);
  myBlankHist->GetYaxis()->SetTitle(yaxisTitle);
  myBlankHist->Draw();

  lo->Draw("same");
  mi->Draw("same");
  hi->Draw("same");
}
