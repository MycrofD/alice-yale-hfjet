#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TMath.h>
#include <TPaveText.h>

#include <iostream>

using namespace std;

Int_t bias = 5;
const char* jetType = "Charged";
const char* radius = "R02";

TString LHC10hPath;
TString LHC11hPath;
TString LHC11hFullPath;

TGraphErrors* PlotSpectraReturnGraph(TString, Int_t, Int_t, Bool_t);
void CompareCharged();
void CompareFull();
void CompareFullScale();
TH1* ChargedScale(TH1* spectrum, Double_t s=1./0.8);

void CompareLHC10h()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  LHC11hPath = Form("OldData/StdDpt_Charged_%s0_Cent0_Bias%d_PtCut0/JetUnfoldPbPbOutput.root",radius,bias);
  LHC11hFullPath = Form("Data/StdDpt_Full_%s0_Cent0_Bias%d_PtCut0/JetUnfoldPbPbOutput.root",radius,bias);

  if (bias==0)
    LHC10hPath = Form("LHC10h_results/Inclusive/%s/UnfoldedJetSpectraWithSystematics%s.root",radius,radius);
  else
    LHC10hPath = Form("LHC10h_results/HadronTriggered/HadronTrigger%d/%s/UnfoldedJetSpectraWithSystematics%s.root",bias,radius,radius);

  if (strcmp(jetType,"Charged")==0)
    CompareCharged();
  else
    CompareFullScale();
}

void CompareCharged()
{
  TFile *LHC10hFile = new TFile(LHC10hPath);
  TFile *LHC11hFile = new TFile(LHC11hPath);

  TGraphErrors      *LHC10hNoSyst = (TGraphErrors*)LHC10hFile->Get("grUnfoldedBestCent0");
  //  TGraphAsymmErrors *LHC10hGraph  = (TGraphAsymmErrors*)LHC10hFile->Get("grUnfoldedSystCent0");
  TH1               *LHC11h       = (TH1*)LHC11hFile->Get("histChi2FinalSpectrum");

  TCanvas *c0 = new TCanvas("c0","c0");
  c0->cd();
  c0->SetLogy();

  LHC11h->GetXaxis()->SetRangeUser(20,110);

  LHC11h->SetMarkerColor(kRed+1);
  LHC11h->SetMarkerStyle(kFullCircle);
  LHC11h->SetMarkerSize(1.3);
  LHC11h->SetLineColor(kRed+1);
  LHC11h->Draw("e1");

  gROOT->LoadMacro("LHC10h_results/plotStoredPrelimsInFrame.C");
  TGraphErrors *LHC10hPrint = PlotSpectraReturnGraph(LHC10hPath, 0, 0, kTRUE);

  LHC11h->Draw("e1 same");

  TLegend *leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->AddEntry(LHC10hPrint, "LHC10h (preliminary)", "p");
  leg->AddEntry(LHC11h, "LHC11h", "pl");
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->Draw();

  TString pdfFileName(LHC11hPath);
  pdfFileName.ReplaceAll("JetUnfoldPbPbOutput.root","CompareLHC10h.pdf");

  c0->Update();
  c0->SaveAs(pdfFileName);

  TH1 *LHC10h = static_cast<TH1*>(LHC11h->Clone("LHC10h"));
  LHC10h->Reset();
  Double_t *y = LHC10hNoSyst->GetY();
  Double_t *x = LHC10hNoSyst->GetX();
  Double_t *ey = LHC10hNoSyst->GetEY();
  Double_t *ex = LHC10hNoSyst->GetEX();
  for (Int_t i = 0; i < LHC10hNoSyst->GetN(); i++) {
    Int_t xbin = LHC10h->GetXaxis()->FindBin(x[i]);
    LHC10h->SetBinContent(xbin, y[i]);
    LHC10h->SetBinError(xbin, ey[i]);
  }
  //LHC10h->Scale(1,"width");

  TH1* ratio = (TH1*)LHC11h->Clone("ratio");
  ratio->Divide(LHC10h);
  ratio->GetYaxis()->SetTitle("LHC11h / LHC10h");
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  ratio->Draw();

  TString pdfRatioFileName(LHC11hPath);
  pdfRatioFileName.ReplaceAll("JetUnfoldPbPbOutput.root","CompareLHC10hRatio.pdf");

  c1->Update();
  c1->SaveAs(pdfRatioFileName);
}

void CompareFull()
{

  cout << "Comparing: " << LHC10hPath << ", " << LHC11hFullPath << endl;

  TFile *LHC10hFile = new TFile(LHC10hPath);
  TFile *LHC11hFile = new TFile(LHC11hFullPath);

  TGraphErrors      *LHC10hNoSyst = (TGraphErrors*)LHC10hFile->Get("grUnfoldedBestCent0");
  //TGraphAsymmErrors *LHC10hGraph  = (TGraphAsymmErrors*)LHC10hFile->Get("grUnfoldedSystCent0");
  TH1               *LHC11h       = (TH1*)LHC11hFile->Get("histSVDFinalSpectrum");

  TH1* ratio = (TH1*)LHC11h->Clone("ratio");

  TH1 *LHC10h = static_cast<TH1*>(LHC11h->Clone("LHC10h"));
  LHC10h->Reset();

  TCanvas *c0 = new TCanvas("c0","c0",600,600);
  c0->cd();
  c0->SetLeftMargin(0.16);
  c0->SetLogy();

  LHC11h->GetYaxis()->SetRangeUser(1e-10,1e-4);

  LHC11h->GetYaxis()->SetTitleOffset(1.8);

  LHC11h->GetXaxis()->SetTitleOffset(1.2);

  Double_t sys = 0.21;

  for (Int_t i = 1; i <= LHC11h->GetNbinsX(); i++) {
    LHC11h->SetBinError(i,LHC11h->GetBinContent(i) * sys);
  }

  LHC11h->SetMarkerColor(kRed+1);
  LHC11h->SetMarkerStyle(kFullCircle);
  LHC11h->SetMarkerSize(1.1);
  LHC11h->SetLineColor(kWhite);
  LHC11h->SetLineStyle(0);
  LHC11h->SetFillColor(kRed-10);
  LHC11h->GetXaxis()->SetRangeUser(30,109.9);
  LHC11h->DrawCopy("e3");

  gROOT->LoadMacro("LHC10h_results/plotStoredPrelimsInFrame.C");
  TGraphErrors *LHC10hPrint = PlotSpectraReturnGraph(LHC10hPath, 0, 0, kTRUE);

  LHC11h->DrawCopy("hist p same");

  TLegend *leg = new TLegend(0.5, 0.68, 0.88, 0.57);
  leg->AddEntry(LHC10hPrint, "Charged jets (QM2012)", "p");
  leg->AddEntry(LHC11h, "Full jets", "pf");
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->Draw();

  TPaveText *text = new TPaveText(0.4,1.0,0.17,0.73,"blNDC");
  text->AddText("0-10% Centrality");
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextSize(0.03);
  text->SetTextAlign(13);
  text->Draw();

  TPaveText *text2 = new TPaveText(0.4,0.9,0.17,0.73,"NDC");
  text2->AddText("Pb-Pb #sqrt{s_{NN}}=2.76 TeV");
  text2->SetBorderSize(0);
  text2->SetFillStyle(0);
  text2->SetTextSize(0.03);
  text2->SetTextAlign(13);
  text2->Draw();

  TPaveText *text3a = new TPaveText(0.92,0.90,0.5,0.84,"NDC");
  text3a->AddText("Anti-k_{T} R=0.2");
  text3a->SetBorderSize(0);
  text3a->SetFillStyle(0);
  text3a->SetTextSize(0.033);
  text3a->SetTextAlign(13);
  text3a->Draw();

  TPaveText *text3 = new TPaveText(0.92,0.84,0.5,0.8,"NDC");
  text3->AddText("p_{T}^{track}   > 0.15 GeV/c");
  text3->SetBorderSize(0);
  text3->SetFillStyle(0);
  text3->SetTextSize(0.028);
  text3->SetTextAlign(13);
  text3->Draw();

  TPaveText *text4 = new TPaveText(0.92,0.8,0.5,0.76,"NDC");
  text4->AddText("p_{T}^{cluster} > 0.30 GeV/c");
  text4->SetBorderSize(0);
  text4->SetFillStyle(0);
  text4->SetTextSize(0.028);
  text4->SetTextAlign(13);
  text4->Draw();
 
  TPaveText *text3b = new TPaveText(0.92,0.76,0.5,0.7,"NDC");
  text3b->AddText("Leading hadron p_{T} > 5 GeV/c");
  text3b->SetBorderSize(0);
  text3b->SetFillStyle(0);
  text3b->SetTextSize(0.028);
  text3b->SetTextAlign(13);
  text3b->Draw();

  c0->Update();
  //c0->SaveAs(Form("%s/Full_%s%s/CompareLHC10h.pdf", path1, radius, LHC11h_paths[type%2]));
  
  Double_t *y = LHC10hNoSyst->GetY();
  Double_t *x = LHC10hNoSyst->GetX();
  Double_t *ey = LHC10hNoSyst->GetEY();
  Double_t *ex = LHC10hNoSyst->GetEX();
  for (Int_t i = 0; i < LHC10hNoSyst->GetN(); i++) {
    Int_t xbin = LHC10h->GetXaxis()->FindBin(x[i]);
    LHC10h->SetBinContent(xbin, LHC10h->GetBinContent(xbin)+y[i]*10);
    LHC10h->SetBinError(xbin, ey[i]);
  }
  LHC10h->Scale(1,"width");
  
  ratio->Divide(LHC10h);
  ratio->GetYaxis()->SetTitle("Full / Charged");
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  ratio->Draw();

  c1->Update();
  //c1->SaveAs(Form("%s/Full_%s%s/CompareLHC10hRatio.pdf", path1, radius, LHC11h_paths[type%2]));
}

void CompareFullScale()
{
  //const char* PythiaName     = Form("PythiaSpectra/Pythia%s.root", LHC11h_paths[type%2]);

  cout << "Comparing: " << LHC10hPath << ", " << LHC11hPath << ", " << LHC11hFullPath << endl;

  TFile *LHC10hFile     = new TFile(LHC10hPath);
  TFile *LHC11hFile     = new TFile(LHC11hPath);
  TFile *LHC11hFullFile = new TFile(LHC11hFullPath);
  //TFile *PythiaFile     = new TFile(PythiaName);

  TGraphErrors      *LHC10hNoSyst = (TGraphErrors*)LHC10hFile->Get("grUnfoldedBestCent0");
  TGraphAsymmErrors *LHC10hGraph  = (TGraphAsymmErrors*)LHC10hFile->Get("grUnfoldedSystCent0");
  TH1               *LHC11h       = 0;
  if (LHC11hFile) 
    LHC11h = (TH1*)LHC11hFile->Get("histSVDFinalSpectrum");
  TH1               *LHC11hFull   = (TH1*)LHC11hFullFile->Get("histSVDFinalSpectrum");
  //TH1               *PythiaRatio  = (TH1*)PythiaFile->Get(Form("Ratio%s",radius));

  TH1 *LHC10h = (TH1*)LHC11hFull->Clone("LHC10h");
  Double_t *y = LHC10hNoSyst->GetY();
  Double_t *x = LHC10hNoSyst->GetX();
  Double_t *ey = LHC10hNoSyst->GetEY();
  for (Int_t i = 0; i < LHC10hNoSyst->GetN(); i++) {
    Int_t xbin = LHC10h->GetXaxis()->FindBin(x[i]);
    LHC10h->SetBinContent(xbin, y[i]);
    LHC10h->SetBinError(xbin, ey[i]);
  }

  
  TH1* LHC10hScaled = ChargedScale(LHC10h);
  LHC10hScaled->SetName("LHC10hScaled");
  LHC10hScaled->SetTitle("LHC10hScaled");

  TH1* LHC11hScaled = 0;

  if (LHC11h) {
    LHC11hScaled = ChargedScale(LHC11h);
    LHC11hScaled->SetName("LHC11hScaled");
    LHC11hScaled->SetTitle("LHC11hScaled");
  }

  /*
  TH1* ratioLHC11h = (TH1*)LHC11hFull->Clone("ratioLHC11h");
  ratioLHC11h->Divide(LHC11h);

  TH1* ratioLHC10h = (TH1*)LHC11hFull->Clone("ratioLHC10h");
  ratioLHC10h->Divide(LHC10h);
  */
  TCanvas *c0 = new TCanvas("c0","c0");
  c0->cd();
  c0->SetLogy();

  LHC11hFull->GetXaxis()->SetRangeUser(30,109.9);
  LHC11hFull->SetMarkerColor(kRed+1);
  LHC11hFull->SetMarkerStyle(kFullCircle);
  LHC11hFull->SetMarkerSize(1.3);
  LHC11hFull->SetLineColor(kRed+1);
  LHC11hFull->Draw("e1");
  LHC11hFull->GetYaxis()->SetRangeUser(1e-10,1e-6);

  if (LHC11hScaled) {
    LHC11hScaled->GetXaxis()->SetRangeUser(30,109.9);
    LHC11hScaled->SetMarkerColor(kRed+1);
    LHC11hScaled->SetMarkerStyle(kOpenCircle);
    LHC11hScaled->SetMarkerSize(1.3);
    LHC11hScaled->SetLineColor(kRed+1);
    LHC11hScaled->Draw("same e1");
  }

  LHC10hScaled->GetXaxis()->SetRangeUser(30,109.9);
  LHC10hScaled->SetMarkerColor(kBlue+1);
  LHC10hScaled->SetMarkerStyle(kOpenCircle);
  LHC10hScaled->SetMarkerSize(1.3);
  LHC10hScaled->SetLineColor(kBlue+1);
  LHC10hScaled->Draw("same e1");
  
  TLegend *leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->AddEntry(LHC11hScaled, "LHC11h charged-jets (scaled)", "pl");
  leg->AddEntry(LHC10hScaled, "LHC10h charged-jets (scaled)", "pl");
  leg->AddEntry(LHC11hFull, "LHC11h full jets", "pl");
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->Draw();

  c0->Update();
  // c0->SaveAs(Form("%s/Full_Meth2_%s%s/CompareLHC10h.pdf", path1, radius, LHC11h_paths[type%2]));
  /*
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  ratioLHC11h->GetYaxis()->SetTitle("Full / Charged");
  ratioLHC11h->GetXaxis()->SetRangeUser(30,109.9);
  ratioLHC11h->SetMarkerColor(kRed+1);
  ratioLHC11h->SetLineColor(kRed+1);
  ratioLHC11h->SetMarkerSize(0.6);
  ratioLHC11h->GetYaxis()->SetRangeUser(2,14);
  ratioLHC11h->Draw();

  ratioLHC10h->GetYaxis()->SetTitle("Full / Charged");
  ratioLHC10h->GetXaxis()->SetRangeUser(30,109.9);
  ratioLHC10h->SetMarkerColor(kBlue+1);
  ratioLHC10h->SetLineColor(kBlue+1);
  ratioLHC10h->SetMarkerSize(0.6);
  ratioLHC10h->GetYaxis()->SetRangeUser(2,14);
  ratioLHC10h->Draw("same");

  /*
  PythiaRatio->GetYaxis()->SetTitle("Full / Charged");
  PythiaRatio->GetXaxis()->SetRangeUser(30,109.9);
  PythiaRatio->SetMarkerColor(kBlack);
  PythiaRatio->SetMarkerSize(0.6);
  PythiaRatio->SetLineColor(kBlack);
  PythiaRatio->GetYaxis()->SetRangeUser(2,14);
  PythiaRatio->Draw("same");
  

  TLegend *legRatio = new TLegend(0.15, 0.67, 0.53, 0.88);
  legRatio->AddEntry(ratioLHC10h, "Full jets / Charged jets (LHC10h)", "pl");
  legRatio->AddEntry(ratioLHC11h, "Full jets / Charged jets (LHC11h)", "pl");
  //legRatio->AddEntry(PythiaRatio, "Full jets / Charged jets (Pythia)", "pl");
  legRatio->SetLineColor(kWhite);
  legRatio->SetFillColor(kWhite);
  legRatio->Draw();

  c1->Update();
  //c1->SaveAs(Form("%s/Full_Meth2_%s%s/CompareLHC10hRatio.pdf", path1, radius, LHC11h_paths[type%2]));
  */
}

TH1* ChargedScale(TH1* spectrum, Double_t s)
{
  TArrayD bins(*(spectrum->GetXaxis()->GetXbins()));
  for (Int_t i = 0; i < bins.GetSize(); i++) {
    bins[i] *= s;
  }

  TH1* result = new TH1D("result","result",bins.GetSize()-1,bins.GetArray());

  for (Int_t i = 1; i <= bins.GetSize(); i++) {
    result->SetBinContent(i, spectrum->GetBinContent(i) / s);
    result->SetBinError(i, spectrum->GetBinError(i) / s);
  }

  return result;
}
