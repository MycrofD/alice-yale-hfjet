#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TList.h>
#include <Riostream.h>
#include <TString.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

Bool_t save = kTRUE;
const char* format = "pdf";

void ComparePtHardBins()
{
  TH1::AddDirectory(kFALSE);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //TFile *file0 = TFile::Open("ThermalModelData_PtHardBin_1-1_Bias5.root");
  //TFile *file1 = TFile::Open("ThermalModelData_PtHardBin_2-10_Bias5.root");

  TFile *file0 = TFile::Open("UnfoldingTestData_bias5_Pt1_2.root");
  TFile *file1 = TFile::Open("UnfoldingTestData_bias5_Pt3_10.root");
  
  TH1* histo0 = static_cast<TH1*>(file0->Get("histInclusiveMeasured"));
  TH1* histo1 = static_cast<TH1*>(file1->Get("histInclusiveMeasured"));
  
  histo0->Rebin(10);
  histo1->Rebin(10);

  histo0->SetLineColor(kRed+2);
  histo0->SetLineWidth(2);

  histo1->SetLineColor(kBlue+2);
  histo1->SetLineWidth(2);

  TH1* sum = static_cast<TH1*>(histo0->Clone("sum"));
  sum->Add(histo1);

  sum->SetLineColor(kBlack);

  TCanvas *canvas = new TCanvas();
  canvas->cd();
  canvas->SetLogy();

  TH1* myBlankHisto = new TH1F("myBlankHisto","myBlankHisto",1000,-250,250);
  myBlankHisto->GetXaxis()->SetRangeUser(-30,120);
  myBlankHisto->GetYaxis()->SetRangeUser(1e-9,1e0);
  myBlankHisto->GetXaxis()->SetTitle("#it{p}_{T,jet}^{det} (GeV/#it{c})");
  myBlankHisto->GetYaxis()->SetTitle("Arb. units");
  myBlankHisto->Draw();

  histo0->Draw("hist same");
  histo1->Draw("hist same");
  sum->Draw("hist same");

  TLegend* leg = new TLegend(0.568966,0.589852,0.863506,0.860465);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(histo0, "5 < #it{p}_{T,hard} < 21 GeV/#it{c}", "l");
  leg->AddEntry(histo1, "#it{p}_{T,hard} > 21 GeV/#it{c}", "l");
  leg->AddEntry(sum, "#it{p}_{T,hard} > 5 GeV/#it{c} (sum)", "l");
  leg->Draw();

  TH1* ratio0 = static_cast<TH1*>(histo0->Clone("ratio"));
  ratio0->Divide(sum);

  TH1* ratio1 = static_cast<TH1*>(histo1->Clone("ratio"));
  ratio1->Divide(sum);

  TCanvas *canvasR = new TCanvas();
  canvasR->cd();
  canvasR->SetGridx();
  canvasR->SetGridy();

  TH1* myBlankHistoR = new TH1F("myBlankHistoR","myBlankHistoR",1000,-250,250);
  myBlankHistoR->GetXaxis()->SetRangeUser(-30,120);
  myBlankHistoR->GetYaxis()->SetRangeUser(0,1.1);
  myBlankHistoR->GetXaxis()->SetTitle("#it{p}_{T,jet}^{det} (GeV/#it{c})");
  myBlankHistoR->GetYaxis()->SetTitle("#frac{#it{N}_{jets}^{#it{p}_{T,hard} <> X GeV/#it{c}}}{#it{N}_{jets}^{#it{p}_{T,hard} > 5 GeV/#it{c}}}");
  myBlankHistoR->Draw();

  ratio0->Draw("hist same");
  ratio1->Draw("hist same");

  if (save) {
    canvas->SaveAs(Form("ComparePtHardBins.%s",format));
    canvasR->SaveAs(Form("ComparePtHardBinsRatio.%s",format));
  }
}
