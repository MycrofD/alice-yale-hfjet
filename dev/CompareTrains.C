#include <TFile.h>
#include <TList.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>

#include <iostream>

using namespace std;

TCanvas *canvas=0;
TCanvas *canvasRatio=0;
TH1* hist1=0;
TH1* hist2=0;
TH1* ratio=0;

TH1* DoCompareTrains(const char* train1, const char* train2, const char* name1, const char* name2, const char* histName="histBayesianFinalSpectrum", const char* opt="", const char* optRatio="");
TH1* Rebin(const char* newname, TH1* histo, Int_t nbinsX, const Double_t *binsX, Int_t keepNorm=0);

Bool_t save = kFALSE;

void CompareBias()
{
  canvasRatio = new TCanvas("canvasRatio_Unbiased-Biased","canvasRatio_Unbiased-Biased");

  canvas=0;
  DoCompareTrains("Data/StdDpt_Full_R020_Cent0_Bias0_PtCut0", "Data/StdDpt_Full_R020_Cent0_Bias5_PtCut0", "Unbiased", "Biased","histSVDFinalSpectrum","","hist");
  canvas->Update();
  return;
  canvas=0;
  DoCompareTrains("Data/StdDpt_Full_R030_Cent0_Bias0_PtCut0", "Data/StdDpt_Full_R030_Cent0_Bias5_PtCut0", "Unbiased", "Biased","histSVDFinalSpectrum","","hist");
  canvas->Update();


  if (save)
    canvas->SaveAs("ComparisonPlots/BiasComparisonBayesian_153-154.pdf");
  TH1 *bayesBias = ratio;
  bayesBias->SetLineColor(kBlue+1);
  canvas=0;
  DoCompareTrains("LT153-154/Full_R020_Unbiased", "LT153-154/Full_R020_Biased", "Unbiased", "Biased","histChi2FinalSpectrum","","same hist");
  canvas->Update();
  if (save)
    canvas->SaveAs("ComparisonPlots/BiasComparisonSVD_153-154.pdf");
  TH1 *svdBias = ratio;
  svdBias->SetLineColor(kRed+1);
  canvas=0;
  DoCompareTrains("LT153-154/Full_R020_Unbiased", "LT153-154/Full_R020_Biased", "Unbiased", "Biased","histChi2FinalSpectrum","","same hist");
  canvas->Update();
  if (save)
    canvas->SaveAs("ComparisonPlots/BiasComparisonChi2_153-154.pdf");
  TH1 *chi2Bias = ratio;
  chi2Bias->SetLineColor(kGreen+2);

  TFile *lhc10h = TFile::Open("LHC10h_results/RatiosBiasedUnbiasedR02.root");
  TH1* lhc10Ratio = (TH1*)lhc10h->Get("grRatioHT0HT5PbPbSyst");
  lhc10Ratio->SetLineColor(kBlack);
  lhc10Ratio->SetMarkerColor(kBlack);
  lhc10Ratio->SetMarkerStyle(kFullCircle);
  lhc10Ratio->SetMarkerSize(0.7);
  lhc10Ratio->Draw("same P");

  TLegend *legBias = new TLegend(0.5, 0.67, 0.88, 0.88);
  legBias->SetLineColor(kWhite);
  legBias->SetFillColor(kWhite);
  legBias->AddEntry(bayesBias, "Bayesian", "l");
  legBias->AddEntry(svdBias, "SVD", "l");
  legBias->AddEntry(chi2Bias, "Chi2", "l");
  legBias->AddEntry(lhc10Ratio, "LHC10h (charged)", "pl");
  legBias->Draw();

  canvasRatio->Update();
  if (save)
    canvasRatio->SaveAs("ComparisonPlots/BiasComparisonRatios_153-154.pdf");
}

void CompareRadii()
{
  canvasRatio = new TCanvas("canvasRatio_R02-R03","canvasRatio_R02-R03");

  canvas=0;
  DoCompareTrains("Data/StdDpt_Full_R030_Cent0_Bias5_PtCut0", "Data/StdDpt_Full_R020_Cent0_Bias5_PtCut0", "R=0.3", "R=0.2","histSVDFinalSpectrum","","hist");
  canvas->Update();

  return;
  if (save)
    canvas->SaveAs("ComparisonPlots/RadiusComparisonBayesian_153-154.pdf");
  TH1 *bayesRadius = ratio;
  bayesRadius->SetLineColor(kBlue+1);
  canvas=0;
  DoCompareTrains("LT153-154/Full_R030_Biased", "LT153-154/Full_R020_Biased", "R=0.3", "R=0.2","histSVDFinalSpectrum","","same hist");
  canvas->Update();
  if (save)
    canvas->SaveAs("ComparisonPlots/RadiusComparisonBayesian_153-154.pdf");
  TH1 *svdRadius = ratio;
  svdRadius->SetLineColor(kRed+1);
  canvas=0;
  DoCompareTrains("LT153-154/Full_R030_Biased", "LT153-154/Full_R020_Biased", "R=0.3", "R=0.2","histChi2FinalSpectrum","","same hist");
  canvas->Update();
  if (save)
    canvas->SaveAs("ComparisonPlots/RadiusComparisonBayesian_153-154.pdf");
  TH1 *chi2Radius = ratio;
  chi2Radius->SetLineColor(kGreen+2);

  TFile *lhc10h = TFile::Open("LHC10h_results/LHC10hANTIKTR02R03RatioBiased5.root");
  TH1* lhc10Ratio = (TH1*)lhc10h->Get("grRatioSystCent0");
  lhc10Ratio->Draw("same P");

  TLegend *legRadius = new TLegend(0.5, 0.67, 0.88, 0.88);
  legRadius->SetLineColor(kWhite);
  legRadius->SetFillColor(kWhite);
  legRadius->AddEntry(bayesRadius, "Bayesian", "l");
  legRadius->AddEntry(svdRadius, "SVD", "l");
  legRadius->AddEntry(chi2Radius, "Chi2", "l");
  legRadius->AddEntry(lhc10Ratio, "LHC10h (charged)", "l");
  legRadius->Draw();

  canvasRatio->Update();
  if (save)
    canvasRatio->SaveAs("ComparisonPlots/RadiusComparisonRatios_153-154.pdf");
}

void CompareTrains()
{
  CompareBias();

  CompareRadii();
 /*
  DoCompareTrains("LT153-154/Full_R020_Biased", "LT153-154/Full_R020_Biased_pTcut1000", "#it{p}_{T,const}^{min}=0.15 GeV/#it{c} (0.30 GeV)", "#it{p}_{T,const}^{min}=1 GeV/#it{c}","histBayesianFinalSpectrum");

  DoCompareTrains("LT153-154/Full_R020_Unbiased", "LT153-154/Full_R020_Unbiased_pTcut1000", "#it{p}_{T,const}^{min}=0.15 GeV/#it{c} (0.30 GeV)", "#it{p}_{T,const}^{min}=1 GeV/#it{c}","histBayesianFinalSpectrum");

  DoCompareTrains("LT153-154/Full_R020_Biased", "LT153-154/Full_R020_Biased_pTcut1000", "#it{p}_{T,const}^{min}=0.15 GeV/#it{c} (0.30 GeV)", "#it{p}_{T,const}^{min}=1 GeV/#it{c}","histSVDFinalSpectrum");

  DoCompareTrains("LT153-154/Full_R020_Unbiased", "LT153-154/Full_R020_Unbiased_pTcut1000", "#it{p}_{T,const}^{min}=0.15 GeV/#it{c} (0.30 GeV)", "#it{p}_{T,const}^{min}=1 GeV/#it{c}","histSVDFinalSpectrum");
 
  DoCompareTrains("LT128-129/Full_R02_Biased", "LT128-129/Full_R02_Biased_Meth1", "Meth.2", "Meth.1");

  DoCompareTrains("LT128-129/Full_R02_Biased", "LT128-129/Full_R02_Biased_Eff10", "5 GeV/c, Eff-5", "5 GeV/c, Eff-10");
  DoCompareTrains("LT128-129/Full_R02_Biased10", "LT128-129/Full_R02_Biased10_Eff10", "10 GeV/c, Eff-5", "10 GeV/c, Eff-10");

  DoCompareTrains("LT128-129/Full_R02_Biased", "LT128-129/Full_R02_Biased_Embedding", "5 GeV/c, random cones", "5 GeV/c, embedding");
  DoCompareTrains("LT128-129/Full_R02_Biased10", "LT128-129/Full_R02_Biased10_Embedding", "10 GeV/c, random cones", "10 GeV/c, embedding");

  DoCompareTrains("LT128-129/Full_R02_Biased", "LT131-132/Full_R02_Biased_HadCorr17", "5 GeV/c, HadCorr=2", "5 GeV/c, HadCorr=1.7");
  DoCompareTrains("LT128-129/Full_R02_Biased10", "LT131-132/Full_R02_Biased10_HadCorr17", "10 GeV/c, HadCorr=2", "10 GeV/c, HadCorr=1.7");

  DoCompareTrains("LT128-129/Full_R02_Biased", "LT135-136/Full_R02_Biased_HadCorr13", "5 GeV/c, HadCorr=2", "5 GeV/c, HadCorr=1.3");
  DoCompareTrains("LT128-129/Full_R02_Biased10", "LT135-136/Full_R02_Biased10_HadCorr13", "10 GeV/c, HadCorr=2", "10 GeV/c, HadCorr=1.3");

  
  TH1 *ratio1 = DoCompareTrains("LT128-129/Full_R02_Biased", "LT128-129/Full_R02_Biased10", "R=0.2, 5 GeV/c hadron triggered, default", "R=0.2, 10 GeV/c hadron triggered, default");
  ratio1->SetName("RatioDefault");
  
  TH1 *ratio2 = DoCompareTrains("LT131-132/Full_R02_Biased_HadCorr17", "LT131-132/Full_R02_Biased10_HadCorr17", "R=0.2, 5 GeV/c hadron triggered, had corr=1.7", "R=0.2, 10 GeV/c hadron triggered, had corr=1.7");
  ratio2->SetName("RatioHadCorr17");

  TH1 *ratio3 = DoCompareTrains("LT135-136/Full_R02_Biased_HadCorr13", "LT135-136/Full_R02_Biased10_HadCorr13", "R=0.2, 5 GeV/c hadron triggered, had corr=1.3", "R=0.2, 10 GeV/c hadron triggered, had corr=1.3");
  ratio3->SetName("RatioHadCorr13");

  TCanvas *c = new TCanvas("CompareRatiosHadCorr","CompareRatiosHadCorr");
  c->cd();
  
  ratio1->SetTitle("Default");
  ratio2->SetTitle("HadCorr = 1.7");
  ratio3->SetTitle("HadCorr = 1.3");

  ratio1->SetLineColor(kRed+1);
  ratio2->SetLineColor(kBlue+1);
  ratio3->SetLineColor(kGreen+2);

  ratio1->SetMarkerColor(kRed+1);
  ratio2->SetMarkerColor(kBlue+1);
  ratio3->SetMarkerColor(kGreen+2);

  ratio1->Draw();  
  ratio2->Draw("same");
  ratio3->Draw("same");

  TLegend *leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->AddEntry(ratio1, ratio1->GetTitle(), "pl");
  leg->AddEntry(ratio2, ratio2->GetTitle(), "pl");
  leg->AddEntry(ratio3, ratio3->GetTitle(), "pl");
  leg->Draw();

  TH1 *ratioRatio17 = static_cast<TH1*>(ratio2->Clone("ratioRatio17"));
  ratioRatio17->Divide(ratio1);

  TH1 *ratioRatio13 = static_cast<TH1*>(ratio3->Clone("ratioRatio13"));
  ratioRatio13->Divide(ratio1);
  
  TCanvas *c2 = new TCanvas("CompareRatiosRatioHadCorr","CompareRatiosRatioHadCorr");
  c2->cd();

  ratioRatio13->Draw();
  ratioRatio17->Draw("same");

  TH1 *ratio4 = DoCompareTrains("LT128-129/Full_R02_Biased_Eff10", "LT128-129/Full_R02_Biased10_Eff10", "R=0.2, 5 GeV/c hadron triggered, eff=-10", "R=0.2, 10 GeV/c hadron triggered, eff=-10");
  ratio4->SetName("RatioEff10");

  TCanvas *c4 = new TCanvas("CompareRatiosEff","CompareRatiosEff");
  c4->cd();

  ratio4->SetTitle("Eff -10%");

  ratio4->SetLineColor(kBlue+1);
  ratio4->SetMarkerColor(kBlue+2);

  ratio1->Draw();  
  ratio4->Draw("same");

  TLegend *leg2 = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg2->SetLineColor(kWhite);
  leg2->SetFillColor(kWhite);
  leg2->AddEntry(ratio1, ratio1->GetTitle(), "pl");
  leg2->AddEntry(ratio4, ratio4->GetTitle(), "pl");
  leg2->Draw();

  TH1 *ratioRatioEff = static_cast<TH1*>(ratio4->Clone("ratioRatioEff"));
  ratioRatioEff->Divide(ratio1);
  
  TCanvas *c5 = new TCanvas("CompareRatiosRatioEff","CompareRatiosRatioEff");
  c5->cd();

  ratioRatioEff->Draw();

  TH1 *ratio5 = DoCompareTrains("LT128-129/Full_R02_Biased_Embedding", "LT128-129/Full_R02_Biased10_Embedding", "R=0.2, 5 GeV/c hadron triggered, embedding", "R=0.2, 10 GeV/c hadron triggered, embedding");
  ratio5->SetName("RatioEmb");

  TCanvas *c7 = new TCanvas("CompareRatiosEmb","CompareRatiosEmb");
  c7->cd();

  ratio5->SetTitle("Embedding");

  ratio5->SetLineColor(kBlue+1);
  ratio5->SetMarkerColor(kBlue+2);

  ratio1->Draw();  
  ratio5->Draw("same");

  TLegend *leg3 = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg3->SetLineColor(kWhite);
  leg3->SetFillColor(kWhite);
  leg3->AddEntry(ratio1, ratio1->GetTitle(), "pl");
  leg3->AddEntry(ratio5, ratio5->GetTitle(), "pl");
  leg3->Draw();

  TH1 *ratioRatioEmb = static_cast<TH1*>(ratio5->Clone("ratioRatioEmb"));
  ratioRatioEmb->Divide(ratio1);
  
  TCanvas *c8 = new TCanvas("CompareRatiosRatioEmb","CompareRatiosRatioEmb");
  c8->cd();

  ratioRatioEmb->Draw();

  
  TH1 *ratioStat = DoCompareTrains("LT128-129/Full_R02_Biased_Set1", "LT128-129/Full_R02_Biased10_Set2", "Set 1 R=0.2, 5 GeV/c hadron triggered, default", "Set 2 R=0.2, 10 GeV/c hadron triggered, default");
  for (Int_t i = ratioStat->GetXaxis()->FindBin(31); i < ratioStat->GetXaxis()->FindBin(31)+9; i++) {
    cout << ratioStat->GetBinError(i) << ", ";
  }
  cout << endl;
    */
}

TH1* DoCompareTrains(const char* train1, const char* train2, const char* name1, const char* name2, const char* histName, const char* opt, const char* optRatio)
{
  gStyle->SetOptStat(1101);
  gStyle->SetOptTitle(0);

  const char *train1FileName = Form("%s/JetUnfoldPbPbOutput.root",train1);
  const char *train2FileName = Form("%s/JetUnfoldPbPbOutput.root",train2);
  
  TFile *train1File = new TFile(train1FileName);
  TFile *train2File = new TFile(train2FileName);

  TH1* train1Hist = (TH1*)train1File->Get(histName);
  if (!train1Hist) {
    cout << "Unable to find histogram " <<  histName << " in file " << train1File->GetName() << endl;
    return 0;
  }
  train1Hist->SetName(name1);
  train1Hist->GetXaxis()->SetRangeUser(30,119.9);
  train1Hist->SetStats(kFALSE);

  TH1* temp = (TH1*)train2File->Get(histName);
  if (!temp) {
    cout << "Unable to find histogram " <<  histName << " in file " << train2File->GetName() << endl;
    return 0;
  }

  TH1* train2Hist = Rebin(name2,temp,train1Hist->GetNbinsX(),train1Hist->GetXaxis()->GetXbins()->GetArray(),1) ;
  train2Hist->GetXaxis()->SetRangeUser(30,119.9);
  train2Hist->SetStats(kFALSE);

  TString train2safe(train2);
  train2safe.ReplaceAll("/","-");

  if (!canvas)
    canvas = new TCanvas(Form("canvas%s-%s-%s",histName,train1,train2),Form("canvas%s-%s",train1Hist->GetName(),train2Hist->GetName()));
  canvas->cd();
  canvas->SetLogy();
  train1Hist->SetLineColor(kBlue+1);
  train1Hist->SetMarkerColor(kBlue+1);
  train1Hist->SetMarkerSize(0.8);
  train1Hist->SetMarkerStyle(kFullTriangleDown);
  train1Hist->Draw(opt);

  train2Hist->SetLineColor(kRed+1);
  train2Hist->SetMarkerColor(kRed+1);
  train2Hist->SetMarkerSize(0.8);
  train2Hist->SetMarkerStyle(kFullTriangleUp);
  train2Hist->Draw("same");

  Double_t min = TMath::Min(train1Hist->GetMinimum(), train2Hist->GetMinimum());
  Double_t max = TMath::Max(train1Hist->GetMaximum(), train2Hist->GetMaximum());

  train1Hist->SetMinimum(min/4);
  train1Hist->SetMaximum(max*4);

  TLegend *leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetLineColor(kWhite);
  leg->SetFillColor(kWhite);
  leg->AddEntry(train1Hist, train1Hist->GetName(), "pl");
  leg->AddEntry(train2Hist, train2Hist->GetName(), "pl");
  leg->Draw();

  Double_t x1 = 0.13;
  Double_t x2 = 0.4;
  Double_t y1 = 0.87;
  Double_t y2 = 0.64;
  TPaveText *pave = new TPaveText(x1,y1,x2,y2,"blNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextSize(0.035);
  pave->SetTextAlign(13);
  pave->AddText(histName);
  pave->Draw();

  if (!canvasRatio) 
    canvasRatio = new TCanvas(Form("canvasRatio%s-%s-%s",histName,train1,train2),Form("canvasRatio%s-%s",train1Hist->GetName(),train2Hist->GetName()));
  canvasRatio->cd();

  ratio = (TH1*)train2Hist->Clone("ratio");
  ratio->SetStats(kFALSE);
  ratio->Divide(train1Hist);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerSize(0.6);
  ratio->SetMarkerStyle(kFullCircle);
  ratio->GetYaxis()->SetTitle(Form("%s / %s", train2Hist->GetName(), train1Hist->GetName()));
  ratio->GetYaxis()->SetRangeUser(0.4,1.8);
  ratio->Draw(optRatio);

  //pave->Draw();

  return ratio;
}

TH1* Rebin(const char* newname, TH1* histo, Int_t nbinsX, const Double_t *binsX, Int_t keepNorm)
{
  TH1* newhisto = new TH1D(newname,newname,nbinsX,binsX);

  for (Int_t x = 0; x <= histo->GetNbinsX(); x++) {
    Double_t xlow    = histo->GetXaxis()->GetBinLowEdge(x)+0.0001;
    Double_t xup     = histo->GetXaxis()->GetBinUpEdge(x)-0.0001;
    Int_t    xlowBin = newhisto->GetXaxis()->FindBin(xlow);
    Int_t    xupBin  = newhisto->GetXaxis()->FindBin(xup);
    if (xlowBin != xupBin) {
      std::cout <<  "Warning, unable to exact rebin!" << std::endl;
    }
    Double_t content = histo->GetBinContent(x);
    Double_t err     = histo->GetBinError(x);
    if (keepNorm) {
      content *= xup - xlow;
      err     *= xup - xlow;
      newhisto->SetBinError(xlowBin, newhisto->GetBinError(xlowBin)+err);
    }
    newhisto->SetBinContent(xlowBin, newhisto->GetBinContent(xlowBin)+content);
  }

  if (keepNorm) {
    for (Int_t x = 0; x <= newhisto->GetNbinsX(); x++) {
      newhisto->SetBinContent(x, newhisto->GetBinContent(x) / (newhisto->GetXaxis()->GetBinUpEdge(x) - newhisto->GetXaxis()->GetBinLowEdge(x)));
      newhisto->SetBinError(x, newhisto->GetBinError(x) / (newhisto->GetXaxis()->GetBinUpEdge(x) - newhisto->GetXaxis()->GetBinLowEdge(x)));
    }
  }
  else {
    newhisto->SetEntries(histo->GetEntries());
  }
  
  return newhisto;
}
