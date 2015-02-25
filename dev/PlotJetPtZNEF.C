#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFile.h>
#include <TList.h>
#include <TLegend.h>
#include <Riostream.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>

const char* spectrumName = "JetsCorrPtVs%s_Bias%d_Cent%d";

const char* varLabel = "#it{z}_{lead}";
const char* var = "Z";

//const char* varLabel = "NEF";
//const char* var = "NEF";

Int_t bias = 0;
Int_t cent = 0;

Int_t kJetType = 0;

Bool_t doSave = kTRUE;
Bool_t drawLabels = kTRUE;
const char* fileFormat = "pdf";

const Int_t nbins = 26;
Double_t    bins[nbins+1] = {-250,-50,-40,-30,-20,-10,-5,0,5,10,15,20,25,30,35,40,45,50,55,60,70,80,100,120,150,200,250};

#define R02

#ifdef R02

const char* train = "374"; // R=0.2, charged leading hadron
//const char* train = "488"; // R=0.2, charged+neutral leading hadron
//const char* train = "491"; // R=0.2, neutral leading hadron
const char* radiusLabel = "Anti-#it{k}_{T}  #it{R} = 0.2";
const char* etaLabel = "|#eta| < 0.5";
Double_t minX = 20;
Double_t maxX = 150;
const char* fileLabel = "R020";
const char* titleLabel = "";

#else 
#ifdef R03

const char* train = "454"; // R=0.3
const char* radiusLabel = "Anti-#it{k}_{T}  #it{R} = 0.3";
const char* etaLabel = "|#eta| < 0.4";
Double_t minX = -60;
Double_t maxX = 120;
const char* fileLabel = "R030";
const char* titleLabel = "";

#else
#ifdef R03_1GEV

const char* train = "458"; // R=0.3
const char* radiusLabel = "#it{R} = 0.3";
const char* etaLabel = "|#eta| < 0.4";
Double_t minX = -60;
Double_t maxX = 140;
const char* fileLabel = "R030_1GeV";
const char* titleLabel = "Det. level const. cut = 1 GeV/#it{c}";

#endif
#endif
#endif

const Int_t NCENT = 2;

Style_t markers[11] = { kFullCircle, kFullSquare, kFullTriangleUp, kFullCircle, kFullCircle, kFullCircle, kFullCircle, kFullCircle, kFullCircle, kFullCircle, kFullCircle };
Color_t colors[11] = { kBlack, kBlue+2, kRed+2, kGreen+2, kMagenta+2, kOrange+2, kCyan+2, kYellow+2, kAzure+2, kViolet+2, kTeal+2 };

void PlotJetPtZNEF()
{
  TH1::AddDirectory(kFALSE);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TString hname(Form(spectrumName, var, bias, cent));
  
  TString fileName(Form("$JETRESULTS/LEGO_TRAIN_%s_JetType%d.root",train,kJetType));

  TFile *file = TFile::Open(fileName);
  TH2* hist = static_cast<TH2*>(file->Get(hname));

  Double_t minPt[3] = {20, 50, 80};
  Double_t maxPt[3] = {30, 60, 140};
  TH1* histos[3] = {0};

  TString canvasName(Form("Jet%sBias%dCent%d",var,bias,cent));

  TCanvas *c = new TCanvas(canvasName, canvasName, 700, 700);
  c->SetLeftMargin(0.15);
  c->SetLogy();
  c->cd();

  TH1* myBlankHisto = new TH1C("myBlankHisto","myBlankHisto",1000,0,1.0);
  myBlankHisto->GetXaxis()->SetRangeUser(0,1.0);
  myBlankHisto->GetXaxis()->SetNdivisions(505,kFALSE);
  myBlankHisto->GetYaxis()->SetRangeUser(1e-2,50);
  myBlankHisto->GetXaxis()->SetTitle(varLabel);
  myBlankHisto->GetYaxis()->SetTitle("Probability density");
  myBlankHisto->GetYaxis()->SetTitleSize(0.045);
  myBlankHisto->GetXaxis()->SetTitleSize(0.045);
  myBlankHisto->GetYaxis()->SetLabelSize(0.050);
  myBlankHisto->GetXaxis()->SetLabelSize(0.050);
  myBlankHisto->GetYaxis()->SetTitleOffset(1.5);
  myBlankHisto->Draw("hist");

  TLegend *leg1 = new TLegend(0.562,0.449,0.961,0.810);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetNColumns(1);
  leg1->SetTextAlign(12);
  leg1->SetTextFont(43);
  leg1->SetTextSize(19);

  for (Int_t i = 0; i < 3; i++) {
    histos[i] = hist->ProjectionY(Form("Jet%sBias%dCent%dPt%.0f",var,bias,cent,(minPt[i]+maxPt[i])/2), hist->GetXaxis()->FindBin(minPt[i]), hist->GetXaxis()->FindBin(maxPt[i]));
    histos[i]->Rebin(10);
    histos[i]->Sumw2();
    histos[i]->Scale(1./histos[i]->Integral(), "width");

    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetMarkerStyle(markers[i]);
    histos[i]->SetMarkerSize(1.2);
    histos[i]->Draw("same");

    leg1->AddEntry(histos[i], Form("%.0f < #it{p}_{T,jet}^{det} < %.0f GeV/#it{c}", minPt[i], maxPt[i]), "pe");
    leg1->AddEntry((TObject*)0, Form("RMS = %.2f", histos[i]->GetRMS()), "");
    leg1->AddEntry((TObject*)0, Form("Mean = %.2f", histos[i]->GetMean()), "");
  }

  if (drawLabels) {
    TPaveText *alice0 = new TPaveText(0.20,0.82,0.90,0.94,"brNDC");
    alice0->SetBorderSize(0);
    alice0->SetFillStyle(0);
    alice0->SetTextAlign(13);
    alice0->SetTextFont(43);
    alice0->SetTextSize(23);
    alice0->AddText("ALICE  Preliminary  Pb-Pb   #sqrt{#it{s}_{NN}} = 2.76 TeV");
    alice0->Draw();

    //TPaveText *description0 = new TPaveText(0.180,0.686,0.409,0.825,"brNDC");
    //TPaveText *description0 = new TPaveText(0.647,0.263,0.876,0.403,"brNDC");
    TPaveText *description0 = new TPaveText(0.335,0.686,0.563,0.828,"brNDC");
    description0->SetBorderSize(0);
    description0->SetFillStyle(0);
    description0->SetTextAlign(12);
    description0->SetTextFont(43);
    description0->SetTextSize(19);
    if (cent==0) description0->AddText("0 - 10% Centrality");
    else description0->AddText("10 - 30% Centrality");
    description0->AddText(radiusLabel);
    description0->AddText(etaLabel);
    description0->AddText("Detector level");
    description0->Draw();
  }

  leg1->Draw();

  if (doSave) c->SaveAs(Form("%s.%s",canvasName.Data(), fileFormat));
}
