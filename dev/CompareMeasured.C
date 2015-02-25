void CompareMeasured()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);

  TFile *fileA = TFile::Open("./Data4/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");
  TFile *fileB = TFile::Open("./Data4/StdDpt_Full_R020_Cent1_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");

  fileA->ls();

  TH1 *histA = (TH1*)fileA->Get("histInclusiveMeasured");
  TH1 *histB = (TH1*)fileB->Get("histInclusiveMeasured");

  histA->Rebin(10);
  histB->Rebin(10);

  histA->Sumw2();
  histB->Sumw2();

  TH1* histeventsA = (TH1*)fileA->Get("histEvents");
  TH1* histeventsB = (TH1*)fileB->Get("histEvents");

  Int_t eventsA = histeventsA->GetEntries();
  Int_t eventsB = histeventsB->GetEntries();

  histA->Scale(1./eventsA);
  histB->Scale(1./eventsB);

  TH1* ratio = (TH1*)histA->Clone("ratio");
  ratio->Divide(histB);

  histA->SetLineColor(kBlack);
  histA->SetMarkerColor(kBlack);
  histA->SetMarkerStyle(kFullCircle);
  histA->SetMarkerSize(0.8);

  histB->SetLineColor(kRed);
  histB->SetMarkerColor(kRed);
  histB->SetMarkerStyle(kFullTriangleUp);
  histB->SetMarkerSize(0.8);

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);

  TCanvas *canvas = new TCanvas();
  canvas->cd();
  canvas->SetLogy();

  canvas->SetLeftMargin(0.12);

  histA->GetXaxis()->SetRangeUser(-40,120);
  histB->GetXaxis()->SetRangeUser(-40,120);
  histA->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN}{dp_{T}} (GeV/c)^{-1}");
  histA->GetYaxis()->SetTitleOffset(1.5);
  histA->GetXaxis()->SetTitleOffset(1.2);

  ratio->GetXaxis()->SetRangeUser(-40,120);
  ratio->GetYaxis()->SetTitle("ratio");

  histA->Draw();
  histB->Draw("same");

  TLegend *leg = new TLegend(0.55,0.60,0.90,0.90);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(13);
  leg->AddEntry(histA, "0-10% centrality", "p");
  leg->AddEntry(histB, "10-30% centrality", "p");
  leg->Draw();

  TCanvas *canvasRatio = new TCanvas();
  canvasRatio->cd();

  ratio->Draw("hist");
}
