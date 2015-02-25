void CompareEfficiency()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TFile *file1 = TFile::Open("Data5/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20_MaxTrackPt1000/JetUnfoldPbPbOutput.root");
  TFile *file2 = TFile::Open("Data5/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  TList *list1 = (TList*)file1->Get("ResponseMatrices");
  TList *list2 = (TList*)file2->Get("ResponseMatrices");

  TH1 *eff1 = (TH1*)list1->FindObject("histDetectorEfficiency");
  TH1 *eff2 = (TH1*)list2->FindObject("histDetectorEfficiency");

  TCanvas *c = new TCanvas();
  c->cd();

  eff1->GetXaxis()->SetRangeUser(0,120);
  eff1->GetYaxis()->SetRangeUser(0.5,1.2);

  eff1->SetLineColor(kBlack);
  eff2->SetLineColor(kRed);

  eff1->Draw("hist");
  eff2->Draw("hist same");

  TLegend *leg  = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->AddEntry(eff1,"#it{p}_{T,charged}^{leading} < 50 GeV/c","pl");
  leg->AddEntry(eff2,"#it{p}_{T,charged}^{leading} < 100 GeV/c","pl");
  leg->Draw();
}
