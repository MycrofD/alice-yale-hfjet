void CompareFoldDeltaPt()
{
  TH1::AddDirectory(kFALSE);

  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetLegendBorderSize(0);

  TFile *RCfile = TFile::Open("./Data5/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");

  TH1* reco = (TH1*)RCfile->Get("histDetectorMCReconstructed");
  reco->Rebin(10);

  TH1* RCdeltaPt = (TH1*)RCfile->Get("histBackgroundResponse");
  RCdeltaPt->SetName("RCdeltaPt");
  RCdeltaPt->Rebin(10);

  TFile *SPEfile = TFile::Open("./Data5/SPEDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");

  TH1* SPEdeltaPt = (TH1*)SPEfile->Get("histBackgroundResponse");
  SPEdeltaPt->SetName("SPEdeltaPt");
  SPEdeltaPt->Rebin(10);

  TH1* RCfold = FoldDeltaPt(reco, RCdeltaPt);
  RCfold->SetName("RCfold");
  TH1* SPEfold = FoldDeltaPt(reco, SPEdeltaPt);
  SPEfold->SetName("SPEfold");

  TH1* ratio = (TH1*)RCfold->Clone("ratio");
  ratio->Divide(SPEfold);
  ratio->GetYaxis()->SetTitle("RC/SPE");
  
  RCfold->SetLineColor(kRed);
  RCfold->SetMarkerColor(kRed);
  RCfold->SetMarkerStyle(kFullSquare);
  RCfold->SetMarkerSize(0.8);

  SPEfold->SetLineColor(kBlue);
  SPEfold->SetMarkerColor(kBlue);
  SPEfold->SetMarkerStyle(kFullCircle);
  SPEfold->SetMarkerSize(0.8);

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);

  TCanvas *c = new TCanvas();
  c->cd();
  c->SetLogy();
  RCfold->Draw("hist");
  SPEfold->Draw("same hist");

  TLegend *leg = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->AddEntry(RCfold, "Random Cones");
  leg->AddEntry(SPEfold, "Single Particle Embedding");

  TCanvas *cr = new TCanvas();
  cr->cd();
  cr->SetGridy();
  cr->SetGridx();
  ratio->Draw("hist");
}

TH1* FoldDeltaPt(TH1* input, TH1* deltaPt)
{
  deltaPt->Sumw2();
  deltaPt->Scale(1./deltaPt->GetEntries());

  TH1* folded = new TH1D("folded","folded;p_{T} (GeV/c);counts",100,-250,250);
  
  Printf("Nbins input = %d", input->GetNbinsX());  
  Printf("Nbins deltaPt = %d", deltaPt->GetNbinsX());

  for (Int_t i = 1; i <= input->GetNbinsX(); i++) {
    for (Int_t j = 1; j <= deltaPt->GetNbinsX(); j++) {
      Int_t ifolded = i + j - 1;
      if (ifolded < 1 || ifolded > 100) continue;
      folded->SetBinContent(ifolded, folded->GetBinContent(ifolded) + input->GetBinContent(i) * deltaPt->GetBinContent(j));
      if (ifolded == 1) {
	Printf("ifolded = %d, i = %d, j = %d",ifolded,i,j);
      }
      Double_t prevErr = folded->GetBinError(ifolded)*folded->GetBinError(ifolded);
      Double_t newRelErr = (input->GetBinError(i)*input->GetBinError(i)) / (input->GetBinContent(i)*input->GetBinContent(i)) + 
	(deltaPt->GetBinError(j)*deltaPt->GetBinError(j)) / (deltaPt->GetBinContent(j)*deltaPt->GetBinContent(j));
      Double_t newErr = prevErr + newRelErr * input->GetBinContent(i) * deltaPt->GetBinContent(j);
      folded->SetBinError(ifolded, TMath::Sqrt(newErr));
    }
  }

  return folded;
}
