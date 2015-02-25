TH1* Rebin(const char* newname, const TH1* histo, Int_t nbinsX, const Double_t *binsX, Int_t keepNorm);
TH1* Rebin(const char* newname, const TH1* histo, const TH1* templatehist, Int_t keepNorm);

Bool_t save = kFALSE;

void Compare()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //TFile *newfile = TFile::Open("./Data/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data/PythiaDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./ThermalToy/StdDpt_C_Full_R020_Bias0_PtHardTrigger/JetUnfoldOutput.root");
  //TFile *oldfile = TFile::Open("./ThermalToy/PythiaDpt_C_Full_R020_Bias0_PtHardTrigger/JetUnfoldOutput.root");

  //TFile *newfile = TFile::Open("./Data3_FullStatistics/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldOutputBinConf2.root");
  //TFile *oldfile = TFile::Open("./Data3_FullStatistics/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldOutputBinConf2.root");

  //TFile *newfile = TFile::Open("./Data4/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");
  //TFile *oldfile = TFile::Open("./Data4/StdDpt_Full_R020_Cent1_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");

  //TFile *newfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent1_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./Data4/StdDpt_Full_R020_Cent1_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data4/SPEDpt_Full_R020_Cent1_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = newfile;

  //TFile *oldfile = TFile::Open("./Data3_FullStatistics/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldOutputBinConf2.root");
  //TFile *newfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./Data3/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data4/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./Data5/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20_MaxTrackPt50/JetUnfoldPbPbInput.root");
  //TFile *oldfile = TFile::Open("./Data5/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbInput.root");

  //TFile *newfile = TFile::Open("./Data7bis/StdDpt_Full_R020_Cent0_Bias0_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  //TFile *oldfile = TFile::Open("./Data10/StdDpt_Full_R020_Cent0_Bias0_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  TFile *newfile = TFile::Open("./Data7quinquies/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  TFile *oldfile = TFile::Open("./Data7/StdDpt_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  //TFile *newfile = TFile::Open("./UnfoldingTestData_bias0_sat_test.root");
  //TFile *oldfile = TFile::Open("./UnfoldingTestData_bias0_test.root");

  TList* newlist = (TList*)newfile->Get("Spectra");
  TList* oldlist = (TList*)oldfile->Get("Spectra");

  //TH1* newhist = (TH1*)newlist->FindObject("histInclusiveSVD_Reg5_PriorInterTruthNonFit0Unfolded");
  //TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveSVD_Reg5_PriorInterTruthNonFit0Unfolded");

  //TH1* newhist = (TH1*)newlist->FindObject("histInclusiveSVD_Reg5_PriorMeasured0Unfolded");
  //TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveChi2_Reg2.00E-01_PriorMeasured0Unfolded");
  //TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveSVD_Reg5_PriorMeasured0Unfolded");
  //TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveSVD_Reg5_PriorMeasured-1Unfolded");

  TH1* newhist_temp = (TH1*)newlist->FindObject("histInclusiveSVD_Reg6_PriorMeasured0Unfolded");
  TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveSVD_Reg5_PriorMeasured0Unfolded");
  TH1* newhist = 0;
  newhist = Rebin("newhist", newhist_temp, oldhist,1);

  //TH1* newhist = (TH1*)newlist->FindObject("histInclusiveMeasured");
  //TH1* oldhist = (TH1*)oldlist->FindObject("histInclusiveMeasured");

  //TH1* newhist = (TH1*)newfile->Get("histInclusiveMeasured");
  //TH1* oldhist = (TH1*)oldfile->Get("histInclusiveMeasured");

  TH1* newhistevents = (TH1*)newfile->Get("histEvents");
  TH1* oldhistevents = (TH1*)oldfile->Get("histEvents");

  Int_t newevents = newhistevents->GetEntries();
  Int_t oldevents = oldhistevents->GetEntries();

  TCanvas *c = new TCanvas();
  c->cd();
  c->SetLogy();

  newhist->Scale(1./newevents);
  oldhist->Scale(1./oldevents);
  
  //newhist->Rebin(10);
  //oldhist->Rebin(10);

  //newhist->Sumw2();
  //oldhist->Sumw2();

  //newhist->GetYaxis()->SetTitle("N_{jets}");

  Double_t min = -50;
  Double_t max = 139.99;

  newhist->GetXaxis()->SetRangeUser(min,max);
  newhist->SetLineColor(kBlack);
  newhist->SetMarkerColor(kBlack);
  newhist->SetMarkerStyle(kFullCircle);

  oldhist->GetXaxis()->SetRangeUser(min,max);
  oldhist->SetLineColor(kRed);
  oldhist->SetMarkerColor(kRed);
  oldhist->SetMarkerStyle(kFullTriangleUp);

  newhist->Draw();
  oldhist->Draw("same");

  TLegend *leg  = new TLegend(0.5, 0.67, 0.88, 0.88);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  //leg->AddEntry(newhist,"#it{p}_{T,charged}^{leading} < 50 GeV/c","pl");
  //leg->AddEntry(oldhist,"#it{p}_{T,charged}^{leading} < 100 GeV/c","pl");
  //leg->AddEntry(newhist,"\"Mirror\" #delta#it{p}_{T}","pl");
  //leg->AddEntry(oldhist,"Nominal","pl");
  leg->AddEntry(newhist,"New","pl");
  leg->AddEntry(oldhist,"Old","pl");
  leg->Draw();

  TH1* ratio = (TH1*)newhist->Clone("ratio");
  ratio->GetXaxis()->SetRangeUser(min,max);
  ratio->Divide(oldhist);
  ratio->SetLineColor(kBlue+1);
  ratio->SetLineWidth(3);
  ratio->GetYaxis()->SetTitle("Ratio New / Old");

  TCanvas *c1 = new TCanvas();
  c1->cd();
  ratio->Draw("hist");

  TString plotName("NewTrackingEfficiency");

  TString spectraName(plotName);
  spectraName += ".pdf";

  TString ratioName(plotName);
  ratioName += "Ratio.pdf";

  if (save) {
    c->SaveAs(spectraName);
    c1->SaveAs(ratioName);
  }
}


//____________________________________________
TH1* Rebin(const char* newname, const TH1* histo, const TH1* templatehist, Int_t keepNorm)
{
  Int_t nbins = templatehist->GetNbinsX();
  Double_t *bins = new Double_t[nbins+1];
  for (Int_t i = 0; i <= nbins; i++) 
    bins[i] = templatehist->GetXaxis()->GetBinUpEdge(i);

  return Rebin(newname, histo, nbins, bins, keepNorm);
}

//____________________________________________
TH1* Rebin(const char* newname, const TH1* histo, Int_t nbinsX, const Double_t *binsX, Int_t keepNorm)
{
  TH1* newhisto = new TH1D(newname,newname,nbinsX,binsX);

  for (Int_t x = 0; x <= histo->GetNbinsX(); x++) {
    Double_t xlow    = histo->GetXaxis()->GetBinLowEdge(x);
    Double_t xup     = histo->GetXaxis()->GetBinUpEdge(x);
    Int_t    xlowBin = newhisto->GetXaxis()->FindBin(xlow+0.0001);
    Int_t    xupBin  = newhisto->GetXaxis()->FindBin(xup-0.0001);
    if (xlowBin != xupBin) {
      Printf("Warning, unable to exact rebin! Old name: %s, new name: %s", histo->GetName(), newname);
    }
    //Printf("Input bin: %.2f - %.2f --> Output bin: %.2f - %.2f", histo->GetXaxis()->GetBinLowEdge(x), histo->GetXaxis()->GetBinUpEdge(x), 
    //	   newhisto->GetXaxis()->GetBinLowEdge(xlowBin), newhisto->GetXaxis()->GetBinUpEdge(xlowBin));
    Double_t content = histo->GetBinContent(x);
    Double_t err     = histo->GetBinError(x);
    if (keepNorm) {
      content *= xup - xlow;
      err     *= xup - xlow;
    }
    newhisto->SetBinContent(xlowBin, newhisto->GetBinContent(xlowBin)+content);
    newhisto->SetBinError(xlowBin, TMath::Sqrt(newhisto->GetBinError(xlowBin)*newhisto->GetBinError(xlowBin)+err*err));
    //Printf("New content %f",newhisto->GetBinContent(xlowBin));
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
