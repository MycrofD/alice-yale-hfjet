void CompareMatrices()
{
  TH1::AddDirectory(kFALSE);

  TFile *file1 = TFile::Open("PythiaToy/PythiaOwnDpt_A_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");
  TFile *file2 = TFile::Open("PythiaToy/PythiaDpt_A_Full_R020_Cent0_Bias5_PtCut0_HadCorr20/JetUnfoldPbPbOutput.root");

  TH2* matrix1 = (TH2*)file1->Get("histDetectorResponseTruth_0");
  TH2* matrix2 = (TH2*)file2->Get("histDetectorResponseTruth_0");

  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1);
  
  c1->cd(1);
  gPad->SetLogz();
  matrix1->Draw("colz");

  c1->cd(2);
  gPad->SetLogz();
  matrix2->Draw("colz");

  TH2* diff = (TH2*)matrix1->Clone("diff");
  diff->Add(matrix2,-1);
  diff->Divide(matrix1);
  
  TCanvas *c2 = new TCanvas();
  c2->cd();
  diff->Draw("colz");
}
