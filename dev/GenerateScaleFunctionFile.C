#include <TFile.h>
#include <TF1.h>

void GenerateScaleFunctionFile(const char* fname = "LHC11h_ScaleFactorFunctions.root")
{
  TFile *file = TFile::Open(fname, "recreate");
  file->cd();

  TF1* sfunc = 0;

  // used for train n. 374
  sfunc = new TF1("LHC11h_HadCorr20_ClustersV2_old", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.76458);
  sfunc->SetParameter(1, -0.0111656);
  sfunc->SetParameter(0, 0.000107296);
  sfunc->Write();

  // fit from full train n. 374
  sfunc = new TF1("LHC11h_HadCorr20_ClustersV2", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.78067);
  sfunc->SetParameter(1, -0.0116032);
  sfunc->SetParameter(0, 0.000114951);
  sfunc->Write();

  // fit from test run of train n. 387
  sfunc = new TF1("LHC11h_HadCorr17_ClustersV2", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.8171);
  sfunc->SetParameter(1, -0.0150141);
  sfunc->SetParameter(0, 0.000193964);
  sfunc->Write();

  // fit from test run of train n. 389
  sfunc = new TF1("LHC11h_HadCorr13_ClustersV2", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.90496);
  sfunc->SetParameter(1, -0.016978);
  sfunc->SetParameter(0, 0.000215038);
  sfunc->Write();

  // fit from full train n. 401
  sfunc = new TF1("LHC11h_HadCorr00_ClustersV2", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 2.0418);
  sfunc->SetParameter(1, -0.0157904);
  sfunc->SetParameter(0, 0.000148058);
  sfunc->Write();

  // used for train n. 408 (fit test train n.408)
  sfunc = new TF1("LHC11h_HadCorr20_Clusters3x3_old", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.85043);
  sfunc->SetParameter(1, -0.0116461);
  sfunc->SetParameter(0, 9.65308e-05);
  sfunc->Write();

  // fit from full train n.408
  sfunc = new TF1("LHC11h_HadCorr20_Clusters3x3", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.8528);
  sfunc->SetParameter(1, -0.0121015);
  sfunc->SetParameter(0, 0.000109308);
  sfunc->Write();

  // fit from full train n. 407
  sfunc = new TF1("LHC11h_HadCorr00_Clusters3x3", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 2.13631);
  sfunc->SetParameter(1, -0.0165245);
  sfunc->SetParameter(0, 0.000142767);
  sfunc->Write();

  // fit from test train n. 457
  sfunc = new TF1("LHC11h_HadCorr20_ClustersV2_ConstCut1GeV", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.63498);
  sfunc->SetParameter(1, -0.0208189);
  sfunc->SetParameter(0, 0.000303435);
  sfunc->Write();

  // GA
  sfunc = new TF1("LHC11h_HadCorr20_ClustersV2_GA", "[0]*x*x+[1]*x+[2]", -1, 100);
  sfunc->SetParameter(2, 1.81208);
  sfunc->SetParameter(1, -0.0105506);
  sfunc->SetParameter(0, 0.000145519);
  sfunc->Write();

  file->Close();
  delete file;
  file = 0;
}
