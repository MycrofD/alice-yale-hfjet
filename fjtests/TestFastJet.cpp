#include <iostream>
#include <vector>
#include <cmath>

#include <fastjet/config.h>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/AreaDefinition.hh>

#include "TestFastJet.h"

int main()
{
  std::cout << "Test for FastJet\n";

  FJWrapper fj1("fj1");
  LoadEvent(fj1);
  fj1.FindJets();

  FJWrapper fj2("fj2");
  LoadEvent(fj2);
  fj2.FindJets();

  fj1.PrintJets();
  fj1.CompareResults(fj2);
}

FJWrapper::~FJWrapper()
{
  if (seqArea) delete seqArea;
}


FJWrapper::FJWrapper(const char* _name) :
  name(_name),
  strategy(fastjet::Best),
  algor(fastjet::antikt_algorithm),
  scheme(fastjet::pt_scheme),
  areaType(fastjet::active_area_explicit_ghosts),
  nGhostRepeats(1),
  ghostArea(0.005),
  maxRap(1.),
  R(0.2),
  gridScatter(1.0),
  ktScatter(0.1),
  meanGhostKt(1e-100),
  seqArea(0),
  inputVectors(),
  incJets()
{

}

void FJWrapper::CompareResults(const FJWrapper& fjw)
{
  printf("Printing jets that are different in %s compared to %s (this)\n", fjw.name.c_str(), name.c_str());
  printf("Found %lu jets in %s, compared to %lu in %s (this)\n", fjw.incJets.size(), fjw.name.c_str(), incJets.size(), name.c_str());
  for (unsigned int i = 0, j = 0; i < incJets.size() && j < fjw.incJets.size(); ) {

    while (incJets[i].perp() < 1e-3 && i < incJets.size()) i++;
    if (i == incJets.size()) break;

    while (fjw.incJets[j].perp() < 1e-3 && j < fjw.incJets.size()) j++;
    if (j == fjw.incJets.size()) break;

    if (std::abs(incJets[i].px() - fjw.incJets[j].px()) > 1e-3 ||
        std::abs(incJets[i].py() - fjw.incJets[j].py()) > 1e-3 ||
        std::abs(incJets[i].pz() - fjw.incJets[j].pz()) > 1e-3) {
      printf("i=%u, px=%.3f, py=%.3f, pz=%.3f, m=%.3f, a=%.3f\n", i, fjw.incJets[i].px(), fjw.incJets[i].py(), fjw.incJets[i].pz(), fjw.incJets[i].m(), fjw.incJets[i].area());
    }
    i++;
    j++;
  }
}

void FJWrapper::FindJets()
{
  fastjet::GhostedAreaSpec ghostedAreaSpec(maxRap, nGhostRepeats, ghostArea, gridScatter, ktScatter, meanGhostKt);
  fastjet::AreaDefinition areaDef(ghostedAreaSpec, areaType);
  fastjet::JetDefinition jetDef(algor, R, scheme, strategy);
  seqArea = new fastjet::ClusterSequenceArea(inputVectors, jetDef, areaDef);
  incJets = seqArea->inclusive_jets(0.0);
}

void FJWrapper::AddInputVector(double px, double py, double pz, double m)
{
  fastjet::PseudoJet inVec(px, py, pz, std::sqrt(m*m+px*px+py*py+pz*pz));
  inputVectors.push_back(inVec);
}

void FJWrapper::PrintJets() const
{
  std::cout << name << std::endl;
  std::cout << "Found " << incJets.size() << " jets" << std::endl;
  unsigned int i = 0;
  for (std::vector<fastjet::PseudoJet>::const_iterator it = incJets.begin(); it != incJets.end(); it++) {
    if ((*it).perp() < 1e-3) continue;
    printf("i=%u, px=%.3f, py=%.3f, pz=%.3f, m=%.3f, a=%.3f\n", i, (*it).px(), (*it).py(), (*it).pz(), (*it).m(), (*it).area());
    i++;
  }
}

void LoadEvent(FJWrapper& fj)
{
  fj.AddInputVector(-1.718,-0.196,1.201,0.139);
  fj.AddInputVector(0.390,-1.113,-0.852,0.139);
  fj.AddInputVector(-1.597,-0.183,1.429,0.139);
  fj.AddInputVector(0.974,-1.116,0.489,0.139);
  fj.AddInputVector(-1.376,-0.199,1.046,0.139);
  fj.AddInputVector(0.972,1.579,1.854,0.139);
  fj.AddInputVector(-1.008,1.203,-1.365,0.139);
  fj.AddInputVector(-1.429,-0.487,0.587,0.139);
  fj.AddInputVector(0.632,0.918,0.474,0.139);
  fj.AddInputVector(1.985,0.393,1.904,0.139);
  fj.AddInputVector(-0.781,0.662,0.954,0.139);
  fj.AddInputVector(0.551,-1.230,-0.788,0.139);
  fj.AddInputVector(-1.385,1.668,1.536,0.139);
  fj.AddInputVector(-0.542,-0.949,-0.529,0.139);
  fj.AddInputVector(1.261,0.034,0.761,0.139);
  fj.AddInputVector(0.347,2.216,1.538,0.139);
  fj.AddInputVector(-1.270,0.160,1.020,0.139);
  fj.AddInputVector(-1.481,-0.876,-1.001,0.139);
  fj.AddInputVector(-1.208,-0.845,-0.969,0.139);
  fj.AddInputVector(-0.810,-1.282,-1.228,0.139);
  fj.AddInputVector(-1.724,2.193,0.530,0.139);
  fj.AddInputVector(-0.939,0.823,0.243,0.139);
  fj.AddInputVector(-1.159,1.677,0.418,0.139);
  fj.AddInputVector(2.140,0.483,1.807,0.139);
  fj.AddInputVector(2.452,0.545,2.103,0.139);
  fj.AddInputVector(1.450,0.153,0.584,0.139);
  fj.AddInputVector(0.766,0.697,0.352,0.139);
  fj.AddInputVector(1.269,0.207,1.007,0.139);
  fj.AddInputVector(-1.478,0.175,0.318,0.139);
  fj.AddInputVector(0.906,-0.946,0.968,0.139);
  fj.AddInputVector(0.825,1.085,0.877,0.139);
  fj.AddInputVector(-0.825,1.319,0.227,0.139);
  fj.AddInputVector(-1.318,-0.374,0.659,0.139);
  fj.AddInputVector(-0.653,1.225,0.319,0.139);
  fj.AddInputVector(0.656,1.171,1.242,0.139);
  fj.AddInputVector(1.005,-0.819,0.916,0.139);
  fj.AddInputVector(-0.797,-1.528,-1.553,0.139);
  fj.AddInputVector(0.811,1.986,0.472,0.139);
  fj.AddInputVector(-1.912,-0.127,0.579,0.139);
  fj.AddInputVector(-2.107,0.268,1.864,0.139);
  fj.AddInputVector(-0.877,-0.512,0.943,0.139);
  fj.AddInputVector(-1.114,-0.637,-1.027,0.139);
  fj.AddInputVector(1.311,-1.127,-0.539,0.139);
  fj.AddInputVector(1.271,1.530,1.097,0.139);
  fj.AddInputVector(-1.018,1.363,0.460,0.139);
  fj.AddInputVector(-0.564,-1.270,-0.995,0.139);
  fj.AddInputVector(-0.443,-1.001,-0.908,0.139);
  fj.AddInputVector(-1.263,0.741,-0.745,0.139);
  fj.AddInputVector(1.099,-1.380,-1.228,0.139);
  fj.AddInputVector(-1.136,-0.543,-0.399,0.139);
  fj.AddInputVector(0.972,-1.356,1.150,0.139);
  fj.AddInputVector(1.983,-0.044,-1.448,0.139);
  fj.AddInputVector(1.062,0.952,0.094,0.139);
  fj.AddInputVector(-0.045,1.486,0.562,0.139);
  fj.AddInputVector(-1.247,-0.047,-0.646,0.139);
  fj.AddInputVector(1.543,0.593,-1.057,0.139);
  fj.AddInputVector(0.064,1.040,-0.492,0.139);
  fj.AddInputVector(-0.809,-0.644,-0.018,0.139);
  fj.AddInputVector(0.846,1.029,0.070,0.139);
  fj.AddInputVector(0.585,1.133,-0.953,0.139);
  fj.AddInputVector(1.973,-0.017,0.205,0.139);
  fj.AddInputVector(1.145,-2.125,-1.797,0.139);
  fj.AddInputVector(-0.553,1.068,0.145,0.139);
  fj.AddInputVector(0.478,-0.968,-0.073,0.139);
  fj.AddInputVector(1.520,-1.361,0.789,0.139);
  fj.AddInputVector(-0.961,0.562,-0.157,0.139);
  fj.AddInputVector(1.426,0.412,1.496,0.139);
  fj.AddInputVector(-1.065,1.450,-0.545,0.139);
  fj.AddInputVector(-1.219,1.287,0.026,0.139);
  fj.AddInputVector(0.693,-1.906,-0.892,0.139);
  fj.AddInputVector(2.757,1.510,-0.495,0.139);
  fj.AddInputVector(1.313,0.357,0.495,0.139);
  fj.AddInputVector(-1.420,0.463,1.294,0.139);
  fj.AddInputVector(-3.394,-0.625,1.926,0.139);
  fj.AddInputVector(-0.801,0.600,0.381,0.139);
  fj.AddInputVector(-1.879,1.096,-1.820,0.139);
  fj.AddInputVector(-1.027,-0.064,-0.716,0.139);
  fj.AddInputVector(-1.517,0.390,1.183,0.139);
  fj.AddInputVector(-1.796,-0.032,1.277,0.139);
  fj.AddInputVector(-0.943,0.392,0.831,0.139);
  fj.AddInputVector(0.842,1.361,-0.745,0.139);
  fj.AddInputVector(-1.030,-0.348,-0.416,0.139);
  fj.AddInputVector(-0.443,1.319,-0.930,0.139);
  fj.AddInputVector(-0.755,-2.368,-1.729,0.139);
  fj.AddInputVector(-1.324,-0.438,1.400,0.139);
  fj.AddInputVector(-1.090,0.148,0.906,0.139);
  fj.AddInputVector(-1.371,0.502,1.264,0.139);
  fj.AddInputVector(1.856,0.272,1.301,0.139);
  fj.AddInputVector(-2.183,0.513,-0.876,0.139);
  fj.AddInputVector(0.784,0.938,0.691,0.139);
  fj.AddInputVector(-0.485,1.444,1.265,0.139);
  fj.AddInputVector(-0.965,0.263,0.551,0.139);
  fj.AddInputVector(0.621,-0.920,-0.746,0.139);
  fj.AddInputVector(1.045,0.371,1.051,0.139);
  fj.AddInputVector(-0.523,1.428,-0.349,0.139);
  fj.AddInputVector(0.846,0.727,0.902,0.139);
  fj.AddInputVector(0.909,1.159,1.004,0.139);
  fj.AddInputVector(1.576,-3.255,-3.009,0.139);
  fj.AddInputVector(0.581,1.355,-1.209,0.139);
  fj.AddInputVector(-0.536,-1.237,0.639,0.139);
  fj.AddInputVector(0.669,-1.145,0.575,0.139);
  fj.AddInputVector(-1.349,1.525,-0.178,0.139);
  fj.AddInputVector(-0.309,-1.022,0.636,0.139);
  fj.AddInputVector(-0.102,2.178,-1.994,0.139);
  fj.AddInputVector(-1.059,-0.301,0.367,0.139);
  fj.AddInputVector(1.173,0.579,1.030,0.139);
  fj.AddInputVector(1.102,-1.230,1.011,0.139);
  fj.AddInputVector(1.755,0.137,-0.717,0.139);
  fj.AddInputVector(1.445,1.505,-0.949,0.139);
  fj.AddInputVector(1.015,0.746,0.618,0.139);
  fj.AddInputVector(-9.770,-0.144,2.189,0.139);
  fj.AddInputVector(2.115,-0.031,-0.431,0.139);
  fj.AddInputVector(-0.510,0.951,0.235,0.139);
  fj.AddInputVector(0.491,-0.950,0.550,0.139);
  fj.AddInputVector(1.159,0.835,1.325,0.139);
  fj.AddInputVector(-1.141,-0.492,0.298,0.139);
  fj.AddInputVector(0.421,-1.123,1.104,0.139);
  fj.AddInputVector(-1.714,-1.095,-0.279,0.139);
  fj.AddInputVector(0.699,1.279,-0.442,0.139);
  fj.AddInputVector(1.103,0.439,-0.980,0.139);
  fj.AddInputVector(0.489,0.933,0.878,0.139);
  fj.AddInputVector(-0.081,1.270,0.649,0.139);
  fj.AddInputVector(0.965,0.305,0.491,0.139);
  fj.AddInputVector(0.154,1.601,0.528,0.139);
  fj.AddInputVector(-1.269,0.168,1.282,0.139);
  fj.AddInputVector(-0.567,0.965,0.844,0.139);
  fj.AddInputVector(-0.331,1.074,0.304,0.139);
  fj.AddInputVector(-2.230,-0.379,1.025,0.139);
  fj.AddInputVector(0.891,-0.733,0.562,0.139);
  fj.AddInputVector(0.625,1.366,0.127,0.139);
  fj.AddInputVector(-0.304,1.251,0.819,0.139);
  fj.AddInputVector(-0.276,-0.988,0.946,0.139);
  fj.AddInputVector(-1.327,1.412,-0.701,0.139);
  fj.AddInputVector(-0.079,1.692,-0.515,0.139);
  fj.AddInputVector(-1.347,-0.075,0.598,0.139);
  fj.AddInputVector(1.634,1.022,-1.648,0.139);
  fj.AddInputVector(-3.017,-0.688,0.640,0.139);
  fj.AddInputVector(-1.669,-0.115,1.586,0.139);
  fj.AddInputVector(1.208,0.126,0.487,0.139);
  fj.AddInputVector(0.229,1.484,-1.089,0.139);
  fj.AddInputVector(-1.141,-0.212,-0.602,0.139);
  fj.AddInputVector(-0.792,-2.082,-0.173,0.139);
  fj.AddInputVector(-1.150,0.872,0.209,0.139);
  fj.AddInputVector(0.522,0.892,0.292,0.139);
  fj.AddInputVector(-2.079,1.823,-1.736,0.139);
  fj.AddInputVector(-1.145,0.517,-0.146,0.139);
  fj.AddInputVector(0.397,1.088,0.432,0.139);
  fj.AddInputVector(-1.531,-1.080,0.166,0.139);
  fj.AddInputVector(-1.008,-0.848,0.581,0.139);
  fj.AddInputVector(-1.150,-1.888,-1.850,0.139);
  fj.AddInputVector(-1.118,-0.785,-0.377,0.139);
  fj.AddInputVector(0.382,0.926,0.142,0.139);
  fj.AddInputVector(0.535,-0.903,0.355,0.139);
  fj.AddInputVector(-0.646,1.106,0.591,0.139);
  fj.AddInputVector(0.817,-0.679,0.486,0.139);
  fj.AddInputVector(-0.481,1.190,-1.231,0.139);
  fj.AddInputVector(-0.991,0.496,0.255,0.139);
  fj.AddInputVector(0.066,1.504,0.917,0.139);
  fj.AddInputVector(-0.483,-1.015,0.115,0.139);
  fj.AddInputVector(-1.071,0.873,1.397,0.139);
  fj.AddInputVector(-0.703,0.969,0.750,0.139);
  fj.AddInputVector(-1.425,-0.139,-0.187,0.139);
  fj.AddInputVector(0.911,0.508,-0.585,0.139);
  fj.AddInputVector(-1.939,-0.698,0.388,0.139);
  fj.AddInputVector(-0.851,-0.592,0.044,0.139);
  fj.AddInputVector(-1.244,-0.984,-1.002,0.139);
  fj.AddInputVector(-0.493,0.913,0.226,0.139);
  fj.AddInputVector(0.801,-1.471,-0.288,0.139);
  fj.AddInputVector(1.106,-1.020,-0.692,0.139);
  fj.AddInputVector(0.784,-0.892,-0.965,0.139);
  fj.AddInputVector(-0.400,-1.137,-0.590,0.139);
  fj.AddInputVector(0.572,0.853,-0.384,0.139);
  fj.AddInputVector(1.473,0.761,-0.540,0.139);
  fj.AddInputVector(-1.260,-0.053,0.487,0.139);
  fj.AddInputVector(-1.234,0.649,0.015,0.139);
  fj.AddInputVector(0.662,-1.619,-0.239,0.139);
  fj.AddInputVector(0.877,1.567,0.551,0.139);
  fj.AddInputVector(-2.010,-0.172,-1.013,0.139);
  fj.AddInputVector(1.081,0.280,0.127,0.139);
  fj.AddInputVector(1.131,0.457,-0.221,0.139);
  fj.AddInputVector(0.533,1.738,-0.149,0.139);
  fj.AddInputVector(-0.831,0.781,0.165,0.139);
  fj.AddInputVector(-0.370,-1.000,-0.206,0.139);
  fj.AddInputVector(0.122,1.571,-0.236,0.139);
  fj.AddInputVector(1.136,0.772,1.318,0.139);
  fj.AddInputVector(0.384,1.119,0.724,0.139);
  fj.AddInputVector(-0.894,0.724,-0.450,0.139);
  fj.AddInputVector(-0.081,1.565,0.854,0.139);
  fj.AddInputVector(-0.506,-1.000,0.222,0.139);
  fj.AddInputVector(0.921,0.559,0.314,0.139);
  fj.AddInputVector(-0.529,-1.272,0.796,0.139);
  fj.AddInputVector(0.533,-0.902,0.210,0.139);
  fj.AddInputVector(0.425,0.940,-0.531,0.139);
  fj.AddInputVector(1.722,-1.840,-0.257,0.139);
  fj.AddInputVector(-1.016,-0.858,0.278,0.139);
  fj.AddInputVector(-0.664,1.224,0.544,0.139);
  fj.AddInputVector(-0.312,-1.011,-0.616,0.139);
  fj.AddInputVector(-0.854,1.242,0.106,0.139);
  fj.AddInputVector(0.553,-0.974,-0.356,0.139);
  fj.AddInputVector(0.343,1.152,0.588,0.139);
  fj.AddInputVector(1.227,0.140,-1.025,0.139);
  fj.AddInputVector(-0.207,1.081,0.713,0.139);
  fj.AddInputVector(0.938,0.625,0.777,0.139);
  fj.AddInputVector(1.401,0.145,-1.010,0.139);
  fj.AddInputVector(1.773,0.196,0.401,0.139);
  fj.AddInputVector(0.872,0.825,1.131,0.139);
  fj.AddInputVector(-0.649,0.839,0.398,0.139);
  fj.AddInputVector(-1.355,-0.421,-0.019,0.139);
  fj.AddInputVector(0.724,0.919,-0.863,0.139);
  fj.AddInputVector(-0.011,1.063,0.449,0.139);
  fj.AddInputVector(-0.385,-1.014,0.662,0.139);
  fj.AddInputVector(0.098,1.210,0.562,0.139);
  fj.AddInputVector(0.422,-0.934,-0.087,0.139);
  fj.AddInputVector(1.083,0.153,0.686,0.139);
  fj.AddInputVector(-1.267,0.164,0.029,0.139);
  fj.AddInputVector(-1.130,0.755,0.945,0.139);
  fj.AddInputVector(1.273,0.544,0.341,0.139);
  fj.AddInputVector(2.161,-1.881,0.193,0.139);
  fj.AddInputVector(-0.784,0.743,0.506,0.139);
  fj.AddInputVector(-0.908,-0.615,-0.394,0.139);
  fj.AddInputVector(1.735,0.110,-1.267,0.139);
  fj.AddInputVector(-1.176,0.771,0.287,0.139);
  fj.AddInputVector(0.740,-1.496,1.570,0.139);
  fj.AddInputVector(0.764,1.285,-0.059,0.139);
  fj.AddInputVector(0.391,-1.022,-0.534,0.139);
  fj.AddInputVector(-0.313,-1.036,1.057,0.139);
  fj.AddInputVector(-4.212,-0.689,1.284,0.139);
  fj.AddInputVector(-0.158,-1.201,-0.988,0.139);
  fj.AddInputVector(-0.136,-1.261,0.433,0.139);
  fj.AddInputVector(2.288,-1.645,2.103,0.139);
  fj.AddInputVector(1.418,-0.856,0.726,0.139);
  fj.AddInputVector(-0.279,-1.783,-1.175,0.139);
  fj.AddInputVector(1.508,-0.969,0.551,0.139);
  fj.AddInputVector(0.108,-1.322,1.296,0.139);
  fj.AddInputVector(1.660,-0.356,-1.138,0.139);
  fj.AddInputVector(1.118,-0.608,-1.116,0.139);
  fj.AddInputVector(-0.053,-1.288,-0.284,0.139);
  fj.AddInputVector(-0.890,-0.912,0.306,0.139);
  fj.AddInputVector(-0.412,0.916,0.810,0.139);
  fj.AddInputVector(1.516,-0.941,1.793,0.139);
  fj.AddInputVector(-1.208,-1.527,1.309,0.139);
  fj.AddInputVector(-0.056,-2.600,-0.980,0.139);
  fj.AddInputVector(2.389,-1.379,-2.317,0.139);
  fj.AddInputVector(5.719,-0.703,1.567,0.139);
  fj.AddInputVector(-0.090,-1.660,1.014,0.139);
  fj.AddInputVector(1.417,-0.043,-0.078,0.139);
  fj.AddInputVector(-0.313,-1.290,1.021,0.139);
  fj.AddInputVector(1.531,-1.019,-1.006,0.139);
  fj.AddInputVector(2.651,-0.478,0.529,0.139);
  fj.AddInputVector(19.032,-0.859,4.776,0.139);
  fj.AddInputVector(2.798,-0.208,-0.794,0.139);
  fj.AddInputVector(-1.595,-0.268,-0.963,0.139);
  fj.AddInputVector(0.168,-1.529,0.865,0.139);
  fj.AddInputVector(0.135,-1.222,0.437,0.139);
  fj.AddInputVector(1.003,-0.597,1.060,0.139);
  fj.AddInputVector(1.810,-0.947,1.291,0.139);
  fj.AddInputVector(1.556,-0.912,0.258,0.139);
  fj.AddInputVector(1.538,-0.610,1.148,0.139);
  fj.AddInputVector(-0.066,-1.081,0.501,0.139);
  fj.AddInputVector(1.064,-0.173,-0.301,0.139);
  fj.AddInputVector(-0.282,-0.962,0.872,0.139);
  fj.AddInputVector(-0.957,-1.037,-1.062,0.139);
  fj.AddInputVector(0.964,-0.574,-0.629,0.139);
  fj.AddInputVector(-1.887,-1.798,0.629,0.139);
  fj.AddInputVector(2.816,-2.077,-1.047,0.139);
  fj.AddInputVector(1.169,-0.216,-0.437,0.139);
  fj.AddInputVector(1.458,-0.507,0.771,0.139);
  fj.AddInputVector(1.269,-0.727,0.482,0.139);
  fj.AddInputVector(-0.948,-1.114,-0.758,0.139);
  fj.AddInputVector(0.229,-1.123,0.961,0.139);
  fj.AddInputVector(-0.025,-1.210,0.553,0.139);
  fj.AddInputVector(2.335,-0.357,1.559,0.139);
  fj.AddInputVector(-0.039,-1.167,0.059,0.139);
  fj.AddInputVector(1.027,-0.259,-0.439,0.139);
  fj.AddInputVector(-0.859,0.549,-0.042,0.139);
  fj.AddInputVector(1.761,-0.702,-1.581,0.139);
  fj.AddInputVector(0.154,-1.312,0.759,0.139);
  fj.AddInputVector(1.988,-0.581,-0.479,0.139);
  fj.AddInputVector(1.102,-0.524,-0.323,0.139);
  fj.AddInputVector(-0.860,-1.185,0.515,0.139);
  fj.AddInputVector(0.122,-1.618,-0.672,0.139);
  fj.AddInputVector(0.037,-1.067,0.139,0.139);
  fj.AddInputVector(1.081,-0.153,-0.508,0.139);
  fj.AddInputVector(1.016,-0.079,0.242,0.139);
  fj.AddInputVector(1.213,-0.495,-0.414,0.139);
  fj.AddInputVector(0.883,-0.473,-0.182,0.139);
  fj.AddInputVector(0.344,-1.087,0.578,0.139);
  fj.AddInputVector(-0.695,-0.839,-0.636,0.139);
  fj.AddInputVector(-0.234,-1.299,0.729,0.139);
  fj.AddInputVector(1.205,-0.466,-0.333,0.139);
  fj.AddInputVector(0.099,-1.083,-0.092,0.139);
  fj.AddInputVector(1.411,-0.488,-0.147,0.139);
  fj.AddInputVector(0.173,-1.021,-0.605,0.139);
  fj.AddInputVector(1.145,-0.731,-0.049,0.139);
  fj.AddInputVector(-0.291,-1.294,-0.888,0.139);
  fj.AddInputVector(1.148,-0.398,-0.235,0.139);
  fj.AddInputVector(1.231,-0.050,0.584,0.139);
  fj.AddInputVector(1.805,-0.450,0.904,0.139);
  fj.AddInputVector(1.614,-0.908,0.465,0.139);
  fj.AddInputVector(1.166,-0.152,0.891,0.139);
  fj.AddInputVector(0.957,-0.381,1.031,0.139);
  fj.AddInputVector(1.051,-0.442,0.871,0.139);
  fj.AddInputVector(-0.321,-1.676,-0.592,0.139);
  fj.AddInputVector(1.004,-0.590,-0.133,0.139);
  fj.AddInputVector(1.250,-0.721,-0.286,0.139);
  fj.AddInputVector(-0.035,-1.025,-0.318,0.139);
  fj.AddInputVector(0.490,-1.678,-0.189,0.139);
  fj.AddInputVector(1.400,-0.373,-0.305,0.139);
  fj.AddInputVector(-0.359,-1.667,-0.202,0.139);
  fj.AddInputVector(-0.702,-0.913,0.331,0.139);
  fj.AddInputVector(-0.409,0.961,0.490,0.139);
  fj.AddInputVector(-0.658,-0.885,0.777,0.139);
  fj.AddInputVector(1.382,-0.187,0.457,0.139);
  fj.AddInputVector(0.885,-0.608,0.894,0.139);
  fj.AddInputVector(0.960,-0.501,0.711,0.139);
  fj.AddInputVector(1.409,-0.108,0.384,0.139);
  fj.AddInputVector(1.424,-0.674,0.932,0.139);
  fj.AddInputVector(0.049,-1.148,1.169,0.139);
  fj.AddInputVector(0.294,-1.013,0.315,0.139);
  fj.AddInputVector(1.021,-0.622,-0.960,0.139);
  fj.AddInputVector(1.185,-0.509,0.617,0.139);
  fj.AddInputVector(0.077,-1.202,-1.028,0.139);
  fj.AddInputVector(-0.099,-1.046,-0.070,0.139);
  fj.AddInputVector(-0.616,-0.820,-0.092,0.139);
  fj.AddInputVector(-1.209,-1.158,-1.372,0.139);
  fj.AddInputVector(2.178,-0.623,-1.991,0.139);
  fj.AddInputVector(0.921,-0.435,-0.884,0.139);
  fj.AddInputVector(0.430,-1.824,1.534,0.139);
  fj.AddInputVector(0.085,-1.239,0.635,0.139);
  fj.AddInputVector(0.119,-1.207,-1.128,0.139);
  fj.AddInputVector(2.471,-0.205,-1.269,0.139);
  fj.AddInputVector(1.261,-0.195,0.312,0.139);
  fj.AddInputVector(1.760,-0.767,0.844,0.139);
  fj.AddInputVector(-1.091,-1.477,-1.537,0.139);
  fj.AddInputVector(-0.212,-2.040,-1.248,0.139);
}
