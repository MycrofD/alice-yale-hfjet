#include <vector>
#include <string>

class FJWrapper {
public:
  FJWrapper(const char* _name);
  ~FJWrapper();

  void AddInputVector(double px, double py, double pz, double m);
  void FindJets();
  void PrintJets() const;
  const std::vector<fastjet::PseudoJet>& GetJets() const { return incJets; }
  const std::string& GetName() const { return name; }
  void CompareResults(const FJWrapper& fjw);

private:
  std::string                     name;
  fastjet::Strategy               strategy;
  fastjet::JetAlgorithm           algor;
  fastjet::RecombinationScheme    scheme;
  fastjet::AreaType               areaType;
  int                             nGhostRepeats;
  double                          ghostArea;
  double                          maxRap;
  double                          R;
  double                          gridScatter;
  double                          ktScatter;
  double                          meanGhostKt;
  fastjet::ClusterSequenceArea   *seqArea;
  std::vector<fastjet::PseudoJet> inputVectors;
  std::vector<fastjet::PseudoJet> incJets;
};

void LoadEvent(FJWrapper& fj);
