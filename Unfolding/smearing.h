#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TVectorD.h>
#include <TF1.h>
#include <vector>
#include <TLatex.h>
#include <TVirtualFitter.h>
#include <TFitResult.h>


class Smearing{
public:
  Smearing(std::vector<TH1*> outputs);
  TH1* GetResult();
  std::vector<TF1*> GetFits();
  std::vector<TH1D*> GetVariations();
private:
  TH1* result;
  std::vector<TH1D*> variations;
  std::vector<TF1*> fits;

};
