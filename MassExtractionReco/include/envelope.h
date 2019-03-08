#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include <time.h>
#include <vector>


using namespace std;

TH1F* GetEnvelopeUp(std::vector<TH1F*> Variations);
TH1F* GetEnvelopeDown(std::vector<TH1F*> Variations);
TH1F* GetError(TH1F* Central, TH1F* Envelope);
TH1F* GetStatError(TH1F* Data);
TH2D* GetCovMatrix(TH1F* Error);
TH1F* GetBackgroundEnvelopeUp(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other);
TH1F* GetBackgroundEnvelopeDown(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other);
TH1F* GetSymmetricError(TH1F* up, TH1F* down);
TH1F* DummyConstVariation(TH1F* Central, double factor);
