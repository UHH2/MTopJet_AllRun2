#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1.h>
#include <TDOMParser.h>
#include <TXMLDocument.h>
#include <time.h>
#include <vector>


using namespace std;


void fill_hist(TTree *, TString);
void fill_pdf(TTree *);


TFile *outputFile;

// variables to store gen or rec info
Double_t massRec; // variables
Bool_t passed_measurement_rec;
Double_t rec_weight;
Double_t gen_weight;

const vector<Double_t> bins = {150, 160, 170, 180, 190, 200};
// const vector<Double_t> bins = {135, 155, 175, 195, 215};
