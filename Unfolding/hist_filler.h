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
#include "TUnfoldBinningXML.h"
#include <time.h>
#include <vector>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>


using namespace std;


void fill_data(TTree *);
void fill_template(TTree *, TString);
void fill_background(TTree *, TString);
void fill_matrix(TTree *, TString);
void fill_modelsys(TTree *, TString);
void fill_pdf(TTree *);
void fill_breweight(TTree *, TString);
void fill_additionaltemplates(TString channel);

TH1* CreateNewSample(vector<TF1*>, double, TH1*, TString);
void PlotBinHists(vector<TGraphErrors*>, vector<TF1*>, TString);



bool ttweight;
TFile *outputFile;
TFile *outputFile_bins;
bool fast;
int fast_factor;

// axis steering
const char *axissteer = "pt[C]";
const char *axissteer_purity = "pt[C]";


// binning schemes
TUnfoldBinning *binning_rec;
TUnfoldBinning *binning_gen;
TUnfoldBinning *binning_gen20;

const TUnfoldBinning *measurement_rec;
const TUnfoldBinning *measurement_gen;
const TUnfoldBinning *measurement_gen20;
const TUnfoldBinning *ptmigration_rec;
const TUnfoldBinning *ptmigration_gen;
const TUnfoldBinning *subptmigration_rec;
const TUnfoldBinning *subptmigration_gen;
const TUnfoldBinning *massmigration_rec;
const TUnfoldBinning *massmigration_gen;
const TUnfoldBinning *leptonptmigration_rec;
const TUnfoldBinning *leptonptmigration_gen;
const TUnfoldBinning *btagmigration_rec;

// variables to store gen or rec info
Double_t massRec, ptRec, massGen, ptGen; // variables
Bool_t passed_measurement_rec, passed_measurement_gen;
Bool_t is_TTbar;
Bool_t passed_ptmigration_rec, passed_ptmigration_gen;
Bool_t passed_subptmigration_rec, passed_subptmigration_gen;
Bool_t passed_leptonptmigration_rec, passed_leptonptmigration_gen;
Bool_t passed_massmigration_rec, passed_massmigration_gen;
Bool_t passed_btagmigration_rec;
Double_t gen_weight, rec_weight, gen_ttfactor; //weights
Double_t gen_2width_factor, gen_4width_factor, gen_8width_factor; //factor to change width

// b quark pt
vector<double> *pt_b = new vector<double>(2);

//weights
Double_t w_central;
Double_t w_nogen;
Double_t w_norec;
Double_t w_correction;
Double_t w_side;

Double_t w_bgr_rec;
Double_t w_all_rec;
Double_t w_sig_rec;
Double_t w_gen;
Double_t w_rec;

Double_t w_gen_w2;
Double_t w_gen_w4;
Double_t w_gen_w8;

// store mtop truth seperately for additional templates
vector<TH1*> h_mtop_truth;
