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


using namespace std;


void fill_data(TTree *);
void fill_template(TTree *, TString);
void fill_background(TTree *);
void fill_matrix(TTree *);
void create_hists();



bool ttweight;
TFile *outputFile;
TFile *outputFile_bins;
bool fast;
int fast_factor;

// binning schemes
TUnfoldBinning *binning_rec;
TUnfoldBinning *binning_gen;
const TUnfoldBinning *measurement_rec;
const TUnfoldBinning *measurement_gen;
const TUnfoldBinning *ptmigration_rec;
const TUnfoldBinning *ptmigration_gen;
const TUnfoldBinning *massmigration_rec;
const TUnfoldBinning *massmigration_gen;
const TUnfoldBinning *btagmigration_rec;


// Hists
TH1 *h_purity_samebin;
TH1 *h_purity_samebin_pt;
TH1 *h_purity_all;

TH1 *h_data;
TH1 *h_mc_rec;
TH1 *h_mc_bgr;
TH1 *h_mc_sig;
TH1 *h_mc_gen;
TH1 *h_mc_truth;

TH2 *h_mc_matrix;

TH1 *h_mc_mtop1665_truth;
TH1 *h_mc_mtop1695_truth;
TH1 *h_mc_mtop1715_truth;
TH1 *h_mc_mtop1735_truth;
TH1 *h_mc_mtop1755_truth;
TH1 *h_mc_mtop1785_truth;
TH1 *h_mc_mtop1665_rec;
TH1 *h_mc_mtop1695_rec;
TH1 *h_mc_mtop1715_rec;
TH1 *h_mc_mtop1735_rec;
TH1 *h_mc_mtop1755_rec;
TH1 *h_mc_mtop1785_rec;


// variables to store gen or rec info
Double_t massRec, ptRec, massGen, ptGen; // variables
Bool_t passed_measurement_rec, passed_measurement_gen, is_TTbar, passed_ptmigration_rec, passed_massmigration_rec, passed_btagmigration_rec, passed_ptmigration_gen, passed_massmigration_gen; // bools
Double_t gen_weight, rec_weight, gen_ttfactor; //weights
Double_t gen_2width_factor, gen_4width_factor, gen_8width_factor; //factor to change width



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





