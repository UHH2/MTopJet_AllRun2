#pragma once
#include <TCanvas.h>
#include <TColor.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TVectorD.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <tuple>
#include <vector>

// #################################################################################################
// Color code ######################################################################################
#define RESET      "\033[0m"
#define BLACK      "\033[30m"      /* Black */
#define RED        "\033[31m"      /* Red */
#define GREEN      "\033[32m"      /* Green */
#define YELLOW     "\033[33m"      /* Yellow */
#define BLUE       "\033[34m"      /* Blue */
#define MAGENTA    "\033[35m"      /* Magenta */
#define CYAN       "\033[36m"      /* Cyan */
#define WHITE      "\033[37m"      /* White */
#define DGRAY      "\033[90m"      /* Dark Gray */

// Background ########################################################################################
#define BACK_DGRAY "\033[100"      /* Dark Gray */
#define BACK_CYAN   "\033[46"       /* Cyan */

using namespace std;

/*
██    ██ ███████  ██████ ████████  ██████  ██████
██    ██ ██      ██         ██    ██    ██ ██   ██
██    ██ █████   ██         ██    ██    ██ ██████
.██  ██  ██      ██         ██    ██    ██ ██   ██
. ████   ███████  ██████    ██     ██████  ██   ██
*/

// -------------------------------------------------------------------------------------------------------
vector<double> square_vector(vector<double> v){
  vector<double> vector_square;
  for(unsigned int i=0; i<v.size(); i++) vector_square.push_back(pow(v[i], 2.));
  return vector_square;
}

// -------------------------------------------------------------------------------------------------------
double GetMaxValue(vector<double> v){
  return *max_element(v.begin(), v.end());
}
double GetMinValue(vector<double> v){
  return *min_element(v.begin(), v.end());
}

// -------------------------------------------------------------------------------------------------------


/*
.██████  ██████  ██    ██ ████████
██      ██    ██ ██    ██    ██
██      ██    ██ ██    ██    ██
██      ██    ██ ██    ██    ██
.██████  ██████   ██████     ██
*/

void print_seperater(){
  cout << "" << endl;
  cout << "*************************************************************************************************" << endl;
  cout << "" << endl;
}

/*
.██████  ██████  ███    ██ ██    ██ ███████ ██████  ████████
██      ██    ██ ████   ██ ██    ██ ██      ██   ██    ██
██      ██    ██ ██ ██  ██ ██    ██ █████   ██████     ██
██      ██    ██ ██  ██ ██  ██  ██  ██      ██   ██    ██
.██████  ██████  ██   ████   ████   ███████ ██   ██    ██
*/

bool stob(std::string s, bool throw_on_error = true){
  auto result = false;    // failure to assert is false
  std::istringstream is(s);
  // first try simple integer conversion
  is >> result;
  if (is.fail()){
    // simple integer failed; try boolean
    is.clear();
    is >> std::boolalpha >> result;
  }

  if (is.fail() && throw_on_error){
    throw std::invalid_argument(s.append(" is not convertable to bool"));
  }
  return result;
}

/*
████████  █████  ██████  ██      ███████
   ██    ██   ██ ██   ██ ██      ██
   ██    ███████ ██████  ██      █████
   ██    ██   ██ ██   ██ ██      ██
   ██    ██   ██ ██████  ███████ ███████
*/

/* Copy: https://stackoverflow.com/questions/14861018/center-text-in-fixed-width-field-with-stream-manipulators-in-c */
template<typename charT, typename traits = std::char_traits<charT> >
class center_helper {
    std::basic_string<charT, traits> str_;
public:
    center_helper(std::basic_string<charT, traits> str) : str_(str) {}
    template<typename a, typename b>
    friend std::basic_ostream<a, b>& operator<<(std::basic_ostream<a, b>& s, const center_helper<a, b>& c);
};

template<typename charT, typename traits = std::char_traits<charT> >
center_helper<charT, traits> centered(std::basic_string<charT, traits> str) {
    return center_helper<charT, traits>(str);
}

// redeclare for std::string directly so we can support anything that implicitly converts to std::string
center_helper<std::string::value_type, std::string::traits_type> centered(const std::string& str) {
    return center_helper<std::string::value_type, std::string::traits_type>(str);
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>& operator<<(std::basic_ostream<charT, traits>& s, const center_helper<charT, traits>& c) {
    std::streamsize w = s.width();
    if (w > c.str_.length()) {
        std::streamsize left = (w + c.str_.length()) / 2;
        s.width(left);
        s << c.str_;
        s.width(w - left);
        s << "";
    } else {
        s << c.str_;
    }
    return s;
}

/*
███████  ██████  ██      ██████  ███████ ██████
██      ██    ██ ██      ██   ██ ██      ██   ██
█████   ██    ██ ██      ██   ██ █████   ██████
██      ██    ██ ██      ██   ██ ██      ██   ██
██       ██████  ███████ ██████  ███████ ██   ██
*/

void creat_folder(TString path, TString name){mkdir(path+"/"+name,0777);}
void creat_folder(TString path){mkdir(path,0777);}

TString creat_folder_and_path(TString path, TString name){
  mkdir(path+"/"+name,0777);
  return path += "/"+name;
}

TString creat_folder_and_path(TString path){
  mkdir(path,0777);
  return path;
}

/*
.██████  ██████  ███    ███ ██████   █████  ██████  ███████
██      ██    ██ ████  ████ ██   ██ ██   ██ ██   ██ ██
██      ██    ██ ██ ████ ██ ██████  ███████ ██████  █████
██      ██    ██ ██  ██  ██ ██      ██   ██ ██   ██ ██
.██████  ██████  ██      ██ ██      ██   ██ ██   ██ ███████
*/

double choose_greater_number(double x1, double x2){
  double x;
  if(x1<x2) x=x2;
  else      x=x1;
  return x;
}

/*
██████   ██████  ██    ██ ███    ██ ██████
██   ██ ██    ██ ██    ██ ████   ██ ██   ██
██████  ██    ██ ██    ██ ██ ██  ██ ██   ██
██   ██ ██    ██ ██    ██ ██  ██ ██ ██   ██
██   ██  ██████   ██████  ██   ████ ██████
*/

int Round(double wert)
{
	return (int) (wert + ((wert < 0)? - 0.5 : 0.5));
}

double Round(double wert, int nachkommastellen)
{
	return ((int) ((wert * pow(10.0, nachkommastellen)) + ((wert < 0)? - 0.5 : 0.5))) / pow(10.0, nachkommastellen);
}
