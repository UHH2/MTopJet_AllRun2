#include "../include/CentralInclude.h"

using namespace std;


/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/

void ConvertUncertTables(TH1*, string, bool);
void ConvertMassContributions(string);
void ConvertCovToTable(TH2*, string, string, bool);


int main(int argc, char* argv[]){

  string directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/Data/combine/";


  TFile *file = new TFile("/nfs/dust/cms/user/schwarzd/CMSSW_8_0_24_patch1/src/UHH2/MTopJet/Unfolding/Results_data_combine.root");
  TH1* unfold = (TH1*)file->Get("Unfold_XS_totuncert");
  TH1* unfold_norm = (TH1*)file->Get("Unfold_norm_totuncert");
  TH2* cov = (TH2*)file->Get("CovTotalXS");
  TH2* cov_norm = (TH2*)file->Get("CovTotalNorm");

  ConvertUncertTables(unfold, directory, false);
  ConvertUncertTables(unfold_norm, directory, true);
  ConvertCovToTable(cov, directory, "COV_TOTAL", false);
  ConvertCovToTable(cov_norm, directory, "COV_TOTAL_NORM", true);
  ConvertMassContributions(directory);


  return 0;
}

//------------------------------------------------------------------------------

void ConvertUncertTables(TH1* unfold, string directory, bool norm){
  vector<string> filenames;
  if(norm){
    filenames.push_back("SYS_EXP_COMPARISION_NORM");
    filenames.push_back("SYS_MODEL_COMPARISION_NORM");
    cout << "converting " << "Uncertainties Table" << endl;
  }
  else{
    filenames.push_back("SYS_EXP_COMPARISION");
    filenames.push_back("SYS_MODEL_COMPARISION");
    cout << "converting " << "Norm Uncertainties Table" << endl;
  }

  // setup output file
  string outname;
  if(norm) outname = "LATEX_Uncertainties_norm";
  else     outname = "LATEX_Uncertainties";
  std::ofstream outs(directory+outname+".txt");
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());

  // read input file
  string uncert;
  double bin1, bin2, bin3, bin4, bin5;
  cout << "\\begin{tabular}{ l | r r r r r  }" << endl;
  cout << "Range in \\mjet [GeV] & 112--132 & 132--152 & 152--172 & 172--192 & 192--232 \\\\" << endl;
  cout << "\\hline" << endl;
  if(norm) cout << "Integrated normalised cross section";
  else     cout << "Integrated cross section [fb]";
  int nbins = unfold->GetXaxis()->GetNbins();
  for(int bin=1; bin<=nbins; bin++){
    cout << " & ";
    double width = unfold->GetBinWidth(bin);
    double entry = unfold->GetBinContent(bin);
    if(norm) cout << std::fixed << std::setprecision(2); // 2
    else     cout << std::fixed << std::setprecision(0); // 0
    cout << entry*width;
  }
  cout << "\\\\" << endl;
  bool firstloop = true;
  for(auto filename: filenames){
    std::ifstream infile(directory+filename+".txt");
    while (infile >> uncert >> bin1 >> bin2 >> bin3 >> bin4 >> bin5){
      bool indent = true;
      if(uncert == "total") continue;
      if(uncert == "stat" && firstloop){
        uncert = "Statistical uncertainty [\\%]";
        cout << "\\hline" << endl;
        indent = false;
      }
      if(uncert == "stat" && !firstloop) continue;
      if(uncert == "exp"){
        uncert = "Experimental uncertainty [\\%]";
        cout << "\\hline" << endl;
        indent = false;
      }
      if(uncert == "model"){
        uncert = "Model uncertainty [\\%]";
        cout << "\\hline" << endl;
        indent = false;
      }
      if(uncert == "jec") uncert = "Jet energy scale";
      if(uncert == "jer") uncert = "Jet energy resolution";
      if(uncert == "cor") uncert = "XCone jet correction";
      if(uncert == "hdamp") uncert = "$h_\\textrm{damp}$";
      if(uncert == "mass") uncert = "Choice of \\mtop";
      if(uncert == "UEtune") uncert = "UE tune";
      if(uncert == "b-tagging") uncert = "b tag";
      if(uncert == "pile-up") uncert = "Pileup";
      if(uncert == "scale") uncert = "Scale";
      if(uncert == "pdf") uncert = "PDF";
      if(uncert == "MCstat") uncert = "MC stat.";
      if(indent) cout << "\\hspace{3mm}";
      cout << std::fixed << std::setprecision(1);
      cout << uncert;
      if(bin1 < 1.0) cout << " & $< 1$";
      else           cout << " & " << bin1;
      if(bin2 < 1.0) cout << " & $< 1$";
      else           cout << " & " << bin2;
      if(bin3 < 1.0) cout << " & $< 1$";
      else           cout << " & " << bin3;
      if(bin4 < 1.0) cout << " & $< 1$";
      else           cout << " & " << bin4;
      if(bin5 < 1.0) cout << " & $< 1$";
      else           cout << " & " << bin5;
      cout << "\\\\"<<endl;
    }
    firstloop = false;
  }
  cout << "\\hline" << endl;
  cout << "Total Uncertainty [\\%]";
  for(int bin=1; bin<=nbins; bin++){
    cout << " & ";
    double uncert = unfold->GetBinError(bin);
    double entry = unfold->GetBinContent(bin);
    double percent = 100*uncert/entry;
    cout << std::fixed << std::setprecision(1);
    cout << percent;
  }
  cout << "\\\\" << endl;
  cout << "\\end{tabular}" << endl;
  std::cout.rdbuf(coutbuf);
}

void ConvertMassContributions(string directory){
  string filename = "Mass_contribs";
  cout << "converting " << directory+filename+".txt" << endl;
  // setup output file
  std::ofstream outs(directory+"LATEX_"+filename+".txt");
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());
  // read input file
  std::ifstream infile(directory+filename+".txt");
  string name;
  double uncert;
  cout << "\\begin{tabular}{ l | r }" << endl;
  cout << "source & uncertainty [GeV] \\\\" << endl;
  cout << "\\hline" << endl;
  while (infile >> name >> uncert){
    bool indent = true;
    if(name == "total"){
      name = "Total uncertainty";
      cout << "\\hline" << endl;
      indent = false;
    }
    if(name == "stat"){
      name = "Statistical uncertainty";
      cout << "\\hline" << endl;
      indent = false;
    }
    if(name == "exp"){
      name = "Experimental uncertainty";
      cout << "\\hline" << endl;
      indent = false;
    }
    if(name == "model"){
      name = "Model uncertainty";
      cout << "\\hline" << endl;
      indent = false;
    }
    if(name == "theo"){
      name = "Theoretical uncertainty";
      cout << "\\hline" << endl;
      indent = false;
    }
    if(name == "jec") name = "Jet energy scale";
    if(name == "jer") name = "Jet energy resolution";
    if(name == "cor") name = "XCone jet correction";
    if(name == "hdamp") name = "$h_\\textrm{damp}$";
    if(name == "mass") name = "Choice of \\mtop";
    if(name == "UEtune") name = "UE tune";
    if(name == "b-tagging") name = "b tag";
    if(name == "pile-up") name = "Pileup";
    if(name == "scale") name = "Scale";
    if(name == "pdf") name = "PDF";
    if(name == "MCstat") name = "MC stat.";
    if(indent) cout << "\\hspace{3mm}";
    cout << std::fixed << std::setprecision(2);
    if(uncert < 0.01) cout << name << " & $< 0.01$ \\\\"<<endl;
    else              cout << name << " & " << uncert << " \\\\"<<endl;
  }
  cout << "\\end{tabular}" << endl;
  std::cout.rdbuf(coutbuf);
}

void ConvertCovToTable(TH2* cov, string directory, string filename, bool norm){
  cout << "converting " << filename << endl;
  // setup output file
  std::ofstream outs(directory+"LATEX_"+filename+".txt");
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());

  cout << "\\begin{tabular}{ c | r r r r r }" << endl;
  cout << "Bin & 1 & 2 & 3 & 4 & 5 \\\\" << endl;
  cout << "\\hline" << endl;
  int nbins = cov->GetXaxis()->GetNbins();
  for(int ybin=1; ybin<=nbins; ybin++){
    cout << ybin;
    for(int skip = 1; skip<ybin; skip++) cout << " & ";
    for(int xbin=ybin; xbin<=nbins; xbin++){
      cout << " & ";
      double width = cov->GetXaxis()->GetBinWidth(xbin);
      double entry = cov->GetBinContent(xbin, ybin);
      if(norm) entry *= 10000;
      cout << std::fixed << std::setprecision(2);
      cout << width*width*entry;
    }
    cout << "\\\\" << endl;
  }
  cout << "\\end{tabular}" << endl;
  std::cout.rdbuf(coutbuf);
}
