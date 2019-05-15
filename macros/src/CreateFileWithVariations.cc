
#include "../include/CentralInclude.h"


// compile with:
// g++ -o CreateFileWithVariations CreateFileWithVariations.cc `root-config --cflags --glibs`

using namespace std;

int main(int argc, char* argv[]){

  // select which Hist to draw
  TString histname = "XCone_cor_SF/M_jet1_";
  TString histname_for_file = "M_jet1";

  TString dirname = "/nfs/dust/cms/user/schwarzd/MTopJet/PostSelection/combine/";
  TString name_out = "MjetSYS";
  vector<TString> processes = {"DATA", "TTbar", "WJets", "SingleTop", "other"};

  TFile * file_out = new TFile(dirname+name_out+".root","RECREATE");

  // first write nominal hist (all Processes)
  for(unsigned int i=0; i<processes.size(); i++){
    TString filename = dirname;
    filename += "uhh2.AnalysisModuleRunner.";
    if(processes[i] == "DATA"){
      filename += "DATA.DATA.root";
    }
    else{
      filename += "MC.";
      filename += processes[i];
      filename += ".root";
    }
    TFile *file = new TFile(filename);
    TH1F* hist = (TH1F*)file->Get(histname);
    file_out->cd();
    hist->Write(histname_for_file+"__"+processes[i]);
  }


  vector<TString> sys = {"BKG", "JEC", "JER", "COR", "MUID", "MUTR", "ELID", "ELTR", "ELRECO", "PU", "BTAGudsg", "BTAGbc", "PDF", "ISR", "FSR", "HDAMP", "GENERATOR", "SCALE"};
  // vector<TString> sys = {"BKG", "JEC", "JER", "COR", "MUID", "MUTR", "ELID", "ELTR", "ELRECO", "PU", "BTAGudsg", "BTAGbc"};
  vector<TString> shift = {"EnvelopeUp", "EnvelopeDown"};

  TString filename = "SYS_";
  TString process = "TTbar";

  for(unsigned int i=0; i<sys.size(); i++){
    for(unsigned int j=0; j<shift.size(); j++){
      TFile *file = new TFile(dirname + filename + sys[i] + ".root");
      TH1F* hist = (TH1F*)file->Get(shift[j]+"/hist");
      TString plusminus;
      if(shift[j] == "EnvelopeUp")        plusminus = "plus";
      else if(shift[j] == "EnvelopeDown") plusminus = "minus";
      file_out->cd();
      hist->Write(histname_for_file+"__"+"TTbar"+"__"+sys[i]+"__"+plusminus);
    }
  }
  file_out->Close();

  return 0;
}
