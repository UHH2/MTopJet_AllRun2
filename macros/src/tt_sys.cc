#include "../include/CentralInclude.h"

// compile with:
// g++ -o tt_sys tt_sys.cc `root-config --cflags --glibs`

using namespace std;

TH1F* AddUncertaintyHists(vector<TH1F*> hists);
TH1F* GetLargestError(TH1F* up, TH1F* down);
TH1F* ConvertToRelative(TH1F* central, TH1F* sys);



int main(int argc, char* argv[]){

  double xmin = 100;
  double xmax = 250;

  TFile* STAT_F = new TFile(dir+"STAT.root");

  TFile* BTAG_F = new TFile(dir+"SYS_BTAG.root");
  TFile* COR_F = new TFile(dir+"SYS_COR.root");
  TFile* JEC_F = new TFile(dir+"SYS_JEC.root");
  TFile* JER_F = new TFile(dir+"SYS_JER.root");
  TFile* MUID_F = new TFile(dir+"SYS_MUID.root");
  TFile* MUTR_F = new TFile(dir+"SYS_MUTR.root");
  TFile* PU_F = new TFile(dir+"SYS_PU.root");
  TFile* SCALE_F = new TFile(dir+"SYS_SCALE.root");
  TFile* BKG_F = new TFile(dir+"SYS_BKG.root");
  TFile *GENERATOR_F = new TFile(dir+"SYS_GENERATOR.root");
  TFile *SHOWER_F = new TFile(dir+"SYS_SHOWER.root");
  TFile *PDF_F = new TFile(dir+"SYS_PDF.root");

  TFile * CENTRAL_F = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");


  TH1F* BTAG_UP = (TH1F*)BTAG_F->Get("ErrorUp/hist");
  TH1F* BTAG_DOWN = (TH1F*)BTAG_F->Get("ErrorDown/hist");
  TH1F* COR_UP = (TH1F*)COR_F->Get("ErrorUp/hist");
  TH1F* COR_DOWN = (TH1F*)COR_F->Get("ErrorDown/hist");
  TH1F* JEC_UP = (TH1F*)JEC_F->Get("ErrorUp/hist");
  TH1F* JEC_DOWN = (TH1F*)JEC_F->Get("ErrorDown/hist");
  TH1F* JER_UP = (TH1F*)JER_F->Get("ErrorUp/hist");
  TH1F* JER_DOWN = (TH1F*)JER_F->Get("ErrorDown/hist");
  TH1F* MUID_UP = (TH1F*)MUID_F->Get("ErrorUp/hist");
  TH1F* MUID_DOWN = (TH1F*)MUID_F->Get("ErrorDown/hist");
  TH1F* MUTR_UP = (TH1F*)MUTR_F->Get("ErrorUp/hist");
  TH1F* MUTR_DOWN = (TH1F*)MUTR_F->Get("ErrorDown/hist");
  TH1F* PU_UP = (TH1F*)PU_F->Get("ErrorUp/hist");
  TH1F* PU_DOWN = (TH1F*)PU_F->Get("ErrorDown/hist");
  TH1F* SCALE_UP = (TH1F*)SCALE_F->Get("ErrorUp/hist");
  TH1F* SCALE_DOWN = (TH1F*)SCALE_F->Get("ErrorDown/hist");
  TH1F* BKG_UP = (TH1F*)BKG_F->Get("ErrorUp/hist");
  TH1F* BKG_DOWN = (TH1F*)BKG_F->Get("ErrorDown/hist");
  TH1F* GENERATOR_UP = (TH1F*)GENERATOR_F->Get("ErrorUp/hist");
  TH1F* GENERATOR_DOWN = (TH1F*)GENERATOR_F->Get("ErrorDown/hist");
  TH1F* SHOWER_UP = (TH1F*)SHOWER_F->Get("ErrorUp/hist");
  TH1F* SHOWER_DOWN = (TH1F*)SHOWER_F->Get("ErrorDown/hist");
  TH1F* PDF_UP = (TH1F*)PDF_F->Get("ErrorUp/hist");
  TH1F* PDF_DOWN = (TH1F*)PDF_F->Get("ErrorDown/hist");

  TH1F* CENTRAL = (TH1F*)CENTRAL_F->Get("XCone_cor/M_jet1_B");

  // only look at largest error
  TH1F* STAT = (TH1F*)STAT_F->Get("Down/hist");
  TH1F* BTAG = GetLargestError(BTAG_UP, BTAG_DOWN);
  TH1F* COR = GetLargestError(COR_UP, COR_DOWN);
  TH1F* JEC = GetLargestError(JEC_UP, JEC_DOWN);
  TH1F* JER = GetLargestError(JER_UP, JER_DOWN);
  TH1F* MUID = GetLargestError(MUID_UP, MUID_DOWN);
  TH1F* MUTR = GetLargestError(MUTR_UP, MUTR_DOWN);
  TH1F* PU = GetLargestError(PU_UP, PU_DOWN);
  TH1F* SCALE = GetLargestError(SCALE_UP, SCALE_DOWN);
  TH1F* BKG = GetLargestError(BKG_UP, BKG_DOWN);
  TH1F* Generator = GetLargestError(GENERATOR_UP, GENERATOR_DOWN);
  TH1F* Shower = GetLargestError(SHOWER_UP, SHOWER_DOWN);
  TH1F* PDF = GetLargestError(PDF_UP, PDF_DOWN);


  // sum up all uncert. on muons and JEC
  vector<TH1F*> JEC_v;
  JEC_v.push_back(JEC);
  JEC_v.push_back(JER);
  JEC_v.push_back(COR);
  TH1F* JETSYS = AddUncertaintyHists(JEC_v);

  vector<TH1F*> Mu_v;
  Mu_v.push_back(MUID);
  Mu_v.push_back(MUTR);
  TH1F* MUSYS = AddUncertaintyHists(Mu_v);


  vector<TH1F*> Exp_v;
  Exp_v.push_back(JETSYS);
  Exp_v.push_back(MUSYS);
  Exp_v.push_back(BTAG);
  Exp_v.push_back(PU);
  Exp_v.push_back(BKG);
  Exp_v.push_back(STAT);
  TH1F* EXPSYS = AddUncertaintyHists(Exp_v);

  vector<TH1F*> Model_v;
  Model_v.push_back(SCALE);
  Model_v.push_back(Generator);
  Model_v.push_back(Shower);
  Model_v.push_back(PDF);
  Model_v.push_back(STAT);
  TH1F* MODELSYS = AddUncertaintyHists(Model_v);


  TH1F* JEC_rel = ConvertToRelative(CENTRAL, JEC);
  TH1F* JER_rel = ConvertToRelative(CENTRAL, JER);
  TH1F* COR_rel = ConvertToRelative(CENTRAL, COR);
  TH1F* JETSYS_rel = ConvertToRelative(CENTRAL, JETSYS);
  TH1F* Muon_rel = ConvertToRelative(CENTRAL, MUSYS);
  TH1F* BTag_rel = ConvertToRelative(CENTRAL, BTAG);
  TH1F* PU_rel = ConvertToRelative(CENTRAL, PU);
  TH1F* Scale_rel = ConvertToRelative(CENTRAL, SCALE);
  TH1F* BKG_rel = ConvertToRelative(CENTRAL, BKG);
  TH1F* Generator_rel = ConvertToRelative(CENTRAL, Generator);
  TH1F* Shower_rel = ConvertToRelative(CENTRAL, Shower);
  TH1F* PDF_rel = ConvertToRelative(CENTRAL, PDF);

  TH1F* Model_rel = ConvertToRelative(CENTRAL, MODELSYS);
  TH1F* Exp_rel = ConvertToRelative(CENTRAL, EXPSYS);
  TH1F* STAT_rel = ConvertToRelative(CENTRAL, STAT);



  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);

  /*
  ██████  ██████   █████  ██     ██          ██ ███████ ████████     ███████ ██    ██ ███████
  ██   ██ ██   ██ ██   ██ ██     ██          ██ ██         ██        ██       ██  ██  ██
  ██   ██ ██████  ███████ ██  █  ██          ██ █████      ██        ███████   ████   ███████
  ██   ██ ██   ██ ██   ██ ██ ███ ██     ██   ██ ██         ██             ██    ██         ██
  ██████  ██   ██ ██   ██  ███ ███       █████  ███████    ██        ███████    ██    ███████
  */


  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  TH1F* plots[] = {JETSYS_rel, JEC_rel, JER_rel, COR_rel};
  // STYLES
  int LineWidth[] = {2, 4, 4, 4};
  int MarkerStyle[] = {-1, 0, 0, 0};
  int FillStyle[] = {3144, 0, 0, 0};
  Color_t col[] = {13, kRed-4, kAzure+7, 1};
  int j=0;
  for(auto i:plots){
    i->SetTitle("");
    i->GetXaxis()->SetTitle("Leading-Jet Mass [GeV]");
    i->GetYaxis()->SetTitle("relative uncertainty [%]");
    i->GetXaxis()->SetRangeUser(xmin, xmax);
    i->GetYaxis()->SetRangeUser(0, 50);
    i->GetXaxis()->SetNdivisions(505);
    i->GetYaxis()->SetNdivisions(505);
    i->GetYaxis()->SetTitleOffset(1.5);
    i->SetLineWidth(LineWidth[j]);
    i->SetFillStyle(FillStyle[j]);
    i->SetLineColor(col[j]);
    i->SetFillColor(col[j]);
    i->SetMarkerStyle(MarkerStyle[j]);
    j++;
  }
  // DRAW
  JETSYS_rel->Draw("HIST");
  JEC_rel->Draw("PX SAME");
  JER_rel->Draw("PX SAME");
  COR_rel->Draw("PX SAME");
  // LEGEND
  TLegend *leg = new TLegend(0.63,0.6,0.88,0.88);
  leg->AddEntry(JETSYS_rel,"jet sys combined","f");
  leg->AddEntry(JEC_rel,"JEC","l");
  leg->AddEntry(JER_rel,"JER","l");
  leg->AddEntry(COR_rel,"add. correction","l");
  leg->Draw();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/SYS_Hists/JetSYS.pdf");

  /*
  ██████  ██████   █████  ██     ██     ███████ ██   ██ ██████      ███████ ██    ██ ███████
  ██   ██ ██   ██ ██   ██ ██     ██     ██       ██ ██  ██   ██     ██       ██  ██  ██
  ██   ██ ██████  ███████ ██  █  ██     █████     ███   ██████      ███████   ████   ███████
  ██   ██ ██   ██ ██   ██ ██ ███ ██     ██       ██ ██  ██               ██    ██         ██
  ██████  ██   ██ ██   ██  ███ ███      ███████ ██   ██ ██          ███████    ██    ███████
  */



  TCanvas *B = new TCanvas("B", "B", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  Exp_rel->SetTitle(" ");
  Exp_rel->GetXaxis()->SetTitle("Leading-Jet Mass [GeV]");
  Exp_rel->GetYaxis()->SetTitle("relative uncertainty [%]");
  Exp_rel->GetXaxis()->SetRangeUser(xmin, xmax);
  Exp_rel->GetYaxis()->SetRangeUser(0, 55);
  Exp_rel->GetXaxis()->SetNdivisions(505);
  Exp_rel->GetYaxis()->SetNdivisions(505);
  Exp_rel->GetYaxis()->SetTitleOffset(1.5);

  for(auto i:{STAT_rel, JETSYS_rel, Muon_rel, BTag_rel, PU_rel, BKG_rel}){
    i->SetLineWidth(4);
    i->SetMarkerStyle(0);
  }

  Exp_rel->SetLineColor(13);
  Exp_rel->SetFillColor(13);
  Exp_rel->SetFillStyle(3144);
  STAT_rel->SetLineColor(kBlack);
  JETSYS_rel->SetLineColor(kBlue);
  Muon_rel->SetLineColor(kAzure+7);
  BTag_rel->SetLineColor(kRed-4);
  PU_rel->SetLineColor(798);
  BKG_rel->SetLineColor(kPink+9);

  Exp_rel->Draw("HIST");
  JETSYS_rel->Draw("PX SAME");
  Muon_rel->Draw("PX SAME");
  BTag_rel->Draw("PX SAME");
  PU_rel->Draw("PX SAME");
  STAT_rel->Draw("PX SAME");
  BKG_rel->Draw("PX SAME");

  TLegend *leg2 = new TLegend(0.63,0.6,0.88,0.88);
  leg2->AddEntry(Exp_rel,"stat #oplus sys","f");
  leg2->AddEntry(STAT_rel,"statistical","l");
  leg2->AddEntry(JETSYS_rel,"jet","l");
  leg2->AddEntry(Muon_rel,"muon","l");
  leg2->AddEntry(BTag_rel,"b tag","l");
  leg2->AddEntry(PU_rel,"pile-up","l");
  leg2->AddEntry(BKG_rel,"background","l");

  leg2->Draw();

  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/SYS_Hists/ExpSYS.pdf");

  /*
  ██████  ██████   █████  ██     ██     ███    ███  ██████  ██████  ███████ ██          ███████ ██    ██ ███████
  ██   ██ ██   ██ ██   ██ ██     ██     ████  ████ ██    ██ ██   ██ ██      ██          ██       ██  ██  ██
  ██   ██ ██████  ███████ ██  █  ██     ██ ████ ██ ██    ██ ██   ██ █████   ██          ███████   ████   ███████
  ██   ██ ██   ██ ██   ██ ██ ███ ██     ██  ██  ██ ██    ██ ██   ██ ██      ██               ██    ██         ██
  ██████  ██   ██ ██   ██  ███ ███      ██      ██  ██████  ██████  ███████ ███████     ███████    ██    ███████
  */


  TCanvas *C = new TCanvas("C", "C", 600, 600);
  gPad->SetLeftMargin(0.1);
  TGaxis::SetMaxDigits(3);
  Model_rel->SetTitle(" ");
  Model_rel->GetXaxis()->SetTitle("Leading-Jet Mass [GeV]");
  Model_rel->GetYaxis()->SetTitle("relative uncertainty [%]");
  Model_rel->GetXaxis()->SetRangeUser(xmin, xmax);
  Model_rel->GetYaxis()->SetRangeUser(0, 100);
  Model_rel->GetXaxis()->SetNdivisions(505);
  Model_rel->GetYaxis()->SetNdivisions(505);
  Model_rel->GetYaxis()->SetTitleOffset(1.5);

  for(auto i:{STAT_rel, Scale_rel, Generator_rel, Shower_rel, PDF_rel}){
    i->SetLineWidth(4);
    i->SetMarkerStyle(0);
  }

  Model_rel->SetLineColor(13);
  Model_rel->SetFillColor(13);
  Model_rel->SetFillStyle(3144);

  Scale_rel->SetLineColor(kAzure+7);
  Generator_rel->SetLineColor(kGreen);
  Shower_rel->SetLineColor(kRed);
  PDF_rel->SetLineColor(798);

  STAT_rel->SetLineColor(1);


  Model_rel->Draw("HIST");
  STAT_rel->Draw("PX SAME");
  Scale_rel->Draw("PX SAME");
  Generator_rel->Draw("PX SAME");
  Shower_rel->Draw("PX SAME");
  PDF_rel->Draw("PX SAME");

  TLegend *leg3 = new TLegend(0.63,0.6,0.88,0.88);
  leg3->AddEntry(Model_rel,"stat #oplus model sys","f");
  leg3->AddEntry(STAT_rel,"statistical","l");
  leg3->AddEntry(Scale_rel,"#mu_{f} and #mu_{r} scales","l");
  leg3->AddEntry(PDF_rel,"PDF","l");
  leg3->AddEntry(Generator_rel,"MC Generator","l");
  leg3->AddEntry(Shower_rel,"Shower Model","l");

  leg3->Draw();

  C->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/SYS_Hists/ModelSYS.pdf");



  return 0;



}

/*
██   ██ ███████ ██      ██████  ███████ ██████      ███████ ██    ██ ███    ██  ██████ ████████ ██  ██████  ███    ██ ███████
██   ██ ██      ██      ██   ██ ██      ██   ██     ██      ██    ██ ████   ██ ██         ██    ██ ██    ██ ████   ██ ██
███████ █████   ██      ██████  █████   ██████      █████   ██    ██ ██ ██  ██ ██         ██    ██ ██    ██ ██ ██  ██ ███████
██   ██ ██      ██      ██      ██      ██   ██     ██      ██    ██ ██  ██ ██ ██         ██    ██ ██    ██ ██  ██ ██      ██
██   ██ ███████ ███████ ██      ███████ ██   ██     ██       ██████  ██   ████  ██████    ██    ██  ██████  ██   ████ ███████
*/

TH1F* AddUncertaintyHists(vector<TH1F*> hists){
  int nbins = hists[0]->GetSize() - 2;
  TH1F* hist = (TH1F*)hists[0]->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double cont = 0;
    for(unsigned int i=0; i<hists.size(); i++){
      cont += pow(hists[i]->GetBinContent(bin), 2);
    }
    double sum = sqrt( cont );
    hist->SetBinContent(bin, sum);
  }
  return hist;
}

TH1F* GetLargestError(TH1F* up, TH1F* down){
  int nbins = up->GetSize() - 2;
  TH1F* hist = (TH1F*)up->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double UPcont = abs(up->GetBinContent(bin));
    double DOWNcont = abs(down->GetBinContent(bin));
    if(UPcont > DOWNcont) hist->SetBinContent(bin, UPcont);
    else                  hist->SetBinContent(bin, DOWNcont);
  }
  return hist;
}


TH1F* ConvertToRelative(TH1F* central, TH1F* sys){
  int nbins = central->GetSize() - 2;
  TH1F* hist = (TH1F*)sys->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double cont_central = central->GetBinContent(bin);
    double cont_sys = sys->GetBinContent(bin);
    double percent;
    if(cont_central == 0) percent = 0;
    else percent = 100*cont_sys/cont_central;
    hist->SetBinContent(bin, percent);
    hist->SetBinError(bin, 0.000000001); // just for plotting reasons
  }
  return hist;
}
