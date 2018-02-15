#include <vector>

void scale_diff(){


  TFile *up = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_up.root");
  TFile *down = new TFile("uhh2.AnalysisModuleRunner.MC.TTbar_down.root");


  std::vector<TH1F*> up_hists, down_hists;
  std::vector<TString> hist_names;
 
  hist_names.push_back("XCone_cor/M_jet1");
  hist_names.push_back("XCone_cor/pt_jet1");
  hist_names.push_back("XCone_cor/pt_jet2");
  hist_names.push_back("EventHists/HT");
  hist_names.push_back("EventHists/MET");
  hist_names.push_back("MuonHists/pt");
  hist_names.push_back("MuonHists/eta");
  hist_names.push_back("MuonHists/phi");


  for(unsigned int i=0; i<hist_names.size(); i++){
    up_hists.push_back((TH1F*)up->Get(hist_names[i]));
    down_hists.push_back((TH1F*)down->Get(hist_names[i]));
  }

  for(unsigned int i=0; i<up_hists.size(); i++){
    up_hists[i]->Scale(1/up_hists[i]->Integral());
    down_hists[i]->Scale(1/down_hists[i]->Integral());
    up_hists[i]->SetLineWidth(3);
    down_hists[i]->SetLineWidth(3);
    down_hists[i]->SetLineColor(kRed);
  }


  TCanvas *c = new TCanvas();
  leg = new TLegend(0.65,0.85,0.85,0.65);
  leg->AddEntry(up_hists[0],"up","l");
  leg->AddEntry(down_hists[0],"down","l");
  c->Divide(3,3);

  for(unsigned int i=0; i<up_hists.size(); i++){
    c->cd(i+1);
    up_hists[i]->Draw("E");
    down_hists[i]->Draw("SAME E");
    leg->Draw();
    gPad->RedrawAxis();
  }


  return;
}
