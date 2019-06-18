#include "plotter.h"


plotter::plotter(TString dir){
  directory = dir;
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);
}

/*
███    ███  █████  ████████ ██████  ██ ██   ██
████  ████ ██   ██    ██    ██   ██ ██  ██ ██
██ ████ ██ ███████    ██    ██████  ██   ███
██  ██  ██ ██   ██    ██    ██   ██ ██  ██ ██
██      ██ ██   ██    ██    ██   ██ ██ ██   ██
*/



void plotter::draw_matrix(TH2* hist_, TString file_name, bool zlog, bool is_migration){
  TH1* hist = (TH1*) hist_->Clone("hist");

  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  if(zlog) gPad->SetLogz();
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(1.5);
  if(is_migration){
    hist->GetXaxis()->SetTitle("generator binning");
    hist->GetYaxis()->SetTitle("detector binning");
  }
  else{
    hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
    hist->GetYaxis()->SetTitle("m_{jet} [GeV]");
  }

  hist->Draw("COLZ");
  hist->Draw("BOX SAME");
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██    ██ ████████ ██████  ██    ██ ████████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██
██    ██ ██    ██    ██    ██████  ██    ██    ██
██    ██ ██    ██    ██    ██      ██    ██    ██
██████   ██████     ██    ██       ██████     ██
*/

void plotter::draw_output(TH1* output_, TH1D* truth_, bool norm, TString file_name){
  // std::vector<double> sys = get_sys_errors();
  // TH1* output_sys = add_error_bar(output, sys);

  TH1* output = (TH1*) output_->Clone("output");
  TH1D* truth = (TH1D*) truth_->Clone("truth");

  TCanvas *c = new TCanvas("c","",600,600);
  double ymax;
  gPad->SetLeftMargin(0.19);

  if(truth->GetMaximum() > output->GetMaximum()) ymax = 1.5 * truth->GetMaximum();
  else ymax = 1.5 * output->GetMaximum();
  // TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.5);
  output->GetXaxis()->SetTitleOffset(0.9);
  output->GetYaxis()->SetTitleSize(0.05);
  output->GetXaxis()->SetTitleSize(0.05);
  output->GetYaxis()->SetNdivisions(505);
  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);
  output->Draw("E1 SAME");
  gStyle->SetEndErrorSize(5);
  truth->SetLineWidth(3);
  truth->SetLineColor(kRed);
  truth->SetLineStyle(2);
  truth->Draw("HIST SAME");
  TLegend *l;
  if(truth->GetSize() > 20) l=new TLegend(0.2,0.6,0.4,0.85);
  else                       l=new TLegend(0.5,0.6,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"Data","ple");
  l->AddEntry(truth,"MC particle level","pl");
  l->SetTextSize(0.04);
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██    ██ ████████ ██████  ██    ██ ████████     ██████   █████  ████████  █████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██   ██ ██   ██    ██    ██   ██
██    ██ ██    ██    ██    ██████  ██    ██    ██        ██   ██ ███████    ██    ███████
██    ██ ██    ██    ██    ██      ██    ██    ██        ██   ██ ██   ██    ██    ██   ██
██████   ██████     ██    ██       ██████     ██        ██████  ██   ██    ██    ██   ██
*/



void plotter::draw_output_data(TH1* output_, TH1* stat_, std::vector<TH1D*> truth_, std::vector<TString> legnames, bool norm, TString file_name){

  TH1* output = (TH1*) output_->Clone("output");
  TH1* stat = (TH1*) stat_->Clone("stat");

  std::vector<TH1D*> truth, truth2;
  for(auto t: truth_){
    truth.push_back( (TH1D*) t->Clone() );
    truth2.push_back( (TH1D*) t->Clone() );
  }

  double max = output->GetMaximum();
  for(auto t: truth){
    if(t->GetMaximum() > max) max = t->GetMaximum();
  }
  double ymax = 1.5 * max;

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.19);
  // TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("#frac{d#sigma}{dm_{jet}} [#frac{fb}{GeV}]");
  output->GetYaxis()->SetTitleOffset(1.5);
  output->GetXaxis()->SetTitleOffset(0.9);
  output->GetYaxis()->SetTitleSize(0.05);
  output->GetXaxis()->SetTitleSize(0.05);
  output->GetYaxis()->SetNdivisions(505);
  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);
  output->Draw("E1");
  stat->SetLineColor(kBlack);
  stat->SetMarkerColor(kBlack);
  stat->SetMarkerStyle(8);
  stat->SetMarkerSize(1);
  gStyle->SetEndErrorSize(5);
  Int_t color[] = {TColor::GetColor("#0059b3"), TColor::GetColor("#e67300")};
  Int_t fillcolor[] = {TColor::GetColor("#99ccff"), TColor::GetColor("#ff8000")};
  // Color_t color[] = {kAzure+2, kOrange+7};
  // Color_t fillcolor[] = {kAzure+10, 798};
  Int_t style[] = {1, 2};
  Int_t fillstyle[] = {3144, 3153};

  for(unsigned int i=0; i<truth.size(); i++){
    truth[i]->SetLineWidth(3);
    truth[i]->SetLineColor(color[i]);
    truth[i]->SetLineStyle(style[i]);
    truth[i]->SetMarkerStyle(0);
    truth[i]->SetFillColor(fillcolor[i]);
    // truth[i]->SetFillStyle(fillstyle[i]);
    if(i==0) truth[i]->Draw("E2 SAME");
    truth2[i]->SetLineWidth(3);
    truth2[i]->SetLineColor(color[i]);
    truth2[i]->SetLineStyle(style[i]);
    truth2[i]->Draw("HIST SAME");
  }
  stat->Draw("E1 SAME");
  output->Draw("E1 SAME");

  TLegend *l=new TLegend(0.56,0.65,0.78,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"Data","ple");
  for(unsigned int i=0; i<truth.size(); i++){
    // l->AddEntry(truth[i],legnames[i],"l");
    if(i==0) l->AddEntry(truth[i],legnames[i],"fl");
    else     l->AddEntry(truth[i],legnames[i],"l");
  }
  l->SetTextSize(0.04);
  l->Draw();
  CMSLabel();
  LumiInfo();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}


/*
██████  ██    ██ ████████ ██████  ██    ██ ████████     ███████ ███    ███ ███████  █████  ██████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██      ████  ████ ██      ██   ██ ██   ██
██    ██ ██    ██    ██    ██████  ██    ██    ██        ███████ ██ ████ ██ █████   ███████ ██████
██    ██ ██    ██    ██    ██      ██    ██    ██             ██ ██  ██  ██ ██      ██   ██ ██   ██
██████   ██████     ██    ██       ██████     ██        ███████ ██      ██ ███████ ██   ██ ██   ██
*/

void plotter::draw_output_smear(std::vector<TH1*> output_, TH1D* truth_, TString file_name){
  // std::vector<double> sys = get_sys_errors();
  // TH1* output_sys = add_error_bar(output, sys);

  std::vector<TH1*> outputs;
  for(int i=0; i<output_.size(); i++){
    TH1* hist = (TH1*) output_[i]->Clone();
    outputs.push_back(hist);
  }
  TH1D* truth = (TH1D*) truth_->Clone("truth");

  TCanvas *c = new TCanvas("c","",600,600);
  double ymax;
  gPad->SetLeftMargin(0.19);

  if(truth->GetMaximum() > outputs[0]->GetMaximum()) ymax = 1.5 * truth->GetMaximum();
  else ymax = 1.5 * outputs[0]->GetMaximum();
  // TGaxis::SetMaxDigits(3);
  for(auto output: outputs){
    output->SetTitle(" ");
    output->GetYaxis()->SetRangeUser(0., ymax);
    output->GetXaxis()->SetTitle("m_{jet} [GeV]");
    output->GetYaxis()->SetTitle("events");
    output->GetYaxis()->SetTitleOffset(1.5);
    output->GetXaxis()->SetTitleOffset(0.9);
    output->GetYaxis()->SetTitleSize(0.05);
    output->GetXaxis()->SetTitleSize(0.05);
    output->GetYaxis()->SetNdivisions(505);
    output->SetLineColor(kBlack);
    output->SetMarkerColor(kBlack);
    output->SetMarkerStyle(8);
    output->SetMarkerSize(1);
    output->Draw("E1 SAME");
  }
  gStyle->SetEndErrorSize(5);
  truth->SetLineWidth(3);
  truth->SetLineColor(kRed);
  truth->SetLineStyle(2);
  truth->Draw("HIST SAME");
  TLegend *l;
  if(truth->GetSize() > 20) l=new TLegend(0.2,0.6,0.4,0.85);
  else                       l=new TLegend(0.5,0.6,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(outputs[0],"data unfolded","pl");
  l->AddEntry(truth,"MC particle level","pl");
  l->SetTextSize(0.04);
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}
/*
██████  ██    ██ ████████ ██████  ██    ██ ████████     ███████ ████████  █████  ████████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██         ██    ██   ██    ██
██    ██ ██    ██    ██    ██████  ██    ██    ██        ███████    ██    ███████    ██
██    ██ ██    ██    ██    ██      ██    ██    ██             ██    ██    ██   ██    ██
██████   ██████     ██    ██       ██████     ██        ███████    ██    ██   ██    ██
*/


void plotter::draw_output_stat(TH1* output_, TH1* stat_, TH1D* truth_, bool norm, TString file_name){
  // std::vector<double> sys = get_sys_errors();
  // TH1* output_sys = add_error_bar(output, sys);

  TH1* output = (TH1*) output_->Clone("output");
  TH1* stat = (TH1*) stat_->Clone("stat");
  TH1D* truth = (TH1D*) truth_->Clone("truth");

  TCanvas *c = new TCanvas("c","",600,600);
  double ymax;
  gPad->SetLeftMargin(0.19);

  if(truth->GetMaximum() > output->GetMaximum()) ymax = 1.5 * truth->GetMaximum();
  else ymax = 1.5 * output->GetMaximum();
  // TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.5);
  output->GetXaxis()->SetTitleOffset(0.9);
  output->GetYaxis()->SetTitleSize(0.05);
  output->GetXaxis()->SetTitleSize(0.05);
  output->GetYaxis()->SetNdivisions(505);
  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);
  output->Draw("E1");
  stat->SetLineColor(kBlack);
  stat->SetMarkerColor(kBlack);
  stat->SetMarkerStyle(8);
  stat->SetMarkerSize(1);
  gStyle->SetEndErrorSize(5);
  truth->SetLineWidth(3);
  truth->SetLineColor(kRed);
  truth->SetLineStyle(2);
  truth->Draw("HIST SAME");
  stat->Draw("E1 SAME");
  output->Draw("E1 SAME");

  TLegend *l=new TLegend(0.5,0.65,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"Data","ple");
  l->AddEntry(truth,"MC particle level","pl");
  l->SetTextSize(0.04);
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██    ██ ████████ ██████  ██    ██ ████████     ███    ███  █████  ███████ ███████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ████  ████ ██   ██ ██      ██
██    ██ ██    ██    ██    ██████  ██    ██    ██        ██ ████ ██ ███████ ███████ ███████
██    ██ ██    ██    ██    ██      ██    ██    ██        ██  ██  ██ ██   ██      ██      ██
██████   ██████     ██    ██       ██████     ██        ██      ██ ██   ██ ███████ ███████
*/

void plotter::draw_output_mass(TH1* output_,  TH1* stat_, std::vector<TH1D*> mtop_templates_, std::vector<bool> show, bool norm, TString file_name){

  TH1* output = (TH1*) output_->Clone("output");
  TH1* stat = (TH1*) stat_->Clone("stat");

  std::vector<TH1D*> mtop_templates, mtop_templates2;
  for(unsigned int i = 0; i < mtop_templates_.size(); i++){
    mtop_templates.push_back((TH1D*) mtop_templates_[i]->Clone(""));
    mtop_templates2.push_back((TH1D*) mtop_templates_[i]->Clone(""));
  }

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.19);

  double max = output->GetMaximum();
  for(unsigned int i = 0; i < mtop_templates.size(); i++){
    if(show[i]){
      double max_temp = mtop_templates[i]->GetMaximum();
      if(max_temp > max) max = max_temp;
    }
  }
  double ymax = 1.5 * max;

  // TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.5);
  output->GetXaxis()->SetTitleOffset(0.9);
  output->GetYaxis()->SetTitleSize(0.05);
  output->GetXaxis()->SetTitleSize(0.05);
  output->GetYaxis()->SetNdivisions(505);
  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);
  output->Draw("E1 SAME");
  stat->SetLineColor(kBlack);
  stat->SetMarkerColor(kBlack);
  stat->SetMarkerStyle(8);
  stat->SetMarkerSize(1);
  gStyle->SetEndErrorSize(5);


  Color_t lcols[] = {kRed-4, kRed-4, kRed-4, TColor::GetColor("#0059b3"), 14, 14, 14};
  Color_t fcols[] = {kRed-10, kRed-10, kRed-10, TColor::GetColor("#99ccff"), 16, 16, 16};
  Color_t style[] = {3, 3, 3, 1, 2, 2, 2};


  for(unsigned int i = 0; i < mtop_templates.size(); i++){
    mtop_templates[i]->SetLineColor(lcols[i]);
    mtop_templates[i]->SetLineStyle(style[i]);
    mtop_templates[i]->SetFillColor(fcols[i]);
    // mtop_templates[i]->SetFillStyle(3144);
    mtop_templates[i]->SetLineWidth(3);
    mtop_templates2[i]->SetLineColor(lcols[i]);
    mtop_templates2[i]->SetLineStyle(style[i]);
    mtop_templates2[i]->SetLineWidth(3);
    if(show[i]){
      mtop_templates[i]->Draw("E2 SAME");
      mtop_templates2[i]->Draw("HIST SAME");
    }
  }
  mtop_templates[3]->Draw("E2 SAME");
  mtop_templates2[3]->Draw("HIST SAME");


  stat->Draw("E1 SAME");
  output->Draw("E1 SAME"); // draw again to set markers in front
  TLegend *l=new TLegend(0.56,0.65,0.78,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"Data","ple");
  if(show[0]) l->AddEntry(mtop_templates[0],"m_{t} = 166.5 GeV","fl");
  if(show[1]) l->AddEntry(mtop_templates[1],"m_{t} = 169.5 GeV","fl");
  if(show[2]) l->AddEntry(mtop_templates[2],"m_{t} = 171.5 GeV","fl");
  if(show[3]) l->AddEntry(mtop_templates[3],"m_{t} = 172.5 GeV","fl");
  if(show[4]) l->AddEntry(mtop_templates[4],"m_{t} = 173.5 GeV","fl");
  if(show[5]) l->AddEntry(mtop_templates[5],"m_{t} = 175.5 GeV","fl");
  if(show[6]) l->AddEntry(mtop_templates[6],"m_{t} = 178.5 GeV","fl");
  l->SetTextSize(0.04);
  l->Draw();
  CMSLabel();
  LumiInfo();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}


/*
████████ ██   ██ ███████  ██████
   ██    ██   ██ ██      ██    ██
   ██    ███████ █████   ██    ██
   ██    ██   ██ ██      ██    ██
   ██    ██   ██ ███████  ██████
*/

void plotter::draw_theo(TH1D* truth_, std::vector<TH1D*> variations_, bool norm, TString file_name){

  TH1D* truth = (TH1D*) truth_->Clone("truth");

  std::vector<TH1D*> variations;
  for(auto var: variations_){
    TH1D* temp = (TH1D*) var->Clone();
    variations.push_back(temp);
  }
  if(norm){
    truth->Scale(1/truth->Integral());
    for(auto var: variations) var->Scale(1/var->Integral());
  }

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.19);
  double ymax = 1.5 * truth->GetMaximum();
  // TGaxis::SetMaxDigits(3);
  truth->SetTitle(" ");
  truth->GetYaxis()->SetRangeUser(0., ymax);
  truth->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) truth->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else truth->GetYaxis()->SetTitle("events");
  truth->GetYaxis()->SetTitleOffset(1.5);
  truth->GetXaxis()->SetTitleOffset(0.9);
  truth->GetYaxis()->SetTitleSize(0.05);
  truth->GetXaxis()->SetTitleSize(0.05);
  truth->GetYaxis()->SetNdivisions(505);
  truth->SetLineColor(kAzure+1);
  truth->SetLineWidth(3);
  truth->Draw("HIST");
  for(auto var: variations){
    var->SetLineWidth(3);
    var->SetLineColor(kRed);
    var->Draw("HIST SAME");
  }
  truth->Draw("HIST SAME");
  TLegend *l=new TLegend(0.5,0.65,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(truth,"nominal sample","pl");
  l->AddEntry(variations[0],"variations","pl");
  l->SetTextSize(0.04);
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██           ██████ ██    ██ ██████  ██    ██ ███████
██          ██      ██    ██ ██   ██ ██    ██ ██
██          ██      ██    ██ ██████  ██    ██ █████
██          ██      ██    ██ ██   ██  ██  ██  ██
███████      ██████  ██████  ██   ██   ████   ███████
*/


void plotter::draw_lcurve(TGraph *lcurve, double x1, double y1,  double x2, double y2, TString file_name){
  TCanvas *c = new TCanvas("c","",600,600);
  lcurve->SetLineColor(kCyan+2);
  lcurve->SetLineWidth(2);
  lcurve->SetTitle(" ;L_{X};L_{Y}");
  lcurve->Draw("AL");
  lcurve->GetXaxis()->SetRangeUser(1.0,  3.0);
  lcurve->GetYaxis()->SetRangeUser(0.0, 10.0);
  lcurve->Draw("AL");
  c->Update();
  TMarker *p1=new TMarker(x1,y1,20);
  p1->SetMarkerSize(1.3);
  p1->SetMarkerColor(kRed);
  p1->Draw();
  TMarker *p2=new TMarker(x2,y2,20);
  p2->SetMarkerSize(1.3);
  p2->SetMarkerColor(1);
  p2->Draw();
  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(p1,"L-Curve scan","pl");
  l->AddEntry(p2,"rho scan","pl");
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██   ██  ██████  ████████  █████  ██    ██
██   ██ ██   ██ ██    ██    ██    ██   ██ ██    ██
██████  ███████ ██    ██    ██    ███████ ██    ██
██   ██ ██   ██ ██    ██    ██    ██   ██ ██    ██
██   ██ ██   ██  ██████     ██    ██   ██  ██████
*/



void plotter::draw_rhotau(TSpline *rhotau, double tau, TString file_name){
  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  double logtau = TMath::Log10(tau);
  std::vector<double> Xpoint = {logtau};
  std::vector<double> Ypoint = {rhotau->Eval(logtau)};
  TGraph * point = new TGraph(1, &Xpoint[0], &Ypoint[0]);
  point->SetTitle(" ");
  point->GetXaxis()->SetTitle("log #tau");
  point->GetYaxis()->SetTitle("#rho(log #tau)");
  point->GetYaxis()->SetTitleOffset(1.5);
  point->GetXaxis()->SetNdivisions(505);
  point->GetYaxis()->SetNdivisions(505);
  point->SetMarkerStyle(20);
  point->SetMarkerSize(1.5);
  point->SetLineColor(1);
  point->Draw("AP");
  point->GetXaxis()->SetLimits(-4.0, -1.5);
  double ymin = Ypoint[0]*0.95;
  double ymax = 1.15;
  point->GetYaxis()->SetRangeUser(ymin, ymax);
  point->Draw("AP");
  c->Update();

  // gPad->SetLogx();

  rhotau->SetLineColor(kCyan+2);
  rhotau->SetLineWidth(2);
  rhotau->Draw("L SAME");
  point->Draw("P SAME");

  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(point,"final #tau value","pl");
  l->AddEntry(rhotau,"#tau scan","l");
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}


/*
██████  ██████   ██████       ██ ███████  ██████ ████████ ██  ██████  ███    ██
██   ██ ██   ██ ██    ██      ██ ██      ██         ██    ██ ██    ██ ████   ██
██████  ██████  ██    ██      ██ █████   ██         ██    ██ ██    ██ ██ ██  ██
██      ██   ██ ██    ██ ██   ██ ██      ██         ██    ██ ██    ██ ██  ██ ██
██      ██   ██  ██████   █████  ███████  ██████    ██    ██  ██████  ██   ████
*/


void plotter::draw_projection(TH1D* proj_, TH1D* compare_, TString file_name ){
  TH1D* proj = (TH1D*) proj_->Clone("proj");
  TH1D* compare = (TH1D*) compare_->Clone("compare");

  // check if distribution agrees with Matrix projection
  // since there are some numerical effects involved in the projection,
  // the warning is only printed if both distributions do not agree within 1%
  int nbins = proj->GetSize() - 2;
  for(int i=1; i<= nbins; i++){
    double nproj = proj->GetBinContent(i);
    double ncomp = compare->GetBinContent(i);
    if(nproj < ncomp * 0.99 && nproj > ncomp * 1.01){
      std::cout << "Projection and Distribution in " << file_name << " do not agree in bin " << i << std::endl;
    }
  }

  TCanvas *c= new TCanvas("Projection Gen","",600,600);
  gPad->SetLeftMargin(0.15);
  proj->SetLineColor(kAzure+7);
  proj->Draw("HIST");
  compare->SetLineColor(kRed);
  compare->SetLineStyle(7);
  compare->Draw("HIST SAME");
  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(proj,"Projection from matrix","l");
  l->AddEntry(compare,"MC Gen","pl");
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ███████  ██████
██   ██ ██      ██
██████  █████   ██
██   ██ ██      ██
██   ██ ███████  ██████
*/

void plotter::draw_rec(TH1D* data_, TH1D* sig_, TH1D* bgr_, TString file_name){
  TH1D* data = (TH1D*) data_->Clone("data");
  TH1D* sig = (TH1D*) sig_->Clone("sig");
  TH1D* bgr = (TH1D*) bgr_->Clone("bgr");

  TCanvas *c= new TCanvas("c","",1200,600);
  gPad->SetLeftMargin(0.15);
  sig->Add(bgr, 1.);
  sig->SetTitle(" ");
  sig->GetYaxis()->SetRangeUser(0., 250);
  sig->GetXaxis()->SetTitle("detector binning");
  sig->GetYaxis()->SetTitle("events");
  sig->GetYaxis()->SetTitleOffset(1.5);
  sig->GetYaxis()->SetNdivisions(505);
  sig->SetFillColor(810);
  sig->SetLineColor(810);
  sig->Draw("HIST");
  bgr->SetFillColor(kGray);
  bgr->SetLineColor(kBlack);
  bgr->SetFillStyle(1001);
  bgr->Draw("HIST SAME");
  data->SetLineColor(kBlack);
  data->SetLineColor(kBlack);
  data->SetLineStyle(1);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  data->Draw("E SAME");
  TLegend *l=new TLegend(0.2,0.7,0.4,0.88);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(data,"Data","pl");
  l->AddEntry(sig,"t#bar{t}","f");
  l->AddEntry(bgr,"Background","f");
  l->Draw();
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}
/*
██ ██████      ██   ██ ██ ███████ ████████
███ ██   ██     ██   ██ ██ ██         ██
██ ██   ██     ███████ ██ ███████    ██
██ ██   ██     ██   ██ ██      ██    ██
██ ██████      ██   ██ ██ ███████    ██
*/


void plotter::draw_1D_hist(TH1D* hist_, TString file_name){
  TH1D* hist = (TH1D*) hist_->Clone("hist");
  TCanvas *c= new TCanvas("Particle Level","",600,600);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
  hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetFillColor(810);
  hist->SetLineColor(810);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██ ███    ██ ██   ██ ██ ███████ ████████ ███████
██   ██ ██ ████   ██ ██   ██ ██ ██         ██    ██
██████  ██ ██ ██  ██ ███████ ██ ███████    ██    ███████
██   ██ ██ ██  ██ ██ ██   ██ ██      ██    ██         ██
██████  ██ ██   ████ ██   ██ ██ ███████    ██    ███████
*/


void plotter::draw_binhist(TGraph* hist_, TString file_name, double xmin, double xmax){
  TGraph* hist = (TGraph*) hist_->Clone();
  TCanvas *c= new TCanvas(" ","",600,600);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle(" ");
  // hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
  // hist->GetYaxis()->SetTitle("events");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetLineColor(kBlack);
  hist->SetMarkerColor(kBlack);
  hist->SetMarkerStyle(20);
  hist->Draw("AP");
  TF1*f1 = new TF1("f1","pol1",0,500);
  hist->Fit("f1","R");
  TF1* fit = hist->GetFunction("f1");
  fit->SetLineColor(kRed);
  fit->Draw("SAME");
  double y1 = fit->Eval(xmin);
  double y2 = fit->Eval(xmax);
  TLine* l1 = new TLine(169, y1, 176, y1);
  TLine* l2 = new TLine(169, y2, 176, y2);
  l1->SetLineWidth(3);
  l1->SetLineColor(13);
  l1->SetLineStyle(2);
  l2->SetLineWidth(3);
  l2->SetLineColor(13);
  l2->SetLineStyle(2);
  l1->Draw("SAME");
  l2->Draw("SAME");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ███████ ██   ████████  █████
██   ██ ██      ██      ██    ██   ██
██   ██ █████   ██      ██    ███████
██   ██ ██      ██      ██    ██   ██
██████  ███████ ███████ ██    ██   ██
*/


void plotter::draw_delta(TH1* hist_, TString file_name){
  TH1* hist = (TH1*) hist_->Clone("hist");
  TCanvas *c= new TCanvas("Particle Level","",600,600);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle(file_name);
  hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
  hist->GetYaxis()->SetTitle("#Delta events");
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetFillColor(810);
  hist->SetLineColor(810);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

void plotter::draw_delta_rel(TH1* hist_, TH1* result_, TString file_name){
  TH1* hist = (TH1*) hist_->Clone("hist");
  TH1* result = (TH1*) result_->Clone("result");
  int Nbins = result->GetXaxis()->GetNbins();
  std::ofstream outs(directory+file_name+".txt");
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());
  for(unsigned int i=1; i<=Nbins; i++){
    double percent = fabs(100*hist->GetBinContent(i)/result->GetBinContent(i));
    hist->SetBinContent(i, percent);
    std::cout << "bin " << i << ": " << percent << "%" << std::endl;
  }
  std::cout.rdbuf(coutbuf);
  TCanvas *c= new TCanvas("Particle Level","",600,600);
  gPad->SetLeftMargin(0.15);
  hist->SetTitle(file_name);
  hist->GetXaxis()->SetTitle("m_{jet} [GeV]");
  hist->GetYaxis()->SetTitle("relative uncertainty [%]");
  hist->GetYaxis()->SetRangeUser(0, 50);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetNdivisions(505);
  hist->SetFillColor(810);
  hist->SetLineColor(810);
  hist->Draw("HIST");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

void plotter::draw_delta_comparison( TH1* total_, TH1* stat_, std::vector<TH1*> MODEL_DELTA, std::vector<TString> UncertNames, TString file_name){
  TH1* total = (TH1*) total_->Clone();
  TH1* stat = (TH1*) stat_->Clone();
  std::vector<TH1*> delta;
  for(unsigned int i=0; i<MODEL_DELTA.size(); i++){
    delta.push_back( (TH1*) MODEL_DELTA[i]->Clone() );
  }
  // print out values in txt
  std::ofstream outs(directory+file_name+".txt");
  int Nbins = total->GetXaxis()->GetNbins();
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());
  std::cout << "total ";
  for(unsigned int bin=1; bin<=Nbins; bin++){
    double percent = total->GetBinContent(bin);
    std::cout << std::fixed << std::setprecision(2) << percent << " ";
  }
  std::cout << std::endl;
  std::cout << "stat ";
  for(unsigned int bin=1; bin<=Nbins; bin++){
    double percent = stat->GetBinContent(bin);
    std::cout << std::fixed << std::setprecision(2) << percent << " ";
  }
  std::cout << std::endl;
  for(unsigned int i=0; i<delta.size(); i++){
    std::cout << UncertNames[i] << " ";
    for(unsigned int bin=1; bin<=Nbins; bin++){
      double percent = delta[i]->GetBinContent(bin);
      std::cout << std::fixed << std::setprecision(2) << percent << " ";
    }
    std::cout << std::endl;
  }
  std::cout.rdbuf(coutbuf);
  //
  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  total->SetTitle("");
  total->GetXaxis()->SetTitle("m_{jet} [GeV]");
  total->GetYaxis()->SetTitle("relative uncertainty [%]");
  total->GetYaxis()->SetTitleOffset(1.5);
  total->GetYaxis()->SetNdivisions(505);
  total->GetYaxis()->SetRangeUser(0, 100);
  total->SetFillColor(13);
  total->SetFillStyle(3144);
  total->SetLineColor(13);
  total->SetMarkerStyle(-1);
  total->Draw("HIST");
  stat->SetLineColor(kBlack);
  stat->SetLineWidth(4);
  stat->SetMarkerStyle(0);
  stat->Draw("B SAME");



  Color_t col[] = {kRed-4, kAzure+7, kGreen, 798, kBlue, kOrange-3, kMagenta, 13, kAzure, 14, kRed+5, kGreen-8};
  int i=0;
  for(auto hist: delta){
    hist->SetLineColor(col[i]);
    hist->SetLineWidth(4);
    hist->SetMarkerStyle(0);
    hist->Draw("B SAME");
    i++;
  }

  // LEGEND
  TLegend *leg = new TLegend(0.4,0.6,0.88,0.88);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->AddEntry(total, "stat #oplus exp. sys", "f");
  leg->AddEntry(stat, "stat", "l");
  for(unsigned int i=0; i<delta.size(); i++){
    if      (UncertNames[i] == "stat")      leg->AddEntry(delta[i],"statistics","l");
    else if (UncertNames[i] == "b-tagging") leg->AddEntry(delta[i],"b tagging","l");
    else if (UncertNames[i] == "pile-up")   leg->AddEntry(delta[i],"pileup","l");
    else if (UncertNames[i] == "jec")       leg->AddEntry(delta[i],"jet energy scale","l");
    else if (UncertNames[i] == "jer")       leg->AddEntry(delta[i],"jet energy resolution","l");
    else if (UncertNames[i] == "cor")       leg->AddEntry(delta[i],"XCone jet correction","l");
    else if (UncertNames[i] == "MuTrigger") leg->AddEntry(delta[i],"muon trigger","l");
    else if (UncertNames[i] == "MuID")      leg->AddEntry(delta[i],"muon ID","l");
    else if (UncertNames[i] == "ElTrigger") leg->AddEntry(delta[i],"electron trigger","l");
    else if (UncertNames[i] == "ElID")      leg->AddEntry(delta[i],"electron ID","l");
    else if (UncertNames[i] == "ElReco")    leg->AddEntry(delta[i],"electron reconstruction","l");
    else                                    leg->AddEntry(delta[i],UncertNames[i],"l");
  }
  leg->Draw();

  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

void plotter::draw_delta_comparison_model( TH1* total_, TH1* stat_, std::vector<TH1*> MODEL_DELTA, std::vector<TString> UncertNames, TString file_name){
  TH1* total = (TH1*) total_->Clone();
  TH1* stat = (TH1*) stat_->Clone();
  std::vector<TH1*> delta;
  for(unsigned int i=0; i<MODEL_DELTA.size(); i++){
    delta.push_back( (TH1*) MODEL_DELTA[i]->Clone() );
  }
  // print out values in txt
  std::ofstream outs(directory+file_name+".txt");
  int Nbins = total->GetXaxis()->GetNbins();
  auto coutbuf = std::cout.rdbuf(outs.rdbuf());
  std::cout << "total ";
  for(unsigned int bin=1; bin<=Nbins; bin++){
    double percent = total->GetBinContent(bin);
    std::cout << std::fixed << std::setprecision(2) << percent << " ";
  }
  std::cout << std::endl;
  std::cout << "stat ";
  for(unsigned int bin=1; bin<=Nbins; bin++){
    double percent = stat->GetBinContent(bin);
    std::cout << std::fixed << std::setprecision(2) << percent << " ";
  }
  std::cout << std::endl;
  for(unsigned int i=0; i<delta.size(); i++){
    std::cout << UncertNames[i] << " ";
    for(unsigned int bin=1; bin<=Nbins; bin++){
      double percent = delta[i]->GetBinContent(bin);
      std::cout << std::fixed << std::setprecision(2) << percent << " ";
    }
    std::cout << std::endl;
  }
  std::cout.rdbuf(coutbuf);
  //
  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  total->SetTitle("");
  total->GetXaxis()->SetTitle("m_{jet} [GeV]");
  total->GetYaxis()->SetTitle("relative uncertainty [%]");
  total->GetYaxis()->SetTitleOffset(1.5);
  total->GetYaxis()->SetNdivisions(505);
  total->GetYaxis()->SetRangeUser(0, 100);
  total->SetFillColor(13);
  total->SetFillStyle(3144);
  total->SetLineColor(13);
  total->SetMarkerStyle(-1);
  total->Draw("HIST");
  stat->SetLineColor(kBlack);
  stat->SetLineWidth(4);
  stat->SetMarkerStyle(0);
  stat->Draw("B SAME");



  Color_t col[] = {kRed-4, kAzure+7, kGreen, 798, kBlue, kOrange-3, kMagenta, 13, kAzure, 14, kRed+5, kGreen-8};
  for(unsigned int i=0; i<delta.size(); i++){
    delta[i]->SetLineColor(col[i]);
    delta[i]->SetLineWidth(4);
    delta[i]->SetMarkerStyle(0);
    delta[i]->Draw("B SAME");
  }

  // LEGEND
  TLegend *leg = new TLegend(0.4,0.6,0.88,0.88);
  leg->SetFillStyle(0);
  leg->SetNColumns(2);
  leg->AddEntry(total, "stat #oplus model sys", "f");
  leg->AddEntry(stat, "stat", "l");
  for(unsigned int i=0; i<delta.size(); i++){
    if      (UncertNames[i] == "mass")      leg->AddEntry(delta[i],"choice of m_{t}","l");
    else if (UncertNames[i] == "stat")      leg->AddEntry(delta[i],"statistics","l");
    else if (UncertNames[i] == "hdamp")     leg->AddEntry(delta[i],"h_{damp}","l");
    else                                    leg->AddEntry(delta[i],UncertNames[i],"l");
  }
  leg->Draw();

  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}



/*
██████  ██    ██ ████████ ██████  ██    ██ ████████     ██████  ███████ ███████ ██    ██ ██████   ██████
██    ██ ██    ██    ██    ██   ██ ██    ██    ██        ██   ██ ██      ██      ██    ██ ██   ██ ██    ██
██    ██ ██    ██    ██    ██████  ██    ██    ██        ██████  ███████ █████   ██    ██ ██   ██ ██    ██
██    ██ ██    ██    ██    ██      ██    ██    ██        ██           ██ ██      ██    ██ ██   ██ ██    ██
██████   ██████     ██    ██       ██████     ██        ██      ███████ ███████  ██████  ██████   ██████
*/



void plotter::draw_output_pseudo(TH1* output_, TH1D* pseudotruth_, TH1D* mctruth_, bool norm, TString file_name){

  TH1* output = (TH1*) output_->Clone("output");
  TH1D* pseudotruth = (TH1D*) pseudotruth_->Clone("pseudotruth");
  TH1D* mctruth = (TH1D*) mctruth_->Clone("mctruth");

  double ymax_temp = 0;
  if(pseudotruth->GetMaximum() > ymax_temp) ymax_temp = pseudotruth->GetMaximum();
  if(mctruth->GetMaximum() > ymax_temp) ymax_temp = mctruth->GetMaximum();
  if(output->GetMaximum() > ymax_temp) ymax_temp = output->GetMaximum();
  double ymax = 1.5 * ymax_temp;


  pseudotruth->SetTitle(" ");
  pseudotruth->GetYaxis()->SetRangeUser(0., ymax);
  pseudotruth->GetXaxis()->SetTitle("m_{jet} [GeV]");
  if(norm) pseudotruth->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else     pseudotruth->GetYaxis()->SetTitle("events");
  pseudotruth->GetYaxis()->SetTitleOffset(1.1);
  pseudotruth->GetXaxis()->SetTitleOffset(0.9);
  pseudotruth->GetYaxis()->SetTitleSize(0.05);
  pseudotruth->GetXaxis()->SetTitleSize(0.05);
  pseudotruth->GetYaxis()->SetNdivisions(505);

  pseudotruth->SetLineWidth(4);
  pseudotruth->SetLineColor(kRed);
  mctruth->SetLineWidth(3);
  mctruth->SetLineStyle(2);
  mctruth->SetLineColor(kBlue);

  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);

  TCanvas *c= new TCanvas("Particle Level","",600,600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  pseudotruth->Draw("HIST SAME");
  mctruth->Draw("HIST SAME");
  output->Draw("E1 SAME");
  TLegend *l;
  if(mctruth->GetSize() > 20) l=new TLegend(0.2,0.6,0.4,0.85);
  else                       l=new TLegend(0.55,0.6,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"pseudo data","pl");
  l->AddEntry(pseudotruth,"pseudo data truth","pl");
  l->AddEntry(mctruth,"MC truth","pl");
  l->SetTextSize(0.04);

  l->Draw();
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████  ██    ██ ██████  ██ ████████ ██    ██
██   ██ ██    ██ ██   ██ ██    ██     ██  ██
██████  ██    ██ ██████  ██    ██      ████
██      ██    ██ ██   ██ ██    ██       ██
██       ██████  ██   ██ ██    ██       ██
*/

void plotter::draw_purity(TH1D* numerator_, TH1D* denominator_, TString file_name){
  TH1D* numerator = (TH1D*) numerator_->Clone("numerator");
  TH1D* denominator = (TH1D*) denominator_->Clone("denominator");

  TGraphAsymmErrors* purity = new TGraphAsymmErrors(numerator, denominator,"cl=0.683 b(1,1) mode");

  purity->SetTitle(" ");
  purity->GetXaxis()->SetTitle("m_{rec}");
  purity->GetYaxis()->SetTitle("purity");
  purity->GetYaxis()->SetRangeUser(0,1);

  purity->GetXaxis()->SetTitleSize(0.05);
  purity->GetYaxis()->SetTitleSize(0.05);
  purity->GetXaxis()->SetTitleOffset(0.9);
  purity->GetYaxis()->SetTitleOffset(0.8);
  purity->GetXaxis()->SetNdivisions(505);
  purity->GetYaxis()->SetNdivisions(505);

  purity->SetMarkerStyle(20);
  purity->SetMarkerSize(0.8);
  purity->SetLineColor(1);

  TCanvas *c= new TCanvas("Purity","",600,600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  purity->Draw("AP");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

void plotter::draw_stability(TH1D* numerator_, TH1D* denominator_, TString file_name){
  TH1D* numerator = (TH1D*) numerator_->Clone("numerator");
  TH1D* denominator = (TH1D*) denominator_->Clone("denominator");


  TGraphAsymmErrors* stability = new TGraphAsymmErrors(numerator, denominator,"cl=0.683 b(1,1) mode");

  stability->SetTitle(" ");
  stability->GetXaxis()->SetTitle("m_{gen}");
  stability->GetYaxis()->SetTitle("stability");
  stability->GetYaxis()->SetRangeUser(0,1);

  stability->GetXaxis()->SetTitleSize(0.05);
  stability->GetYaxis()->SetTitleSize(0.05);
  stability->GetXaxis()->SetTitleOffset(0.9);
  stability->GetYaxis()->SetTitleOffset(0.8);
  stability->GetXaxis()->SetNdivisions(505);
  stability->GetYaxis()->SetNdivisions(505);

  stability->SetMarkerStyle(20);
  stability->SetMarkerSize(0.8);
  stability->SetLineColor(1);

  TCanvas *c= new TCanvas("Purity","",600,600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  stability->Draw("AP");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██████ ██   ██ ██     ██████
██      ██   ██ ██          ██
██      ███████ ██      █████
██      ██   ██ ██     ██
██████ ██   ██ ██     ███████
*/

void plotter::draw_chi2(TF1 * fit_, std::vector<double> masses_, std::vector<double> chi2_, double mass, double uncert, TString file_name){
  TF1 * fit = (TF1*)fit_->Clone("fit");
  TVectorD masses(masses_.size());
  TVectorD chi2(chi2_.size());
  for(int i=0; i<masses_.size(); i++) masses[i] = masses_[i];
  for(int i=0; i<chi2_.size(); i++) chi2[i] = chi2_[i];

  TGraph* chi_hist = new TGraph(masses,chi2);
  TCanvas *c = new TCanvas("Chi2", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  // TGaxis::SetMaxDigits(3);
  chi_hist->SetTitle(" ");
  chi_hist->GetXaxis()->SetTitle("m_{t} [GeV]");
  chi_hist->GetYaxis()->SetTitle("#chi^{2}");
  chi_hist->GetYaxis()->SetTitleOffset(1.1);
  chi_hist->GetXaxis()->SetTitleOffset(0.9);
  chi_hist->GetYaxis()->SetTitleSize(0.05);
  chi_hist->GetXaxis()->SetTitleSize(0.05);
  chi_hist->GetXaxis()->SetNdivisions(505);
  chi_hist->GetYaxis()->SetNdivisions(505);
  chi_hist->SetMarkerStyle(20);
  chi_hist->SetMarkerSize(1.5);
  chi_hist->SetLineColor(1);
  chi_hist->Draw("AP");
  fit->Draw("SAME");

  // write extracted mass value into plot
  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextFont(43);
  text.SetTextSize(18);
  char mass_text[32];
  sprintf(mass_text, "%.5g", mass);
  char uncert_text[32];
  if(uncert < 1) sprintf(uncert_text, "%.3g", uncert);
  else           sprintf(uncert_text, "%.4g", uncert);
  TString masstext = "m_{t} = ";
  masstext += mass_text;
  masstext += " #pm ";
  masstext += uncert_text;
  text.DrawLatex(.4,.6, masstext);
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
  return;
}

/*
██████  ██  █████  ███████
██   ██ ██ ██   ██ ██
██████  ██ ███████ ███████
██   ██ ██ ██   ██      ██
██████  ██ ██   ██ ███████
*/

void plotter::draw_bias(TH1* output_, TH1D* truth_, TH1* bias_, TString file_name){
  TH1D* bias = (TH1D*) bias_->Clone("bias");
  TH1D* truth = (TH1D*) truth_->Clone("truth");
  TH1D* output = (TH1D*) output_->Clone("output");

  double ymax_temp = 0;
  if(truth->GetMaximum() > ymax_temp) ymax_temp = truth->GetMaximum();
  if(bias->GetMaximum() > ymax_temp) ymax_temp = bias->GetMaximum();
  if(output->GetMaximum() > ymax_temp) ymax_temp = output->GetMaximum();
  double ymax = 1.5 * ymax_temp;

  TCanvas *c= new TCanvas("Bias Distribution","",600,600);
  // TGaxis::SetMaxDigits(3);
  gPad->SetLeftMargin(0.15);
  truth->SetTitle(file_name);
  truth->GetYaxis()->SetRangeUser(0., ymax);
  truth->GetXaxis()->SetTitle("m_{jet} [GeV]");
  truth->GetYaxis()->SetTitle("events");
  truth->GetYaxis()->SetTitleOffset(1.1);
  truth->GetXaxis()->SetTitleOffset(0.9);
  truth->GetYaxis()->SetTitleSize(0.05);
  truth->GetXaxis()->SetTitleSize(0.05);
  truth->GetYaxis()->SetNdivisions(505);

  truth->SetLineWidth(4);
  truth->SetLineColor(kRed);
  bias->SetLineWidth(4);
  bias->SetLineColor(kAzure+7);
  output->SetLineColor(kBlack);
  output->SetMarkerColor(kBlack);
  output->SetMarkerStyle(8);
  output->SetMarkerSize(1);

  truth->Draw("HIST SAME");
  bias->Draw("HIST SAME");
  output->Draw("E1 SAME");
  TLegend *l=new TLegend(0.55,0.6,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"unfolded result","pl");
  l->AddEntry(truth,"truth","pl");
  l->AddEntry(bias,"bias distribution","pl");
  l->SetTextSize(0.04);
  l->Draw();
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
███████ ███    ███ ███████  █████  ██████
██      ████  ████ ██      ██   ██ ██   ██
███████ ██ ████ ██ █████   ███████ ██████
██ ██  ██  ██ ██      ██   ██ ██   ██
███████ ██      ██ ███████ ██   ██ ██   ██
*/



void plotter::draw_smearFit(TH1D* variation, TF1* fit_, TString file_name){
  TH1D* var = (TH1D*) variation->Clone();
  TF1* fit = (TF1*) fit_->Clone();
  TCanvas *c= new TCanvas("","",600,600);
  gPad->SetLeftMargin(0.15);
  double ymax = var->GetMaximum() * 1.5;
  var->GetYaxis()->SetRangeUser(0., ymax);
  var->GetXaxis()->SetTitle("bin content");
  var->GetYaxis()->SetTitle("events");
  var->GetYaxis()->SetTitleOffset(1.5);
  var->GetYaxis()->SetNdivisions(505);
  var->SetFillColor(kGray);
  var->SetLineColor(kBlack);
  var->Draw("HIST");
  fit->SetLineColor(kRed);
  fit->SetLineWidth(4);
  fit->Draw("SAME");
  gPad->RedrawAxis();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

/*
██   ██ ███████ ██      ██████  ███████ ██████
██   ██ ██      ██      ██   ██ ██      ██   ██
███████ █████   ██      ██████  █████   ██████
██   ██ ██      ██      ██      ██      ██   ██
██   ██ ███████ ███████ ██      ███████ ██   ██
*/


TH1* plotter::add_error_bar(TH1* hist, std::vector<double> errors){
  float stat_err;
  for(int i=1; i<= errors.size(); i++){
    stat_err = hist->GetBinError(i);
    hist->SetBinError(i, sqrt(errors[i-1] * errors[i-1] + stat_err * stat_err));
  }
}

TH1D* get_difference(TH1D* hist1_, TH1D* hist2_){
  TH1D* hist1 = (TH1D*) hist1_->Clone("hist1");
  TH1D* hist2 = (TH1D*) hist2_->Clone("hist2");
  int Nbins = hist1->GetSize() - 2;
  TH1D* h_diff;
  for(int i = 1; i<=Nbins; i++){
    double diff = std::abs(hist1->GetBinContent(i) - hist2->GetBinContent(i));
    h_diff->Fill(i, diff);
  }
  return h_diff;
}


void plotter::CMSLabel(){
  TString cmstext = "CMS";
  TLatex *text = new TLatex(3.5, 24, cmstext);
  text->SetNDC();
  text->SetTextAlign(13);
  text->SetX(0.25);
  text->SetTextFont(62);
  text->SetTextSize(0.07);
  text->SetY(0.83);
  text->Draw();
}

void plotter::LumiInfo(){
  double lumi = 35.9;
  TString infotext = TString::Format("%3.1f fb^{-1} (13 TeV)", lumi);
  TLatex *text1 = new TLatex(3.5, 24, infotext);
  text1->SetNDC();
  text1->SetTextAlign(33);
  text1->SetTextFont(42);
  text1->SetTextSize(0.04);
  text1->SetX(0.90);
  text1->SetY(0.945);
  text1->Draw();
}
