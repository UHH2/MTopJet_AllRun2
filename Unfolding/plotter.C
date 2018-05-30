#include "plotter.h"


plotter::plotter(TString dir){
  directory = dir;
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
}

/*
███    ███  █████  ████████ ██████  ██ ██   ██
████  ████ ██   ██    ██    ██   ██ ██  ██ ██
██ ████ ██ ███████    ██    ██████  ██   ███
██  ██  ██ ██   ██    ██    ██   ██ ██  ██ ██
██      ██ ██   ██    ██    ██   ██ ██ ██   ██
*/



void plotter::draw_matrix(TH2* hist_, TString file_name, bool zlog){
  TH1* hist = (TH1*) hist_->Clone("hist");

  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  if(zlog) gPad->SetLogz();
  hist->SetTitle(" ");
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetTitle("generator binning");
  hist->GetYaxis()->SetTitle("detector binning");
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
  gPad->SetLeftMargin(0.15);

  if(truth->GetMaximum() > output->GetMaximum()) ymax = 1.1 * truth->GetMaximum();
  else ymax = 1.1 * output->GetMaximum();
  TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.1);
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
  TLegend *l=new TLegend(0.5,0.65,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"data unfolded","pl");
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
  gPad->SetLeftMargin(0.15);

  if(truth->GetMaximum() > output->GetMaximum()) ymax = 1.1 * truth->GetMaximum();
  else ymax = 1.1 * output->GetMaximum();
  TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.1);
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
  l->AddEntry(output,"data unfolded","pl");
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

void plotter::draw_output_mass(TH1* output_, std::vector<TH1D*> mtop_templates_, std::vector<bool> show, bool norm, TString file_name){

  TH1* output = (TH1*) output_->Clone("output");
  std::vector<TH1D*> mtop_templates;
  for(unsigned int i = 0; i < mtop_templates_.size(); i++){
    mtop_templates.push_back((TH1D*) mtop_templates_[i]->Clone(""));
  }

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);

  double max = output->GetMaximum();
  for(unsigned int i = 0; i < mtop_templates.size(); i++){
    if(show[i]){
      double max_temp = mtop_templates[i]->GetMaximum();
      if(max_temp > max) max = max_temp;
    }
  }
  double ymax = 1.1 * max;

  TGaxis::SetMaxDigits(3);
  output->SetTitle(" ");
  output->GetYaxis()->SetRangeUser(0., ymax);
  output->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  if(norm) output->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else output->GetYaxis()->SetTitle("events");
  output->GetYaxis()->SetTitleOffset(1.1);
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

  mtop_templates[0]->SetLineColor(kRed);
  mtop_templates[1]->SetLineColor(kGreen);
  mtop_templates[2]->SetLineColor(kOrange);
  mtop_templates[3]->SetLineColor(13);
  mtop_templates[4]->SetLineColor(kBlue);
  mtop_templates[5]->SetLineColor(kMagenta);
  mtop_templates[6]->SetLineColor(kAzure+7);

  for(unsigned int i = 0; i < mtop_templates.size(); i++){
    mtop_templates[i]->SetLineWidth(3);
    if(show[i]) mtop_templates[i]->Draw("HIST SAME");
  }
  output->Draw("E1 SAME"); // draw again to set markers in front
  TLegend *l=new TLegend(0.5,0.65,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(output,"data unfolded","pl");
  if(show[0]) l->AddEntry(mtop_templates[0],"m_{top} = 166.5 GeV","pl");
  if(show[1]) l->AddEntry(mtop_templates[1],"m_{top} = 169.5 GeV","pl");
  if(show[2]) l->AddEntry(mtop_templates[2],"m_{top} = 171.5 GeV","pl");
  if(show[3]) l->AddEntry(mtop_templates[3],"m_{top} = 172.5 GeV","pl");
  if(show[4]) l->AddEntry(mtop_templates[4],"m_{top} = 173.5 GeV","pl");
  if(show[5]) l->AddEntry(mtop_templates[5],"m_{top} = 175.5 GeV","pl");
  if(show[6]) l->AddEntry(mtop_templates[6],"m_{top} = 178.5 GeV","pl");
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


void plotter::draw_lcurve(TGraph *lcurve, double x1, double y1, TString file_name){
  TCanvas *c = new TCanvas("c","",600,600);
  lcurve->SetLineColor(kCyan+2);
  lcurve->SetLineWidth(2);
  lcurve->SetTitle("L curve;L_{X};L_{Y}");
  lcurve->Draw("AL");
  lcurve->GetXaxis()->SetRangeUser(1.9,2.2);
  lcurve->Draw("AL");
  c->Update();
  TMarker *p1=new TMarker(x1,y1,20);
  p1->SetMarkerSize(1.3);
  p1->SetMarkerColor(kRed);
  p1->Draw();
  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(p1,"L-Curve scan","pl");
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

  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  sig->Add(bgr, 1.);
  sig->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  sig->GetYaxis()->SetTitle("events");
  sig->GetYaxis()->SetTitleOffset(1.5);
  sig->GetYaxis()->SetNdivisions(505);
  sig->SetFillColor(810);
  sig->SetLineColor(810);
  sig->Draw("HIST SAME");
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
  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
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
  hist->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
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
  hist->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
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

void plotter::draw_delta_comparison( TH1* total_, std::vector<TH1*> MODEL_DELTA, std::vector<TString> UncertNames, TString file_name){
  TH1* total = (TH1*) total_->Clone();
  std::vector<TH1*> delta;
  for(unsigned int i=0; i<MODEL_DELTA.size(); i++){
    delta.push_back( (TH1*) MODEL_DELTA[i]->Clone() );
  }

  TCanvas *c= new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15);
  total->SetTitle("");
  total->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  total->GetYaxis()->SetTitle("relative uncertainty [%]");
  total->GetYaxis()->SetTitleOffset(1.5);
  total->GetYaxis()->SetNdivisions(505);
  total->GetYaxis()->SetRangeUser(0, 1.5*total->GetMaximum());
  total->SetFillColor(13);
  total->SetFillStyle(3144);
  total->SetLineColor(13);
  total->SetMarkerStyle(-1);
  total->Draw("HIST");

  Color_t col[] = {kRed-4, kAzure+7, kGreen, 798, kBlue, kOrange-3, 12};
  int i=0;
  for(auto hist: delta){
    gPad->SetLeftMargin(0.15);
    hist->SetLineColor(col[i]);
    hist->SetLineWidth(4);
    hist->SetMarkerStyle(0);
    hist->Draw("PX SAME");
    i++;
  }

  // LEGEND
  TLegend *leg = new TLegend(0.63,0.6,0.88,0.88);
  leg->AddEntry(total, "exp. sys combined", "f");
  for(unsigned int i=0; i<delta.size(); i++) leg->AddEntry(delta[i],UncertNames[i],"l");
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
  double ymax = 1.1 * ymax_temp;


  pseudotruth->SetTitle(" ");
  pseudotruth->GetYaxis()->SetRangeUser(0., ymax);
  pseudotruth->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
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
  TGaxis::SetMaxDigits(3);
  pseudotruth->Draw("HIST SAME");
  mctruth->Draw("HIST SAME");
  output->Draw("E1 SAME");
  TLegend *l=new TLegend(0.5,0.6,0.85,0.85);
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

  TH1D* purity = numerator; // just to set correct binning
  purity->Divide(numerator, denominator, 1., 1., "B");

  purity->SetTitle(" ");
  purity->GetXaxis()->SetTitle("m_{gen}");
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
  TGaxis::SetMaxDigits(3);
  purity->Draw("E1");
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

void plotter::draw_chi2(TF1 * fit_, std::vector<double> masses_, std::vector<double> chi2_, TString file_name){
  TF1 * fit = (TF1*)fit_->Clone("fit");
  TVectorD masses(masses_.size());
  TVectorD chi2(chi2_.size());
  for(int i=0; i<masses_.size(); i++) masses[i] = masses_[i];
  for(int i=0; i<chi2_.size(); i++) chi2[i] = chi2_[i];

  TGraph* chi_hist = new TGraph(masses,chi2);
  TCanvas *c = new TCanvas("Chi2", "", 600, 600);
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  chi_hist->SetTitle(" ");
  chi_hist->GetXaxis()->SetTitle("m_{top} [GeV]");
  chi_hist->GetYaxis()->SetTitle("#chi^{2}");
  chi_hist->GetYaxis()->SetTitleOffset(1.5);
  chi_hist->GetXaxis()->SetNdivisions(505);
  chi_hist->GetYaxis()->SetNdivisions(505);
  chi_hist->SetMarkerStyle(20);
  chi_hist->SetMarkerSize(1.5);
  chi_hist->SetLineColor(1);
  chi_hist->Draw("AP");
  fit->Draw("SAME");
  c->SaveAs(directory + file_name + ".pdf");
  return;
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
