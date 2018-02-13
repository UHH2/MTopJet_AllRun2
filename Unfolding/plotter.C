#include "plotter.h"


plotter::plotter(TString dir){
  directory = dir;
  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
}

void plotter::draw_matrix(TH2* hist, TString file_name, bool zlog){
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

void plotter::draw_output(TH1* output, TH1D* truth, bool norm, TString file_name){
  // std::vector<double> sys = get_sys_errors();
  // TH1* output_sys = add_error_bar(output, sys);

  TCanvas *c = new TCanvas("c","",600,600);
  double ymax;
  gPad->SetLeftMargin(0.15); 
  if(norm){
    double norm_out = 1/output->Integral();
    double norm_truth = 1/truth->Integral();
    output->Scale(norm_out, "width");
    truth->Scale(norm_truth, "width");
  }
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

void plotter::draw_output_mass(TH1* output, std::vector<TH1D*> mtop_templates, std::vector<bool> show, bool norm, TString file_name){
  // std::vector<double> sys = get_sys_errors();
  // TH1* output_sys = add_error_bar(output, sys);

  TCanvas *c = new TCanvas("c","",600,600);
  gPad->SetLeftMargin(0.15); 
  if(norm){
    double norm_out = 1/output->Integral();
    output->Scale(norm_out, "width");
    for(unsigned int i = 0; i < mtop_templates.size(); i++){
      double norm_truth = 1/mtop_templates[i]->Integral();
      mtop_templates[i]->Scale(norm_truth, "width");
    }
  }
  double ymax;
  ymax = 1.1 * output->GetMaximum();
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
  mtop_templates[0]->SetLineColor(kBlue);
  mtop_templates[1]->SetLineColor(kRed);
  mtop_templates[2]->SetLineColor(kGreen-1);
  mtop_templates[3]->SetLineColor(kAzure+7);
  mtop_templates[4]->SetLineColor(kRed-2);
  mtop_templates[5]->SetLineColor(798);

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
  if(show[3]) l->AddEntry(mtop_templates[3],"m_{top} = 173.5 GeV","pl");
  if(show[4]) l->AddEntry(mtop_templates[4],"m_{top} = 175.5 GeV","pl");
  if(show[5]) l->AddEntry(mtop_templates[5],"m_{top} = 178.5 GeV","pl");
  l->SetTextSize(0.04);
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}

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

void plotter::draw_projection(TH1D* proj, TH1D* compare, TString file_name ){ 
  TCanvas *c= new TCanvas("Projection Gen","",600,600);
  gPad->SetLeftMargin(0.15); 
  proj->Draw("HIST");
  compare->SetMarkerColor(kBlack);
  compare->SetMarkerStyle(8);
  compare->SetMarkerSize(1); 
  compare->Draw("E1 SAME");
  TLegend *l=new TLegend(0.55,0.65,0.85,0.8);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->AddEntry(proj,"Projection from matrix","l");
  l->AddEntry(compare,"MC Gen","pl");
  l->Draw();
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}


void plotter::draw_rec(TH1D* data, TH1D* sig, TH1D* bgr, TString file_name){
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

  void plotter::draw_1D_hist(TH1D* hist, TString file_name){
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



void plotter::draw_output_pseudo(std::vector<TH1*> output,std::vector<TH1D*> truth, std::vector<TH1D*> mc_truth, bool norm, TString file_name){

  if(norm){
    for(unsigned int i=0; i<truth.size(); i++){
      double norm_truth = 1/truth[i]->Integral();
      truth[i]->Scale(norm_truth, "width");
    }
    for(unsigned int i=0; i<mc_truth.size(); i++){
      double norm_mctruth = 1/mc_truth[i]->Integral();
      mc_truth[i]->Scale(norm_mctruth, "width");
    }
    for(unsigned int i=0; i<output.size(); i++){
      double norm_out = 1/output[i]->Integral();
      output[i]->Scale(norm_out, "width");
    }
  }

  double ymax_temp = 0;
  for(unsigned int i=0; i<truth.size(); i++){
    if(truth[i]->GetMaximum() > ymax_temp) ymax_temp = truth[i]->GetMaximum();
  }
  for(unsigned int i=0; i<mc_truth.size(); i++){
    if(mc_truth[i]->GetMaximum() > ymax_temp) ymax_temp = mc_truth[i]->GetMaximum();
  }
  for(unsigned int i=0; i<output.size(); i++){
    if(output[i]->GetMaximum() > ymax_temp) ymax_temp = output[i]->GetMaximum();
  }
  double ymax = 1.1 * ymax_temp;


  truth[0]->SetTitle(" ");
  truth[0]->GetYaxis()->SetRangeUser(0., ymax);
  truth[0]->GetXaxis()->SetTitle("Leading Jet Mass [GeV]");
  if(norm) truth[0]->GetYaxis()->SetTitle("#frac{1}{#sigma} #frac{d#sigma}{dm_{jet}} [#frac{1}{GeV}]");
  else     truth[0]->GetYaxis()->SetTitle("events");
  truth[0]->GetYaxis()->SetTitleOffset(1.1);
  truth[0]->GetXaxis()->SetTitleOffset(0.9);
  truth[0]->GetYaxis()->SetTitleSize(0.05);
  truth[0]->GetXaxis()->SetTitleSize(0.05);
  truth[0]->GetYaxis()->SetNdivisions(505);

  for(unsigned int i=0; i<truth.size(); i++){
    truth[i]->SetLineWidth(4);
    truth[i]->SetLineColor(kRed);
    mc_truth[i]->SetLineWidth(3);
    mc_truth[i]->SetLineStyle(2);
    mc_truth[i]->SetLineColor(kBlue);
  }
  int N_bins = output[0]->GetXaxis()->FindBin(5000) - 1;

  std::vector<TGraphErrors*> output_new;

  double x[N_bins];
  double y[N_bins];
  double err[N_bins];

  for(unsigned int i=0; i<output.size(); i++){
    TAxis* xaxis = output[i]->GetXaxis();
    for(int j=1; j<=N_bins; j++){
      x[j-1] = xaxis->GetBinCenter(j) + (- 0.4 + i*0.2666667) * xaxis->GetBinWidth(j);
      y[j-1] = output[i]->GetBinContent(j);
      err[j-1] = output[i]->GetBinError(j);
    }
    TGraphErrors* output_temp = new TGraphErrors(N_bins, x, y, 0, err);
    output_new.push_back(output_temp);
  }

  for(unsigned int i=0; i<output_new.size(); i++){
    output_new[i]->SetMarkerColor(kBlack);
    output_new[i]->SetMarkerStyle(20 + i);
    if(i==3)output_new[i]->SetMarkerStyle(33);
  }


  TCanvas *c= new TCanvas("Particle Level","",600,600);
  gPad->SetLeftMargin(0.15); 
  gStyle->SetErrorX(0);
  TGaxis::SetMaxDigits(3);   
  for(unsigned int i=0; i<truth.size(); i++) truth[i]->Draw("HIST SAME");
  for(unsigned int i=0; i<mc_truth.size(); i++) mc_truth[i]->Draw("HIST SAME");
  for(unsigned int i=0; i<output_new.size(); i++) output_new[i]->Draw("P");
  TLegend *l=new TLegend(0.5,0.6,0.85,0.85);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  for(unsigned int i=0; i<output.size(); i++){
    TString label = "pseudo data " + std::to_string(i+1);
    l->AddEntry(output_new[i],label,"pl");
  }
  l->AddEntry(truth[0],"pseudo data truth","pl");
  l->AddEntry(mc_truth[0],"pseudo mc truth","pl");
  l->SetTextSize(0.04);

  l->Draw();
  gPad->RedrawAxis(); 
  c->SaveAs(directory + file_name + ".pdf");
  delete c;
}



void plotter::draw_purity(TH1D* numerator, TH1D* denominator, TString file_name){

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


TH1* plotter::add_error_bar(TH1* hist, std::vector<double> errors){
  float stat_err;
  for(int i=1; i<= errors.size(); i++){
    stat_err = hist->GetBinError(i);
    hist->SetBinError(i, sqrt(errors[i-1] * errors[i-1] + stat_err * stat_err));
  }
}

TH1D* get_difference(TH1D* hist1, TH1D* hist2){
  int Nbins = hist1->GetSize() - 2;
  TH1D* h_diff;
  for(int i = 1; i<=Nbins; i++){
    double diff = std::abs(hist1->GetBinContent(i) - hist2->GetBinContent(i));
    h_diff->Fill(i, diff);
  }
  return h_diff;
}
