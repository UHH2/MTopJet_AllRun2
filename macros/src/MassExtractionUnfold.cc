#include "../include/CentralInclude.h"

using namespace std;

double xmin = 0;
double xmax = 0;

// class to transform fit parameters into uncorrelated ones
class TransformParameters{
public:
  TransformParameters(TFitResultPtr fitfunction, TString oldformula, int n_param);
  TString getNewFunction();
private:
  TString newFormula;
};
////

// class to vary correction
class Variation{
public:
  Variation(vector<double> params, vector<double> errors_up, vector<double> errors_down, TString function);
  vector<TF1*> getVariedFits() { return VariedFits; };
  TGraph *CalculateEnvelope(TF1* central, string updown);

private:
  vector<TF1*> VariedFits;
};
////

vector<double> CalculateParameterUncertainty(TGraphErrors* data_, TString function, TFitResultPtr fitresult, vector<double> parameters, TString updown);

/*
███    ███  █████  ██ ███    ██
████  ████ ██   ██ ██ ████   ██
██ ████ ██ ███████ ██ ██ ██  ██
██  ██  ██ ██   ██ ██ ██  ██ ██
██      ██ ██   ██ ██ ██   ████
*/


int main(int argc, char* argv[]){

  // true masses
  vector<double> truth = {169.5, 171.5, 172.5, 173.5, 175.5};


  // read txt files and get measured masses with uncert.
  vector<double> measured_stat;
  vector<double> error_stat;
  vector<double> measured;
  vector<double> error;
  string directory = "/afs/desy.de/user/s/schwarzd/Plots/Unfolding/";
  vector<string> massdir = {"Pseudo1695", "Pseudo1715", "Pseudo1", "Pseudo1735", "Pseudo1755"};
  for(auto mdir: massdir){
    std::ifstream infile(directory+mdir+"/Mass.txt");
    double m, e, m2, e2;
    while (infile >> m >> e >> m2 >> e2){
      cout << "stat only:        " << m << " +- " << e << endl;
      cout << "with mass uncert: " << m2 << " +- " << e2 << endl;
      measured_stat.push_back(m);
      error_stat.push_back(e);
      measured.push_back(m2);
      error.push_back(e2);
    }
  }

  TGraphErrors* masses_stat = new TGraphErrors(5, &truth[0], &measured_stat[0], 0, &error_stat[0]);
  TGraphErrors* masses = new TGraphErrors(5, &truth[0], &measured[0], 0, &error[0]);

  // diagonal line where measured = generated
  xmin = 167;
  xmax = 178;
  double ymin = xmin;
  double ymax = xmax;

  TLine *diag_line = new TLine(xmin, ymin, xmax, ymax);
  diag_line->SetLineColor(kRed);
  diag_line->SetLineWidth(3);
  diag_line->SetLineStyle(7);

  // fit line for mass calibration
  // TString fitformula = "[0] + [1]*x";
  // TF1* fit = new TF1("fit", fitformula, xmin, xmax);
  // TFitResultPtr fitresult = calib_gen->Fit(fit, "SR0Q");


  // now transform parameters to get uncorrelated ones
  // TransformParameters Trafo(fitresult, fitformula, 2);
  // TString fitformula2 = Trafo.getNewFunction();
  // TF1* fit2 = new TF1("fit2", fitformula2, xmin, xmax);

  // Use Parameters from previous transformation as start value
  // fit2->SetParameter(0, fitresult->Parameter(0));
  // fit2->SetParameter(1, fitresult->Parameter(1));
  // TFitResultPtr fitresult2 = calib_gen->Fit(fit2, "SR0Q");

  // Do Transformation Again
  // this has to be done because the numerical method does not deliver a
  // diagonal covariance matrix after one try
  // TransformParameters Trafo2(fitresult2, fitformula2, 2);
  // TString fitformula3 = Trafo2.getNewFunction();
  // TF1* fit3 = new TF1("fit3", fitformula3, xmin, xmax);
  // Use Parameters from previous transformation as start value
  // fit3->SetParameter(0, fitresult2->Parameter(0));
  // fit3->SetParameter(1, fitresult2->Parameter(1));
  // TFitResultPtr fitresult3 = calib_gen->Fit(fit3, "SR0Q");

  // vector<double> params, errors, errors_up, errors_down;
  // params.push_back(fitresult3->Parameter(0));
  // params.push_back(fitresult3->Parameter(1));
  // errors.push_back(fitresult3->ParError(0));
  // errors.push_back(fitresult3->ParError(1));
  // errors_up   = CalculateParameterUncertainty(calib_gen, fitformula3, fitresult3, params, "up");
  // errors_down = CalculateParameterUncertainty(calib_gen, fitformula3, fitresult3, params, "down");
  //
  //
  // Variation Variations(params, errors_up, errors_down, fitformula3);
  // TGraph* Area = Variations.CalculateEnvelope(fit3, "area");
  // TGraph* up = Variations.CalculateEnvelope(fit3, "up");
  // TGraph* down = Variations.CalculateEnvelope(fit3, "down");

  /*
  ██████  ██       ██████  ████████
  ██   ██ ██      ██    ██    ██
  ██████  ██      ██    ██    ██
  ██      ██      ██    ██    ██
  ██      ███████  ██████     ██
  */


  gStyle->SetOptStat(kFALSE);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetLegendBorderSize(0);

  TCanvas *A = new TCanvas("A", "A", 600, 600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  masses_stat->SetTitle(" ");
  masses_stat->GetXaxis()->SetLimits(xmin, xmax);
  masses_stat->GetHistogram()->SetMinimum(ymin);
  masses_stat->GetHistogram()->SetMaximum(ymax);
  masses_stat->GetXaxis()->SetTitle("m_{top} truth [GeV]");
  masses_stat->GetYaxis()->SetTitle("m_{top} measured [GeV]");
  masses_stat->GetYaxis()->SetTitleOffset(1.5);
  masses_stat->GetXaxis()->SetNdivisions(505);
  masses_stat->GetYaxis()->SetNdivisions(505);
  masses_stat->SetMarkerStyle(20);
  masses_stat->SetMarkerSize(1.3);
  masses_stat->SetLineColor(1);
  masses_stat->Draw("AP");
  // Area->SetFillColor(16);
  // Area->Draw("f SAME");
  // fit3->SetLineColor(kAzure+7);
  // fit3->SetLineWidth(3);
  // fit3->SetLineStyle(1);
  diag_line->Draw("SAME");
  // fit3->Draw("SAME");
  masses_stat->Draw("P SAME");
  // up->Draw("l SAME");
  // down->Draw("l SAME");
  TLegend *leg = new TLegend(0.6,0.15,0.88,0.38);
  leg->AddEntry(masses_stat, "extracted m_{top}^{MC} (stat only)", "l");
  leg->AddEntry(diag_line, "perfect measurement", "l");
  // leg->AddEntry(fit3, "calibration fit", "l");
  // leg->AddEntry(Area, "fit uncertainty", "f");
  leg->Draw();
  gPad->RedrawAxis();
  A->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/MassCalibratrion_stat.pdf");

  TCanvas *B = new TCanvas("B", "B", 600, 600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.15);
  TGaxis::SetMaxDigits(3);
  masses->SetTitle(" ");
  masses->GetXaxis()->SetLimits(xmin, xmax);
  masses->GetHistogram()->SetMinimum(ymin);
  masses->GetHistogram()->SetMaximum(ymax);
  masses->GetXaxis()->SetTitle("m_{top} truth [GeV]");
  masses->GetYaxis()->SetTitle("m_{top} measured [GeV]");
  masses->GetYaxis()->SetTitleOffset(1.5);
  masses->GetXaxis()->SetNdivisions(505);
  masses->GetYaxis()->SetNdivisions(505);
  masses->SetMarkerStyle(20);
  masses->SetMarkerSize(1.3);
  masses->SetLineColor(1);
  masses->Draw("AP");
  // Area->SetFillColor(16);
  // Area->Draw("f SAME");
  // fit3->SetLineColor(kAzure+7);
  // fit3->SetLineWidth(3);
  // fit3->SetLineStyle(1);
  diag_line->Draw("SAME");
  // fit3->Draw("SAME");
  masses->Draw("P SAME");
  // up->Draw("l SAME");
  // down->Draw("l SAME");
  TLegend *leg2 = new TLegend(0.6,0.15,0.88,0.38);
  leg2->AddEntry(masses, "extracted m_{top}^{MC}", "l");
  leg2->AddEntry(diag_line, "perfect measurement", "l");
  // leg->AddEntry(fit3, "calibration fit", "l");
  // leg->AddEntry(Area, "fit uncertainty", "f");
  leg2->Draw();
  gPad->RedrawAxis();
  B->SaveAs("/afs/desy.de/user/s/schwarzd/Plots/Unfolding/MassCalibratrion.pdf");

}



/*
████████ ██████   █████  ███    ██ ███████ ███████  ██████  ██████  ███    ███
   ██    ██   ██ ██   ██ ████   ██ ██      ██      ██    ██ ██   ██ ████  ████
   ██    ██████  ███████ ██ ██  ██ ███████ █████   ██    ██ ██████  ██ ████ ██
   ██    ██   ██ ██   ██ ██  ██ ██      ██ ██      ██    ██ ██   ██ ██  ██  ██
   ██    ██   ██ ██   ██ ██   ████ ███████ ██       ██████  ██   ██ ██      ██
*/



TransformParameters::TransformParameters(TFitResultPtr fitresult, TString oldformula, int n_param){
  //
  // To get uncorrelated Parameters from a fit, a matrix equation has to be solved.
  // From the regular fit, one gets a set of Parameters P = (a,b,c,...).
  // Now you have to find a Transformation O^T Cov(P) O = D, where D is a Diagonal Matrix with Eigenvalues
  // of Cov(P) on its diagonal.
  // The new (uncorrelated) Parameters are then P' = O^T P
  //

  // get Covariance Matrix
  TMatrixDSym oldCov = fitresult->GetCovarianceMatrix();
  // extract Eigenvectors
  const TMatrixDSymEigen oldEigen(oldCov);
  // get Transformation Matrix O and O^T
  const TMatrixD TrafoMatrix = oldEigen.GetEigenVectors();
  TMatrixD TranspMatrix(TrafoMatrix);
  TranspMatrix.T();

  // get old parameters
  vector<double> oldParams;
  for(int i=0; i<n_param; i++){
    oldParams.emplace_back(fitresult->Parameter(i));
  }
  const TVectorD V_oldParams(n_param,&oldParams[0]);

  // calculate new parameters in new basis
  TDecompSVD ParameterSVD(TrafoMatrix);
  Bool_t ok;
  const TVectorD V_newParams= ParameterSVD.Solve(V_oldParams, ok);

  // transform back new parameters into old ones
  // a = a (a',b',c',...)
  TDecompSVD ParameterSVD2(TranspMatrix);
  Bool_t ok2;
  const TVectorD V_old_from_newParams= ParameterSVD2.Solve(V_newParams, ok2);

  //parametrize old parameters by using new ones with Trafo Matrix
  TMatrixD M_Dummy(TrafoMatrix);

  // Fill vector of vectors: 1 vector for each old parameter
  // containing the values of the new parameters expressing the old one.
  vector<vector<double>> TransformationValues;
  for(int i=0; i<n_param; i++){
    vector<double> singleParam;
    for(int j=0; j<n_param; j++){
      singleParam.emplace_back(TMatrixDRow(M_Dummy,i)[j]);
    }
    TransformationValues.emplace_back(singleParam);
    for(unsigned int j=0; j<singleParam.size(); j++){
      singleParam.pop_back();
    }
  }

  // now write new fit function with new parameters
  oldformula;
  for(int i=0; i<n_param; i++){
    TString idx = "[";
    idx += i;
    idx += "]";
    TString new_expr = "( [0] *(";
    new_expr += TransformationValues[i][0];
    new_expr += ")";
    for(int j=1; j<n_param; j++){
      TString newidx = "";
      newidx += j;
      new_expr += " + ["+newidx+"] *(";
      new_expr += TransformationValues[i][j];
      new_expr += ")";
    }
    new_expr += ")";
    oldformula.ReplaceAll(idx,new_expr);
  }

  newFormula = oldformula;
}

// Some Get Functions
TString TransformParameters::getNewFunction(){
  return newFormula;
}

/*
██    ██  █████  ██████  ██  █████  ████████ ██  ██████  ███    ██
██    ██ ██   ██ ██   ██ ██ ██   ██    ██    ██ ██    ██ ████   ██
██    ██ ███████ ██████  ██ ███████    ██    ██ ██    ██ ██ ██  ██
 ██  ██  ██   ██ ██   ██ ██ ██   ██    ██    ██ ██    ██ ██  ██ ██
  ████   ██   ██ ██   ██ ██ ██   ██    ██    ██  ██████  ██   ████
*/



Variation::Variation(vector<double> params, vector<double> errors_up, vector<double> errors_down, TString function){
  vector<double> pup;
  vector<double> pdown;

  int Nparams = params.size();

  for(unsigned int i=0; i<Nparams; i++){
    pup.push_back(params[i] + sqrt(errors_up[i]*errors_up[i]));
    pdown.push_back(params[i] - sqrt(errors_down[i]*errors_down[i]));
  }

  // Only one parameter of the function is either set to 'up' or 'down'.
  // Therefore, for 4 parameters there are 4 functions for 'up' and 4 functions for 'down' variations.
  int Nvariations = Nparams*2;

  // now get a function for every possible variation
  // first, iterate through 'up' variations (first half of variations)
  // then, iterate through 'down' variations (second half)
  // the i==j statement makes sure to iterate the variations through every parameter
  for(unsigned int i=0; i<Nvariations; i++){
    TF1* f = new TF1("f",function,xmin,xmax);
    for(unsigned int j=0; j<Nparams; j++){
      if(i<Nparams){
        if(i==j) f->SetParameter(j, pup[j]);
        else     f->SetParameter(j, params[j]);
      }
      else{
        if(i==j) f->SetParameter(j, pdown[j]);
        else     f->SetParameter(j, params[j]);
      }
    }
    VariedFits.push_back(f);
  }
}


TGraph *Variation::CalculateEnvelope(TF1* central, string updown){

  double c;
  double v;
  double delta;
  const int steps = 100;
  double stepsize = (xmax - xmin)/steps;
  double x;
  TVectorD up(steps);
  TVectorD down(steps);
  TVectorD xvalues(steps);

  for(int i=0; i<steps; i++){
    x = xmin + stepsize * i;
    c = central->Eval(x);
    double down2 = 0;
    double up2 = 0;
    for(unsigned int j=0; j<VariedFits.size(); j++){
      v = VariedFits[j]->Eval(x);
      delta = c-v;
      if(delta < 0) delta *= -1;
      if(v < c) down2 += pow(delta,2);
      if(v >= c) up2 +=  pow(delta,2);
    }
    up[i] = c + sqrt(up2);
    down[i] = c - sqrt(down2);
    xvalues[i] = x;
  }

  TGraph *upgraph = new TGraph(xvalues, up);
  TGraph *downgraph = new TGraph(xvalues, down);

  TGraph *area = new TGraph(2*steps);
  for (int i=0;i<steps;i++) {
    area->SetPoint(i,xvalues[i],up[i]);
    area->SetPoint(steps+i,xvalues[steps-i-1],down[steps-i-1]);
  }

  if(updown == "up") return upgraph;
  if(updown == "down") return downgraph;
  else return area;

}

/*
██    ██  █████  ██████  ██    ██     ███████ ██████  ██████   ██████  ██████  ███████
██    ██ ██   ██ ██   ██  ██  ██      ██      ██   ██ ██   ██ ██    ██ ██   ██ ██
██    ██ ███████ ██████    ████       █████   ██████  ██████  ██    ██ ██████  ███████
 ██  ██  ██   ██ ██   ██    ██        ██      ██   ██ ██   ██ ██    ██ ██   ██      ██
  ████   ██   ██ ██   ██    ██        ███████ ██   ██ ██   ██  ██████  ██   ██ ███████
*/


// vary one parameter until the original chi2 grows by 1
// the delta to the original parameter is considered the 1sigma uncertainty
vector<double> CalculateParameterUncertainty(TGraphErrors* data_, TString function, TFitResultPtr fitresult, vector<double> parameters, TString updown){
  TGraphErrors* data = (TGraphErrors*) data_->Clone();
  TF1 *fit = new TF1("fit", function, xmin, xmax);
  double originalChi2 = fitresult->Chi2();
  vector<double> uncertainties;
  for(unsigned int i=0; i<parameters.size(); i++){
    int i_step = 0;
    double stepsize = 0.00001;
    double delta = 0;
    double diff = 0;
    bool foundUncert = false;
    double newChi2 = originalChi2;
    while(!foundUncert){
      i_step++;
      if(updown == "up")        delta =  i_step*sqrt(parameters[i]*parameters[i])*stepsize;
      else if(updown == "down") delta = -i_step*sqrt(parameters[i]*parameters[i])*stepsize;
      else throw std::invalid_argument("Select up or down in CalculateParameterUncertainty");
      for(unsigned int j=0; j<parameters.size(); j++){
        if(i != j) fit->FixParameter(j, parameters[j]);
        else       fit->FixParameter(j, parameters[j]+delta);
      }
      TFitResultPtr newresult = data->Fit(fit,"QSR");
      newChi2 = newresult->Chi2();
      diff = newChi2 - originalChi2;
      if(diff >= 1) foundUncert = true;
      if(i_step > 50000){
        cout << "had to break parameter uncertainty calculation after 50000 steps!!!" << endl;
        foundUncert = true;
      }
    }
    cout << "It took " << i_step << " steps to find uncertainty of parameter " << i
    << " resulting in a Chi2 difference of " << diff << endl;
    // cout << "Chi2 before = " << originalChi2 << " Chi2 after = " << newChi2 << endl;
    // cout << "Parameter before = " << parameters[i] << " Parameter after = " << parameters[i]+delta << endl;
    uncertainties.push_back(delta);
  }
  return uncertainties;
}
