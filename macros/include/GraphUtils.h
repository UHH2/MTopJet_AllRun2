#pragma once
#include <TCanvas.h>
#include <TColor.h>
#include <TDecompSVD.h>
#include <TDecompLU.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TGraph2DErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrix.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TFitResult.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualFitter.h>
#include <TVectorD.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#include <sstream>
#include <string>
#include <tuple>

#define RESET   "\033[0m"
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define BLUE    "\033[34m"      /* Blue */

using namespace std;
using namespace std::chrono;

/*
██████  ██████
.    ██ ██   ██
.█████  ██   ██
██      ██   ██
███████ ██████
*/

// #################################################################################################
// Numeric solution for Z values of 2D plot ########################################################

vector<vector<double>> FindXY(TF2 *function, double zfix, double xmin, double xmax, double ymin, double ymax, int steps = 1000, double accuracy = 0, bool count=false){
  vector<double> points;
  vector<vector<double>> all_points;
  double dx = (xmax-xmin)/steps;
  double dy = (ymax-ymin)/steps;
  double x, y, z;
  double n_wide = 100;
  double wide_steps = steps/n_wide;

  // First step is to decrease the area to increase the speed of the algorithm #######################
  // x_limit -----------------------------------------------------------------------------------------
  bool reached_x_limit_up   = false;
  bool reached_x_limit_down = false;
  int step_limit_x_up   = 0;
  int step_limit_x_down = 0;

  for(int xpar=0; xpar<n_wide+1; xpar++){
    x=xmin+wide_steps*xpar*dx;
    if(x>xmax) break;
    else{
      for(int ypar=0; ypar<n_wide+1; ypar++){
        if(reached_x_limit_down && reached_x_limit_up && ypar==0) break;
        y=ymin+wide_steps*ypar*dy;
        if(y>ymax) break;
        z=function->Eval(x, y);

        if((z+accuracy)>zfix){
          if(reached_x_limit_down) reached_x_limit_up = true;
          if((!reached_x_limit_down) && ypar==n_wide) step_limit_x_down = wide_steps*xpar;
          if(reached_x_limit_down && reached_x_limit_up && ypar==n_wide)    step_limit_x_up   = wide_steps*xpar;
        }
        else{
          if(!reached_x_limit_down) reached_x_limit_down=true;
          if(reached_x_limit_down)  reached_x_limit_up   = false;
          break;
        }
      }
    }
  }
  // cout << "\nup: " << step_limit_x_up << "   down: " << step_limit_x_down << endl;

  // y_limit -----------------------------------------------------------------------------------------
  bool reached_y_limit_up   = false;
  bool reached_y_limit_down = false;
  int step_limit_y_up   = 0;
  int step_limit_y_down = 0;

  for(int ypar=0; ypar<n_wide+1; ypar++){
    y=ymin+wide_steps*ypar*dy;
    if(y>ymax) break;
    else{
      for(int xpar=0; xpar<n_wide+1; xpar++){
        if(reached_y_limit_down && reached_y_limit_up && xpar==0) break;
        x=xmin+wide_steps*xpar*dx;
        if(x>xmax) break;
        z=function->Eval(x, y);

        if((z+accuracy)>zfix){
          if(reached_y_limit_down) reached_y_limit_up = true;
          if(!reached_y_limit_down && xpar==n_wide)  step_limit_y_down = wide_steps*ypar;
          if(reached_y_limit_down  && reached_y_limit_up && xpar==n_wide)  step_limit_y_up   = wide_steps*ypar;
        }
        else{
          if(!reached_y_limit_down) reached_y_limit_down = true;
          if(reached_y_limit_down)  reached_y_limit_up   = false;
          break;
        }
      }
    }
  }
  // cout << "up: " << step_limit_y_up << "   down: " << step_limit_y_down << endl;

  /* First step is to decrease the area to increase the speed of the algorithm*/
  // For 2000 -> calculation time decreases!
  // 2463774ms   0
  // 1221039ms  10
  //  994403ms  20
  //  897007ms  50
  //  827621ms 100

  if(count){
    cout << "\nxup: " << step_limit_x_up << " | xdown: " << step_limit_x_down << endl;
    cout << "yup: " << step_limit_y_up << " | ydown: " << step_limit_y_down << endl;
  }

  auto start = high_resolution_clock::now(); // Calculation time - start
  for(int xpar=step_limit_x_down; xpar<step_limit_x_up; xpar++){
    if(count&&(xpar%1000==0)){
      auto stop = high_resolution_clock::now();  // Calculation time - stop
      auto duration = duration_cast<milliseconds>(stop - start);
      cout << "Passed x = " << xpar << " (" << duration.count()/1000 << "s)" << endl;
    }
    int count_zfix_reached=0;
    x=xmin+xpar*dx;
    for(int ypar=step_limit_y_down; ypar<step_limit_y_up; ypar++){
      y=ymin+ypar*dy;
      z=function->Eval(x, y);
      if(abs(z-zfix)<accuracy){
        count_zfix_reached++;
        points.push_back(x);
        points.push_back(y);
        points.push_back(z);
        all_points.push_back(points);
        points={};
      }
      // if(count_zfix_reached==2) break;
    }
  }

  return all_points;
}


/*
██████   ██████  ██ ███    ██ ████████
██   ██ ██    ██ ██ ████   ██    ██
██████  ██    ██ ██ ██ ██  ██    ██
██      ██    ██ ██ ██  ██ ██    ██
██       ██████  ██ ██   ████    ██
*/

double ReturnIndex_High(vector<double> vec){
  unsigned int d = vec.size();
  bool isMax=false;
  double index;
  for(unsigned int point=0; point<d; point++){
    for(unsigned int point1=0; point1<d; point1++){
      if(point!=point1){
        bool bigger = vec[point]>vec[point1];
        if(!bigger){
          isMax=false;
          break; // to avoid point that are bigger then the last point
        }
        else isMax=true;
      }
    }
    if(isMax) index=point;
  }
  return index;
}


double ReturnIndex_Low(vector<double> vec){
  unsigned int d = vec.size();
  bool isMin=false;
  double index;
  for(unsigned int point=0; point<d; point++){
    for(unsigned int point1=0; point1<d; point1++){
      if(point!=point1){
        bool smaller = vec[point]<vec[point1];
        if(!smaller){
          isMin=false;
          break; // to avoid point that are bigger then the last point
        }
        else isMin=true;
      }
    }
    if(isMin) index=point;
  }
  return index;
}
