#include "sys_handler.h"

void sys_handler::write_syserror(std::string name, std::vector<TH1*> unfolded_varied, std::vector<TH1D*> truth){

  unsigned int N_bins = truth[0]->GetSize() - 2; //do not count underflow/overflow
  std::vector<double> errors;
  for(unsigned int j = 1; j <= N_bins; j++){
    double diff = 0;
    double diff_temp = 0;
    for(unsigned int i = 0; i < unfolded_varied.size(); i++){
      diff_temp = abs(unfolded_varied[i]->GetBinContent(j) - truth[i]->GetBinContent(j));
      if(diff_temp >= diff) diff = diff_temp;
    }
    errors.push_back(diff);
  }

  std::ofstream sys_file;
  sys_file.open(name);
  for(unsigned int i = 0; i < errors.size(); i++){
    sys_file << std::fixed << std::setprecision(10) << i+1 << " " 
	     << errors[i] << std::endl;
  }
  sys_file.close();
  return;
}

std::vector<double> sys_handler::get_sys_errors(){
  std::ifstream file("SYS_TEST.txt");
  int bin;
  double err;
  std::vector<double> errors;
  while(file >> bin >> err){
    errors.push_back(err);
  }
  file.close();
  return errors;
}
