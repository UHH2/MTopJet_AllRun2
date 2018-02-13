#include <iostream>
#include <fstream>
#include <vector>
#include <TH1.h>


class sys_handler{

 public:
  void write_syserror(std::string, std::vector<TH1*>, std::vector<TH1D*>);
  std::vector<double> get_sys_errors();

 private:

};
