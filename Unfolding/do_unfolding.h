#include "unfold.h"
#include "plotter.h"
#include "sys_handler.h"
#include "chi2fit.h"
#include "chi2sys.h"
#include "normalise.h"
#include "modelsys.h"
#include "smearing.h"

TString directory;
TString input_file;
TString output_file;
TString binning_xml;


class Tools{
 public:
  TH2D* GetDiagonalCovMatrix(TH1D* hist);
};
