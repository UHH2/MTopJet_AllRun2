#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/MTopJet/include/ModuleBASE.h>

using namespace std;

class SplitSampleModule : public ModuleBASE {

 public:
  explicit SplitSampleModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

 protected:
  int counter;
  int divide_by;
  int sample_nr;

};

SplitSampleModule::SplitSampleModule(uhh2::Context& ctx){
  divide_by = 4; // this variable selects in how many samples you divide the original sample
  sample_nr = 1; // select which of the samples you keep
  counter = 1;
}

bool SplitSampleModule::process(uhh2::Event& event){
  bool keep_event = false;
  if(counter == sample_nr) keep_event = true;

  if(counter == divide_by) counter = 1;
  else counter++;

  return keep_event;
}

UHH2_REGISTER_ANALYSIS_MODULE(SplitSampleModule)
