#include "../include/CreatHists.h"
#include "../include/CentralInclude.h"
#include "../include/HistogramUtils.h"
#include "../include/GraphUtils.h"
#include "../include/Utils.h"
#include "../include/tdrstyle_all.h"

using namespace std;

TString year;

int main(int argc, char* argv[]){

  TH1F* color_100 = new TH1F("colors_100", "colors_100", 100, 0, 100);
  TH1F* color_200 = new TH1F("colors_200", "colors_200", 100, 100, 200);
  TH1F* color_300 = new TH1F("colors_300", "colors_300", 100, 200, 300);
  TH1F* color_400 = new TH1F("colors_400", "colors_400", 100, 300, 400);
  TH1F* color_500 = new TH1F("colors_500", "colors_500", 100, 400, 500);
  TH1F* color_600 = new TH1F("colors_600", "colors_600", 100, 500, 600);
  TH1F* color_700 = new TH1F("colors_700", "colors_700", 100, 600, 700);
  TH1F* color_800 = new TH1F("colors_800", "colors_800", 100, 700, 800);
  TH1F* color_900 = new TH1F("colors_900", "colors_900", 100, 800, 900);
  TH1F* color;
  int icolor = 0;
  int color_max = 5;
  while(icolor<color_max){
    for(unsigned int i=0; i<color->GetNbinsX(); i++){
      color_100->GetXaxis()->SetRangeUser(0, 100);
      color_100->GetYaxis()->SetRangeUser(0, 1.2);
      color_200->GetXaxis()->SetRangeUser(0, 100);
      color_200->GetYaxis()->SetRangeUser(0, 1.2);
      color_300->GetXaxis()->SetRangeUser(0, 100);
      color_300->GetYaxis()->SetRangeUser(0, 1.2);
      color_400->GetXaxis()->SetRangeUser(0, 100);
      color_400->GetYaxis()->SetRangeUser(0, 1.2);
      color_500->GetXaxis()->SetRangeUser(0, 100);
      color_500->GetYaxis()->SetRangeUser(0, 1.2);
      color_600->GetXaxis()->SetRangeUser(0, 100);
      color_600->GetYaxis()->SetRangeUser(0, 1.2);
      color_700->GetXaxis()->SetRangeUser(0, 100);
      color_700->GetYaxis()->SetRangeUser(0, 1.2);
      color_800->GetXaxis()->SetRangeUser(0, 100);
      color_800->GetYaxis()->SetRangeUser(0, 1.2);
      color_900->GetXaxis()->SetRangeUser(0, 100);
      color_900->GetYaxis()->SetRangeUser(0, 1.2);
    }
  }
}
