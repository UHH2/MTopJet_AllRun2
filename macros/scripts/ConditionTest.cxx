#include <iostream>
#include <memory>

using namespace std;

int ConditionTest(TString year){
  // TString year = argv[1];
  TString test = year.EqualTo("2016")?"2016":"2017";
  cout << test << endl;

  test =
  year.EqualTo("2016")?"2016"
  :year.EqualTo("2017")?"2017"
  :year.EqualTo("2018")?"2018"
  :"faild";
  cout << test << endl;
  return true;
}
