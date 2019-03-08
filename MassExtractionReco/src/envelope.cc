#include "../include/envelope.h"


using namespace std;

int main(int argc, char* argv[])
{
  TFile *file = new TFile("Histograms.root");

  std::vector<TH1F*> ScaleVariations,BTagVariations,PUVariations;
  std::vector<TH1F*>  MuIdVariations, MuTrVariations;
  std::vector<TH1F*>  JECVariations, JERVariations, CORVariations;
  std::vector<TH1F*> GeneratorVariations, ShowerVariations, PDFVariations;

  std::vector<TString> ScaleNames = {"SCALE_upup", "SCALE_upnone", "SCALE_noneup", "SCALE_downdown", "SCALE_downnone", "SCALE_nonedown"};
  std::vector<TString> BTagNames = {"btagbcup", "btagbcdown", "btagudsgup", "btagudsgdown",};
  std::vector<TString> PUNames = {"puup", "pudown"};
  std::vector<TString> MuIdNames = {"muidup", "muiddown"};
  std::vector<TString> MuTrNames = {"mutrup", "mutrdown"};
  std::vector<TString> JECNames = {"jecup", "jecdown"};
  std::vector<TString> JERNames = {"jerup", "jerdown"};
  std::vector<TString> CORNames = {"corup", "cordown"};
  std::vector<TString> GeneratorNames = {"Generator"};
  std::vector<TString> ShowerNames = {"Shower"};
  std::vector<TString> PDFNames;
  for(unsigned int i=0; i<100; i++){
    TString name = "PDF_";
    name += i;
    PDFNames.push_back(name);
  }




  for(auto name: ScaleNames) ScaleVariations.push_back((TH1F*)file->Get(name));
  for(auto name: BTagNames) BTagVariations.push_back((TH1F*)file->Get(name));
  for(auto name: PUNames) PUVariations.push_back((TH1F*)file->Get(name));
  for(auto name: MuIdNames) MuIdVariations.push_back((TH1F*)file->Get(name));
  for(auto name: MuTrNames) MuTrVariations.push_back((TH1F*)file->Get(name));
  for(auto name: JECNames) JECVariations.push_back((TH1F*)file->Get(name));
  for(auto name: JERNames) JERVariations.push_back((TH1F*)file->Get(name));
  for(auto name: CORNames) CORVariations.push_back((TH1F*)file->Get(name));
  for(auto name: GeneratorNames) GeneratorVariations.push_back((TH1F*)file->Get(name));
  for(auto name: ShowerNames) ShowerVariations.push_back((TH1F*)file->Get(name));
  for(auto name: PDFNames) PDFVariations.push_back((TH1F*)file->Get(name));

  TH1F* Central = (TH1F*)file->Get("tt");
  TH1F* WJets =   (TH1F*)file->Get("WJets");
  TH1F* ST =      (TH1F*)file->Get("SingleTop");
  TH1F* other =   (TH1F*)file->Get("other");
  TH1F* Data =    (TH1F*)file->Get("data");
  TH1F* StatError = GetStatError(Data);

  std::vector< TH1F* > EnvelopeUp;
  std::vector< TH1F* > EnvelopeDown;
  std::vector< TH1F* > ErrorUp;
  std::vector< TH1F* > ErrorDown;
  std::vector< TH1F* > ErrorSym;
  std::vector< TH2D* > CovMatrix;

  TH1F* BKGEnvelopeUp = GetBackgroundEnvelopeUp(Central, WJets, ST, other);
  TH1F* BKGErrorUp = GetError(Central, BKGEnvelopeUp);
  TH1F* BKGEnvelopeDown = GetBackgroundEnvelopeDown(Central, WJets, ST, other);
  TH1F* BKGErrorDown = GetError(Central, BKGEnvelopeDown);
  TH1F* BKGErrorSym = GetSymmetricError(BKGErrorUp, BKGErrorDown);
  EnvelopeUp.push_back(BKGEnvelopeUp);
  EnvelopeDown.push_back(BKGEnvelopeDown);
  ErrorUp.push_back(BKGErrorUp);
  ErrorDown.push_back(BKGErrorDown);
  ErrorSym.push_back(BKGErrorSym);
  CovMatrix.push_back(GetCovMatrix(BKGErrorSym));

  TH1F* ScaleEnvelopeUp = GetEnvelopeUp(ScaleVariations);
  TH1F* ScaleErrorUp = GetError(Central, ScaleEnvelopeUp);
  TH1F* ScaleEnvelopeDown = GetEnvelopeDown(ScaleVariations);
  TH1F* ScaleErrorDown = GetError(Central, ScaleEnvelopeDown);
  TH1F* ScaleErrorSym = GetSymmetricError(ScaleErrorUp, ScaleErrorDown);
  EnvelopeUp.push_back(ScaleEnvelopeUp);
  EnvelopeDown.push_back(ScaleEnvelopeDown);
  ErrorUp.push_back(ScaleErrorUp);
  ErrorDown.push_back(ScaleErrorDown);
  ErrorSym.push_back(ScaleErrorSym);
  CovMatrix.push_back(GetCovMatrix(ScaleErrorSym));

  TH1F* BTagEnvelopeUp = GetEnvelopeUp(BTagVariations);
  TH1F* BTagErrorUp = GetError(Central, BTagEnvelopeUp);
  TH1F* BTagEnvelopeDown = GetEnvelopeDown(BTagVariations);
  TH1F* BTagErrorDown = GetError(Central, BTagEnvelopeDown);
  TH1F* BTagErrorSym = GetSymmetricError(BTagErrorUp, BTagErrorDown);
  EnvelopeUp.push_back(BTagEnvelopeUp);
  EnvelopeDown.push_back(BTagEnvelopeDown);
  ErrorUp.push_back(BTagErrorUp);
  ErrorDown.push_back(BTagErrorDown);
  ErrorSym.push_back(BTagErrorSym);
  CovMatrix.push_back(GetCovMatrix(BTagErrorSym));

  TH1F* PUEnvelopeUp = GetEnvelopeUp(PUVariations);
  TH1F* PUErrorUp = GetError(Central, PUEnvelopeUp);
  TH1F* PUEnvelopeDown = GetEnvelopeDown(PUVariations);
  TH1F* PUErrorDown = GetError(Central, PUEnvelopeDown);
  TH1F* PUErrorSym = GetSymmetricError(PUErrorUp, PUErrorDown);
  EnvelopeUp.push_back(PUEnvelopeUp);
  EnvelopeDown.push_back(PUEnvelopeDown);
  ErrorUp.push_back(PUErrorUp);
  ErrorDown.push_back(PUErrorDown);
  ErrorSym.push_back(PUErrorSym);
  CovMatrix.push_back(GetCovMatrix(PUErrorSym));

  TH1F* MuIdEnvelopeUp = GetEnvelopeUp(MuIdVariations);
  TH1F* MuIdErrorUp = GetError(Central, MuIdEnvelopeUp);
  TH1F* MuIdEnvelopeDown = GetEnvelopeDown(MuIdVariations);
  TH1F* MuIdErrorDown = GetError(Central, MuIdEnvelopeDown);
  TH1F* MuIdErrorSym = GetSymmetricError(MuIdErrorUp, MuIdErrorDown);
  EnvelopeUp.push_back(MuIdEnvelopeUp);
  EnvelopeDown.push_back(MuIdEnvelopeDown);
  ErrorUp.push_back(MuIdErrorUp);
  ErrorDown.push_back(MuIdErrorDown);
  ErrorSym.push_back(MuIdErrorSym);
  CovMatrix.push_back(GetCovMatrix(MuIdErrorSym));

  TH1F* MuTrEnvelopeUp = GetEnvelopeUp(MuTrVariations);
  TH1F* MuTrErrorUp = GetError(Central, MuTrEnvelopeUp);
  TH1F* MuTrEnvelopeDown = GetEnvelopeDown(MuTrVariations);
  TH1F* MuTrErrorDown = GetError(Central, MuTrEnvelopeDown);
  TH1F* MuTrErrorSym = GetSymmetricError(MuTrErrorUp, MuTrErrorDown);
  EnvelopeUp.push_back(MuTrEnvelopeUp);
  EnvelopeDown.push_back(MuTrEnvelopeDown);
  ErrorUp.push_back(MuTrErrorUp);
  ErrorDown.push_back(MuTrErrorDown);
  ErrorSym.push_back(MuTrErrorSym);
  CovMatrix.push_back(GetCovMatrix(MuTrErrorSym));

  TH1F* JECEnvelopeUp = GetEnvelopeUp(JECVariations);
  TH1F* JECErrorUp = GetError(Central, JECEnvelopeUp);
  TH1F* JECEnvelopeDown = GetEnvelopeDown(JECVariations);
  TH1F* JECErrorDown = GetError(Central, JECEnvelopeDown);
  TH1F* JECErrorSym = GetSymmetricError(JECErrorUp, JECErrorDown);
  EnvelopeUp.push_back(JECEnvelopeUp);
  EnvelopeDown.push_back(JECEnvelopeDown);
  ErrorUp.push_back(JECErrorUp);
  ErrorDown.push_back(JECErrorDown);
  ErrorSym.push_back(JECErrorSym);
  CovMatrix.push_back(GetCovMatrix(JECErrorSym));

  TH1F* JEREnvelopeUp = GetEnvelopeUp(JERVariations);
  TH1F* JERErrorUp = GetError(Central, JEREnvelopeUp);
  TH1F* JEREnvelopeDown = GetEnvelopeDown(JERVariations);
  TH1F* JERErrorDown = GetError(Central, JEREnvelopeDown);
  TH1F* JERErrorSym = GetSymmetricError(JERErrorUp, JERErrorDown);
  EnvelopeUp.push_back(JEREnvelopeUp);
  EnvelopeDown.push_back(JEREnvelopeDown);
  ErrorUp.push_back(JERErrorUp);
  ErrorDown.push_back(JERErrorDown);
  ErrorSym.push_back(JERErrorSym);
  CovMatrix.push_back(GetCovMatrix(JERErrorSym));

  TH1F* COREnvelopeUp = GetEnvelopeUp(CORVariations);
  TH1F* CORErrorUp = GetError(Central, COREnvelopeUp);
  TH1F* COREnvelopeDown = GetEnvelopeDown(CORVariations);
  TH1F* CORErrorDown = GetError(Central, COREnvelopeDown);
  TH1F* CORErrorSym = GetSymmetricError(CORErrorUp, CORErrorDown);
  EnvelopeUp.push_back(COREnvelopeUp);
  EnvelopeDown.push_back(COREnvelopeDown);
  ErrorUp.push_back(CORErrorUp);
  ErrorDown.push_back(CORErrorDown);
  ErrorSym.push_back(CORErrorSym);
  CovMatrix.push_back(GetCovMatrix(CORErrorSym));

  TH1F* GeneratorEnvelopeUp = GetEnvelopeUp(GeneratorVariations);
  TH1F* GeneratorErrorUp = GetError(Central, GeneratorEnvelopeUp);
  TH1F* GeneratorEnvelopeDown = GetEnvelopeDown(GeneratorVariations);
  TH1F* GeneratorErrorDown = GetError(Central, GeneratorEnvelopeDown);
  TH1F* GeneratorErrorSym = GetSymmetricError(GeneratorErrorUp, GeneratorErrorDown);
  EnvelopeUp.push_back(GeneratorEnvelopeUp);
  EnvelopeDown.push_back(GeneratorEnvelopeDown);
  ErrorUp.push_back(GeneratorErrorUp);
  ErrorDown.push_back(GeneratorErrorDown);
  ErrorSym.push_back(GeneratorErrorSym);
  CovMatrix.push_back(GetCovMatrix(GeneratorErrorSym));

  TH1F* ShowerEnvelopeUp = GetEnvelopeUp(ShowerVariations);
  TH1F* ShowerErrorUp = GetError(Central, ShowerEnvelopeUp);
  TH1F* ShowerEnvelopeDown = GetEnvelopeDown(ShowerVariations);
  TH1F* ShowerErrorDown = GetError(Central, ShowerEnvelopeDown);
  TH1F* ShowerErrorSym = GetSymmetricError(ShowerErrorUp, ShowerErrorDown);
  EnvelopeUp.push_back(ShowerEnvelopeUp);
  EnvelopeDown.push_back(ShowerEnvelopeDown);
  ErrorUp.push_back(ShowerErrorUp);
  ErrorDown.push_back(ShowerErrorDown);
  ErrorSym.push_back(ShowerErrorSym);
  CovMatrix.push_back(GetCovMatrix(ShowerErrorSym));

  TH1F* PDFEnvelopeUp = GetEnvelopeUp(PDFVariations);
  TH1F* PDFErrorUp = GetError(Central, PDFEnvelopeUp);
  TH1F* PDFEnvelopeDown = GetEnvelopeDown(PDFVariations);
  TH1F* PDFErrorDown = GetError(Central, PDFEnvelopeDown);
  TH1F* PDFErrorSym = GetSymmetricError(PDFErrorUp, PDFErrorDown);
  EnvelopeUp.push_back(PDFEnvelopeUp);
  EnvelopeDown.push_back(PDFEnvelopeDown);
  ErrorUp.push_back(PDFErrorUp);
  ErrorDown.push_back(PDFErrorDown);
  ErrorSym.push_back(PDFErrorSym);
  CovMatrix.push_back(GetCovMatrix(PDFErrorSym));

  // create a dummy sys with constant shift (for testing)
  TH1F* DummyUp = DummyConstVariation(Central, 2);
  TH1F* DummyDown = DummyConstVariation(Central, 1);
  TH1F* DummyErrorUp = GetError(Central, DummyUp);
  TH1F* DummyErrorDown = GetError(Central, DummyDown);
  TH1F* DummyErrorSym = GetSymmetricError(DummyErrorUp, DummyErrorDown);
  TH2D* DummyCov = GetCovMatrix(DummyErrorSym);

  TFile* Dummy = new TFile("SYS_DUMMY.root","RECREATE");
  Dummy->mkdir("CovMatrix");
  Dummy->mkdir("ErrorUp");
  Dummy->mkdir("ErrorDown");
  Dummy->mkdir("ErrorSym");
  Dummy->mkdir("EnvelopeUp");
  Dummy->mkdir("EnvelopeDown");
  Dummy->cd("CovMatrix");
  DummyCov->Write();
  Dummy->cd("ErrorUp");
  DummyErrorUp->Write();
  Dummy->cd("ErrorDown");
  DummyErrorDown->Write();
  Dummy->cd("ErrorSym");
  DummyErrorSym->Write();
  Dummy->cd("EnvelopeUp");
  DummyUp->Write();
  Dummy->cd("EnvelopeDown");
  DummyDown->Write();
  Dummy->Close();
  ////

  TFile *Stat, *BKG;

  Stat = new TFile("STAT.root","RECREATE");
  Stat->mkdir("Up");
  Stat->mkdir("Down");
  Stat->cd("Up");
  StatError->Write();
  Stat->cd("Down");
  StatError->Write();
  Stat->Close();


  std::vector<TFile*> SysFiles;
  SysFiles.push_back(new TFile("SYS_BKG.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_SCALE.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_BTAG.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_PU.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_MUID.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_MUTR.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_JEC.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_JER.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_COR.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_GENERATOR.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_SHOWER.root","RECREATE"));
  SysFiles.push_back(new TFile("SYS_PDF.root","RECREATE"));

  if( SysFiles.size() != EnvelopeUp.size() ) cout << "Not same number of Files and Envelope Up Hists!" << endl;
  if( SysFiles.size() != EnvelopeDown.size() ) cout << "Not same number of Files and Envelope Down Hists!" << endl;
  if( SysFiles.size() != ErrorUp.size() ) cout << "Not same number of Files and Error Up Hists!" << endl;
  if( SysFiles.size() != ErrorDown.size() ) cout << "Not same number of Files and Error Down Hists!" << endl;
  if( SysFiles.size() != ErrorSym.size() ) cout << "Not same number of Files and Error Sym Hists!" << endl;
  if( SysFiles.size() != CovMatrix.size() ) cout << "Not same number of Files and Cov Matrices!" << endl;

  std::cout << "Write Files..." << std::endl;
  for(unsigned int i=0; i<SysFiles.size(); i++){
    SysFiles[i]->mkdir("ErrorUp");
    SysFiles[i]->mkdir("ErrorDown");
    SysFiles[i]->mkdir("ErrorSym");
    SysFiles[i]->mkdir("EnvelopeUp");
    SysFiles[i]->mkdir("EnvelopeDown");
    SysFiles[i]->mkdir("CovMatrix");
    SysFiles[i]->cd("ErrorUp");
    ErrorUp[i]->Write();
    SysFiles[i]->cd("ErrorDown");
    ErrorDown[i]->Write();
    SysFiles[i]->cd("ErrorSym");
    ErrorSym[i]->Write();
    SysFiles[i]->cd("EnvelopeUp");
    EnvelopeUp[i]->Write();
    SysFiles[i]->cd("EnvelopeDown");
    EnvelopeDown[i]->Write();
    SysFiles[i]->cd("CovMatrix");
    CovMatrix[i]->Write();
    SysFiles[i]->Close();
  }



  std::cout<<"Finished!" << std::endl;
  return 0;

}


TH1F* GetEnvelopeUp(std::vector<TH1F*> Variations){
  int nbins = Variations[0]->GetSize() - 2;
  TH1F* hist = (TH1F*)Variations[0]->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double max = 0;
    for(int j=0; j<Variations.size(); j++){
      double cont = Variations[j]->GetBinContent(bin);
      if(cont > max) max = cont;
    }
    hist->SetBinContent(bin, max);
  }
  return hist;
}

TH1F* GetEnvelopeDown(std::vector<TH1F*> Variations){
  int nbins = Variations[0]->GetSize() - 2;
  TH1F* hist = (TH1F*)Variations[0]->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double min = 100000000000000;
    for(int j=0; j<Variations.size(); j++){
      double cont = Variations[j]->GetBinContent(bin);
      if(cont < min) min = cont;
    }
    hist->SetBinContent(bin, min);
  }
  return hist;
}

TH1F* GetError(TH1F* Central, TH1F* Envelope){
  int nbins = Central->GetSize() - 2;
  TH1F* hist = (TH1F*)Central->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double diff = Envelope->GetBinContent(bin) - Central->GetBinContent(bin);
    hist->SetBinContent(bin, diff);
  }
  return hist;
}

TH1F* GetStatError(TH1F* Data){
  int nbins = Data->GetSize() - 2;
  TH1F* hist = (TH1F*)Data->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double err = Data->GetBinError(bin);
    hist->SetBinContent(bin, err);
  }
  return hist;
}

TH2D* GetCovMatrix(TH1F* Error){
  int nbins = Error->GetSize() - 2;
  TH2D* cov = new TH2D("cov", "cov", nbins, 1, nbins+1, nbins, 1, nbins+1);
  for(int i=1; i <= nbins; i++){
    for(int j=1; j <= nbins; j++){
      double CovEntry = Error->GetBinContent(i) * Error->GetBinContent(j);
      cov->Fill(i,j,CovEntry);
    }
  }
  return cov;
}

TH1F* DummyConstVariation(TH1F* Central, double factor){
  int nbins = Central->GetSize() - 2;
  TH1F* hist = (TH1F*)Central->Clone("hist");
  for(int i=1; i <= nbins; i++){
    double entry = factor * Central->GetBinContent(i);
    hist->SetBinContent(i, entry);
  }
  return hist;
}

TH1F* GetBackgroundEnvelopeUp(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other){
  int nbins = Central->GetSize() - 2;
  TH1F* hist = (TH1F*)Central->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double nominal = Central->GetBinContent(bin);
    double bkg = WJets->GetBinContent(bin) * 0.19;
    bkg += ST->GetBinContent(bin) * 0.23;
    bkg += other->GetBinContent(bin) * 1.00;
    hist->SetBinContent(bin, nominal+bkg);
  }
  return hist;
}

TH1F* GetBackgroundEnvelopeDown(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other){
  int nbins = Central->GetSize() - 2;
  TH1F* hist = (TH1F*)Central->Clone("hist");
  hist->Reset();
  for(int bin=1; bin<=nbins; bin++){
    double nominal = Central->GetBinContent(bin);
    double bkg = WJets->GetBinContent(bin) * 0.19;
    bkg += ST->GetBinContent(bin) * 0.23;
    bkg += other->GetBinContent(bin) * 1.00;
    hist->SetBinContent(bin, nominal-bkg);
  }
  return hist;
}

TH1F* GetSymmetricError(TH1F* up, TH1F* down){
  // calculate absolute values from up and down Variations
  // but the signe matters because of anti-correlations
  // thus, take the largest absolute value, with the sign of the up variation
  // one could also take the sign of the down variation
  int nbins = up->GetSize() - 2;
  TH1F* sym = (TH1F*)up->Clone("sym");
  sym->Reset();
  for(int i=1; i<=nbins; i++){
    double a = up->GetBinContent(i);
    double b = down->GetBinContent(i);
    double sign;
    if(a != 0) sign = a/abs(a);
    else sign = 1.;
    double err;
    if(abs(a) > abs(b)) err = abs(a);
    else                err = abs(b);
    double xvalue = up->GetBinCenter(i);
    double value = err*sign;
    sym->Fill(xvalue, err*sign);
  }
  return sym;
}
