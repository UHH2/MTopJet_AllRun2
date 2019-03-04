#include "../include/CentralInclude.h"

// compile with:
// g++ -o tt_envelope tt_envelope.cc `root-config --cflags --glibs`

using namespace std;

TH1F* GetEnvelopeUp(std::vector<TH1F*> Variations);
TH1F* GetEnvelopeDown(std::vector<TH1F*> Variations);
TH1F* GetError(TH1F* Central, TH1F* Envelope);
TH1F* GetStatError(TH1F* Data);
TH2D* GetCovMatrix(TH1F* Error);
TH1F* GetBackgroundEnvelopeUp(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other);
TH1F* GetBackgroundEnvelopeDown(TH1F* Central, TH1F* WJets, TH1F* ST, TH1F* other);
TH1F* GetSymmetricError(TH1F* up, TH1F* down);
TH1F* DummyConstVariation(TH1F* Central, double factor);

int main(int argc, char* argv[])
{

  std::vector<TFile*> ScaleVariations_f, BTagVariations_f, PUVariations_f;
  std::vector<TFile*> MuIdVariations_f, MuTrVariations_f;
  std::vector<TFile*> JECVariations_f, JERVariations_f, CORVariations_f;
  std::vector<TFile*> GeneratorVariations_f, ShowerVariations_f, PDFVariations_f;

  ScaleVariations_f.push_back(new TFile(dir+"SCALE_upup/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  ScaleVariations_f.push_back(new TFile(dir+"SCALE_upnone/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  ScaleVariations_f.push_back(new TFile(dir+"SCALE_noneup/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  ScaleVariations_f.push_back(new TFile(dir+"SCALE_downdown/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  ScaleVariations_f.push_back(new TFile(dir+"SCALE_downnone/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  ScaleVariations_f.push_back(new TFile(dir+"SCALE_nonedown/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  BTagVariations_f.push_back(new TFile(dir+"BTAG_bc_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  BTagVariations_f.push_back(new TFile(dir+"BTAG_bc_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  BTagVariations_f.push_back(new TFile(dir+"BTAG_udsg_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  BTagVariations_f.push_back(new TFile(dir+"BTAG_udsg_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  PUVariations_f.push_back(new TFile(dir+"PU_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  PUVariations_f.push_back(new TFile(dir+"PU_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  MuIdVariations_f.push_back(new TFile(dir+"MUID_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  MuIdVariations_f.push_back(new TFile(dir+"MUID_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  MuTrVariations_f.push_back(new TFile(dir+"MUTR_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  MuTrVariations_f.push_back(new TFile(dir+"MUTR_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  JECVariations_f.push_back(new TFile(dir+"JEC_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  JECVariations_f.push_back(new TFile(dir+"JEC_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  JERVariations_f.push_back(new TFile(dir+"JER_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  JERVariations_f.push_back(new TFile(dir+"JER_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  CORVariations_f.push_back(new TFile(dir+"COR_up/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  CORVariations_f.push_back(new TFile(dir+"COR_down/uhh2.AnalysisModuleRunner.MC.TTbar.root"));
  GeneratorVariations_f.push_back(new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_amcatnlo-pythia.root"));
  ShowerVariations_f.push_back(new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar_powheg-herwig.root"));
  PDFVariations_f.push_back(new TFile(dir+"PDF_up/PDF_Variations.root"));
  PDFVariations_f.push_back(new TFile(dir+"PDF_down/PDF_Variations.root"));

  TFile* WJets_f = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.WJets.root");
  TFile* SingleTop_f = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.SingleTop.root");
  TFile* Other_f = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.other.root");

  TFile * Central_f = new TFile(dir+"uhh2.AnalysisModuleRunner.MC.TTbar.root");
  TFile * Data_f = new TFile(dir+"uhh2.AnalysisModuleRunner.DATA.DATA.root");


  TString histdir, histname;
  histdir = "XCone_cor/";
  histname = "M_jet1";

  std::vector<TH1F*> ScaleVariations,BTagVariations,PUVariations;
  std::vector<TH1F*>  MuIdVariations, MuTrVariations;
  std::vector<TH1F*>  JECVariations, JERVariations, CORVariations;
  std::vector<TH1F*> GeneratorVariations, ShowerVariations, PDFVariations;

  for(unsigned int i=0; i<ScaleVariations_f.size(); i++){
    ScaleVariations.push_back((TH1F*)ScaleVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<BTagVariations_f.size(); i++){
    BTagVariations.push_back((TH1F*)BTagVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<PUVariations_f.size(); i++){
    PUVariations.push_back((TH1F*)PUVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<MuIdVariations_f.size(); i++){
    MuIdVariations.push_back((TH1F*)MuIdVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<MuTrVariations_f.size(); i++){
    MuTrVariations.push_back((TH1F*)MuTrVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<JECVariations_f.size(); i++){
    JECVariations.push_back((TH1F*)JECVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<JERVariations_f.size(); i++){
    JERVariations.push_back((TH1F*)JERVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<CORVariations_f.size(); i++){
    CORVariations.push_back((TH1F*)CORVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<GeneratorVariations_f.size(); i++){
    GeneratorVariations.push_back((TH1F*)GeneratorVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<ShowerVariations_f.size(); i++){
    ShowerVariations.push_back((TH1F*)ShowerVariations_f[i]->Get(histdir + histname));
  }
  for(unsigned int i=0; i<PDFVariations_f.size(); i++){
    // pdf variations have a slightly different directory structure
    PDFVariations.push_back((TH1F*)PDFVariations_f[i]->Get(histname));
  }

  TH1F* Central = (TH1F*)Central_f->Get(histdir + histname);

  TH1F* WJets = (TH1F*)WJets_f->Get(histdir + histname);
  TH1F* ST = (TH1F*)SingleTop_f->Get(histdir + histname);
  TH1F* other = (TH1F*)Other_f->Get(histdir + histname);



  TH1F* Data = (TH1F*)Data_f->Get(histdir + histname);
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

  TFile* Dummy = new TFile(dir+"SYS_DUMMY.root","RECREATE");
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

  Stat = new TFile(dir+"STAT.root","RECREATE");
  Stat->mkdir("Up");
  Stat->mkdir("Down");
  Stat->cd("Up");
  StatError->Write();
  Stat->cd("Down");
  StatError->Write();
  Stat->Close();


  std::vector<TFile*> SysFiles;
  SysFiles.push_back(new TFile(dir+"SYS_BKG.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_SCALE.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_BTAG.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_PU.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_MUID.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_MUTR.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_JEC.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_JER.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_COR.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_GENERATOR.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_SHOWER.root","RECREATE"));
  SysFiles.push_back(new TFile(dir+"SYS_PDF.root","RECREATE"));

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
