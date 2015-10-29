#include <iostream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "TH1S.h"
//#include "TH2S.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"

void plot(float offset, float gain){
  TFile *splines = new TFile("/user/wimmer/analysis/gui/ELoss/1.5mm.root");
  
  TSpline3* proton = (TSpline3*)splines->Get("proton_dee");
  vector<double> de;
  vector<double> er;
  de.resize(proton->GetNp());
  er.resize(proton->GetNp());
  for(int i=0;i<proton->GetNp();i++){
    proton->GetKnot(i, er[i], de[i]);
    if(er[i]+de[i]>84.6 && er[i]+de[i]<86.4)
      cout << "etot " << er[i]+de[i] << " de " << de[i] << " er " << er[i] << endl;
  }

  
  TFile *spectra = new TFile("/projects/e06006/calibrated/runs_031_035.root");
  //TFile *spectra = new TFile("/projects/e06006/calibrated/test34.root");
  
  TH2S *spec = (TH2S*)spectra->Get("f_csi_5_3");
  spec->Draw("colz");
  //proton->Draw("same");
  cout << proton->GetNp() << endl;
  const int n = proton->GetNp();
  double x[200];
  double y[200];
  for(UShort_t i=0;i<proton->GetNp();i++){
    proton->GetKnot(i, x[i], y[i]);
    x[i] -= offset;
    x[i] /=gain;
  }
  TGraph* graph = new TGraph(proton->GetNp(), x, y);
  TSpline3* spline = new TSpline3("pro",graph);
  spline->Draw("same");
}
