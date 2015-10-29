#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include "TEnv.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"

#include "CommandLineInterface.hh"
#include "Kinematics.hh"
#include "Nucleus.hh"
#include "Reconstruction.hh"
#ifndef rad2deg
#define rad2deg                       180./(TMath::Pi())
#endif
#ifndef deg2rad
#define deg2rad                       (TMath::Pi())/180.
#endif
using namespace std;
TSpline3* theory[4];
TSpline3* theoryres;
int s1;
int s[2];
double fit(Double_t *x, Double_t *par){
  //cout << " s1 " << s1 << endl;
  //cout << " points " << theory[s1]->GetNp() << endl;
  if(par[0]<=0)
    return 1e6;
  double theta[theory[s1]->GetNp()];
  double sigma[theory[s1]->GetNp()];
  for(int i=0;i<theory[s1]->GetNp();i++){      
    theory[s1]->GetKnot(i, theta[i], sigma[i]);      
    //cout << i << "\t" << theta[i] << "\t" << sigma[i] << endl;
    sigma[i]*=par[0];
  }
  theoryres = new TSpline3("fittedspline",theta,sigma,theory[s1]->GetNp());
  //return 0;
  //cout << "par[0] " << par[0] <<endl;
  return theoryres->Eval(x[0]);
}
double fit2(Double_t *x, Double_t *par){
  if(par[0]<=0)
    return 1e6;
  double theta[2][theory[s[0]]->GetNp()];
  double sigma[2][theory[s[0]]->GetNp()];
  double sigmares[theory[s[0]]->GetNp()];
  for(int i=0;i<theory[s[0]]->GetNp();i++){      
    sigmares[i] = 0;
    theory[s[0]]->GetKnot(i, theta[0][i], sigma[0][i]);      
    //cout << i << "\t" << theta[i] << "\t" << sigma[i] << endl;
    sigma[0][i]*=par[0];
    sigmares[i] += sigma[0][i];
  }
  for(int i=0;i<theory[s[1]]->GetNp();i++){      
    theory[s[1]]->GetKnot(i, theta[1][i], sigma[1][i]);      
    //cout << i << "\t" << theta[i] << "\t" << sigma[i] << endl;
    sigma[1][i]*=par[1];
    sigmares[i] += sigma[1][i];
  }
  theoryres = new TSpline3("fittedspline",theta[0],sigmares,theory[s[0]]->GetNp());
  //return 0;
  //cout << "par[0] " << par[0] <<endl;
  return theoryres->Eval(x[0]);
}
TSpline3* scaled(Double_t factor, TSpline3* fitted){
  double theta[fitted->GetNp()];
  double sigma[fitted->GetNp()];
  for(int i=0;i<fitted->GetNp();i++){      
    fitted->GetKnot(i, theta[i], sigma[i]);      
    sigma[i]*=factor;
  }
  TSpline3 *res = new TSpline3("result",theta,sigma,fitted->GetNp());  
  return res;
}
TSpline3* scaled2(Double_t factor, TSpline3* fitted, Double_t factor2, TSpline3* fitted2){
  double theta[fitted->GetNp()];
  double sigma[fitted->GetNp()];
  double sigmares[fitted->GetNp()];
  for(int i=0;i<fitted->GetNp();i++){      
    sigmares[i] = 0;
    //cout << sigmares[i] << "\t" ;
    fitted->GetKnot(i, theta[i], sigma[i]);      
    sigma[i]*=factor;
    //cout <<" + " <<sigma[i] << " = ";
    sigmares[i] += sigma[i];
    //cout << sigmares[i] << "\t" ;
    fitted2->GetKnot(i, theta[i], sigma[i]);      
    sigma[i]*=factor2;
    //cout <<" + " <<sigma[i] << " = ";
    sigmares[i] += sigma[i];
    //cout << sigmares[i] << endl ;
  }
  TSpline3 *res = new TSpline3("result",theta,sigmares,fitted->GetNp());  
  return res;
}
int main(int argc, char* argv[]){
  char* OutputFile = NULL;
  char* SettFile = NULL;
  int gs = 0;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "settings file", &SettFile);
  interface->Add("-g", "ground state (if neon set higher than 2, if magnesium 0 gs 1 ex, if >10 assym errors only one fit contrib)", &gs);
  interface->CheckFlags(argc, argv);
  if(OutputFile == NULL){
    cerr<<"You have to provide the output file!"<<endl;
    exit(1);
  }
  if(SettFile == NULL){
    cerr<<"You have to provide the settings file!"<<endl;
    exit(1);
  }
  TEnv *sett = new TEnv((char*)SettFile);
  const char* InputFile = NULL;
  const char* TheoryFile = NULL;

  TColor *orange = gROOT->GetColor(5);
  orange->SetRGB(1.0, 0.612, 0.002); 
  
  TColor *green = gROOT->GetColor(3);
  green->SetRGB(0.15, 0.7, 0.15);
  gStyle->SetOptTitle(0); 

  bool axis = true;
  TAxis *xa, *ya;
  double xval[2]={0,180};
  double yval[2]={10,5e5};
  TGraph *dummy;

  int nrofdata;
  int nrofsplines;
  int lim[2];
  InputFile = sett->GetValue("InputFile","/home/kwimmer/fresco/30Mg_3H/backwarderrors/gsdata.dat");
  TheoryFile = sett->GetValue("TheoryFile","/home/kwimmer/fresco/30Mg_3H/configuration/mgall1d2p.root");
  nrofdata = sett->GetValue("NrofData",0);
  nrofsplines = sett->GetValue("NrofSplines",0);
  lim[0] = sett->GetValue("FitLowLim",0);
  lim[1] = sett->GetValue("FitHighLim",180);
  double start[nrofsplines];
  for(int i=0;i<nrofsplines;i++){
    start[i] = sett->GetValue(Form("StartVal.%d",i),1.);
  }
  cout << "InputFile " << InputFile << endl;
  cout << "TheoryFile " << TheoryFile << endl;
  cout << "nrofdata " << nrofdata << endl;
  cout << "nrofsplines " << nrofsplines << endl;
  TFile *output = new TFile(Form("%s.root",OutputFile),"recreate");
  //TFile *infile = new TFile(InputFile);
  TFile *theofile = new TFile(TheoryFile);

  ifstream datat(InputFile);
  double tan[nrofdata],tanl[nrofdata],tanh[nrofdata],tsi[nrofdata],tsil[nrofdata],tsih[nrofdata];
  for(int i=0; i<nrofdata; i++){
    if(gs<2){
      datat >> tan[i] >> tanl[i] >> tanh[i] >> tsi[i] >> tsil[i]>> tsih[i];
      if(i>18)//???????????????????????????????
	tsih[i]*=2;
    }
    else if(gs<10){
      datat >> tan[i] >> tanl[i]>> tsi[i] >> tsil[i];
      tanh[i] = tanl[i];
      if(gs<4){
	if(i>18)
	  tsih[i]*=2;
      }
      tsih[i] = tsil[i];
    }
    else{
      datat >> tan[i] >> tanl[i] >> tanh[i] >> tsi[i] >> tsil[i]>> tsih[i];
      if(i>18)//???????????????????????????????
	tsih[i]*=2;
    }
  }
  TGraphAsymmErrors *tdata = new TGraphAsymmErrors(nrofdata,tan,tsi,tanl,tanh,tsil,tsih);
  //tdata->Draw("P");



  TSpline3* theorysp[nrofsplines];
  for(int i=0; i<nrofsplines; i++){
    theorysp[i] = (TSpline3*)theofile->Get(Form("spline_%d",i));
    theory[i] = theorysp[i];
  }
  for(int i=0; i<4; i++){
    if(gs==1)
      theory[i] = theorysp[i];
    else if(gs==0)
      theory[i] = theorysp[i+4];      
  }
  double startpar[nrofsplines];
  TF1* fitted[nrofsplines];
  TCanvas *c = new TCanvas("c","c",0,0,400,300);
  TF1* fittedp;
  TF1* fittedf;
  TF1* fitteds;
  TF1* fittedd;
  TF1* fitted2;
  //fitted2->SetParLimits(1,0.83,1.55);
  //fitted->Write(Form("fit_cm_%d",i),TObject::kOverwrite);

  c->cd(0);
  TPad *pad = new TPad("pad","",0,0,1,1,-1,-1,0);
  pad->Draw();
  pad->cd();
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.1);
  pad->SetBottomMargin(0.12);
  pad->SetTopMargin(0.015);
  pad->SetLogy(1);
  //yval[0] = cm->GetMinimum()*0.5;
  //yval[1] = cm->GetMaximum()*2;
  double x0,x1;
  if(axis){
    tdata->ComputeRange(x0,yval[0],x1,yval[1]);
    //cout << x0 << "\t" << x1 << "\t" << yval[0] << "\t" << yval[1] << endl;
    yval[0]/=2.;
    yval[1]*=2.;
  }
  dummy = new TGraph(2,xval,yval);
  xa = dummy->GetXaxis();
  xa->SetRangeUser(0.,180.);
  xa->SetTitle("#vartheta_{cm} [#circ]");
  xa->SetTitleOffset(1.0);
  xa->SetTitleSize(0.05);
  xa->SetLabelOffset(0.008);
  xa->SetLabelColor(1);
  xa->SetLabelSize(0.05);
  xa->SetTickLength(0.03);
    
  ya = dummy->GetYaxis();
  ya->SetTitle("d#sigma/d#Omega [a.u.]");
  ya->SetTitleOffset(0.9);
  ya->SetTitleSize(0.05);
  ya->SetLabelOffset(0.005);
  ya->SetLabelColor(1);
  ya->SetLabelSize(0.05);
  ya->SetTickLength(0.03);
  dummy->Draw("AP");
  tdata->Draw("P");
  double spec[4];
  double chi[4];
  double spece[4];
  TSpline3 *result[3];
  if(gs<2){
    //fitted->Draw("same");
    fittedp = new TF1("fittedp",fit,0,180,1);
    int sp=0;
    int sf=0;
    if(gs==1){
      cout << " (p3/2)^2 " << endl;
      s1 = 2;
      sp = 2;
    }
    else{
      cout << " (s1/2)^2 " << endl;
      s1 = 0;
      sp = 0;
    }
    fittedp->SetParameter(0,1);
    tdata->Fit(fittedp,"RN");
    spec[0] = fittedp->GetParameter(0);
    
    fittedf = new TF1("fittedf",fit,0,180,1);
    if(gs==1){
      cout << " (f7/2)^2 " << endl;
      s1 = 3;
      sf = 3;
    }
    else{
      cout << " (d3/2)^2 " << endl;
      s1 = 1;
      sf = 1;
    }
    fittedf->SetParameter(0,10);
    tdata->Fit(fittedf,"RN");
    spec[1] = fittedf->GetParameter(0);
    
    fitted2 = new TF1("fitted2",fit2,0,180,2);
    if(gs==1){
      cout << " a(p3/2)^2 + b(f7/2)^2 " << endl;
      s[0] = 2;
      s[1] = 3;
    }
    else{
      cout << " a(s1/2)^2 + b(d3/2)^2 " << endl;
      s[0] = 0;
      s[1] = 1;
    }
    fitted2->SetParameter(0,0.16);
    fitted2->SetParameter(1,1.89);
    fitted2->SetParLimits(0,0.15,0.17);
    tdata->Fit(fitted2,"RN");
    
    spec[2] = fitted2->GetParameter(0);
    spec[3] = fitted2->GetParameter(1);
    int k=0;
    for(int i=0;i<4;i++){
      if(gs==0)
	k=i+4;
      else
	k=i;
      theorysp[k]->SetLineColor(1+i);
      theorysp[k]->SetLineWidth(1);    
      theorysp[k]->SetLineStyle(2);    
      theorysp[k]->Draw("same");
      cout << spec[i] << "\t";
    }
    cout << endl;
    //theorysp[i]->Write(Form("fit_range_cm_%d",i),TObject::kOverwrite);
    //theorysp[i]->SetLineColor(3);
    //theorysp[i]->SetLineWidth(1); 
    //theorysp[i]->Draw("same");
    //TLatex *txt = new TLatex(20,yval[1]/10,Form("%s",Histos[0][i]->GetName()));
    //txt->SetTextSize(0.07);
    //txt->Draw();
    //}
    if(gs==0){
      result[0] = scaled2(fitted2->GetParameter(0), theorysp[s[0]+4],fitted2->GetParameter(1), theorysp[s[1]+4]);
      result[1] = scaled(fittedp->GetParameter(0), theorysp[0+4]);
      result[2] = scaled(fittedf->GetParameter(0), theorysp[1+4]);
    }
    else{
      result[0] = scaled2(fitted2->GetParameter(0), theorysp[s[0]],fitted2->GetParameter(1), theorysp[s[1]]);
      result[1] = scaled(fittedp->GetParameter(0), theorysp[2]);
      result[2] = scaled(fittedf->GetParameter(0), theorysp[3]);
    }
    
    result[0]->SetLineColor(5);
    result[0]->Draw("same");
    
    result[1]->SetLineColor(sp+1);
    result[1]->Draw("same");
    
    result[2]->SetLineColor(sf+1);
    result[2]->Draw("same");
    
    /*
      cout << "fittedp->GetParameter(0) "<< fittedp->GetParameter(0) << endl;
      fittedp->SetParameter(0,spec[0]);
      fittedp->SetLineColor(2);
      fittedp->SetLineWidth(1);
      fittedp->Draw("same");
      
      fittedf->SetParameter(0,spec[1]);
      fittedf->SetLineColor(4);
      fittedf->SetLineWidth(1);
      fittedf->Draw("same");
      
      fitted2->SetParameter(0,spec[2]);
      fitted2->SetParameter(1,spec[3]);
      fitted2->SetLineColor(3);
      fitted2->SetLineWidth(1);
      fitted2->Draw("same");
    */
    output->cd();
    result[0]->Write("result0",TObject::kOverwrite);
    result[1]->Write("result1",TObject::kOverwrite);
    result[2]->Write("result2",TObject::kOverwrite);
  
    output->Close();
  }
  else{
    for(int i=0;i<nrofsplines;i++){
      fitted[i] = new TF1(Form("fitted%d",i),fit,lim[0],lim[1],1);
      s1=i;
      fitted[i]->SetParameter(0,start[i]);
      tdata->Fit(fitted[i],"RN");
      chi[i] = fitted[i]->GetChisquare();
      spec[i] = fitted[i]->GetParameter(0);
      spece[i] = fitted[i]->GetParError(0);
      result[i] = scaled(fitted[i]->GetParameter(0), theorysp[i]);
      theorysp[i]->SetLineStyle(7);
      theorysp[i]->Draw("same");
      result[i]->Draw("same");
    }
    cout << "spectroscopic factors " << endl;
    for(int i=0;i<nrofsplines;i++)
      cout << spec[i] << " +- " << spece[i] << "\t";
    cout << endl;
    cout << "chi " << endl;
    for(int i=0;i<nrofsplines;i++)
      cout << chi[i] << "\t";
    cout << endl;
  }
  output->Close();
  c->SaveAs(Form("%s.ps",OutputFile));
}
