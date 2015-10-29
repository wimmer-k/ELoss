#include <iostream>
#include <iomanip>

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"

#include "CommandLineInterface.hh"
#include "Compound.hh"
#include "Reconstruction.hh"

using namespace std;

int main(int argc, char* argv[]){
  char* OutputFile = NULL;
  double stepsize = 0.2; //in MeV
  double maxenrec = 100.; //in MeV

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-s", "stepsize", &stepsize);
  interface->Add("-m", "maxenrec", &maxenrec);
  interface->CheckFlags(argc, argv);

  if(OutputFile == NULL){
    cerr<<"You have to provide at least the output file!"<<endl;
    exit(1);
  }
  TFile* outfile = new TFile(OutputFile,"recreate");
  if(outfile->IsZombie()){
    return 4;
  }
  char* massfile = (char*)"/user/wimmer/analysis/eloss/mass.dat";

  Nucleus *proton = new Nucleus(1,0,massfile);
  Nucleus *deuteron = new Nucleus(1,1,massfile);
  Nucleus *triton = new Nucleus(1,2,massfile);
  Nucleus *alpha = new Nucleus(2,2,massfile);
  Nucleus *he3 = new Nucleus(2,1,massfile);

  Nucleus *silicon = new Nucleus(14,14,massfile);
  Compound *detector = new Compound(silicon);

  Nucleus *caesium = new Nucleus(55,78,massfile);
  Nucleus *iodine = new Nucleus(53,74,massfile);

  Compound *csi = new Compound(caesium,1,iodine,1);


  Reconstruction *protondeltaE = new Reconstruction(proton, detector);
  Reconstruction *deuterondeltaE = new Reconstruction(deuteron, detector);
  Reconstruction *tritondeltaE = new Reconstruction(triton, detector);
  Reconstruction *alphadeltaE = new Reconstruction(alpha, detector);
  Reconstruction *he3deltaE = new Reconstruction(he3, detector);

  Reconstruction *protonE = new Reconstruction(proton, csi);
  Reconstruction *deuteronE = new Reconstruction(deuteron, csi);
  Reconstruction *tritonE = new Reconstruction(triton, csi);
  Reconstruction *alphaE = new Reconstruction(alpha, csi);
  Reconstruction *he3E = new Reconstruction(he3, csi);

  TSpline3* protondee;
  protondeltaE->SetTargetThickness(1500*2.33/10.);
  protondee= protondeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  protondee->SetName(Form("proton_dee"));
  protondee->Write("",TObject::kOverwrite);

  TSpline3* protonmin;
  protondeltaE->SetTargetThickness(65*2.33/10.);
  protonmin= protondeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  protonmin->SetName(Form("proton_min"));
  protonmin->Write("",TObject::kOverwrite);

  TSpline3* protonmax;
  protondeltaE->SetTargetThickness(1500*2.33/10./cos(20* TMath::Pi()/180.));
  protonmax= protondeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  protonmax->SetName(Form("proton_max"));
  protonmax->Write("",TObject::kOverwrite);


  TSpline3* deuterondee;
  deuterondeltaE->SetTargetThickness(1500*2.33/10.);
  deuterondee= deuterondeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  deuterondee->SetName(Form("deuteron_dee"));
  deuterondee->Write("",TObject::kOverwrite);
  TSpline3* tritondee;
  tritondeltaE->SetTargetThickness(1500*2.33/10.);
  tritondee= tritondeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  tritondee->SetName(Form("triton_dee"));
  tritondee->Write("",TObject::kOverwrite);
  TSpline3* alphadee;
  alphadeltaE->SetTargetThickness(1500*2.33/10.);
  alphadee= alphadeltaE->Energy2EnergyAfter(maxenrec,stepsize);
  alphadee->SetName(Form("alpha_dee"));
  alphadee->Write("",TObject::kOverwrite);
  TSpline3* he3dee;
  he3deltaE->SetTargetThickness(1500*2.33/10.);
  he3dee= he3deltaE->Energy2EnergyAfter(maxenrec,stepsize);
  he3dee->SetName(Form("he3_dee"));
  he3dee->Write("",TObject::kOverwrite);

  protonE->SetTargetThickness(40000*4.51/10.);
  deuteronE->SetTargetThickness(40000*4.51/10.);
  tritonE->SetTargetThickness(40000*4.51/10.);
  alphaE->SetTargetThickness(40000*4.51/10.);
  he3E->SetTargetThickness(40000*4.51/10.);
  
  //cout << protondee->GetNp() << " " << protonerest->GetNp() << endl;
  vector<double> de;
  vector<double> er;
  vector<double> de2;
  vector<double> er2;
  TGraph* graph;
  int stop;

  //proton
  de.resize(protondee->GetNp());
  er.resize(protondee->GetNp());
  de2.resize(protondee->GetNp());
  er2.resize(protondee->GetNp());
  stop = 0;
  for(int i=0;i<protondee->GetNp();i++){
    protondee->GetKnot(i, er[i], de[i]);
    er2[i] = protonE->EnergyAfter(er[i],-5);
    //cout << de[i] << " " << er[i] << "   " << er2[i] << endl;
    if(er2[i]<0.01)
      stop++;
    else
      er[i] = er[i]-er2[i];
  }
  delete graph;
  graph = new TGraph(stop, &er[0], &de[0]);
  TSpline3* protonstop = new TSpline3("protonstop",graph);
  protonstop->SetName(Form("protonstop"));
  protonstop->Write("",TObject::kOverwrite);
  for(int i=0;i<protonstop->GetNp();i++)
    protonstop->GetKnot(i, er[protonstop->GetNp()-i-1], de[protonstop->GetNp()-i-1]);
  delete graph;
  graph = new TGraph(stop, &de[0], &er[0]);
  protonstop = new TSpline3("protonstop_inv",graph);
  protonstop->SetName(Form("protonstop_inv"));
  protonstop->Write("",TObject::kOverwrite);
  
  if(protondee->GetNp()>stop){
    delete graph;
    graph = new TGraph(protondee->GetNp()-stop, &er[stop], &de[stop]);
    TSpline3* protonpunch = new TSpline3("protonpunch",graph);
    protonpunch->SetName(Form("protonpunch"));
    protonpunch->Write("",TObject::kOverwrite);
    
    for(int i=0;i<protonpunch->GetNp();i++)
      protonpunch->GetKnot(i, er[protonpunch->GetNp()-i-1], de[protonpunch->GetNp()-i-1]);
    delete graph;
    graph = new TGraph(protonpunch->GetNp(), &de[0], &er[0]);
    protonpunch = new TSpline3("protonpunch_inv",graph);
    protonpunch->SetName(Form("protonpunch_inv"));
    protonpunch->Write("",TObject::kOverwrite);
    
  }
  //deuteron
  de.resize(deuterondee->GetNp());
  er.resize(deuterondee->GetNp());
  de2.resize(deuterondee->GetNp());
  er2.resize(deuterondee->GetNp());
  stop = 0;
  for(int i=0;i<deuterondee->GetNp();i++){
    deuterondee->GetKnot(i, er[i], de[i]);
    er2[i] = deuteronE->EnergyAfter(er[i],-5);
    //cout << de[i] << " " << er[i] << "   " << er2[i] << endl;
    if(er2[i]<0.01)
      stop++;
    else
      er[i] = er[i]-er2[i];
  }
  delete graph;
  graph = new TGraph(stop, &er[0], &de[0]);
  TSpline3* deuteronstop = new TSpline3("deuteronstop",graph);
  deuteronstop->SetName(Form("deuteronstop"));
  deuteronstop->Write("",TObject::kOverwrite);
  for(int i=0;i<deuteronstop->GetNp();i++)
    deuteronstop->GetKnot(i, er[deuteronstop->GetNp()-i-1], de[deuteronstop->GetNp()-i-1]);
  delete graph;
  graph = new TGraph(stop, &de[0], &er[0]);
  deuteronstop = new TSpline3("deuteronstop_inv",graph);
  deuteronstop->SetName(Form("deuteronstop_inv"));
  deuteronstop->Write("",TObject::kOverwrite);

  if(deuterondee->GetNp()>stop){
    delete graph;
    graph = new TGraph(deuterondee->GetNp()-stop, &er[stop], &de[stop]);
    TSpline3* deuteronpunch = new TSpline3("deuteronpunch",graph);
    deuteronpunch->SetName(Form("deuteronpunch"));
    deuteronpunch->Write("",TObject::kOverwrite);

    for(int i=0;i<deuteronpunch->GetNp();i++)
      deuteronpunch->GetKnot(i, er[deuteronpunch->GetNp()-i-1], de[deuteronpunch->GetNp()-i-1]);
    delete graph;
    graph = new TGraph(deuteronpunch->GetNp(), &de[0], &er[0]);
    deuteronpunch = new TSpline3("deuteronpunch_inv",graph);
    deuteronpunch->SetName(Form("deuteronpunch_inv"));
    deuteronpunch->Write("",TObject::kOverwrite);
  }
  //triton
  de.resize(tritondee->GetNp());
  er.resize(tritondee->GetNp());
  de2.resize(tritondee->GetNp());
  er2.resize(tritondee->GetNp());
  stop = 0;
  for(int i=0;i<tritondee->GetNp();i++){
    tritondee->GetKnot(i, er[i], de[i]);
    er2[i] = tritonE->EnergyAfter(er[i],-5);
    //cout << de[i] << " " << er[i] << "   " << er2[i] << endl;
    if(er2[i]<0.01)
      stop++;
    else
      er[i] = er[i]-er2[i];
  }
  delete graph;
  graph = new TGraph(stop, &er[0], &de[0]);
  TSpline3* tritonstop = new TSpline3("tritonstop",graph);
  tritonstop->SetName(Form("tritonstop"));
  tritonstop->Write("",TObject::kOverwrite);
  for(int i=0;i<tritonstop->GetNp();i++)
    tritonstop->GetKnot(i, er[tritonstop->GetNp()-i-1], de[tritonstop->GetNp()-i-1]);
  delete graph;
  graph = new TGraph(stop, &de[0], &er[0]);
  tritonstop = new TSpline3("tritonstop_inv",graph);
  tritonstop->SetName(Form("tritonstop_inv"));
  tritonstop->Write("",TObject::kOverwrite);  
  
  if(tritondee->GetNp()>stop){
    delete graph;
    graph = new TGraph(tritondee->GetNp()-stop, &er[stop], &de[stop]);
    TSpline3* tritonpunch = new TSpline3("tritonpunch",graph);
    tritonpunch->SetName(Form("tritonpunch"));
    tritonpunch->Write("",TObject::kOverwrite);

    for(int i=0;i<tritonpunch->GetNp();i++)
      tritonpunch->GetKnot(i, er[tritonpunch->GetNp()-i-1], de[tritonpunch->GetNp()-i-1]);
    delete graph;
    graph = new TGraph(tritonpunch->GetNp(), &de[0], &er[0]);
    tritonpunch = new TSpline3("tritonpunch_inv",graph);
    tritonpunch->SetName(Form("tritonpunch_inv"));
    tritonpunch->Write("",TObject::kOverwrite);
  }
  //alpha
  de.resize(alphadee->GetNp());
  er.resize(alphadee->GetNp());
  de2.resize(alphadee->GetNp());
  er2.resize(alphadee->GetNp());
  stop = 0;
  for(int i=0;i<alphadee->GetNp();i++){
    alphadee->GetKnot(i, er[i], de[i]);
    er2[i] = alphaE->EnergyAfter(er[i],-5);
    //cout << de[i] << " " << er[i] << "   " << er2[i] << endl;
    if(er2[i]<0.01)
      stop++;
    else
      er[i] = er[i]-er2[i];
  }
  delete graph;
  graph = new TGraph(stop, &er[0], &de[0]);
  TSpline3* alphastop = new TSpline3("alphastop",graph);
  alphastop->SetName(Form("alphastop"));
  alphastop->Write("",TObject::kOverwrite);
  for(int i=0;i<alphastop->GetNp();i++)
    alphastop->GetKnot(i, er[alphastop->GetNp()-i-1], de[alphastop->GetNp()-i-1]);
  delete graph;
  graph = new TGraph(stop, &de[0], &er[0]);
  alphastop = new TSpline3("alphastop_inv",graph);
  alphastop->SetName(Form("alphastop_inv"));
  alphastop->Write("",TObject::kOverwrite);  
  
  if(alphadee->GetNp()>stop){
    delete graph;
    graph = new TGraph(alphadee->GetNp()-stop, &er[stop], &de[stop]);
    TSpline3* alphapunch = new TSpline3("alphapunch",graph);
    alphapunch->SetName(Form("alphapunch"));
    alphapunch->Write("",TObject::kOverwrite);

    for(int i=0;i<alphapunch->GetNp();i++)
      alphapunch->GetKnot(i, er[alphapunch->GetNp()-i-1], de[alphapunch->GetNp()-i-1]);
    delete graph;
    graph = new TGraph(alphapunch->GetNp(), &de[0], &er[0]);
    alphapunch = new TSpline3("alphapunch_inv",graph);
    alphapunch->SetName(Form("alphapunch_inv"));
    alphapunch->Write("",TObject::kOverwrite);
  }
  //he3
  de.resize(he3dee->GetNp());
  er.resize(he3dee->GetNp());
  de2.resize(he3dee->GetNp());
  er2.resize(he3dee->GetNp());
  stop = 0;
  for(int i=0;i<he3dee->GetNp();i++){
    he3dee->GetKnot(i, er[i], de[i]);
    er2[i] = he3E->EnergyAfter(er[i],-5);
    //cout << de[i] << " " << er[i] << "   " << er2[i] << endl;
    if(er2[i]<0.01)
      stop++;
    else
      er[i] = er[i]-er2[i];
  }
  delete graph;
  graph = new TGraph(stop, &er[0], &de[0]);
  TSpline3* he3stop = new TSpline3("he3stop",graph);
  he3stop->SetName(Form("he3stop"));
  he3stop->Write("",TObject::kOverwrite);
  for(int i=0;i<he3stop->GetNp();i++)
    he3stop->GetKnot(i, er[he3stop->GetNp()-i-1], de[he3stop->GetNp()-i-1]);
  delete graph;
  graph = new TGraph(stop, &de[0], &er[0]);
  he3stop = new TSpline3("he3stop_inv",graph);
  he3stop->SetName(Form("he3stop_inv"));
  he3stop->Write("",TObject::kOverwrite);  
  if(he3dee->GetNp()>stop){
    delete graph;
    graph = new TGraph(he3dee->GetNp()-stop, &er[stop], &de[stop]);
    TSpline3* he3punch = new TSpline3("he3punch",graph);
    he3punch->SetName(Form("he3punch"));
    he3punch->Write("",TObject::kOverwrite);

    for(int i=0;i<he3punch->GetNp();i++)
      he3punch->GetKnot(i, er[he3punch->GetNp()-i-1], de[he3punch->GetNp()-i-1]);
    delete graph;
    graph = new TGraph(he3punch->GetNp(), &de[0], &er[0]);
    he3punch = new TSpline3("he3punch_inv",graph);
    he3punch->SetName(Form("he3punch_inv"));
    he3punch->Write("",TObject::kOverwrite);
  }
  //protondeltaE->SetTargetThickness(65*2.33/10.);
  //cout << "400 "<< protondeltaE->EnergyAfter(400,-5) << endl;
  //cout << "300 "<< protondeltaE->EnergyAfter(300,-5) << endl;
  //cout << "200 "<< protondeltaE->EnergyAfter(200,-5) << endl;
  //cout << "100 "<< protondeltaE->EnergyAfter(100,-5) << endl;
  //cout << " 50 "<< protondeltaE->EnergyAfter(50,-5) << endl;
  outfile->Close();
}
