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
#include "Nucleus.hh"
#include "Reconstruction.hh"

using namespace std;

int main(int argc, char* argv[]){
  char* OutputFile = NULL;
  char* targetc = NULL;
  vector<int> target;
  vector<int> projectile;
  vector<double> enpu;
  double thick;

  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-o", "outputfile", &OutputFile);
  interface->Add("-p", "projectile N | Z", &projectile);
  interface->Add("-t", "target N | Z", &target);  
  interface->Add("-tc", "target compund", &targetc);  
  interface->Add("-e", "energy [MeV/u], or energy from, to, nsteps", &enpu);  
  interface->Add("-th", "target thickness [mg/cm^2]", &thick);  

  interface->CheckFlags(argc, argv);

  //TFile* outfile = new TFile(OutputFile,"recreate");
  //if(outfile->IsZombie()){
  //  return 4;
  //}
  char* massFile = (char*)"/home/wimmer/progs/reaction/mass.dat";

  Nucleus *proj;
  Nucleus *targ;

  if(projectile.size() == 2){
    //cout << "projectile Z " << projectile[1] << " N " << projectile[0] << endl;
    proj = new Nucleus(projectile[1],projectile[0], massFile);
  }
  else{
    //cerr<<"no or incorrect Projectile provided!";
    for(int i=0; i<projectile.size(); i++){
      cerr<<projectile[i]<<" ";
    }
    cerr<<endl;
    exit(1);
  }
  Compound *targetmat;
  if(target.size() == 2){
    //cout << "target Z " << target[1] << " N " << target[0] << endl;
    targ = new Nucleus(target[1],target[0], massFile);
    //cout << "target created starting compound" << endl; 
    targetmat = new Compound(targ);
    cout << "calculating energy loss of " <<  proj->GetSymbol() << " in " << thick << " mg/cm^2 " << targ->GetSymbol() << endl;
  }
  else if(targetc!=NULL){
    targetmat = new Compound(targetc);
    cout << "calculating energy loss of " <<  proj->GetSymbol() << " in " << thick << " mg/cm^2 " << targetmat->GetSymbol() << endl;
  }
  else{    
    cerr<<"flag -p provided but not two arguments following: ";
    for(int i=0; i<target.size(); i++){
      cerr<<target[i]<<" ";
    }
    exit(1);
  }
  Reconstruction *eloss = new Reconstruction(proj, targetmat);
  eloss->SetTargetThickness(thick);

  if(enpu.size()==1){
    cout << "beam energy " << enpu[0] << " [MeV/u] = " << enpu[0]*proj->GetA() << " [MeV]" << endl;
    eloss->Print(enpu[0]*proj->GetA(),enpu[0]*proj->GetA(),1);
  }
  else if(enpu.size()==3){
    cout << "beam energy from " << enpu[0] << " [MeV/u] = " << enpu[0]*proj->GetA() << " [MeV] to " << enpu[1] << " [MeV/u] = " << enpu[1]*proj->GetA() << " [MeV] in " << enpu[2] << " steps"<< endl;
    eloss->Print(enpu[0]*proj->GetA(),enpu[1]*proj->GetA(),enpu[2]);
  }
  else{
    cout << "error in energy input " << endl;
  }
  /*
  cout << enpu << " [MeV/u] = " << enpu*proj->GetA() << " [MeV]" << endl;
  cout << " after " << endl;
  cout << eloss->EnergyAfter(enpu*proj->GetA(),-5)/proj->GetA() << " [MeV/u] " << eloss->EnergyAfter(enpu*proj->GetA(),-5) << " [MeV]" << endl;
  //outfile->Close();
  cout << "-----------------------" << endl;
  eloss->Print(enpu*proj->GetA());
  */
  //cout << "-----------------------" << endl;
  //eloss->Print(100*proj->GetA(),150*proj->GetA(),6);
  //cout << "-----------------------" << endl;
  return 0;
}
