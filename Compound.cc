#include "Compound.hh"

using namespace std;

#define debug

static char* massfile=(char*)"/home/wimmer/progs/eloss/mass.dat";
static double amu = 931.494043;
Compound::Compound(char* symbol){
  int length = strlen(symbol);
  if(length == 0){
    cerr<<"error, type Material, excisting PE, DPE, MY, TTI and DTI"<<endl;
    exit(1);
  }
  fSymbol = symbol;
  if(strstr(symbol,"PE")){
    cout << "Polyethylene!" << endl;
    //H4C2
    SetNofElements(2);
    fNuclei = new Nucleus*[2];
    fFrac = new double[2];

    fNuclei[0] = new Nucleus(1,0,massfile);
    fNuclei[1] = new Nucleus(6,6,massfile);

    fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
    fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

    fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
  }
  else if(strstr(symbol,"DPE")){
    cout << "deuterated Polyethylene!" << endl;
    //D4C2
    SetNofElements(2);
    fNuclei = new Nucleus*[2];
    fFrac = new double[2];

    fNuclei[0] = new Nucleus(1,1,massfile);
    fNuclei[1] = new Nucleus(6,6,massfile);

    fFrac[0] = 4.*fNuclei[0]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());
    fFrac[1] = 2.*fNuclei[1]->GetMass()/(4.*fNuclei[0]->GetMass() + 2.*fNuclei[1]->GetMass());

    fMass = fNuclei[0]->GetMass()*4. + fNuclei[1]->GetMass()*2.;
  }
  else if(strstr(symbol,"MY")){
    cout << "Mylar!" << endl;
    //H8C10O4
    SetNofElements(3);
    fNuclei = new Nucleus*[3];
    fFrac = new double[3];

    fNuclei[0] = new Nucleus(1,0,massfile);
    fNuclei[1] = new Nucleus(6,6,massfile);
    fNuclei[2] = new Nucleus(8,8,massfile);

    fFrac[0] = 8.*fNuclei[0]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
    fFrac[1] = 10.*fNuclei[1]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());
    fFrac[2] = 4.*fNuclei[2]->GetMass()/(8.*fNuclei[0]->GetMass() + 10.*fNuclei[1]->GetMass() + 4.*fNuclei[2]->GetMass());

    fMass = fNuclei[0]->GetMass()*8. + fNuclei[1]->GetMass()*10. + fNuclei[2]->GetMass()*4.;


  }
  else if(strstr(symbol,"TTI")){
    cout << "Tritiated Titanium Target!" << endl;
    // ratioTTI = atomic ratio Tritium/Titanium
    if(isalpha(symbol[0])){
      cerr<<"give atomic ratio of Tritium to Titanium!"<< endl;
      exit(1);
    }
    else{
      double ratio = atof(symbol);
      //cout << ratio << endl;
      SetNofElements(2);
      fNuclei = new Nucleus*[2];
      fFrac = new double[2];
      fNuclei[0] = new Nucleus(1,2,massfile);
      fNuclei[1] = new Nucleus(22,26,massfile);
      fFrac[0] = ratio*fNuclei[0]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
      fFrac[1] = fNuclei[1]->GetMass()/(ratio*fNuclei[0]->GetMass()+fNuclei[1]->GetMass());
      
      fMass = fNuclei[0]->GetMass()*ratio + fNuclei[1]->GetMass();
    }
      
  }
  else if(strstr(symbol,"DTI")){
    //cout << "Deuterated Titanium Target!" << endl;
    if(isalpha(symbol[0])){
      cerr<<"give atomic ratio of Deuterium to Titanium!"<< endl;
      exit(1);
    }
    else{
      double ratio = atof(symbol);
      cout << ratio << endl;
      SetNofElements(2);
      fNuclei = new Nucleus*[2];
      fFrac = new double[2];

      fNuclei[0] = new Nucleus(1,1,massfile);
      fFrac[0] = ratio/(1+ratio);
      fNuclei[1] = new Nucleus(22,26,massfile);
      fFrac[1] = 1/(1+ratio);

      fMass = fNuclei[0]->GetMass()*fFrac[0] + fNuclei[1]->GetMass()*fFrac[1];
    }
  }
  else{
    cerr<<"Compound not implemented yet!"<< endl;
    exit(1);
  }
}
Compound::Compound(Nucleus* target){
  SetNofElements(1);
  fNuclei = new Nucleus*[1];
  fFrac = new double[1];
  fNuclei[0] = target;
  fFrac[0] = 1;
  
}
Compound::Compound(Nucleus* n1, double f1, Nucleus* n2, double f2){
  SetNofElements(2);
  fNuclei = new Nucleus*[2];
  fFrac = new double[2];
  fNuclei[0] = n1;
  fNuclei[1] = n2;

  fFrac[0] = f1*fNuclei[0]->GetMass()/(f1*fNuclei[0]->GetMass() + f2*fNuclei[1]->GetMass());
  fFrac[1] = f2*fNuclei[1]->GetMass()/(f1*fNuclei[0]->GetMass() + f2*fNuclei[1]->GetMass());

  fMass = fNuclei[0]->GetMass()*fFrac[0] + fNuclei[1]->GetMass()*fFrac[1];
  //cout << fNuclei[0]->GetSymbol() << " " << fNuclei[0]->GetMass() << " " << fFrac[0] << endl;
  //cout << fNuclei[1]->GetSymbol() << " " << fNuclei[1]->GetMass() << " " << fFrac[1] << endl;
  //cout << fMass << endl;
}
Compound::~Compound(){
  if(fFrac!=NULL)
    delete[] fFrac;
  if(fNuclei!=NULL){
    for(int i=0;i<GetNofElements();i++){
      if(fNuclei[i]!=NULL)
	delete fNuclei[i];
    }
    delete[] fNuclei;
  }
  if(fSymbol!=NULL)
    delete[] fSymbol;
}
Nucleus* Compound::GetNucleus(int i){
  if(i<GetNofElements())
    return fNuclei[i];
  else
    return NULL;
}
double Compound::GetFrac(int i){
  //cout << "i" << i << "GetNofElements()" << GetNofElements() << endl;
  if(i<GetNofElements())
    return fFrac[i];
  else
    return 0;
}
