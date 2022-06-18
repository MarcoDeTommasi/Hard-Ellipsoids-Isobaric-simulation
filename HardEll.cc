#include "HardEll.h"
#include <iostream>
#include "he_overlap.c"

HardEll::HardEll():Molecule(0,0,0){};

HardEll::HardEll(const Vector3D& posizione, const double& a,const double& b, const double& c ):Molecule(posizione){
  saxi_[0]=a;
  saxi_[1]=b;
  saxi_[2]=c;
}

//Importante: devo definire necessariamente la matrice di orientazione chiamando CalcoloU prima di chimare questa altrimenti non ha senso
bool HardEll::EllOverlap(HardEll he){
  double r[3],rhe[3];
  r[0]=posizione_.x();
  r[1]=posizione_.y();
  r[2]=posizione_.z();
  rhe[0]=(he.posizione()).x();
  rhe[1]=(he.posizione()).y();
  rhe[2]=(he.posizione()).z();
  if(check_overlap_pw_c(saxi_,r,OrientMatrix_,he.saxi_,rhe,he.OrientMatrix_)<0.0){
    return true;
  }else{
    return false;
  }
}

void HardEll::SetA(const double& a){
  saxi_[0]=a;
}
void HardEll::SetB(const double& b){
  saxi_[1]=b;
}
void HardEll::SetC(const double& c){
  saxi_[2]=c;
}

Vector3D HardEll::saxi()const{
  Vector3D temp(saxi_[0],saxi_[1],saxi_[2]);
  return temp;
}

const HardEll& HardEll::operator=(const HardEll& rhs){
  saxi_[0]= rhs.saxi_[0];
  saxi_[1]= rhs.saxi_[1];
  saxi_[2]= rhs.saxi_[2];
  x_[0]=rhs.x_[0];
  x_[1]=rhs.x_[1];
  x_[2]=rhs.x_[2];
  posizione_=rhs.posizione_;
  CalcoloU();
  return *this;
}


