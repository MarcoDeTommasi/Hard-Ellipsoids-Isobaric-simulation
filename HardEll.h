#ifndef HardEll_h
#define HardEll_h

#include "Molecule.h"
#include <iostream>


class HardEll: public Molecule{
public:
  HardEll();
  HardEll(const Vector3D& posizione,const double& a=1, const double& b=1, const double& c=2);

  bool EllOverlap(HardEll he);
  
  void SetA(const double& a);
  void SetB(const double& b);
  void SetC(const double& c);
  
  double a()const{return saxi_[0];}
  double b()const{return saxi_[1];}
  double c()const{return saxi_[2];}
  Vector3D saxi()const;
  
 const HardEll& operator=(const HardEll& rhs);

private:
  double saxi_[3];
 
};



#endif
