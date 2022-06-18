#ifndef Particle_h
#define Particle_h

#include "Vector3D.h"
#include <iostream>

using namespace std;

class Particle{
 public:
  Particle(const Vector3D& posizione);
  Particle(const double& x,const double& y,const double& z);
  Particle();

  Vector3D posizione()const{return posizione_;}
  void SetPos(const Vector3D& posizione);
  void RandomPos(const double& Lx, const double& Ly, const double& Lz);
  // double RandMarsenne();
  
  friend ostream& operator<<(ostream& os, const Particle& p);

  double x(){
    return posizione_.x();
  }
  double y(){
    return posizione_.y();
  }
  double z(){
    return posizione_.z();
  }


  void SetX(const double x);
  void SetY(const double y);
  void SetZ(const double z);
  
protected:
  Vector3D posizione_;
  
};
#endif
