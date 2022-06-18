#ifndef Molecule_h
#define Molecule_h

#include "Vector3D.h"
#include "Particle.h"
#include "RandomMarsenne.cc"
//devo inserire unicamente l'orientazione tridimensionale

class Molecule: public Particle{
public:
  
  Molecule();
  Molecule(const Vector3D& pos);
  Molecule(const double& x, const double& y, const double& z);

  void RandomOrient();
  void CalcoloU();
  void PrintOrient() const;
  void SetOrientation(const Vector3D& x, const Vector3D& y, const Vector3D& z);

  Vector3D VectX(){return x_[0];}
  Vector3D VectY(){return x_[1];}
  Vector3D VectZ(){return x_[2];}
  
  
  Vector3D RandomVect();
  
protected:
  Vector3D x_[3];
  double OrientMatrix_[3][3];
 
  
};

#endif
