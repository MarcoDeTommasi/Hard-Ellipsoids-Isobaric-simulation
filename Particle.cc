#include "Particle.h"
#include <math.h>
#include <random>
#include "RandomMarsenne.cc"
using std::cout;
using std::endl;
using std::ostream;
using unidist=uniform_real_distribution<double>;

Particle::Particle( const Vector3D& posizione){
  posizione_=posizione;
}

Particle::Particle(const double& x, const double& y, const double& z ){
  posizione_.setCoord(x,y,z);

}

Particle::Particle(){
  posizione_.setCoord(0,0,0);
}

void Particle::SetPos(const Vector3D& posizione){
  posizione_=posizione;
}

/*double Particle::RandMarsenne(){
    unidist xsi(-1.0,1.0);
    return xsi(mt);
    }*/

void Particle::RandomPos(const double& Lx, const double& Ly, const double& Lz){
  RandomMarsenne gen;
    posizione_.setCoord(gen.RandM01()*Lx, gen.RandM01()*Ly,gen.RandM01()*Lz);
}

void Particle::SetX(const double x){
  posizione_.setX(x);
}

void Particle::SetY(const double y){
  posizione_.setY(y);
}
void Particle::SetZ(const double z){
  posizione_.setZ(z);
}


ostream& operator<<(ostream& os, const Particle& p){
  os<<"Particella nella posizione "<<p.posizione();
  return os;
}
