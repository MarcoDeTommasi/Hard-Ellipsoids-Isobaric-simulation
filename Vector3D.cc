#include <iostream>
#include <math.h>
#include "Vector3D.h"

using std::cout;
using std::endl;

Vector3D::Vector3D(const double& x,const double& y, const double& z){
  x_[0]=x;
  x_[1]=y;
  x_[2]=z;
}
Vector3D::Vector3D(const Vector3D& rhs){
  x_[0]=rhs.x_[0];
  x_[1]=rhs.x_[1];
  x_[2]=rhs.x_[2];
}
void Vector3D::setCoord(const double& x, const double& y,const double& z){
  x_[0]=x;
  x_[1]=y;
  x_[2]=z;
}

void Vector3D::setX(const double& x){
  x_[0]=x;
}

void Vector3D::setY(const double& y){
  x_[1]=y;
}

void Vector3D::setZ(const double& z){
  x_[2]=z;
}


double Vector3D::mod() const {
  return sqrt(x_[0]*x_[0]+x_[1]*x_[1]+x_[2]*x_[2]);
}
double Vector3D::ScalarProduct( const Vector3D& rhs)const{
 return x_[0]*rhs.x_[0]+x_[1]*rhs.x_[1]+x_[2]*rhs.x_[2];
}
Vector3D Vector3D::VectorProduct( const Vector3D& rhs) const{
 return  Vector3D(x_[1]*rhs.x_[2]-x_[2]*rhs.x_[1],x_[2]*rhs.x_[0]-x_[0]*rhs.x_[2],x_[0]*rhs.x_[1]-x_[1]*rhs.x_[0]);
 }
Vector3D Vector3D::operator+(const Vector3D& rhs) const{
  return Vector3D(x_[0]+rhs.x_[0], x_[1]+rhs.x_[1],x_[2]+rhs.x_[2]);
}
Vector3D Vector3D::operator-(const Vector3D& rhs) const{
  return Vector3D(x_[0]-rhs.x_[0], x_[1]-rhs.x_[1],x_[2]-rhs.x_[2]);
}
Vector3D Vector3D::operator*(const double& rhs) const{
  return Vector3D(x_[0]*rhs, x_[1]*rhs,x_[2]*rhs);
}
Vector3D Vector3D::operator/(const double& rhs) const{
  return Vector3D(x_[0]/rhs, x_[1]/rhs,x_[2]/rhs);
}
const Vector3D& Vector3D::operator=(const double& a) {
  x_[0]=a; 
  x_[1]=a;
  x_[2]=a;    
  return *this;
}
const Vector3D& Vector3D::operator=(const Vector3D& rhs) {
  x_[0]=rhs.x_[0];
  x_[1]=rhs.x_[1];
  x_[2]=rhs.x_[2];
  return *this;
}
const Vector3D& Vector3D::operator+=(const Vector3D& rhs) {
  x_[0]+=rhs.x_[0]; 
  x_[1]+=rhs.x_[1];
  x_[2]+=rhs.x_[2];    
  return *this;
}

bool Vector3D::operator!=(const Vector3D& rhs) {
  bool value;
  if(x_[0]==rhs.x_[0] && x_[1]==rhs.x_[1] && x_[2]==rhs.x_[2] ){
    value=false;
  }else{
    value=true;
    }
  return value;
}



std::ostream& operator<<(std::ostream& os, const Vector3D& rhs){
  using namespace std;
  os <<"("<<rhs.x_[0]<<","<<rhs.x_[1]<<","<<rhs.x_[2]<<")";
    return os;
    }

const Vector3D operator*(const double& lhs, const Vector3D& rhs){
  return Vector3D(rhs.x_[0]*lhs , rhs.x_[1]*lhs, rhs.x_[2]*lhs);
}
