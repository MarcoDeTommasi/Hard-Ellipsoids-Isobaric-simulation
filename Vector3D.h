#ifndef Vector3D_h
#define Vector3D_h
#include <ostream>

class Vector3D{
 public:
  Vector3D(const double& x=0, const double& y=0, const double& z=0);
  Vector3D(const Vector3D& vect);
  void setCoord(const double& x,const double& y, const double& z);
  void setX(const double& x);
  void setY(const double& y);
  void setZ(const double& z);
  double mod() const;
  double ScalarProduct(const Vector3D& rhs) const ;
  double x()const{return x_[0];}
  double y()const{return x_[1];}
  double z()const{return x_[2];}
  Vector3D VectorProduct(const Vector3D& rhs) const;
  Vector3D operator+(const Vector3D& rhs) const;
  Vector3D operator-(const Vector3D& rhs) const;
  const Vector3D& operator=(const double& a );
  const Vector3D& operator=(const Vector3D& rhs);
  const Vector3D& operator+=(const Vector3D& rhs);
  bool operator!=(const Vector3D& rhs);
  Vector3D operator*(const double& rhs) const; // vector times scalar
  Vector3D operator/(const double& rhs) const;  // vector divided by scalar 

  friend const Vector3D operator*(const double& lhs,const Vector3D& rhs);
  friend std::ostream& operator<<(std::ostream& os, const Vector3D& rhs);
  

 private:
  double x_[3];
};
#endif
