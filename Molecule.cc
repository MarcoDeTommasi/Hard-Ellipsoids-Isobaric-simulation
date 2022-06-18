#include "Molecule.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <math.h>

using std::endl;
using std::cout;

Molecule::Molecule():Particle(0,0,0){
  
};
Molecule::Molecule(const Vector3D& pos):Particle(pos){}
Molecule::Molecule(const double& x, const double& y, const double& z):Particle(x,y,z){}

//qui impongo che il vettore ottenuto sia normalizzato
Vector3D Molecule::RandomVect(){
  RandomMarsenne gen;
  double  xisq, xi1, xi2, xi;
   double v[3];
  xisq = 1.0;
  while (xisq >= 1.0)
    {
      xi1  =1.0 - 2.0*gen.RandM01();
      xi2  =1-0 - 2.0*gen.RandM01();
      xisq = xi1 * xi1 + xi2 * xi2;
    }

  xi = sqrt(fabs(1.0 - xisq));
  v[0] = 2.0 * xi1 * xi;
  v[1] = 2.0 * xi2 * xi;
  v[2] = 1.0 - 2.0 * xisq;

  Vector3D RandVect(v[0],v[1],v[2]);
  
  return RandVect;
}
//in questa funzione genero i tre vettori che mi danno il sistema di riferimento associato alla molecola, inoltre prendo una matrice (3x3) e vi copio dentro i risultati ottenuti
void Molecule::RandomOrient(){
  x_[0]=RandomVect();
  Vector3D temp(1,0,0), temp1(0,1,0);
  x_[1]=temp;
  if(x_[0]!=temp==false){
    x_[1]=temp1;
  }
  x_[1]=x_[1]-(x_[0].ScalarProduct(x_[1]))/(x_[0].ScalarProduct(x_[0])) *x_[0];
  x_[1]=x_[1]/x_[1].mod();
  x_[2]=x_[0].VectorProduct(x_[1]);
  x_[2]=x_[2]/x_[2].mod();
  //check norma del vettore
  /*for(int i=0;i<3;i++){
    cout<<"check vettore"<< x_[i]<<"modulo : "<< x_[i].mod()<<endl;  
    }*/
  CalcoloU();
}
void Molecule::SetOrientation(const Vector3D& x, const Vector3D& y, const Vector3D& z){
  x_[0]=x;
  x_[1]=y;
  x_[2]=z;
  CalcoloU();
}
void Molecule::CalcoloU(){
 for (int i=0;i<3;i++){
    OrientMatrix_[i][0]=x_[i].x();
    OrientMatrix_[i][1]=x_[i].y();
    OrientMatrix_[i][2]=x_[i].z();
  }
}
void Molecule::PrintOrient() const{
  cout<<"Matrice di Orientazione:"<<endl;
  for(int x=0;x<3;x++)  
    {
      for(int y=0;y<3;y++)  
        {
	  cout<<OrientMatrix_[x][y]<<" ";  
        }
    cout<<endl; 
    } 
}
