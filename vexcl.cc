#ifndef vexcl_cc
#define vexcl_cc


#include "HardEll.h"
#include <vector>
#include <fstream>

using namespace std;

class vexcl{
private:
  vector<double> L;
  long int N;
  double X0;
  HardEll part[2];
  
 public:

  vexcl(vector<double>& taglia , const long int Ntrials, const double& X0t ){
    for(int i=0;i<3;i++){
      L.push_back(taglia[i]);
    }
    N=Ntrials;
    X0=X0t;
  }
  
  vexcl(){
    N=10000;
    for(int i=0;i<3;i++){
      L[i]=10;
    }
    X0=2;
  }
  

  void setLength(const vector<double>& taglia){
    for(int i=0;i<3;i++){
      L[i]=taglia[i];
    }
  }
  void setTrials(long int Nt){
    N = Nt;
  }
  
  void setX0(const double X0t){
    X0=X0t;
  }
  
 double calc_vexcl(){
  double m=0;
  Vector3D posPart1(L[0]/2,L[1]/2,L[2]/2);
  part[0].SetPos(posPart1);
  part[0].RandomOrient();
  part[0].SetA(X0);
  part[0].SetB(1);
  part[0].SetC(1);
  part[1].SetA(X0);
  part[1].SetB(1);
  part[1].SetC(1);
  for(int i=0; i<N; i++){
    part[1].RandomPos(L[0],L[1],L[2]);
    part[1].RandomOrient();
    if(part[0].EllOverlap(part[1])==true){
      m+=1;
    }
  }
  return (double)((m/N)*L[0]*L[1]*L[2]);
} 
  
  
};


#endif
