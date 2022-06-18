#ifndef vexcl_h
#define vexcl_h


#include "HardEll.h"
#include <vector>
#include <fstream>


class vexcl{
  double  L[3];
  long int N;
  double X0;
  
  vexcl(const double* taglia[3], const long int Nelem, const double& X0t);
  vexcl();
  
 public:
  vector<HardEll> part;

  void calc_vexcl();

  void setLength(double Lb);
  void setX0(const double& X0t);
  void setTrials(long int Ntrials);
  
  
 
};

#endif
