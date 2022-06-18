#include "Montecarlo.cc"
#include "vexcl.cc"

int main(){
 
   vector<double> L(3);
    L[0]=20;
    L[1]=20;
    L[2]=20;
  
    long int McSteps=500000;
    double beta=1;
    long int N=250;
    //double pressure=0.015;
    int nbins=200;
    vector<double> hist(nbins,0.0);
    vector<HardEll> he(N);;
    double phi=N/(L[0]*L[1]*L[2])*4/3*M_PI*2;
    double a=2;
    double b=1;
    double c=1;
    bool verbose=true;

    //calcolo volume escluso
    long int trials=100000;
    vexcl v1(L, trials, 2 );
    cout<<"volume escluso dell'ellissoide con Xo="<<a<<" :"<<v1.calc_vexcl()<<endl;

    //calcolo della g(r)
    long int McSteps1=100000;
    double pres=0.015;
    Montecarlo m1(he,L,nbins,2,1,1);
    m1.FillRet();
    m1.MCsteps(McSteps1,true,false,1,0.3,900,1,pres,false, true);
    m1.printgr();
      
    for(int i=2; i<=12;i++){
      double pressure= i*0.0025;
      Montecarlo m2(he,L,nbins,2,1,1);
      m2.FillRet();
      m2.MCsteps(McSteps,true,true,1,0.3,900,beta,pressure,verbose,false);
    
      }
}
