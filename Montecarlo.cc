#ifndef Montecarlo_h
#define Montecarlo_h

#include <iostream>
#include <fstream>
#include <random>
#include "Vector3D.h"
#include "Particle.h"
#include "Molecule.h"
#include "HardEll.h"
#include "RandomMarsenne.cc"

using namespace std;

class Montecarlo{

protected:
  long int accepttot_;
  long int acceptvol_;
  long int acceptvolM_;
  vector<Vector3D> r0_;
private:
  double L_[3];
  vector<HardEll> hardell_;
  vector<double> gr_;
  double a_,b_,c_;

  
public: 
  
  Montecarlo(const vector<HardEll>& he, vector<double>& L,int nbins,double a,double b, double c){
    for(int i=0; i<he.size();i++){
      hardell_.push_back(he[i]);
      r0_.push_back(he[i].posizione());
    }
    L_[0]=L[0];
    L_[1]=L[1];
    L_[2]=L[2];
    a_=a;
    b_=b;
    c_=c;
    accepttot_=0;
    acceptvol_=0;
    acceptvolM_=0;
    for(int i=0; i<nbins;i++){
      gr_.push_back(0.0);
    }
  }
  
  void FillRet(){
    long int k=0;
    vector<double> x(3,0.0);
    HardEll temp;
    Vector3D d(0,0,0),distance(0,0,0);
    temp.SetA(a_);
    temp.SetB(b_);
    temp.SetC(c_);
    for(long int i=0;i<hardell_.size();i++){
      hardell_[i].RandomPos(L_[0],L_[1],L_[2]);
      hardell_[i].RandomOrient();
      hardell_[i].SetA(a_);
      hardell_[i].SetB(b_);
      hardell_[i].SetC(c_);
      while(k<=i-1 && i>0){
	k=0;
	for(long int j=i-1;j>=0;j--){
	  
	  distance=hardell_[i].posizione()-hardell_[j].posizione();

	  if(hardell_[i].EllOverlap(hardell_[j])==true){
	    
	    hardell_[i].RandomPos(L_[0],L_[1],L_[2]);
	    hardell_[i].RandomOrient();
	    
	  }else if(abs(distance.x())>L_[0]/2 || abs(distance.y())>L_[1]/2 || abs(distance.z())>L_[2]/2){
	    
	    d.setX(L_[0]*(int)(distance.x()/(L_[0]/2)));
	    d.setY(L_[1]*(int)(distance.y()/(L_[1]/2)));
	    d.setZ(L_[2]*(int)(distance.z()/(L_[2]/2)));

	    temp.SetPos(hardell_[j].posizione()+d);
	    temp.SetOrientation(hardell_[j].VectX(),hardell_[j].VectY(),hardell_[j].VectZ());
	    
	    if(hardell_[i].EllOverlap(temp)==true){
	      hardell_[i].RandomPos(L_[0],L_[1],L_[2]);
	      hardell_[i].RandomOrient();
	    }else{
	      k++;
	    }
	  }else{
	    k++;
	  }
	}
      }
      r0_[i]=hardell_[i].posizione();
    }   
  }
  
  
  HardEll  traslation(HardEll he,const double delta){
    double x[3];
    RandomMarsenne gen;
    x[0]=he.x()+delta*(gen.RandM01()-0.5);
    x[1]=he.y()+delta*(gen.RandM01()-0.5);
    x[2]=he.z()+delta*(gen.RandM01()-0.5);

    for(int g=0; g<3;g++){
      if(x[g]>L_[g]){
	x[g]=x[g]-L_[g];
      }else if(x[g]<0){
	x[g]=x[g]+L_[g];
      }
    }
    he.SetX(x[0]);
    he.SetY(x[1]);
    he.SetZ(x[2]);
    return he;
  }
  
  HardEll  rotation(HardEll he,const double gamma){
    Vector3D temp=he.RandomVect();
    Molecule mtemp;
    mtemp.SetOrientation(he.VectX(),he.VectY(),he.VectZ());
    Vector3D xprimo=(temp*gamma+mtemp.VectX())/(gamma*temp+mtemp.VectX()).mod();
    Vector3D normal=(mtemp.VectX()).VectorProduct(xprimo);
    normal=normal/normal.mod();
    double alpha=acos(xprimo.ScalarProduct(mtemp.VectX()))/((mtemp.VectX()).mod()*xprimo.mod());

    double y1p=(1-(1-cos(alpha))*(normal.z()*normal.z()+normal.y()*normal.y()))*(mtemp.VectY()).x() +(-sin(alpha)*normal.z()+(1-cos(alpha))*normal.x()*normal.y())*(mtemp.VectY()).y() +(sin(alpha)*normal.y()+(1-cos(alpha))*(normal.z()*normal.x()))*(mtemp.VectY()).z();
    double y2p= (sin(alpha)*normal.z()+(1-cos(alpha))*(normal.x()*normal.y()))*(mtemp.VectY()).x() +(1-(1-cos(alpha))*(normal.z()*normal.z()+normal.x()*normal.x()))*(mtemp.VectY()).y()+(-sin(alpha)*normal.x()+(1-cos(alpha))*normal.z()*normal.y())*((mtemp.VectY()).z());
    double y3p= (-sin(alpha)*normal.y()+(1-cos(alpha))*(normal.z()*normal.x()))*(mtemp.VectY()).x()+(sin(alpha)*normal.x()+(1-cos(alpha))*normal.z()*normal.y())*((mtemp.VectY()).y())+(1-(1-cos(alpha))*(normal.y()*normal.y()+normal.x()*normal.x()))*(mtemp.VectY().z());

    Vector3D yprimo(y1p,y2p,y3p);

    double z1p=(1-(1-cos(alpha))*(normal.z()*normal.z()+normal.y()*normal.y()))*(mtemp.VectZ()).x() +(-sin(alpha)*normal.z()+(1-cos(alpha))*normal.x()*normal.y())*(mtemp.VectZ()).y() +(sin(alpha)*normal.y()+(1-cos(alpha))*(normal.z()*normal.x()))*(mtemp.VectZ()).z();
    double z2p= (sin(alpha)*normal.z()+(1-cos(alpha))*(normal.x()*normal.y()))*(mtemp.VectZ()).x() +(1-(1-cos(alpha))*(normal.z()*normal.z()+normal.x()*normal.x()))*(mtemp.VectZ()).y()+(-sin(alpha)*normal.x()+(1-cos(alpha))*normal.z()*normal.y())*((mtemp.VectZ()).z());
    double z3p= (-sin(alpha)*normal.y()+(1-cos(alpha))*(normal.z()*normal.x()))*(mtemp.VectZ()).x()+(sin(alpha)*normal.x()+(1-cos(alpha))*normal.z()*normal.y())*((mtemp.VectZ()).y())+(1-(1-cos(alpha))*(normal.y()*normal.y()+normal.x()*normal.x()))*(mtemp.VectZ().z());

    Vector3D zprimo(z1p,z2p,z3p);

    he.SetOrientation(xprimo,yprimo,zprimo);
    return he;
  }

  
  HardEll checkMove(HardEll he, long int k, HardEll hepremove){
    long int j=0;
    Vector3D distance(0,0,0);
    while(j<hardell_.size()){
      if(j!=k){
	distance=he.posizione()-hardell_[j].posizione();
	if(he.EllOverlap(hardell_[j])==true){
	j=hardell_.size();
	he=hepremove;
	
	}else if(abs(distance.x())>L_[0]/2 || abs(distance.y())>L_[1]/2 || abs(distance.z())>L_[2]/2){
	  
	  he=checkborder(he,hardell_[j],hepremove,j);
	  
	}
      }
      j++;
    }
    return he;
  }
  
  HardEll checkborder(HardEll he1,HardEll he2, HardEll hepremove,long int& j){
    Vector3D d(0,0,0);
    Vector3D distance=he1.posizione()-he2.posizione();
    HardEll temp=he2;

    d.setX(L_[0]*(int)(distance.x()/(L_[0]/2)));
    d.setY(L_[1]*(int)(distance.y()/(L_[1]/2)));
    d.setZ(L_[2]*(int)(distance.z()/(L_[2]/2)));

    temp.SetPos(he2.posizione()+d);
    
    if(he1.EllOverlap(temp)==true){
      he1=hepremove;
      j=hardell_.size();
    }
    return he1;
  }
  
  void  MCsteps(long int McSteps,bool Rotaz, bool Volume,double delta,double gamma,double deltaVmax,double beta,double pressure,bool verbose, bool RadialFunction){
    RandomMarsenne gen;
    vector<double> rho(McSteps,0.0);
    long int k=0,trasl=0,rot=0,accept=0,acceptrot=0,voltype,trialVol=0;
    double movetype,acceptM=0,trialMrot=0,trialM=0,p=0,trialMvol=0,acceptMrot=0,acceptMvol=0;
    HardEll hepremove;

    for(long int i=0 ; i<McSteps; i++){
      k=(long int)round(gen.RandM01()*(hardell_.size()-1));
      voltype=(long int)round(gen.RandM01()*(hardell_.size()-1));
      movetype=gen.RandM01();
      if(Volume==true && voltype==1){//cosi ho prob 1/N
	volmv(McSteps,beta,pressure,deltaVmax);
	trialVol++;
	trialMvol++;
      }else{
	if(movetype<0.5 || Rotaz==false) {
	  trasl++;
	  trialM++;
	  hepremove=hardell_[k];
	  hardell_[k]=traslation(hardell_[k],delta);
	  hardell_[k]=checkMove(hardell_[k],k,hepremove);
	  if((hardell_[k].posizione()!=hepremove.posizione())==true){
	    accept++;
	    acceptM++;
	  }
	}else{
	  rot++;
	  trialMrot++;
	  hepremove=hardell_[k];
	  hardell_[k]=rotation(hardell_[k],gamma);
	  hardell_[k]=checkMove(hardell_[k],k,hepremove);

	  
	  if((hardell_[k].VectX()!=hepremove.VectX())==true && (hardell_[k].VectY()!=hepremove.VectY())==true && (hardell_[k].VectZ()!=hepremove.VectZ())==true){
	    acceptrot++;
	    acceptMrot++;
	  }
	}
      }
      
      if(i%2000==0 && i!=0 && p==0 && Volume==true){
	if(EquilibrationCriterion(L_[0]/2)==true){
	  if(verbose==true){
	    cout<<"il sistema è equilibrato, fine fase di calibrazione."<<endl;
	    cout<<"msd : "<< msd()<<endl;
	  }
	  p=1;
	  accept=0;
	  trasl=0;
	  acceptrot=0;
	  rot=0;
	  acceptvol_=0;
	  trialVol=0;
	  cout<<"Inizio calcolo della densità e della g(r)"<<endl;
	  cout<<"delta trovato: "<<delta<<endl;
	  cout<<"gamma trovato: "<<gamma<<endl;
	  cout<<"Delta Vmax trovato: "<<deltaVmax<<endl;
	  
	}else{
	  if(verbose==true){
	    cout<<"il sistema non è equilibrato, continuo la fase preliminare.. "<<endl;
	    cout<<"msd:"<<msd()<<endl;
	  }
	  
	  if(acceptM>trialM/5*3){
	    delta=delta*1.1;
	    acceptM=0;
	    trialM=0;
	  }else if(acceptM<trialM/5*2){
	    delta=delta/1.1;	    
	    acceptM=0;
	    trialM=0;
	  }else{
	    cout<<"delta trovato:"<<delta<< " accepted: "<<acceptM<< " trials: "<<trialM<<endl;
	    acceptM=0;
	    trialM=0;
	  }	  
	  if(acceptMrot>trialMrot/5*3){
	    gamma=gamma*1.1;
	    
	    acceptMrot=0;
	    trialMrot=0;
	  }else if(acceptMrot<trialMrot/5*2){
	    gamma=gamma/1.1;	     
	    acceptMrot=0;
	    trialMrot=0;
	  }else{
	    cout<<"gamma trovato:"<<gamma<< " accepted: "<<acceptMrot<< " trials: "<<trialMrot<<endl;
	    acceptMrot=0;
	    trialMrot=0;
	  }	  
	  if(acceptvolM_>trialMvol/5*3 ){
	    deltaVmax=deltaVmax*1.1;
	    acceptvolM_=0;
	    trialMvol=0;
	  }else if(acceptvolM_<trialMvol/5*2){
	    deltaVmax=deltaVmax/1.1;
	    acceptvolM_=0;
	    trialMvol=0;
	  }else{	    
	    cout<<"Delta Vmax trovato:"<<deltaVmax << " accepted: "<<acceptvolM_<< " trials: "<<trialMvol<<endl;
	    acceptvolM_=0;
	    trialMvol=0;	     
	  }
	}
      }else if(p==1 && Volume==true){
	rho[i]=hardell_.size()/(L_[0]*L_[1]*L_[2]);
      }
      if(RadialFunction==true){
	gr(McSteps);
      }
    }
    accepttot_=accept+acceptrot+acceptvol_;
    if(verbose==true){
      cout<<"dati di fine simulazione"<<endl;
      cout<<"traslational: accepted: "<<accept<< " trial: "<< trasl <<endl;
      cout<<"Rotational:   accepted: "<<acceptrot<< " trial: "<< rot <<endl;
      cout<<"volume accepted: " <<acceptvol_ <<" trial: " <<trialVol <<endl;
      cout<<"accept totali : "<<accepttot_<<endl;
    }
    if(Volume==true){
    Printrho(McSteps,rho,pressure);
    }
  }
    
  
  void Printrho(long int McSteps, vector<double>& rho, double pressure){
    ofstream file;
    string dens="density",dat=".dat";
    string pres=to_string(pressure);
    double MeanRho=0;
    file.open(dens+pres+dat);
      cout<<"file :"<<dens+pres+dat<<" creato."<<endl;
      for(int i=0;i<McSteps;i++){
	file<<i<<" "<<rho[i]<<endl;
	MeanRho=MeanRho+rho[i];
      }
      file.close();
      MeanRho=MeanRho/McSteps;
      cout<<"rho media:"<< MeanRho<<endl;
    } 

  
    void PrintConfiguration(string nomefile){
      ofstream file;
      file.open(nomefile);
      for(int i=0;i<hardell_.size();i++){
	file<<hardell_[i].x()<<" "<< hardell_[i].y() << " " << hardell_[i].z()<<endl;
      }
      file.close();
    }

 void volmv(long int McSteps,double beta,double pressure,double deltaVmax){
   //vero solo nel caso di cubo
   double  volIniz=L_[0]*L_[1]*L_[2],newVol=0;
   if(deltaVmax>volIniz){
     deltaVmax=volIniz;
   }
   bool reject=false;
   double arg=0,maximum=0;
   RandomMarsenne gen;
   newVol=volIniz+(gen.RandM01()-0.5)*deltaVmax;
   double newL=pow(newVol,1./3.);;
   for(int i=0; i<hardell_.size();i++){
     hardell_[i].SetPos(hardell_[i].posizione()*(newL/L_[0]));
   }
   for(int j=0;j<hardell_.size();j++){
     for(int k=0;k<j;k++){
       if(hardell_[j].EllOverlap(hardell_[k])==true){
	 //reject the move
	 reject=true;
       }
     }
   }
   arg=-beta*(pressure*(newVol-volIniz))+hardell_.size()*(log(newVol/volIniz));
   if((exp(arg)<1)){
     maximum=exp(arg);
   }else{
     maximum=1;
   }
   double random01=gen.RandM01();
   if(random01>maximum || reject==true){
     //restoring the old configuration
     for(int m=0; m<hardell_.size();m++){
       hardell_[m].SetPos(hardell_[m].posizione()*(L_[0]/newL));
     }
   }else{
     L_[0]=newL;
     L_[1]=newL;
     L_[2]=newL;
     acceptvol_++;
     acceptvolM_++;
   }
 }
  

 void gr(long int number){
   int nbins= gr_.size();
   vector<double>hist(nbins,0.0);
   double rho=hardell_.size()/(L_[0]*L_[1]*L_[2]);
   Vector3D distance(0,0,0),d(0,0,0);
   double width=(sqrt(3)*L_[0]/2)/nbins;
   long int m=0;
   double deltaV=0;
   for(int i=0;i<hardell_.size();i++){
     for(int j=i+1;j<hardell_.size();j++){
	
       distance=hardell_[i].posizione()-hardell_[j].posizione();

       d.setX(L_[0]*(int)(distance.x()/(L_[0]/2)));
       d.setY(L_[0]*(int)(distance.y()/(L_[0]/2)));
       d.setZ(L_[0]*(int)(distance.z()/(L_[0]/2)));

       m=round((distance-d).mod()/width);

       deltaV=4*M_PI*pow((distance-d).mod(),2)*width;
       hist[m]=hist[m]+1/(rho*deltaV*hardell_.size());
     }
   }
   for(int i=0; i<nbins;i++){
     gr_[i]=gr_[i]+hist[i];
   }
 }
  
  void printgr(){
    double width=(sqrt(3)*L_[0]/2)/gr_.size();
    ofstream file;
    cout<<"file 'grclass.dat' creato"<<endl;
    file.open("grclass.dat");
    for(int i=0; i<gr_.size();i++){
      file<<i*width<<" "<< gr_[i]<<endl;
    }
    file.close();
  }
  
  double msd(){
    double meansquare=0;
    for(int i=0;i<hardell_.size();i++){
      meansquare= meansquare+pow((hardell_[i].posizione()-r0_[i]).mod(),2)/hardell_.size();
    }
    return meansquare;
  }
  
  bool EquilibrationCriterion(double sigma){
    if(msd()>pow(sigma,2)){
      return true;
    }else{
      return false;
    }
  }
  
};

#endif
