
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;


double MATH_Integration::integralFunction(double a, double b,double rho,double delta,double k) {

  List res;
  if(k > 0) res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["subdivisions"] = mSubd,_["rho"] =rho,_["delta"] =delta,_["k"]=k);   
  else res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["subdivisions"] = mSubd,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
//   double integ=res["value"];
  return integ;
}



void MATH_Polynom::reduce(double eps){      
      int i=0,dmax=mCoef.size()-1;
      std::vector<double>::iterator it1=mCoef.begin();
      std::vector<double>::iterator itmax;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	if(*it <= eps) *it=0;
	if(*it > 0) {
	  itmax=it; 
	  dmax=i ;
	}

      }

      mCoef=std::vector<double>(it1,itmax+1);
      mDeg=dmax;
    };
    

void MATH_Polynom::square_fft(){
  
  int n=mCoef.size();
  
  std::vector<double> temp=mCoef;
  
  n*=2;
  n--;
   
  temp.resize(n);
  
  double tp=log2((double)(n));
  
  
  int N=floor(tp);
  
  if(N != tp) {
    N++;
    N=pow(2.,N);
    temp.resize(N);
  }
  
  
  
  
  ComplexVector CFFT=(*mFFT)(temp,_["inverse"]=false);
  
  Rcomplex c;
  ComplexVector::iterator it;
  for(it=CFFT.begin() ; it!=CFFT.end(); ++it) {
    c=*it;
    *it=c*c;
  }
  
  CFFT = (*mFFT)(CFFT,_["inverse"]=true);
  
  mCoef.resize(n);
  mDeg=n-1;
  
  it=CFFT.begin();
  for(std::vector<double>::iterator itC=mCoef.begin() ; itC!=mCoef.end() ; ++itC , ++it){
    *itC = (*it).r;
    (*itC)/=N;
  }
    
  
}

