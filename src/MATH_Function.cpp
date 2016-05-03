
#include "MATH_Function.h"

#include <math.h>

using namespace Rcpp;


double MATH_Integration::testintegralFunction(double a, double b,double rho,double delta) {
  std::string name="CLONE_P0_WD";
  setFunctionName(name);
  List res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
  std::cout<<integ<<std::endl;

  name="CLONE_dP0_dr_WD";
  setFunctionName(name);
  res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  integ=as<double>(res["value"]);

  return integ;
}


double MATH_Integration::integralFunction(double a, double b,double rho,double delta,double k) {

  List res;
  if(k > 0) res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta,_["k"]=k);   
  else res=(*mIntegrate)(*mIntegrand, a, b, _["rel.tol"] = mReltol,_["rho"] =rho,_["delta"] =delta);
  double integ=as<double>(res["value"]);
//   double integ=res["value"];
  return integ;
}



void MATH_Polynom::reduce(double eps){      
      int i=0,dmax;
      std::vector<double>::iterator it1=mCoef.begin();
      std::vector<double>::iterator itmax;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
// 	std::cout<<"C["<<i<<"] ="<<*it<<std::endl;
	if(*it <= eps) *it=0;
	if(*it > 0) {
	  itmax=it; 
	  dmax=i ;
	}
// 	if(*it > 0) dmax=i;
      }
//       mDeg=dmax;
      
//       std::cout<<"Initial size of mCoef ="<<mCoef.size()<<std::endl;
//       for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it) {
// 	std::cout<<' '<<*it;
//       } 
//       std::cout<<' '<<std::endl;

//       std::cout<<"Resize mCoef"<<std::endl;
      mCoef=std::vector<double>(it1,itmax+1);
//       std::cout<<"new size of mCoef ="<<mCoef.size()<<std::endl;
//       for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it) {
// 	std::cout<<' '<<*it;
//       }
//       std::cout<<' '<<std::endl;
//       std::cout<<"new deg ="<<mCoef.size()-1<<std::endl;
//       mDeg=mCoef.size()-1;
      mDeg=dmax;
//       mCoef.resize(dmax+1);
    };
    
    
    
// void MATH_Polynom::square_conv(){
//   
//   int n=mCoef.size();
//   
//   n*=2;
//   n--;
//   
//   std::vector<double> temp(n);
//   
//   int i=0,j;
//   for(std::vector<double>::iterator it1=mCoef.begin() ; it1 != mCoef.end() ; ++it1, i++){
//     j=0;
//     for(std::vector<double>::iterator it2=mCoef.begin() ; it2 != mCoef.end() ; ++it2, j++){
//       temp[i+j]+=(*it1)*(*it2);
//     }
//   }
//   
//   setCoef(temp);
//   
// }


void MATH_Polynom::square_fft(){
  
  int n=mCoef.size();
  
  std::vector<double> temp=mCoef;
//   std::vector<double> temp2(n-1);
  
  n*=2;
  n--;
   
  temp.resize(n);
  
  double tp=log2(n);
  
  
  int N=floor(tp);
  
  if(N != tp) {
//     std::cout<<"Rajout de 0 numéro2"<<std::endl;
    N++;
    N=pow(2,N);
//     std::cout<<"N ="<<N<<std::endl;
    temp.resize(N);
//     for(std::vector<double>::iterator it=temp.begin()+n+1 ; it!=temp.end() ; ++it) *it =0;
  }
  
//   std::cout<<"size of temp = "<<temp.size()<<std::endl;
  
  
  
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


// NumericVector MATH_ZeroEquation::unirootFunction(){
//   
//   NumericVector output;
//   NumericVector Int=getInterval();
//   double binf=Int[0],bsup=Int[1];
//   
//   if((*mFunc)(binf)*(*mFunc)(bsup) > 0){
//     output[0]=1 ; output[1] = false;
//   } else {
//     List res=(*mSolver)(*mFunc,_["interval"]=mInterval);
//     output[0]=res["root"];output[1]=true;
//   }
//   return output; 
// }

