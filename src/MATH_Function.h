#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include <Rcpp.h>
using namespace Rcpp;

class MATH_Integration {

private:

  double mReltol;
  
  int mSubd;

  Function* mIntegrate=NULL;

  Function* mIntegrand=NULL;

  List fcts;


protected:


public:

    MATH_Integration() {};
    MATH_Integration(List fns, double reltol, int subd){
      mReltol=reltol;
      mSubd=subd;
      
      if(!mIntegrate) delete mIntegrate;
      mIntegrate=new Function("integrate");
      fcts=fns;

      if(!mIntegrand) delete mIntegrand;
      mIntegrand=new Function("identity");  // allouer


    };

    ~MATH_Integration(){
      if(!mIntegrate) delete mIntegrate;
      if(!mIntegrand) delete mIntegrand;
    };


    void setFunctionName(std::string name){

      *mIntegrand = fcts[name];
    };

  /*
   * Integrals functions
   */
  double integralFunction(double a, double b,double rho, double delta,double k);

};



class MATH_Polynom {
  
private:
  
  void init_all(int deg){
    setDegree(deg);
    init_fft();
  };
  
  void init_fft(){
    if(!mFFT) delete mFFT;
    mFFT= new Function("fft");
  };
  
protected:
   std::vector<double> mCoef;
   int mDeg;
   
   Function* mFFT=NULL;


public:

    MATH_Polynom() {
      init_all(1);
    };
    
    MATH_Polynom(std::vector<double> C){
      setCoef(C);
      init_fft();
    }
    
    ~MATH_Polynom(){
      if(!mFFT) delete mFFT;
    };
    
    
    void setDegree(int deg){
      mCoef.resize(deg+1);
      mCoef[deg]=1;
      mDeg=deg;
    };
    
    void setCoef(std::vector<double> C){
      mCoef=C;
      mDeg=C.size()-1;
    };
    
    
    std::vector<double> getCoef(){
      return mCoef;
    }
    
    int getDegree(){
      return mDeg;
    };
    
    double& operator[](int i) {
      return mCoef[i];
    };
    
    const double& operator[](int i) const {
      return mCoef[i];
    };
    
    MATH_Polynom& operator*=(double f) {
      int i=0;
      for(std::vector<double>::iterator it=mCoef.begin() ; it!=mCoef.end() ; ++it,i++) {
	(*it)*=f;
      }
      return *this;
    };
    
    MATH_Polynom& operator+=(double f) {
      mCoef[0]+=f;
      return *this;
    };
    
    void reduce(double eps);

    void square_fft();
    

};



#endif
