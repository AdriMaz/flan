#ifndef MATH_FUNCTION_H
#define MATH_FUNCTION_H

#include <Rcpp.h>
using namespace Rcpp;

class MATH_Integration {

private:

  double mReltol;
  
  double mSubd;

  Function* mIntegrate;

  Function* mIntegrand;

  List fcts;


protected:


public:

    MATH_Integration() {};
    MATH_Integration(List fns, double reltol, double subd){
      mReltol=reltol;
      mSubd=subd;
      mIntegrate=new Function("integrate");
//       List fns=Environment::global_env().get(".integrands");
      fcts=fns;

      mIntegrand=new Function("identity");  // allouer


    };

    ~MATH_Integration(){};

// List getFns(){return fcts;};

    void setFunctionName(std::string name){
    //   std::cout<<"reltol=<"<<mReltol<<">"<<std::endl;
    //       std::cout<<"Size de fcts="<<fcts.size()<<std::endl;
    //   Function func(fcts[name]);
      *mIntegrand = fcts[name];
    };

  /*
   * Integrals functions
   */

  double testintegralFunction(double a, double b,double rho, double delta);
  double integralFunction(double a, double b,double rho, double delta,double k);

};



class MATH_Polynom {
  
private:
  
  void init_all(int deg){
    setDegree(deg);
    init_fft();
  };
  
  void init_fft(){
    mFFT= new Function("fft");
  };
  
protected:
   std::vector<double> mCoef;
   int mDeg;
   
   Function* mFFT;


public:

    MATH_Polynom() {
      init_all(1);
    };
    
    MATH_Polynom(std::vector<double> C){
      setCoef(C);
      init_fft();
    }
    
    ~MATH_Polynom(){};
    
    
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
// 	std::cout<<"Avt *f : C["<<i<<"] ="<<*it<<std::endl;
	(*it)*=f;
// 	std::cout<<"Après *f : C["<<i<<"] ="<<*it<<std::endl;
      }
      return *this;
    };
    
    MATH_Polynom& operator+=(double f) {
      mCoef[0]+=f;
      return *this;
    };
    
    void reduce(double eps);

    void square_fft();
    
//     void square_conv();


};


// class MATH_ZeroEquation {
// 
// private: 
//   
//   Function* mSolver;
//   
//   Function* mFunc;
//   
//   NumericVector mInterval; 
// 
// public:
//   
//     MATH_ZeroEquation() {};
//     
//     
//     
//     MATH_ZeroEquation(double binf,double bsup){
// //       mReltol=reltol_;
//       mSolver=new Function("uniroot");
//       setInterval(binf,bsup);
// //       List fns=Environment::global_env().get(".integrands");
// //       fcts=fns;
//       mFunc=new Function("identity");  // allouer
//     };
// 
//     ~MATH_ZeroEquation(){};
//     
//     void setInterval(double binf,double bsup){
//       mInterval[0]=binf;
//       mInterval[1]=bsup;
//     };
//     
//     void setFunction(Function func){
//       *mFunc=func;
//     }
//     
//     NumericVector getInterval(){
//       return mInterval;
//     };
//     
//     
//     
//     
//     double unirootFunction();
//     
//     
//     
// };




#endif
