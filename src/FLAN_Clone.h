/**
  * FLAN Software
  *
  * @author 2015-2020 Adrien Mazoyer  <adrien.mazoyer@imag.fr>
  * @see The GNU Public License (GPL)
  */
/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef FLAN_CLONE_H
#define FLAN_CLONE_H

#include "MATH_Function.h"

using namespace Rcpp ;

//#include "FLAN_Function.h"

// Lifetime Distributions class
class FLAN_Dist {

private:

  std::string mName;                // Name of the distribution
  std::vector<double> mParams;      // Parameter(s) of the distribution
  
public:

  FLAN_Dist(){};

  FLAN_Dist(List dist){
    std::string name=dist["name"];
    mName=name;
    if(mName.compare("lnorm") == 0) {
      mParams.resize(2);
      mParams[0]=dist["meanlog"];mParams[1]=dist["sdlog"];
//       init_integrator();
//       init_solver();
    } else if(mName.compare("gamma") == 0) {
      mParams.resize(2);
      mParams[0]=dist["shape"];mParams[1]=dist["scale"];
    } else if(mName.compare("exp") == 0) {
      mParams.resize(1);
      mParams[0]=dist["rate"];
    } else if(mName.compare("dirac") == 0) {
      mParams.resize(1);
      mParams[0]=dist["location"];
    }
  };

  ~FLAN_Dist(){};


  // Rescales the parameters to unit growth rate.
//   void adjustGrowthRate(double death){
// 
//     double m=2*(1-death);
//     
//     if(mName.compare("exp")==0){
//       
//       mParams[0]=2/m;
//       
//     } else if(mName.compare("dirac")==0){
//       
//       mParams[0]=log(m);
//       
//     } else 
//       if(mName.compare("lnorm")==0){
// 
//       double meanlog=mParams[0];
//       double sdlog=mParams[1];
//     
//       double s;
// 
//       // [...] cf source R
// 
//     } else if(mName.compare("gamma")==0){
//       double shape=mParams[0];
//       mParams[1]=pow(m,1./shape)-1;
//     }
// 
//   };

  std::string getDistName() {
    return mName;
  }

  std::vector<double> getDistParams() {
    return mParams;
  }

};


  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_SimClone class to compute sample of clone size
  */

class FLAN_SimClone {

private:

    double mFitness;
    double mDeath;
    FLAN_Dist *mDist;   // Lifetime distribution
    FLAN_Dist *mDistN;  // Lifetime distribution for normal cells (drawing function)

protected:
  static const double DEATH_EPS_SIM;     // Threshold for death

public:
    FLAN_SimClone(){};

    ~FLAN_SimClone(){};

    FLAN_SimClone(double rho,double death, List dist){

      mFitness=rho;
      mDeath=death;
      
      mDist= new FLAN_Dist(dist);
//       mDist->adjustGrowthRate(mDeath);

    };
    
    FLAN_SimClone(double rho,double death, List distn, List distm){

      mFitness=rho;
      mDeath=death;
      
      mDist= new FLAN_Dist(distm);
      mDistN= new FLAN_Dist(distn);
//       mDist->adjustGrowthRate(mDeath);

    };
      // create object to generate samples
    FLAN_SimClone(double rho,double death, FLAN_Dist *dist){

      mFitness=rho;
      mDeath=death;
      mDist=dist;

    };

    // -----------------------
    // Sample computation
    // ------------------------

    NumericVector computeSample(int n);
    
    int splitTimes(double t);
    
//     List splitTimes_draw(double t);

};


  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class for the distribution of a clone size
  */


class FLAN_Clone {

private:
  
//   void init_Zero(){
//     mZero=new MATH_ZeroEquation(0.01,100);
//   };


protected:
    static const double DEATH_EPS_DIST;     // Threshold for death
    double mFitness;    // Realtive fitness
    double mDeath;      // Death probability
    
//     MATH_ZeroEquation* mZero;

    FLAN_Clone(){};

    // create object for GF method
    FLAN_Clone(double death){
      mDeath=death;
      
//       init_Zero();
    };

    // create object to compute distribution
    FLAN_Clone(double rho,double death){
      mFitness=rho;
      mDeath=death;
      
//       init_Zero();
    };


    ~FLAN_Clone(){};


    // Set attributes
    void setFitness(double rho){
      mFitness=rho;
    }
    void setDeath(double death){
      mDeath=death;
    }


    // Get attributes

    double getFitness(){
      return mFitness;
    }
    double getDeath(){
      return mDeath;
    }



public:

    /*! \brief compute the probability Pk[k]=k+1 , k in [0,pkMax[
     * @return true if the computing succeeds otherwise it returns false because of wrong parameter of the distribution
     */
//     virtual List get() = 0;
    virtual NumericVector computeProbability(int m) = 0;

    /*! \brief compute the first derivatives of pk with respect to the parameters
     * derivatives[i] is the 1er order derivatives with respect to parameter i
     * @return true if the computing succeeds otherwise it returns false because of wrong parameter of the distribution
     *
     */
    
    virtual List computeProbability1DerivativeRho(int m) = 0;


    /*/////////////////////////////////////////////////////////////////////////////////////////////////////
     * Function for GF estimates
     */////////////////////////////////////////////////////////////////////////////////////////////////////

//     bool solveCumulativeFunctionEquation(double z1,double z2,
// 					 double death,double y,
// 					 double rho);
//     /*compute function to solve root equation
//      */
//     double computeFunction(int x);


//     double solveGeneratingFunction(double y, double z1,double z2);
    
    virtual double computeGeneratingFunction(double z) {
      std::vector<double> Z (1);
      Z[0]=z;
      return computeGeneratingFunction2(mFitness,Z)[0];
    };
    
    virtual std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z) = 0;
    
    virtual double computeGeneratingFunction1DerivativeRho(double z) = 0;

};



 /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is suppsoed to be exponential
  */

// DEFINE_SPTR(FLAN_ExponentialClone);

class FLAN_ExponentialClone : public FLAN_Clone {
//   SP_OBJECT(FLAN_ExponentialClone);

  protected:



  private:


    MATH_Integration* mIntegrator;

    void init() {
      double flantol=Environment::global_env()[".flantol"];
      double flansubd=Environment::global_env()[".flansubd"];
      List fns=Environment::global_env().get(".integrands");

//       std::cout<<"Size ="<<integrands.size()<<std::endl;
      mIntegrator=new MATH_Integration(fns,flantol,flansubd);
    }

  public:



    FLAN_ExponentialClone():FLAN_Clone() {
      init();
    };
    FLAN_ExponentialClone(double death):FLAN_Clone(death) {
      init();
    };
    FLAN_ExponentialClone(double rho,double death):FLAN_Clone(rho,death) {
      init();
    };
    ~FLAN_ExponentialClone(){};

    /* Compute the probability P[X=k]
     * Stor it in vector mProb
     */
//     List get(){
// 	     return mIntegrator->getFns();
//     };

    NumericVector computeProbability(int m);
    
    List computeProbability1DerivativeRho(int m) ;

//     double computeGeneratingFunction(double z)  ;
    
    std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z);

    double computeGeneratingFunction1DerivativeRho(double z)  ;

};

//
//  /*//////////////////////////////////////////////////////////////////////////////////////////////
//   * FLAN_Clone class when the lifetime model is suppsoed to be constant
//   */


class FLAN_DiracClone : public FLAN_Clone {
//   SP_OBJECT(FLAN_DiracClone);

private:
  
  void init_death(){
      mPol=MATH_Polynom();
    };

protected:
    
    MATH_Polynom mPol;

public:

    FLAN_DiracClone():FLAN_Clone() {};
    FLAN_DiracClone(double death):FLAN_Clone(death) {
      if(death > 0) init_death();
    };
    FLAN_DiracClone(double rho,double death):FLAN_Clone(rho,death) {
      if(death > 0) init_death();
    };
    ~FLAN_DiracClone(){};


    
//     List get(){};

    NumericVector computeProbability(int m);
    
    List computeProbability1DerivativeRho(int m);

//     double computeGeneratingFunction(double z)  ;
    
    std::vector<double> computeGeneratingFunction2(double rho,std::vector<double> Z);

    double computeGeneratingFunction1DerivativeRho(double z)  ;

};
#endif
