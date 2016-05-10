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


#include "FLAN_Clone.h"

const double FLAN_SimClone::DEATH_EPS_SIM=1.e-4;

// std::vector<double> FLAN_SimClone::computeSample(int n) {
NumericVector FLAN_SimClone::computeSample(int n) {

//   std::vector<double> sample(n);
  
//   RNGScope rngScope;

//   std::cout<<"RHO ="<<mFitness<<std::endl;
  std::string name=mDist->getDistName();
  
  NumericVector sample=rexp(n,mFitness);
  

  if(name.compare("dirac") == 0){

    // specialized method
    if (mDeath<DEATH_EPS_SIM) {
        // case without mDeath
      for(NumericVector::iterator it = sample.begin(); it != sample.end(); ++it) {
	  *it = pow(2.,floor(*it/log(2.)));
        }
    } else {
        // case with death
        double t,a=log(2*(1-mDeath));
	for(NumericVector::iterator it = sample.begin(); it != sample.end(); ++it) {
	  double cg=0;
	  int m=1;
	  int nAlives=0,nDeaths=0;
	  t=*it;
            while ((cg<t-a) && (m>0)) {
// 		nDeaths=R::rbinom(m,mDeath);
		nDeaths=(int)(rbinom(1,m,mDeath)[0]);
                // the number of alives among m
                nAlives=m-nDeaths;
                // the new alive cells
                m+=nAlives-nDeaths;
                // current generation
                cg+=a;
            }
//             std::cout<<"cg final ="<<cg<<" | m final ="<<m<<" | Tps de mort (t-a) = "<<t-a<<std::endl;
//             sample[k] =m;
	      *it=m;
        }
    }
   // end if dirac distribution
  } else if(name.compare("exp") == 0) {
//       std::cout<<"Dist = exp"<<std::endl;
      double p,up;
      double geoD,beD;

      if (mDeath<DEATH_EPS_SIM) {
  //         for (int k=0;k<n;k++) {
	for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {
	      p=exp(-(*it));
	      
	      geoD=rgeom(1,p)[0];

	      if(geoD>=0) *it=geoD+1;
	      else if (geoD<0 || geoD!=geoD)*it=-1e-5;
	  }
      } else {
	for(NumericVector::iterator it = sample.begin();it != sample.end(); ++it) {

	    p=exp(-(*it));
	    
	    up=(1-2*mDeath)/(1-mDeath*(1+p));
	    
	    beD=rbinom(1,1,up)[0];
	    
	    if (beD==1){

	      geoD=rgeom(1,p*up)[0];
	      
	      if(geoD>=0) *it=geoD+1;
	      else if (geoD<0 || geoD!=geoD) *it=-1e-55;
	      
	    } else *it=0;
	}
      }
    // end if exponential distribution
  } else {

    // transform the sample

//     std::vector<double> splitTimesList;
    int nc=0;
    for(NumericVector::iterator it = sample.begin();
	  it != sample.end(); ++it) {
	     nc=splitTimes(*it);
// 	     ,splitTimesList);
        //sample[k]=splitTimesList.getSize()-1;
        *it=nc;
    }

  }

  return sample;
}


int FLAN_SimClone::splitTimes(double t) {
// ,std::vector<double>& splitTimesList) {

    std::string name=mDist->getDistName();
    std::vector<double> params=mDist->getDistParams();

    // return the clone size
    int cloneSize=0;

    //current generation
    int g=0;

    // make a sample of size 1
    //division times of current generation
    std::vector<double> st (1);
    if(name.compare("lnorm") == 0) st=as< std::vector<double> >(rlnorm(1,params[0],params[1]));
    
    if(name.compare("gamma") == 0) st=as< std::vector<double> >(rgamma(1,params[0],params[1]));

    // take only the times if before t
    // ng is the number of simulated value less than t
    // number of alive cells
    int ng=0, nActiveAlives;
    std::vector<double> gtf;
    
    // If no cells death
    if(mDeath<DEATH_EPS_SIM){
//       std::cout<<"Initialization..."<<std::endl;
      if(st[0]<t){
// 	splitTimesList.push_back(st[0]);
        ng=1;
      } else {
	// no division before t
	ng=0;
	cloneSize++;
      }
      nActiveAlives=ng;
      int i=0;
//       std::cout<<"Main loop"<<std::endl;
      while ((ng>0) && (ng<1.e+6)) {
	g++;
	i++;
	gtf.resize(2*ng);
	if(name.compare("lnorm") == 0) gtf=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
	if(name.compare("gamma") == 0) gtf=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));
	
	//  split times of daughters with duplicating split times
        for (int i=0;i<ng;i++) {
            gtf[2*i]+=st[i];
            gtf[2*i+1]+=st[i];
        }

        ng*=2;
        // keep only values less than t
        // daughters that still divide
        // keep their split times
	
	nActiveAlives=0;
	st.resize(ng);
        for (int i=0;i<ng;i++) {
            if (gtf[i]<t) {
                // division at  t of an alive cell
                st[nActiveAlives]=gtf[i];
                nActiveAlives++;
            } else {
                // no division before t
                cloneSize++;
            }
        }
//         splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());
        // ng is the size of st
        //dividing cells in next generation
        ng=nActiveAlives;
      }
    // If cells death with prob mDeath
    } else {
    
      NumericVector u=runif(1,0.,1.);
      if((st[0]<t) && (u[0]>mDeath)){
// 	splitTimes.add(st[0]);
        ng=1;
      } else if (st[0]>=t){
	// no division before t
	ng=0;
	cloneSize++;
      } else if(u[0] <= mDeath) {
	// the cell dies out
	ng=0;
      }
      nActiveAlives=ng;
      
      while ((ng>0) && (ng<1.e+6)) {
	g++;
// 	NumericVector gtf;
	gtf.resize(2*ng);
	if(name.compare("lnorm") == 0) gtf=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
	if(name.compare("gamma") == 0) gtf=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));
	
	//  split times of daughters with duplicating split times
        for (int i=0;i<ng;i++) {
            gtf[2*i]+=st[i];
            gtf[2*i+1]+=st[i];
        }

        ng*=2;
	
	u=runif(ng,0.,1.);
	
	nActiveAlives=0;
	st.resize(ng);
	for (int i=0;i<ng;i++) {

            if ((gtf[i]<t) && (u[i]>mDeath)) {
                // division at  t of an alive cell
                st[nActiveAlives]=gtf[i];
                nActiveAlives++;
            } else if (gtf[i]>=t) {
                // no division before t
                cloneSize++;
            }
        }
        // append st to spliTimes list
        //stack split times
//         splitTimesList.append(st);
// 	splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());

        // ng is the size of st
        //dividing cells in next generation
        ng=nActiveAlives;
      }
  }

    return cloneSize;
}


// List FLAN_SimClone::splitTimes_draw(double t){
//   
//     std::string name=mDist->getDistName();
//     std::vector<double> params=mDist->getDistParams();
//     std::vector<double> paramsN=mDistN->getDistParams();
// 
//     // return the clone size
//     int cloneSize;
// 
//     //current generation
//     int g=0;
// 
//     // make a sample of size 1
//     //division times of current generation
//     std::vector<double> st (1);
//     
//     if(name.compare("exp") == 0) st=as< std::vector<double> >(rexp(1,paramsN[0]));
//     if(name.compare("dirac") == 0) st=params[0];
//     if(name.compare("lnorm") == 0) st=as< std::vector<double> >(rlnorm(1,paramsN[0],paramsN[1]));
//     if(name.compare("gamma") == 0) st=as< std::vector<double> >(rgamma(1,paramsN[0],paramsN[1]));
// 
//     // take only the times if before t
//     // ng is the number of simulated value less than t
//     // number of alive cells
//     std::vector<int> ce (1);         // cells of current generation
//     ce[0]=1;
//     std::vector<int> mut (1);
//     mut[0]=0;         //mutants in current generation
//     int nd=0;         // number of dividing cells
//     std::vector<double> bd (1);      // birth date of current generation
//     bd[0]=1;
//     std::vector<double> de (1); // death date of current generation
//     de[0]=bd[0]+st[0];
//     std::vector<int> ord (1);        // location of cells in current generation
//     ord[0]=1;
//     std::vector<bool> div (1);
//     
//     std::vector<int> CE (1);         // All cells
//     std::vector<double> MUT (1);           // All mutants
//     std::vector<double> BD (1);            // All birth dates 
//     std::vector<double> DE (1);            // All death dates
//     std::vector<double> ORD (1);            // All locations
//     
// 
//     // If no cells death
//     if(mDeath<DEATH_EPS_SIM){
//           div[0]=(de[0]<t)?;
//     
// 	  de[0]=(de[0]<t)?de[0]:t; 
// //       std::cout<<"Initialization..."<<std::endl;
//       if(st[0]<t) nd=1;
//       else nd=0;  	// no division before t
//       
//       cloneSize=nd;
//       
//       CE[0]=ce[0] ; MUT[0]=mut[0] ; BD[0]=bd[0];
//       DE[0]=de[0] ; ORD[0]=ord[0];
//       
//       
//       
//       int nActiveAlives=ng;
//       int i=0,nce;
//       std::vector<double> stn;
//       std::vector<double> stm;
// //       std::cout<<"Main loop"<<std::endl;
//       while ((nd>0) && (nx<1.e+4)) {
// 	g++;
// 	i++;
// 	
// 	ce.insert(ce.end(), ce.begin(), ce.end());
// 	
// 	for(int i=0;i<ce.size();i++){
// 	  if(div[i]) ce[i]*=2;
// 	  if(i>=ce.size()/2) ce[i]++;
// 	}
// 	
// 	
// 	
// 	stn.resize(2*ng);
// 	stm.resize(2*ng);
// 	if(name.compare("lnorm") == 0) stn=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
// 	if(name.compare("gamma") == 0) stm=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));
// 	
// 	//  split times of daughters with duplicating split times
//         for (int i=0;i<ng;i++) {
//             gtf[2*i]+=st[i];
//             gtf[2*i+1]+=st[i];
//         }
// 
//         ng*=2;
//         // keep only values less than t
//         // daughters that still divide
//         // keep their split times
// 	
// 	nActiveAlives=0;
// 	st.resize(ng);
//         for (int i=0;i<ng;i++) {
//             if (gtf[i]<t) {
//                 // division at  t of an alive cell
//                 st[nActiveAlives]=gtf[i];
//                 nActiveAlives++;
//             } else {
//                 // no division before t
//                 cloneSize++;
//             }
//         }
// //         splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());
//         // ng is the size of st
//         //dividing cells in next generation
//         ng=nActiveAlives;
//       }
//     // If cells death with prob mDeath
//     } else {
//     
//       NumericVector u=runif(1,0.,1.);
//       if((st[0]<t) && (u[0]>mDeath)){
// // 	splitTimes.add(st[0]);
//         ng=1;
//       } else if (st[0]>=t){
// 	// no division before t
// 	ng=0;
// 	cloneSize++;
//       } else if(u[0] <= mDeath) {
// 	// the cell dies out
// 	ng=0;
//       }
//       nActiveAlives=ng;
//       
//       while ((ng>0) && (ng<1.e+6)) {
// 	g++;
// // 	NumericVector gtf;
// 	gtf.resize(2*ng);
// 	if(name.compare("lnorm") == 0) gtf=as< std::vector<double> >(rlnorm(2*ng,params[0],params[1]));
// 	if(name.compare("gamma") == 0) gtf=as< std::vector<double> >(rgamma(2*ng,params[0],params[1]));
// 	
// 	//  split times of daughters with duplicating split times
//         for (int i=0;i<ng;i++) {
//             gtf[2*i]+=st[i];
//             gtf[2*i+1]+=st[i];
//         }
// 
//         ng*=2;
// 	
// 	u=runif(ng,0.,1.);
// 	
// 	nActiveAlives=0;
// 	st.resize(ng);
// 	for (int i=0;i<ng;i++) {
// 
//             if ((gtf[i]<t) && (u[i]>mDeath)) {
//                 // division at  t of an alive cell
//                 st[nActiveAlives]=gtf[i];
//                 nActiveAlives++;
//             } else if (gtf[i]>=t) {
//                 // no division before t
//                 cloneSize++;
//             }
//         }
//         // append st to spliTimes list
//         //stack split times
// //         splitTimesList.append(st);
// // 	splitTimesList.insert(splitTimesList.end(), st.begin(), st.end());
// 
//         // ng is the size of st
//         //dividing cells in next generation
//         ng=nActiveAlives;
//       }
//   }
// 
//     return cloneSize; 
//   
// }



  /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class for the distribution of a clone size
  */
const double FLAN_Clone::DEATH_EPS_DIST=1.e-4;

 /*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is supposed to be exponential
  */


// create object for GF method


NumericVector FLAN_ExponentialClone::computeProbability(int m){

  std::vector<double> P(m+1);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;
//     if(m == 0) return NumericVector(P[0]);
//     if(m == 0) return P[0];

    if (m > 0){
      int k=1;
      for(std::vector<double>::iterator it=P.begin()+1 ; it!= P.end() ; ++it,k++) *it=mFitness*R::beta(mFitness+1,k);
    }
  } else {
    double d1=mDeath/(1-mDeath);
    int m_max=1000;
//     std::cout<<"function type set"<<std::endl;
//     std::cout<<mIntegrator<<std::endl;
//     mIntegrator->testintegralFunction(0,1,1,0.05);
//     std::cout<<"done1"<<std::endl;
    mIntegrator->setFunctionName("CLONE_P0_WD");
//     std::cout<<"done"<<std::endl;
    //integrate the function in [0,1]
    double I;
    // m=0 probability
    I=mIntegrator->integralFunction(0.,1.,mFitness,d1,0.);
//     std::cout<<"integral computed"<<std::endl;
    P[0]=I*d1*mFitness;

//     sum_pk+=p_k;
//     if (m==0) return NumericVector(P[0]);
    
    if(m > 0){

      double d2=(1.-2.*mDeath)/(1-mDeath);
      d2*=d2;

      //m>0 probability
      int m1=m;
      if (m1>=m_max) m1=m_max;
  //     int k=1;
      mIntegrator->setFunctionName("CLONE_PK_WD");
      std::vector<double>::iterator it=P.begin()+1;
      for (int k=1;k<=m1;k++,++it) {
// 	    I=mIntegrator->integralFunction(0.,k,mFitness,d1,k);
// 	    *it=I*d2*mFitness*pow(k,-mFitness-1.);
	    I=mIntegrator->integralFunction(0.,1,mFitness,d1,k);
	    *it=I*d2*mFitness;
      }

      // equivalent computation
      double a=pow(d2,(1.-mFitness)/2.)*mFitness*R::gammafn(mFitness+1);
  //     k=m1+1;
      for (int k=m1+1;k<=m;k++,++it) {
  //     for(std::vector<double>::iterator it=P.begin()+m1+1 ; it!= P.end() ; ++it,k++){
	*it=a*pow(k,-mFitness-1);
  //         sum_pk+=p_k;
      }
    }
 }

  return NumericVector(P.begin(),P.end());

}

List FLAN_ExponentialClone::computeProbability1DerivativeRho(int m){

  std::vector<double> P(m+1);
  std::vector<double> dP_dr(m+1);
//   double pk;
  
  if(mDeath<DEATH_EPS_DIST){
    
    P[0]=0;
    dP_dr[0]=0;
    if(m == 0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);

//     int k=1;
    int k=1;
    std::vector<double>::iterator itdP=dP_dr.begin()+1 ;
    double dg=R::digamma(mFitness+1);
//     for(int k=1; k<=m ; k++,++itP,++itdP){
    for(std::vector<double>::iterator itP=P.begin()+1 ; itP!= P.end(); ++itP,++itdP,k++){
      *itP=mFitness*R::beta(mFitness+1,k);
      *itdP=(*itP)*(1/mFitness+dg-R::digamma(mFitness+k+1));
//       *it=pk*(1/mFitness+dg-R::digamma(mFitness+k+1));
    }
  } else {
    double d1=mDeath/(1-mDeath);
    
    double I;
    
    int m_max=1000;
    // set the function type
    mIntegrator->setFunctionName("CLONE_P0_WD");
    I=mIntegrator->integralFunction(0.,1.,mFitness,d1,0.);
    P[0]=I*d1*mFitness;
    
    mIntegrator->setFunctionName("CLONE_dP0_dr_WD");
    //integrate the function in [0,1]
    // m=0 probability
    I=mIntegrator->integralFunction(0.,1.,mFitness,d1,0.);
    dP_dr[0]=I*d1*mFitness+P[0]/mFitness;
//     dP_dr[0]=I*d1*mFitness+pk/mFitness;

//     sum_pk+=p_k;
    if (m==0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);

    
    double d2=(1.-2.*mDeath)/(1-mDeath);
    d2*=d2;
    //m>0 probability
    int m1=m;
    if (m1>=m_max) m1=m_max;
//     int k=1;
    std::vector<double>::iterator itP=P.begin()+1 ;
    std::vector<double>::iterator itdP=dP_dr.begin()+1 ;
    for (int k=1;k<=m1;k++,++itP,++itdP) {
	mIntegrator->setFunctionName("CLONE_PK_WD");
// 	I=mIntegrator->integralFunction(0.,k,mFitness,d1,k);
// 	*itP=I*d2*mFitness*pow(k,-mFitness-1.);
	I=mIntegrator->integralFunction(0.,1,mFitness,d1,k);
	*itP=I*d2*mFitness;
	
	mIntegrator->setFunctionName("CLONE_dPK_dr_WD");
//         I=mIntegrator->integralFunction(0.,k,mFitness,d1,k);
//         *itdP=(*itP)/mFitness+I*d2*mFitness*pow(k,-mFitness-1.);
	I=mIntegrator->integralFunction(0.,1,mFitness,d1,k);
        *itdP=(*itP)/mFitness+I*d2*mFitness;
    }

    // equivalent computation
    double gfn=R::gammafn(mFitness+1);
    double rhodg=mFitness*R::digamma(mFitness+1);
    double a=pow(d2,(1-mFitness)/2.);
    double b=a*mFitness*R::gammafn(mFitness+1);
    double c=-0.5*log(d2)*mFitness+1.;
    double kpr;
//     for(std::vector<double>::iterator it=dP_dr.begin()+m1+1 ; it!= dP_dr.end() ; ++it, k++){
    for(int k=m1+1; k<=m ;k++ , ++itP,++itdP){
	
      kpr=pow(k,-mFitness-1);
      *itP=a*kpr;
      
      *itdP=b*kpr*(
            gfn*(c-mFitness*log(k))
            +rhodg);
    }
  }

  return List::create(_["P"]=NumericVector(P.begin(),P.end()),
		      _["dP_dr"]=NumericVector(dP_dr.begin(),dP_dr.end()));

}


std::vector<double> FLAN_ExponentialClone::computeGeneratingFunction2(double rho,std::vector<double> Z){
  
  double eps =1e-8;
  
  std::vector<double> H(Z.size());
  std::vector<double>::iterator itH=H.begin();
  
  for (std::vector<double>::iterator itZ = Z.begin() ; itZ != Z.end() ; ++itZ, ++itH){
    if(fabs(*itZ) <= eps) *itH=0.;
    else if(fabs(1-(*itZ)) <= eps) *itH=1.;
    else {
      double dstar,zstar,I,h;
      
      dstar=mDeath/(1-mDeath);
      zstar=((*itZ)-dstar)/(1-(*itZ));
      
      mIntegrator->setFunctionName("CLONE_PGF");
      I=mIntegrator->integralFunction(0.,1.,rho,zstar,0);
      
      *itH=dstar+I*zstar*(1-dstar)*rho;
      

    }
  }
  
  return H;
    
}


double FLAN_ExponentialClone::computeGeneratingFunction1DerivativeRho(double z){
  
  double eps =1e-8;
  
  if(fabs(z) <= eps) return 0.;
  else if(fabs(1-z) <= eps) return 0.;
  else {
    double dstar,zstar,I1,I2,h;
    
    dstar=mDeath/(1-mDeath);
    zstar=(z-dstar)/(1-z);
    
    mIntegrator->setFunctionName("CLONE_PGF");
    I1=mIntegrator->integralFunction(0.,1.,mFitness,zstar,0);
    
    mIntegrator->setFunctionName("CLONE_PGF_dr");
    I2=mIntegrator->integralFunction(0.,1.,mFitness,zstar,0);
    
    h=zstar*(1-dstar)*(I1+I2*mFitness);
    
    return h;
  }
  
}





/*//////////////////////////////////////////////////////////////////////////////////////////////
  * FLAN_Clone class when the lifetime model is suppsoed to be constant
  */

NumericVector FLAN_DiracClone::computeProbability(int m){

  std::vector<double> P(m+1);

  if(mDeath<DEATH_EPS_DIST){

//     std::cout<<"Size of P="<<P.size()<<std::endl;
    P[0]=0;

//     if(m == 0) return NumericVector(P[0]);
//     if(m == 0) return P[0];
    if(m > 0){
//     for(NumericVector::iterator it=P.begin()+1;it!=P.end(); ++it) *it=0;
      int n=floor(log(m)/log(2));
        //take only the indices index=2^ind for ind in [0,n] where index < m
      for (int k=0;k<=n;k++) P[pow(2,k)]=(1-pow(2,-mFitness))*pow(2,-k*mFitness);
    }

  } else {   
//     for(NumericVector::iterator it=P.begin()+1;it!=P.end(); ++it) *it=0;
        
    // case with death >0
    double eps=1.e-6;
    double deps=1.e-8;
    
    int itmax=25;
  
    // initiate the polynome X
    mPol.setDegree(1);
//     mPol[0]=0;
//     mPol[1]=1;
//     std::cout<<"mPol[0]"<<mPol[0]<<std::endl;
//     std::cout<<"mPol[1]"<<mPol[1]<<std::endl;
    

    // minimum number of iterations
    int i=0;
    double umd=1-mDeath;
    double a=log(2*umd);
    double t=exp(-mFitness*a);
    double ti=1;//t^i
    int Pol_degree=1;
    
    double sumP;
    double sumP_old;
    

    // initialize P
    P[0]=0;P[1]=1;
    sumP=1;
    double err = 1;
    int dmax;
    int k;
    // loop
//         while (((Pdegree<=m) || (fabs(sumPk-sumPk_old)>deps)) && (i<25)) {
    std::vector<double>::iterator itP;
    while (((Pol_degree<=m) || (err>eps)) && (i<itmax)) {
	// next iteration
// 	  std::cout<<"it n° "<<i<<std::endl;
	i++;
	ti*=t;
	sumP_old=sumP;
	
	// update P
	//P:=death+(1-death)*P^2
// 	std::cout<<"Compute square of mPol"<<std::endl;
	mPol.square_fft();
// 	    P*=P;

	mPol*=umd;
	mPol+=mDeath;
	
	mPol.reduce(deps);
	Pol_degree=mPol.getDegree();
	// update pk=sum_{i=1}^{i=number of polynoms computer} Pi[k] t^i
	dmax=(m<Pol_degree)?m:Pol_degree;
	sumP=0;
// 	std::cout<<"dmax ="<<dmax<<std::endl;
	itP=P.begin();
	for (k=0;k<=dmax;k++,++itP) {
	    *itP+=mPol[k]*ti;
	    sumP+=*itP;
	}
//             std::cout<<"pk[0]="<<pk[0]<<std::endl;
      
// 	std::cout << "sum P["<<i<<"]="<<sumP<<" Pol:"<<Pol_degree<<" err="<<(fabs(1-sumP/sumP_old))<<std::endl;
	// next iteration
	err = fabs(1-sumP_old/sumP);
// 	std::cout<<"it n°"<<i<<"err =)"<<err<<std::endl;
// 	err = fabs(1-sumP/sumP_old);
    } //end loop

    // multiply by (1-t)
    for(itP=P.begin() ; itP!=P.end();++itP) (*itP)*=(1-t);
    
  }

  return NumericVector(P.begin(),P.end());

}


List FLAN_DiracClone::computeProbability1DerivativeRho(int m){

  std::vector<double> P(m+1);
  std::vector<double> dP_dr(m+1);
  
//   std::vector<double> P=computeProbability(m);

  if(mDeath<DEATH_EPS_DIST){

    P[0]=0;
    dP_dr[0]=0;
    if(m == 0) return List::create(_["P"]=P[0],_["dP_dr"]=dP_dr[0]);
    
//     NumericVector::iterator itP=P.begin()+1;
//     for(NumericVector::iterator itdP=dP_dr.begin()+1 ; itdP!=dP_dr.end(); ++itdP, ++itP) {
//       *itP=0;
//       *itdP=0;
//     }
//     
    int n=floor(log(m)/log(2));
    double index;

    for(int k=0;k<=n;k++) {
      index=pow(2,k);
      P[index]=(1-pow(2,-mFitness))*pow(2,-k*mFitness);
      dP_dr[index]=log(2)*(pow(2,-mFitness*(k+1))-k*P[index]);
    }

  } else {
     // case with death >0
    double eps=1.e-6;
    double deps=1.e-8;
    
    int itmax=25;
  
    // initiate the polynome X
    mPol.setDegree(1);
//     mPol[0]=0;
//     mPol[1]=1;
    

    // minimum number of iterations
    int i=0;
    double umd=1-mDeath;
    double a=log(2*umd);
    double t=exp(-mFitness*a);
    double umt=1-t;
    double ti=1;//t^i
    double iti=0; // i*t^i
    int Pol_degree=1;
    
    double sumP;
    double sumP_old;
    

    // initialize P
    P[0]=0;P[1]=1;
    dP_dr[0]=0;dP_dr[1]=t-umt;
    sumP=1;
    double err = 1;
    int dmax;
    int k;
    // loop
    std::vector<double>::iterator itP;
    std::vector<double>::iterator itdP;
    while (((Pol_degree<=m) || (err>eps)) && (i<itmax)) {
	// next iteration
	i++;
	ti*=t;
	iti=i*ti;
	sumP_old=sumP;
		  
	// update P
	mPol.square_fft();
	  
	// update P
	mPol*=(1-mDeath);
	mPol+=mDeath;
	mPol.reduce(deps);
      
	Pol_degree=mPol.getDegree();

	dmax=(m<Pol_degree)?m:Pol_degree;
	sumP=0;
	
	itP=P.begin();
	itdP=dP_dr.begin();
	for (k=0;k<=dmax;k++,++itP,++itdP) {
	    
	  *itP+=mPol[k]*ti;
	  sumP+=*itP;
	  
	  *itdP-=mPol[k]*iti;
	}
	  
      err = fabs(1-sumP_old/sumP);
      
    }
    
    itdP=dP_dr.begin();
    for(itP=P.begin() ; itP!=P.end() ; ++itP,++itdP){
      (*itdP)*=umt;
      (*itdP)+=t*(*itP);
      (*itP)*=umt;
    }
      
  }
  
  return List::create(_["P"]=NumericVector(P.begin(),P.end()),
		      _["dP_dr"]=NumericVector(dP_dr.begin(),dP_dr.end()));

}


std::vector<double> FLAN_DiracClone::computeGeneratingFunction2(double rho,std::vector<double> Z){

  
    std::vector<double> H (Z.size());
    std::vector<double>::iterator itH=H.begin();
    double eps=1.e-8;

    for(std::vector<double>::iterator itZ=Z.begin() ; itZ != Z.end() ; ++itZ, ++itH){
      // z=0 return 0
      if (fabs((*itZ))<eps) *itH=0.;

      // z=1 return 1
      if (fabs(1-(*itZ))<eps) *itH=1.;

      // otherwize
      double s=0;
      if (mDeath<DEATH_EPS_DIST) {
	  double a=pow(2,-rho);
	  int n=floor(4-log(fabs(log((*itZ))))/log(2))+1;
  // 	n++;

	  for (int k=0;k<=n;k++) {
	      s+=pow((*itZ),pow(2,k))*pow(a,k);
	  }
	  s*=(1-a);
      } else {
	  double a=log(2*(1-mDeath));
	  double dstar=mDeath/(1-mDeath);
	  int n=floor(-log(eps)/(rho*a))+1;
	  double tp=exp(-rho*a);
	  double tpi=1;
	  double bi=(*itZ);
	  s=(*itZ);
	  for (int i=1;i<=n;i++) {
	      bi=mDeath+(1-mDeath)*bi*bi;
	      tpi*=tp;
	      s+=tpi*bi;
	  }
	  s*=(1-tp);
	  s+=dstar*tpi*tp;
      }
      *itH=s;
    }
    
  return H;

}
//
//
double FLAN_DiracClone::computeGeneratingFunction1DerivativeRho(double z) {

    double eps=1.e-8;
    // z=0 return 0
    if (fabs(z)<eps) return 0;
    // z=1 return 1
    if (fabs(1-z)<eps) return 0;
    // otherwise

    if (mDeath<DEATH_EPS_DIST) {
        double a=pow(2,-mFitness);
        int n=floor(4-log(fabs(log(z)))/log(2))+1;
// 	n++;
        double s1=0,s2=0,t=0;
        for (int k=0;k<=n;k++) {
            t=pow(z,pow(2,k))*pow(a,k);
            s1+=t;
            s2+=k*t;

        }
        return log(2)*(a*s1-(1-a)*s2);
    } else {
        double a=log(2*(1-mDeath));
        double dstar=mDeath/(1-mDeath);
        int n=floor(-log(eps)/(mFitness*a))+1;
// 	n++;
        double tp=exp(-mFitness*a);
        double tpi=1;
        double bi=z;
        double s=z;
        double ds=0;
        for (int i=1;i<=n;i++) {
            bi=mDeath+(1-mDeath)*bi*bi;
            tpi*=tp;
            s+=tpi*bi;
            ds-=i*tpi*bi;
        }
        s=a*(ds*(1-tp)+s*tp);
        // no usefull when n is big
//         s-=dstar*a*tpi*tp*(n+1);
	s+=dstar*tpi*tp;
        return s;
    }
}
//
//
//
