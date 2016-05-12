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

#include "FLAN_Sim.h"

using namespace Rcpp;

    // --------------------------------------------------------------------------------------------------------------
    //////////////////////////////////////////// Samples computation ////////////////////////////////////////////////
    // --------------------------------------------------------------------------------------------------------------

    /* Generates a sample of size n of couples (mutant count ; final count)
     * The final counts are sampled with the log-normal
     * distribution with coefficient of variation cvfn and mean mfn.
     *
     * n: number of experiments
     * mfn: mean of final counts
     * Cfn: variation number of final counts
     * alphapi: mean number of mutations probability of mutation (depends of cvfn)
     * rho: fitness parameter
     * death: probability of death
     */


List FLAN_Sim::computeSamplesMutantsFinalsNumber(int n)  {

   
  RNGScope rngScope;
  
  NumericVector mutantCount(n);

  if (mCvfn>0) {
      double sdLog2=log(1+mCvfn*mCvfn);
      double sdLog=sqrt(sdLog2);
      double meanLog=log(mMfn)-sdLog2/2;
      

//       for(NumericVector::iterator it = finalCount.begin();
// 	  it != finalCount.end(); ++it) {
// 	    *it=R::rlnorm(meanLog,sdLog);
//       }
      
      NumericVector finalCount=rlnorm(n,meanLog,sdLog);

      mutantCount=computeSampleMutantsNumber(n,finalCount);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=finalCount);


  } else {
      mutantCount=computeSampleMutantsNumber(n);
      return List::create(_["mc"]=mutantCount,
			  _["fn"]=R_NilValue);
  }
  //export in R
//   return List::create(_["mc"]=mutantCount,_["fn"]=finalCount);
  
//   return List::create(_["mc"]=NumericVector(mutantCount.begin(),mutantCount.end()),
// 		      _["fn"]=NumericVector(finalCount.begin(),finalCount.end())
//   );


}



    /* Generates a sample of size n of mutant counts
     * n: sample size
     * alpha: mean number of mutations
     * rho: fitness parameter
     * death: probability of death
     */



NumericVector FLAN_Sim::computeSampleMutantsNumber(int n)  {
  
    double s,sj;
    // set the size of the mutants
    NumericVector mutantCount=rpois(n,mMut);
    int mc;
    int i=0,j;
    bool testneg;
    for (NumericVector::iterator it = mutantCount.begin(); it != mutantCount.end(); ++it,i++) {
//       int i=0;i<n;i++) {
        // simulate a poisson number

// 	  mc=(int)(mcs[i]);
	  mc=(int)(*it);

        if (mc>0) {
            //Poisson sum of Yule variables
            NumericVector sample=mClone->computeSample(mc);

            //s = sum_j sample[j]
//             int m=sample.getSize();
            s=0;
	    j=0;
	    testneg=false;
//             for (int j=0;j<m;j++){
	    while(j<mc && !testneg){
	      sj=sample[j];
	      if(sj < 0){
// 		std::cout<<"Negaative number !"<<std::endl;
		testneg=true;
		s=sj;
// 		s=-1e-5;
	      } else {
		s+=sj;
		j++;
	      }
	    }
	    *it=s;
        } else
            *it=0;
    }
    
    return mutantCount;
}


    /* Generates a sample of size n of mutant counts, given a sample of final counts.
     *
     * n: number of experiments
     * pi: mutation probability
     * rho: fitness parameter
     * death: probability of death
     * finalCount:
     */

NumericVector FLAN_Sim::computeSampleMutantsNumber(int n,
					  NumericVector& finalCount)  {

    std::vector <double> mutantCount(n);
//     FLAN_SimClone* clone=new FLAN_SimClone(mFitness,mDeath,mDist);

    // poisson number of mutation

//     NumericVector sample;
    double s,sj;
    // set the size of the mutants
//     mutantCount.resize(n);

//     double mcs;
    double lambda;
    int mc,j;
    bool testneg;
//     int i=0;
    NumericVector::iterator itfn=finalCount.begin();
    for (std::vector <double>::iterator itmc = mutantCount.begin();
	  itmc != mutantCount.end(); ++itmc, ++itfn) {
        // simulate a poisson number
      lambda=mMut*(*itfn);

      mc=(int)(rpois(1,lambda)[0]);
      
      if (mc>0) {
            //Poisson sum of Yule variables
            NumericVector sample=mClone->computeSample(mc);
            //s = sum_j sample[j]
//             int m=sample.getSize();
            s=0;
	    j=0;
	    sj;
	    testneg=false;
//             for (int j=0;j<m;j++){
	    while(j<mc && !testneg){
	      sj=sample[j];
	      if(sj < 0){
		testneg=true;
		s=sj;
// 		s=-1e-5;
	      } else {
		s+=sj;
		j++;
	      }
	    }
	    *itmc=s;
        } else
            *itmc=0;
    }
    
    return NumericVector(mutantCount.begin(),mutantCount.end());

}
