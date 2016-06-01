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
#include "FLAN_MutationModel.h"
    // --------------------------
    // Probability methods
    //---------------------------
    

NumericVector FLAN_MutationModel::computeProbability(int m) {
    
    //vector of clone probabilities
   
    NumericVector P=mClone->computeProbability(m);
    
    return deduceProbability(m,P);
}


NumericVector FLAN_MutationModel::deduceProbability(int m,NumericVector& P) {
    
    // initialize Q
    std::vector<double> Q(m+1);
    
    
    //initial probability q0
    Q[0]=exp(-mMutNumber*(1.-P[0]));
    
    if(m > 0){
      double s=0;
      for (int k=1;k<=m;k++) {
	  // convolution
	  s=0;
  // 	i=1;
	  for (int i=1;i<=k;i++) {
	      s+=i*P[i]*Q[k-i];
	  }
	  Q[k]=(mMutNumber/k)*s;
      }
    }
      
    return NumericVector(Q.begin(),Q.end());
}


List FLAN_MutationModel::computeProbability1DerivativeAlpha(int m
// 							    NumericVector& Q,
// 							    NumericVector& dQ_da
									  ) {
    
                                                 
    
    // compute the probabilities of mClone
    NumericVector P=mClone->computeProbability(m);
    
    return deduceProbability1DerivativeAlpha(m,P);
}


List FLAN_MutationModel::deduceProbability1DerivativeAlpha(int m,
							    NumericVector& P
// 							    NumericVector& Q,
// 							    NumericVector& dQ_da
							  ) {
    
                                                    
    //initial probability Q0
    std::vector<double> Q(m+1);
    std::vector<double> dQ_da(m+1);

    
    
    //initial probability q0
    Q[0]=exp(-mMutNumber*(1.-P[0]));
        // first derivatives of Q with respect to mMutNumber
    dQ_da[0]=-(1-P[0])*Q[0];
   
    if (m == 0) return List::create(_["Q"]=Q[0],
				    _["dQ_da"]=dQ_da[0]
				   );
    
    double s,ds_da; 
    for (int k=1;k<=m;k++) {
        s=0;ds_da=0;
        for (int i=1;i<=k;i++) {

            s+=i*P[i]*Q[k-i];
	    ds_da+=P[i]*Q[k-i];

        }
	Q[k]=(mMutNumber/k)*s;
        dQ_da[k]=ds_da-Q[k];
    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_da"]=NumericVector(dQ_da.begin(),dQ_da.end())
		       );
}


List FLAN_MutationModel::computeProbability1DerivativeRho(int m) {
    
         
    
    // compute the probabilities and their derivative with respect to mFitness of mClone
    List P_dP_dr=mClone->computeProbability1DerivativeRho(m);
    
    NumericVector P=P_dP_dr["P"];
    NumericVector dP_dr=P_dP_dr["dP_dr"];
    
    // compute the probabilities Q

    return deduceProbability1DerivativeRho(m,P,dP_dr);
}




List FLAN_MutationModel::deduceProbability1DerivativeRho(int m,
							  NumericVector& P,
							  NumericVector& dP_dr) {
    
                                                    
    //initial probability Q0

    std::vector<double> Q(m+1);
    std::vector<double> dQ_dr(m+1);

    
    // first derivatives of Q with respect to mFitness
    Q[0]=exp(-mMutNumber*(1.-P[0]));
    dQ_dr[0]=mMutNumber*(dP_dr[0])*Q[0];
    
//     if (m==0) return true;
    if (m==0)  return List::create(_["Q"]=Q,
				   _["dQ_dr"]=dQ_dr
				  );
    
    double s=0,ds_dr=0;

    for (int k=1;k<=m;k++) {
        s=0;ds_dr=0;
        for (int i=1;i<=k;i++) {
            s +=i*P[i]*Q[k-i];
	    ds_dr+=dP_dr[i]*Q[k-i];

        }
        Q[k] =(mMutNumber/k)*s;
	dQ_dr[k]=mMutNumber*ds_dr;


    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_dr"]=NumericVector(dQ_dr.begin(),dQ_dr.end())
			);
}



List FLAN_MutationModel::computeProbability1DerivativesAlphaRho(int m
// 								NumericVector& Q,
// 								NumericVector& dQ_da,
// 								NumericVector& dQ_dr
							       ) {
    
                                                 
    // compute the derivative with respect to mFitness of GY
    
    List P_dP=mClone->computeProbability1DerivativeRho(m);
    
    NumericVector P=P_dP["P"];
    NumericVector dP_dr=P_dP["dP_dr"];
    
    
    return deduceProbability1DerivativesAlphaRho(m,P,dP_dr);
}


List FLAN_MutationModel::deduceProbability1DerivativesAlphaRho(int m,
								NumericVector& P,
								NumericVector& dP_dr
// 								NumericVector& Q,
// 								NumericVector& dQ_da,
// 								NumericVector& dQ_dr
							      ) {
    
                                                    
    //initial probability Q0
    std::vector<double> Q(m+1);
    std::vector<double> dQ_da(m+1);
    std::vector<double> dQ_dr(m+1);

    
    
    Q[0]=exp(-mMutNumber*(1-P[0]));
    // first derivatives of Q with respect to mMutNumber
    dQ_da[0]=-(1-P[0])*Q[0];
    
    // first derivatives of Q with respect to mFitness
    dQ_dr[0]=mMutNumber*(dP_dr[0])*Q[0];

//     if (m==0) return true;
    if (m==0) return List::create(_["Q"]=Q,
				  _["dQ_da"]=dQ_da,
				  _["dQ_dr"]=dQ_dr
				 );
    
    double s,ds_da,ds_dr;

    for (int k=1;k<=m;k++) {
        s=0;ds_da=0;ds_dr=0;
        for (int i=1;i<=k;i++) {
            s +=i*P[i]*Q[k-i];
	    ds_da+=P[i]*Q[k-i];
	    ds_dr+=dP_dr[i]*Q[k-i];

        }
        Q[k]=(mMutNumber/k)*s;
	dQ_da[k]=ds_da-Q[k];
	dQ_dr[k]=mMutNumber*ds_dr;

    }
    return List::create(_["Q"]=NumericVector(Q.begin(),Q.end()),
			_["dQ_da"]=NumericVector(dQ_da.begin(),dQ_da.end()),
			_["dQ_dr"]=NumericVector(dQ_dr.begin(),dQ_dr.end())
	    );
}


  

NumericVector FLAN_MutationModel::computeCumulativeFunction(int m) {
		
    std::vector<double> cumsum(m+1);
        
    NumericVector Q=computeProbability(m);
    
    std::partial_sum(Q.begin(),Q.end(),cumsum.begin(),std::plus<double>());
    
    
    if(!mLT) {
      for(std::vector<double>::iterator it=cumsum.begin();it!=cumsum.end();++it) {
	(*it)*=-1;
	(*it)+=1;
      }
    }
    
    return NumericVector(cumsum.begin(),cumsum.end());
}


// GF method


List FLAN_MutationModel::MutationGFEstimation(){
  
  int n=mSample.size();
  
//   double z3=mTuning["z3"];
  double z3=0.8;
  
  z3=pow(z3,1./mScale);
  
  double g=0;
  
  double alpha,hr,gg,sd_alpha;
  
  for(NumericVector::iterator it=mSample.begin() ; it!=mSample.end() ; ++it) g+=pow(z3,*it);
  
  g/=n;
  
  hr=mClone->computeGeneratingFunction(z3);
  
  alpha=log(g)/(hr-1);

  setMutNumber(alpha);
  
  gg=exp(mMutNumber*(hr-1));
  
  sd_alpha=covariance2(z3,z3)/pow(gg*log(g)/mMutNumber,2);
  sd_alpha/=n;
  sd_alpha=sqrt(sd_alpha);
  
  double tp,pmut,sd_pmut;
  if(mMfn > 0){
    pmut=alpha/mMfn;
    sd_pmut=sd_alpha/mMfn;
    
    if(mCvfn > 0){
      tp=alpha*(1-hr)*pow(mCvfn,2.);
      pmut*=(1+tp/2);
      sd_pmut*=(1+tp);
    }
    
    return(List::create(_["mutprob"]=pmut,_["sd.mutprob"]=sd_pmut));
  } else return(List::create(_["mutations"]=alpha,_["sd.mutations"]=sd_alpha));
}


double FLAN_MutationModel::computeGeneratingFunction(double z) {
    return exp(mMutNumber*(mClone->computeGeneratingFunction(z)-1));
}
// 
// 
double FLAN_MutationModel::covariance2(double z1, double z2) {
    return computeGeneratingFunction(z1*z2)-
        computeGeneratingFunction(z1)*
        computeGeneratingFunction(z2);
}
// 
NumericMatrix FLAN_MutationModel::covariance(double z1, double z2,double z3) {

    
    NumericMatrix M(3,3);
    M(0,0)=covariance2(z1,z1);
    M(0,1)=covariance2(z1,z2);
    M(0,2)=covariance2(z1,z3);

  
    M(1,1)=covariance2(z2,z2);
    M(1,2)=covariance2(z2,z3);

    M(2,2)=covariance2(z3,z3);
    for (int i=0;i<3;i++) {
        for (int j=0;j<i;j++) {
            M(i,j)=M(j,i);
        }
    }

    double ccf_z1=mClone->computeGeneratingFunction(z1);
    double ccf_z2=mClone->computeGeneratingFunction(z2);
    double ccf_z3=mClone->computeGeneratingFunction(z3);

    double dccf_z1=mClone->computeGeneratingFunction1DerivativeRho(z1);
    double dccf_z2=mClone->computeGeneratingFunction1DerivativeRho(z2);
    double dccf_z3=mClone->computeGeneratingFunction1DerivativeRho(z3);
    

    double d01=(ccf_z2-1)*dccf_z1-(ccf_z1-1)*dccf_z2;
    double d11=(ccf_z1-1)*dccf_z2-(ccf_z2-1)*dccf_z1;
    
    
    NumericMatrix CO(3,2);
    CO(0,1)= (ccf_z2-1)/(mMutNumber*d01*computeGeneratingFunction(z1));
    
    CO(1,1)=(ccf_z1-1)/(mMutNumber*d11*computeGeneratingFunction(z2));
    CO(2,1)= 0;
    
    CO(0,0)=(mMutNumber*dccf_z3*CO(0,1))/(1-ccf_z3);
    CO(1,0)=(mMutNumber*dccf_z3*CO(1,1))/(1-ccf_z3);
    CO(2,0)=1./(computeGeneratingFunction(z3)*(ccf_z3-1));
    
    NumericMatrix cov(2,2);
    double d, c;
    for (int i=0;i<2;i++) {
        for (int j=0;j<2;j++) {
            d=0;
            for (int k=0;k<3;k++) {
                c=0;
                for (int s=0;s<3;s++) {
                    c+=M(k,s)*CO(s,j);
                }
                d+=CO(k,i)*c;
            }
            cov(i,j)=d;
        }
    }
    
    return cov;
    
}
// 
// 
//  // Unbiased estimation of pi and its standart deviation if fluctuation of final counts
List FLAN_MutationModel::unbiasPiEstimation(double sd,double z,
					    double mfn,double cvfn)  {
    
    double pm=mMutNumber/mfn;
    double sd_pm=sd/mfn;
    if(cvfn > 0){
      double f=mMutNumber*(1-mClone->computeGeneratingFunction(z))*cvfn*cvfn;
      pm*=1+f/2;
      sd_pm*=1+f;
    }
    
    return List::create(_["mutprob"] = pm, _["sd.mutprob"]=sd_pm);
}
