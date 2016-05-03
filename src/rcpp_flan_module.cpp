#include "rcpp_flan_module.h"

// static void finalizer_of_flan( Flan* ptr ){
//       printf("finalizer sim_vam has been called\n");
// }

RCPP_MODULE(flan_module) {
  class_<FLAN_Sim>("FlanSim")
	.constructor<List>()
	.method("rflan",&FLAN_Sim::computeSamplesMutantsFinalsNumber,"compute sample mutants")
  ;

  class_<FLAN_MutationModel>("FlanMutMod")
	.constructor<List>()
// 	.constructor<double,double,std::string,double,>()
	.method("pflan",&FLAN_MutationModel::computeCumulativeFunction,"compute cumulative function")
	.method("dflan",&FLAN_MutationModel::computeProbability,"compute probability")
	.method("dflanda",&FLAN_MutationModel::computeProbability1DerivativeAlpha,"compute probability derivative wtt alpha")
	.method("dflandr",&FLAN_MutationModel::computeProbability1DerivativeRho,"compute probability derivative wtt rho")
	.method("dflangrad",&FLAN_MutationModel::computeProbability1DerivativesAlphaRho,"compute probability derivative wtt alpha and rho")
	.method("deduce.dflan",&FLAN_MutationModel::deduceProbability,"compute probability")
	.method("deduce.dflanda",&FLAN_MutationModel::deduceProbability1DerivativeAlpha,"compute probability derivative wtt alpha")
// 	.method("pgf",&FLAN_MutationModel::computeGeneratingFunction,"compute generaying function")
	.method("MutationGFEstimation",&FLAN_MutationModel::MutationGFEstimation,"estimate alpha with GF method")
	.method("CovGFEstimation",&FLAN_MutationModel::covariance,"standard deviation of GF method")
// 	.method("getfcts",&FLAN_MutationModel::getFns,"get function")
	.method("unbias.mutprob",&FLAN_MutationModel::unbiasPiEstimation,"unbias mutprob estimation")
  ;
  
  class_<FLAN_ExponentialClone>("FlanExpClone")
	.constructor<double>()
	.constructor<double,double>()
	.method("dclone",&FLAN_ExponentialClone::computeProbability,"compute probability")
// 	.method("dclonedr",&FLAN_ExponentialClone::computeProbability1DerivativeRho,"compute probability")
	.method("pgf",&FLAN_ExponentialClone::computeGeneratingFunction2,"compute generaying function")
  ;
  
  class_<FLAN_DiracClone>("FlanDirClone")
	.constructor<double>()
	.constructor<double,double>()
	.method("dclone",&FLAN_DiracClone::computeProbability,"compute probability")
// 	.method("dclonedr",&FLAN_DiracClone::computeProbability1DerivativeRho,"compute probability")
	.method("pgf",&FLAN_DiracClone::computeGeneratingFunction2,"compute generaying function")
  ;
  
//   class_<FLAN_SimClone>("FlanSimClone")
// 	.constructor<double,double,List>()
// 	.method("rclone",&FLAN_SimClone::computeSample,"compute clone sample")
//   ;

  
//   class_<MATH_Polynom>("MathPol")
//   .constructor< std::vector<double> >()
// //   .constructor<int>()
// //   .method("square_conv",&MATH_Polynom::square_conv,"compute prod of a polynom with itself")
//   .method("square_fft",&MATH_Polynom::square_fft,"compute prod of a polynom with itself")
//   .method("get.coef",&MATH_Polynom::getCoef,"get coef of polynom")
//   .method("get.deg",&MATH_Polynom::getDegree,"get degree of polynom")
//   .method("reduce",&MATH_Polynom::reduce,"get coef of polynom")
//   ;

  //function( "newMaintenancePolicy", &newMaintenancePolicy );

}
