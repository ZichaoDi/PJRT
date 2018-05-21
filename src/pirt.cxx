#include<iostream>
#include "petsctao.h"
#include "geom.hpp"
#include "init.hpp"
#include "solver.hpp"
#include <fstream>


int main(int argc,char **args)
{
  PetscReal tmp_real;
  detector_geometry det;
  Mat A;
  Vec sample, recon_sample, data;
  Tao tao;

  initialize_pirt(argc,args,det,&A,&recon_sample,&data);
  construct_matrix(A,det);

  /*
   * Solve, changing solver at runtime with:
   * ./pirt -tao_type lmvm
   * or
   * ./pirt -tau_type ntr
   * etc. Some require hessians which we have not yet coded
   */

  if (det.synthetic){
    //Generate the 'data' by doing a MatMult
    MatMult(A,recon_sample,data);
    //Make a new array for our recovered sample, and copy
    VecDuplicate(recon_sample,&sample);
    VecCopy(recon_sample,sample);
    solve(A,recon_sample,data);
  } else {
    solve(A,recon_sample,data);
  }


  //Get some information about the vector
  VecMin(recon_sample,NULL,&tmp_real);
  PetscPrintf(PETSC_COMM_WORLD,"Minimum value in solution: %f\n",tmp_real);
  VecMax(recon_sample,NULL,&tmp_real);
  PetscPrintf(PETSC_COMM_WORLD,"Maximum value in solution: %f\n\n",tmp_real);

  print_solution_to_file(recon_sample,det.solFile);
  //Uncomment below to print the full vector
  //VecView(rec_sample,PETSC_VIEWER_STDOUT_SELF);

  if (det.synthetic){
    /*
     * Get the difference of the solution and the ground truth
     * for synthetic data
     */
    //calculate sample = recon_sample + -1 * sample
    VecAXPY(sample,-1.0,recon_sample);
    //get norm
    VecNorm(sample,NORM_2,&tmp_real);
    PetscPrintf(PETSC_COMM_WORLD,"Norm of difference between truth and solution: %f\n",tmp_real);
    VecDestroy(&sample);
  }

  finalize_pirt(&A,&recon_sample,&data);
  return 0;
}
