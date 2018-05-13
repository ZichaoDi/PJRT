
#include <math.h>
#include <stdlib.h>
#include "petsc.h"
#include "petsctao.h"

void calculate_A(Mat);

int main(int argc,char **args){
  PetscInt n_meas,n_var,n_theta,n_tau,ierr;
  PetscInt nid,np;
  Tao tao;
  Mat A;
  Vec x,b;
  num_meas = 10;
  num_var  = 10;
  /* Initialize PETSc and MPI */
  ierr = PetscInitialize(&argc,&args,(char*)0,NULL);
  /* Get core's id */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&nid);
  /* Get number of processors */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&np);

  //  read_b(b);
  n_theta = 10;
  n_tau   = 10;
  n_meas  = n_theta * n_tau;
  n_var   = pow(n_tau,2);


  /* Setup petsc matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);
  // Maybe define allocation
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n_meas,n_var);
  ierr = MatSetFromOptions(A);
  //Set amount of memory needed - never underallocate
  ierr = MatMPIAIJSetPreallocation(A,n_tau,NULL,n_tau,NULL);
  ierr = MatCreateVecs(A,&x,&b);

  calculate_A(A);

  //Setup TAO solver

  /* Create TAO solver with initial solution method */
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao); CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOLMVM); CHKERRQ(ierr);

  /* Set solution vec and an initial guess */
  ierr = VecSet(x,0); CHKERRQ(ierr);
  ierr = TaoSetInitialVector(tao,x); CHKERRQ(ierr);

  /* Set routines for function, gradient, etc */
  /* ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,NULL); CHKERRQ(ierr); */
  /* ierr = TaoSetHessianRoutine(tao,H,H,FormHessian,NULL); CHKERRQ(ierr); */

  /* Check for TAO command line options */
  ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

  /* Solve the application */
  ierr = TaoSolve(tao); CHKERRQ(ierr);

  /* Free data structures */
  ierr = TaoDestroy(&tao); CHKERRQ(ierr);
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

  return;
}

void calculate_A(Mat A){
  PetscInt Istart,Iend,i,j,num_rows,num_cols,ierr;
  PetscScalar mat_val;
  PetscScalar omega[4] = {0.0,1.0,2.0,3.0};
  
  //Get the rows which this core owns
  ierr = MatGetOwnershipRange(A,&Istart,&Iend);
  ierr = MatGetSize(A,&num_rows,&num_cols);
  for (i=Istart;i<Iend;i++){
    for (j=i;j<num_cols;j++){
      mat_val = 1.0;
      ierr = MatSetValue(A,i,j,mat_val,ADD_VALUES);
    }
  }
}

void read_b(Vec *b){
  PetscInt size,ierr;
  //read file
}

