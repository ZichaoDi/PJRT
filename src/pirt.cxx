
#include<iostream>
#include "petsc.h"
#include "geom.hpp"

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  int nid,np;
  int N,numThetan,nTau;
  polygon_t rectangle;
  double tol = 1e-2;
  double omega[4] = {-2, 2, -2, 2};
  double sparsity = 0.1;
  int m[2];
  PetscInt Istart,Iend;
  linestring_t line;
  std::vector<point_t> intersect_pts;
  detector_geometry det;
  Mat A;

  // Initialize PETSc and MPI
  ierr = PetscInitialize(&argc,&args,(char*)0,NULL); CHKERRQ(ierr);
  /* Get core's id */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&nid); CHKERRQ(ierr);
  /* Get number of processors */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&np); CHKERRQ(ierr);

  N = 30;
  numThetan = 60;
  nTau = (int)ceil(sqrt(2*N*N));

  m[0] = N; m[1] = N;

  ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
  ierr = MatSetType(A,MATMPIAIJ); CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nTau*numThetan,m[0]*m[1]); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,sparsity*nTau*numThetan,NULL,sparsity*m[0]*m[1],NULL); CHKERRQ(ierr);
  det = detector_geometry(omega,tol,numThetan,m,nTau);

  ierr = MatGetOwnershipRange(A,&Istart,&Iend); CHKERRQ(ierr);

  for (int i=Istart;i<Iend;i++){
    //i' = (j-1)*numThetan + i where j=tau and i=theta
    int tau_idx   = i/numThetan;
    int theta_idx = i%numThetan;
    setMatrixElements(tau_idx, theta_idx, det, A);
    // intersectionSet(tau_idx, theta_idx, det);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
