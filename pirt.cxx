
#include<iostream>
#include "petsc.h"

int main(int argc,char **args)
{
  PetscInt ierr;
  int nid,np;
  int N,numTheta,nTau;


  N = 30;
  numTheta = 60;
  nTau = ceil(sqrt(2*N*N));


  std::cout<<"Hello World" << std::endl;
  /* Initialize PETSc and MPI */
  ierr = PetscInitialize(&argc,&args,(char*)0,NULL); CHKERRQ(ierr);
  /* Get core's id */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&nid); CHKERRQ(ierr);
  /* Get number of processors */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&np); CHKERRQ(ierr);


  return 0;
}
