
#include<iostream>
#include "petsctao.h"
#include "geom.hpp"
#include <fstream>
void read_synthetic_data(Mat *A,Vec *sample,Vec *answer,int m[2],int *nTau,int numThetan);
PetscErrorCode FormFunctionGradient(Tao tao,Vec X,PetscReal *f, Vec G,void *ptr);
Mat A;
Vec data;
int main(int argc,char **args)
{
  int nid,np;
  int N,numThetan,nTau;
  polygon_t rectangle;
  double tol = 1e-2;
  double omega[4] = {-2, 2, -2, 2};
  double sparsity = 0.05;
  int m[2];
  PetscInt Istart,Iend;
  linestring_t line;
  std::vector<point_t> intersect_pts;
  detector_geometry det;
  int synthetic_data = 1;
  Vec sample, rec_sample;
  Tao tao;

  // Initialize PETSc and MPI
  PetscInitialize(&argc,&args,(char*)0,NULL);
  /* Get core's id */
  MPI_Comm_rank(PETSC_COMM_WORLD,&nid);
  /* Get number of processors */
  MPI_Comm_size(PETSC_COMM_WORLD,&np);

  N = 300;
  numThetan = 160;
  nTau = (int)ceil(sqrt(2*N*N));
  m[0] = N; m[1] = N;

  if (synthetic_data) {
    numThetan = 130;
    read_synthetic_data(&A,&sample,&data,m,&nTau,numThetan);
  } else {
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetType(A,MATMPIAIJ);
    MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nTau*numThetan,m[0]*m[1]);
    MatSetFromOptions(A);
    MatMPIAIJSetPreallocation(A,sparsity*nTau*numThetan,NULL,sparsity*m[0]*m[1],NULL);
  }
  det = detector_geometry(omega,tol,numThetan,m,nTau);

  MatGetOwnershipRange(A,&Istart,&Iend);

  //Set up matrix
  for (int i=Istart;i<Iend;i++){
    //i' = (j-1)*numThetan + i where j=tau and i=theta
    int tau_idx   = i/numThetan;
    int theta_idx = i%numThetan;
    setMatrixElements(tau_idx, theta_idx, det, A);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

  MatMult(A,sample,data);
  VecDuplicate(sample,&rec_sample);

  //Set TAO things
  TaoCreate(PETSC_COMM_SELF,&tao);
  TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,NULL);
  VecSet(rec_sample, 0.0);
  TaoSetInitialVector(tao,rec_sample);

  /*
   * This allows one to set the solver at runtime; i.e.
   * ./pirt -tao_type lmvm
   * or
   * ./pirt -tau_type ntr
   * etc. Some require hessians which we have not yet coded
   */
  TaoSetFromOptions(tao);
  TaoSolve(tao);

  PetscReal norm;
  VecMin(rec_sample,NULL,&norm);
  std::cout << "Min val " << norm << std::endl;
  VecMax(rec_sample,NULL,&norm);
  std::cout << "Max val " << norm << std::endl;

  //Uncomment below to print the full vector
  //VecView(rec_sample,PETSC_VIEWER_STDOUT_SELF);


  //calculate rec_sample = rec_sample + -1 * sample
  VecAXPY(rec_sample,-1.0,sample);
  VecNorm(rec_sample,NORM_2,&norm);
  std::cout << "Diff norm: " << norm << std::endl;
  MatNorm(A,NORM_FROBENIUS,&norm);
  if (nid==0) std::cout << "Norm: " << norm << "\n";
  PetscFinalize();
  return 0;
}

PetscErrorCode FormFunctionGradient(Tao tao,Vec X,PetscReal *f, Vec G,void *ptr){
  Vec tmp;
  PetscReal norm;

  MatCreateVecs(A,NULL,&tmp);
  //tmp = A * X
  MatMult(A,X,tmp);

  // tmp = tmp + -1.0 * data
  VecAXPY(tmp,-1.0,data);
  VecNorm(tmp,NORM_2,&norm);
  *f = 0.5*norm*norm;

  MatMultTranspose(A,tmp,G);

  VecDestroy(&tmp);
  PetscFunctionReturn(0);
}

void read_synthetic_data(Mat *A,Vec *sample,Vec *answer,int m[2],int *nTau, int numThetan){
  std::ifstream inFile;
  PetscScalar Aij;
  PetscInt i,j,iMax,jMax;
  inFile.open("sample.dat");
  if (inFile.fail()){
    std::cout << "Error in read file\n";
    exit(0);
  }
  iMax = -1; jMax = -1;

  while (!inFile.eof()){
    inFile >> i;
    inFile >> j;
    inFile >> Aij;

    if(i>iMax) iMax = i;
    if(j>jMax) jMax = j;
  }
  m[0] = iMax; m[1] = jMax;
  *nTau = (int)ceil(sqrt(2*m[0]*m[1]));

  MatCreate(PETSC_COMM_WORLD,A);
  MatSetType(*A,MATMPIAIJ);
  MatSetSizes(*A,PETSC_DECIDE,PETSC_DECIDE,(*nTau)*numThetan,iMax*jMax);
  MatSetFromOptions(*A);
  MatMPIAIJSetPreallocation(*A,iMax*2,NULL,iMax*2,NULL);

  MatCreateVecs(*A,sample,answer);
  //Reset file to start
  inFile.clear();
  inFile.seekg(0);
  while (!inFile.eof()){
    inFile >> i;
    inFile >> j;
    inFile >> Aij;
    i = i-1; j = j-1;
    VecSetValue(*sample,i*iMax+j,Aij,INSERT_VALUES);
    if(i>iMax) iMax = i;
    if(j>jMax) jMax = j;
  }
  return;
}
