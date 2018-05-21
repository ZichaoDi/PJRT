#include "geom.hpp"
#include "init.hpp"
#include "solver.hpp"
#include <fstream>

int nid,np;
void initialize_pirt(int argc,char **args,detector_geometry& det,Mat *A,
                     Vec *sample_x,Vec *data_b){
  PetscBool synthetic_data = PETSC_TRUE;
  PetscReal tol=1e-2,sparsity=0.1;
  PetscReal omega[4] = {-2, 2, -2, 2};
  PetscInt nGrid=100,numThetan=30,nTau,m[2];
  char inFile[255],solFile[255];

  // Initialize PETSc and MPI
  PetscInitialize(&argc,&args,(char*)0,NULL);
  /* Get core's id */
  MPI_Comm_rank(PETSC_COMM_WORLD,&nid);
  /* Get number of processors */
  MPI_Comm_size(PETSC_COMM_WORLD,&np);

  PetscPrintf(PETSC_COMM_WORLD,"PIRT -- Parallel Iterative Reconstruction Tomography\n");

  nTau = (int)ceil(sqrt(2*nGrid*nGrid));
  PetscOptionsGetInt(NULL,NULL,"-nTheta",&numThetan,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nTau",&nTau,NULL);
  PetscOptionsGetInt(NULL,NULL,"-nGrid",&nGrid,NULL);
  PetscOptionsGetReal(NULL,NULL,"-box_scale",&tol,NULL);
  PetscOptionsGetBool(NULL,NULL,"-synthetic",&synthetic_data,NULL);
  if (synthetic_data){
    strcpy(inFile,"sample.dat");
  }
  PetscOptionsGetString(NULL,NULL,"-inFile",inFile,sizeof(inFile),NULL);

  strcpy(solFile,"solution.dat");
  PetscOptionsGetString(NULL,NULL,"-solFile",solFile,sizeof(solFile),NULL);

  m[0] = nGrid; m[1] = nGrid;

  if (synthetic_data) {
    //Read in a synthetic data file
    PetscPrintf(PETSC_COMM_WORLD,"Reading in synthetic data from %s\n",inFile);
    read_synthetic_data(inFile,A,sample_x,data_b,m,&nTau,numThetan);
  } else {
    //TODO: Add in read real data file
    MatCreate(PETSC_COMM_WORLD,A);
    MatSetType(*A,MATMPIAIJ);
    MatSetSizes(*A,PETSC_DECIDE,PETSC_DECIDE,nTau*numThetan,m[0]*m[1]);
    MatSetFromOptions(*A);
    MatMPIAIJSetPreallocation(*A,sparsity*nTau*numThetan,NULL,sparsity*m[0]*m[1],NULL);
  }

  det = detector_geometry(omega,tol,numThetan,m,nTau,solFile);
  return;
}


void read_synthetic_data(char filename[255],Mat *A,Vec *sample_x,
                         Vec *data_b,PetscInt m[2],PetscInt *nTau, PetscInt numThetan){
  std::ifstream inFile;
  PetscScalar Aij;
  PetscInt i,j,iMax,jMax;
  inFile.open(filename);
  if (inFile.fail()){
    std::cout << "Error in read synthetic data file\n";
    exit(0);
  }
  iMax = -1; jMax = -1;
  //Read through the file once to get iMax, jMax
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

  MatCreateVecs(*A,sample_x,data_b);
  //Reset file to start
  inFile.clear();
  inFile.seekg(0);
  while (!inFile.eof()){
    inFile >> i;
    inFile >> j;
    inFile >> Aij;
    i = i-1; j = j-1;
    VecSetValue(*sample_x,i*iMax+j,Aij,INSERT_VALUES);
    if(i>iMax) iMax = i;
    if(j>jMax) jMax = j;
  }
  VecAssemblyBegin(*sample_x);
  VecAssemblyEnd(*sample_x);
  return;
}

void print_solution_to_file(Vec recon_sample,char filename[255]){
  PetscViewer viewer;

  PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
  PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_SYMMODU);
  VecView(recon_sample,viewer);
  PetscViewerPopFormat(viewer);
  PetscViewerDestroy(&viewer);
  PetscPrintf(PETSC_COMM_WORLD,"Solution printed to file: %s\n",filename);
  
  return;
}

void finalize_pirt(Mat *A,Vec *recon_sample,Vec *data){

  MatDestroy(A);
  VecDestroy(recon_sample);
  VecDestroy(data);
  PetscFinalize();
}
