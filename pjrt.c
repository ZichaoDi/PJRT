
#include <math.h>
#include <stdlib.h>
#include "petsc.h"
#include "petsctao.h"

double calculate_L_ijv(double*,double*,double*,double*, double);

int main(int argc,char **args){
  int i,j,n_tau,n_theta,n_e,n_E,n_x,n_y,n_xy,e,E,v,nid,np;
  int n_total_1,n_total_2,my_start,my_end,i_p,j_p,n_v;
  PetscInt ierr;
  double A_ijveE, M_eE, I_ij, L_ijv,L_ijv_sum;
  double omega[4]={-2,2,-2,2};
  double primary_source[2],primary_detector[2],source[2],detector[2];
  double xbox[4]={omega[0],omega[0],omega[1],omega[1],omega[0]};
  double ybox[4]={omega[2],omega[3],omega[3],omega[2],omega[2]};
  double **M, *theta, *MU_e, **x_sample;
  Tao tao;
  Mat A;
  Vec x;

  n_x=100;
  n_y=100;
  n_tau = floor(sqrt(n_x*n_xy+n_y*n_y));
  n_theta = 73;
  n_xy = n_x * n_y;
  n_E = 2000;
  n_e = 3;

  n_total_1 = n_tau*n_theta*n_E;
  n_total_2 = n_xy*n_e;

  /* Initialize PETSc and MPI */
  ierr = PetscInitialize(&argc,&args,(char*)0,NULL);

  /* Setup petsc matrix */
  ierr = MatCreate(PETSC_COMM_WORLD,&A);
  // Maybe define allocation
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n_total_1,n_total_2);
  ierr = MatSetFromOptions(A);

  ierr = MatMPIAIJSetPreallocation(A,0.05*n_total_1,NULL,0.05*n_total_2,NULL);

  /* Get core's id */
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&nid);
  /* Get number of processors */
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&np);

  // Allocate a 2D array; M is local to each core
  M = (double **)malloc(n_e * sizeof(double *));
  for(e=0;e<n_e;e++){
    M[e] = (double *) malloc(n_E * sizeof(double));
  }

  x_sample = (double **)malloc(n_v * sizeof(double *));
  for(v=0;v<n_v;v++){
    x_sample[e] = (double *) malloc(n_e * sizeof(double));
  }

  theta = (double *) malloc(n_theta * sizeof(double));
  theta = (double *) malloc(n_e * sizeof(double));

  my_start = (n_tau/np)*nid;     //0   10  20   30
  my_end   = (n_tau/np)*(nid+1); //10  20  30   40

  primary_source[0] = 1; primary_source[1] = 2;
  primary_detector[0] = 1; primary_detector[1] = 2;

  L_ijv_sum = 0;

  for(i=my_start;i<my_end;i++){
    for(j=0;j<n_theta;j++){
      source[0] = sin(theta[j])*(primary_source[0]); //Something like this
      source[1] = sin(theta[j])*(primary_source[1]);

      detector[0] = sin(theta[j])*(primary_detector[0]); //Something like this
      detector[1] = sin(theta[j])*(primary_detector[1]);

      I_ij = 0;
      for(v=0;v<n_v;v++){
        L_ijv = calculate_L_ijv(source,detector,xbox,ybox,theta[j]);
        if (L_ijv!=0.0) {
          for (e=0;e<n_e;e++){
            I_ij += x_sample[v][e]*L_ijv*MU_e[e];
            // O_ijve = calculate_O(); //no self-absorption for now

            for (E=0;E<n_E;E++){
              M_eE = M[e][E];

              A_ijveE    = L_ijv*I_ij*M_eE;

              i_p        = (i-1)*n_theta+j+(E-1)*(n_theta*n_tau);
              j_p        = (v-1)*n_e+e;
              //A[i_p,j_p] = A_ijveE;
              if (A_ijveE!=0){
                ierr = MatSetValue(A,i_p,j_p,A_ijveE,INSERT_VALUES);
              }
            }
          }
        L_ijv_sum += L_ijv_sum;
        }
      }
    }
  }

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

double calculate_L_ijv(double* source,double* detector,double* xbox,double* ybox, double theta){
  double L_ijv;
  L_ijv = 1.0;
  return L_ijv;
}

