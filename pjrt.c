
#include <math.h>
#include <stdlib.h>
//#include "petsc.h"

int main(int argc,char **args){
  int i,j,n_tau,n_theta,n_e,n_E,n_x,n_y,n_xy,e,E,v;
  int n_total_1,n_total_2,my_start,my_end,i_p,j_p;
  PetscInt ierr;
  double A_ijveE;
  double omega[4]={-2,2,-2,2};
  double primary_source[2],primary_detector[2],source[2],detector[2];
  double *theta;
  double xbox[4]={omega[0],omega[0],omega[1],omega[1],omega[0]};
  double ybox[4]={omega[2],omega[3],omega[3],omega[2],omega[2]};



  Mat A, x, MU_e;

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

  my_start = (n_tau/np)*nid;     //0   10  20   30
  my_end   = (n_tau/np)*(nid+1); //10  20  30   40

  primary_source[0] = 1; primary_source[1] = 2;
  primary_detector[0] = 1; primary_detector[1] = 2;

  for(i=my_start;i<my_end;i++){
    for(j=0;j<n_theta;j++){
      source[0] = sin(theta[j])*(primary_source[0]); //Something like this
      source[1] = sin(theta[j])*(primary_source[1]);

      detector[0] = sin(theta[j])*(primary_detector[0]); //Something like this
      detector[1] = sin(theta[j])*(primary_detector[1]);

      [sub_v,Lvec] = calculate_L_ijv(source,detector,xbox,ybox,theta[j]);
      int n_subv = sizeof(sub_v);
      I_ij = 0;
      for(v=0;v<n_subv;v++){
        for (e=0;e<n_e;e++){
          I_ij += x[sub_v[v],e]*Lvec[v]*MU_e[e,1];
          // O_ijve = calculate_O(); not self-absorption for now; 

          for (E=0;E<n_E;E++){
            M_eE = M[e,E];

            A_ijveE    = Lvec[v]*I_ij*M_eE;

            i_p        = (i-1)*n_theta+j+(E-1)*(n_theta*n_tau);
            j_p        = (sub_v[v]-1)*n_e+e;
            //A[i_p,j_p] = A_ijveE;
            if (A_ijveE!=0){
              ierr = MatSetValue(A,i_p,j_p,A_ijveE,SET_VALUES);
            }
          }
        }
        L_ijv_sum += L_ijv_sum;
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
  ierr = TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,NULL); CHKERRQ(ierr);
  ierr = TaoSetHessianRoutine(tao,H,H,FormHessian,NULL); CHKERRQ(ierr);

  /* Check for TAO command line options */
  ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

  /* Solve the application */
  ierr = TaoSolve(tao); CHKERRQ(ierr);

  /* Free data structures */
  ierr = TaoDestroy(&tao); CHKERRQ(ierr);
  ierr = VecDestroy(&x); CHKERRQ(ierr);
  ierr = MatDestroy(&A); CHKERRQ(ierr);

}

double calculate_L_ijv(double* source,double* detector,double* xbox,double* ybox, double* theta){


  return L_ijv;
}

