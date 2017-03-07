
#include <math.h>
#include <stdlib.h>

int main(int argc,char **args){
  int i,j,n_tau,n_theta,n_e,n_E,n_x,n_y,n_xy,e,E,v;
  double L_ijv;
  Mat A;

  n_tau = 100;
  n_theta = 1000;
  n_xy = n_x * n_y;
  n_E = 100;
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

      for(v=0;v<n_xy;v++){
        L_ijv = calculate_L_ijv(source,detector,pixel_x,pixel_v);
        I_ijv = 0;
        for (e=0;e<n_e;e++){
          I_ijv += calculate_I_ijv(L_ijv_sum,sample[v,e]);

          O_ijve = calculate_O();

          for (E=0;E<n_E;E++){
            M_eE = M[e,E];

            A_ijveE    = L_ijv*I_ijv*O_ijve*M_eE;

            i_p        = ...;
            j_p        = ...;
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

  TAO_Solve(A,x,b);
}

double calculate_L_ijv(double* source,double* detector,double* pixel_x,double* pixel_y){


  return L_ijv;
}

double calculate_I_ijv(double L_ijv_sum,double this_sample){
  double I_ijv;

  return I_ijv;
}
