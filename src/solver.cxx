#include "solver.hpp"

void construct_matrix(Mat A, detector_geometry det){
  PetscInt Istart,Iend,tau_idx,theta_idx;

  MatGetOwnershipRange(A,&Istart,&Iend);

  //Set up matrix
  for (int i=Istart;i<Iend;i++){
    //i' = (j-1)*numThetan + i where j=tau and i=theta
    tau_idx   = i/det.numThetan;
    theta_idx = i%det.numThetan;
    setMatrixElements(tau_idx, theta_idx, det, A);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
  return;
}

void solve(Mat A, Vec sample, Vec data){
  Tao tao;
  PetscInt n_its;
  PetscReal f_val,gnorm,cnorm,gatol,grtol,gttol;
  TaoConvergedReason reason;
  solver_data solv;

  //Set TAO things
  TaoCreate(PETSC_COMM_SELF,&tao);

  //Sotre A, data for use in gradient / hessian
  solv.A    = A;
  solv.data = data;

  TaoSetObjectiveAndGradientRoutine(tao,FormFunctionGradient,(void *)&solv);

  //Set initial guess as all 0.0
  VecSet(sample, 0.0);
  TaoSetInitialVector(tao,sample);

  //Get options from command line
  TaoSetFromOptions(tao);

  //Solve!
  PetscPrintf(PETSC_COMM_WORLD,"Solving...\n\n\n");
  TaoSolve(tao);

  //Get information about our solution
  TaoGetSolutionStatus(tao, &n_its, &f_val, &gnorm, &cnorm, NULL, &reason);
  TaoGetTolerances(tao,&gatol,&grtol,&gttol);

  if (reason==TAO_CONVERGED_GATOL){
    PetscPrintf(PETSC_COMM_WORLD,"Converged because norm of gradient is less than %e\n",gatol);
  } else if (reason==TAO_CONVERGED_GRTOL){
    PetscPrintf(PETSC_COMM_WORLD,
                "Converged because relative norm of gradient is less than %e\n",grtol);
  } else if (reason==TAO_CONVERGED_GTTOL) {
    PetscPrintf(PETSC_COMM_WORLD,"Converged because norm of gradient reduced by %e\n",gttol);
  } else if (reason==TAO_DIVERGED_MAXITS){
    PetscPrintf(PETSC_COMM_WORLD,"Did not converge! Reached maximum iterations.\n");
    PetscPrintf(PETSC_COMM_WORLD,"    Try increasing the maximum with -tao_max_it <nits>\n");
  } else if (reason==TAO_DIVERGED_MAXFCN){
    PetscPrintf(PETSC_COMM_WORLD,"Did not converge! Reached maximum function evaluations.\n");
    PetscPrintf(PETSC_COMM_WORLD,"    Try increasing the maximum with -tao_max_funcs <nfcn>\n");
  } else {
    PetscPrintf(PETSC_COMM_WORLD,"Other converged reason. Code: %d\n",reason);
  }

  PetscPrintf(PETSC_COMM_WORLD,"Number of iterations: %d\n",n_its);
  PetscPrintf(PETSC_COMM_WORLD,"Final function value: %e\n",f_val);
  PetscPrintf(PETSC_COMM_WORLD,"Final gradient norm: %e\n\n\n",gnorm);

  TaoDestroy(&tao);
  return;
}

PetscErrorCode FormFunctionGradient(Tao tao,Vec X,PetscReal *f, Vec G,void *ptr){
  Vec tmp;
  PetscReal norm;
  Mat A;
  Vec data;
  solver_data *solv = (solver_data *) ptr;
  A = solv->A;
  data = solv->data;
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

