
#include<iostream>
#include "petsc.h"
#include "geom.hpp"

int main(int argc,char **args)
{
  // PetscInt ierr;
  // int nid,np;
  int N,numThetan,nTau;
  polygon_t rectangle;
  double tol = 1e-2;
  double omega[4] = {-2, 2, -2, 2};
  int m[2];
  linestring_t line;
  std::vector<point_t> intersect_pts;
  detector_geometry det;

  N = 30;
  numThetan = 60;
  nTau = (int)ceil(sqrt(2*N*N));

  m[0] = N; m[1] = N;

  det = detector_geometry(omega,tol,numThetan,m,nTau);
  std::cout << "constructed" << std::endl;

  bg::append(line, det.DetKnot0[20]);
  bg::append(line, det.SourceKnot0[20]);
  bg::intersection(line,det.rectangle,intersect_pts);
  std::cout << "2 " << std::endl;
  // std::cout << DetKnot0[20].get<0>() << " " << intersect_pts[i].get<0>() << std::endl;
  for(unsigned int j=0; j<intersect_pts.size(); j++){
    std::cout << intersect_pts[j].get<0>() << " " << intersect_pts[j].get<1>() << std::endl;
  }

  /* Initialize PETSc and MPI */
  // ierr = PetscInitialize(&argc,&args,(char*)0,NULL); CHKERRQ(ierr);
  // /* Get core's id */
  // ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&nid); CHKERRQ(ierr);
  // /* Get number of processors */
  // ierr = MPI_Comm_size(PETSC_COMM_WORLD,&np); CHKERRQ(ierr);


  return 0;
}
