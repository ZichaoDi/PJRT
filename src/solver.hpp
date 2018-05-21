#ifndef SOLVER_H
#define SOLVER_H
#include "geom.hpp"
#include <petsctao.h>

typedef struct {
  Mat A;
  Vec data;
} solver_data;

void construct_matrix(Mat A, detector_geometry det);
void solve(Mat A, Vec sample, Vec data);
PetscErrorCode FormFunctionGradient(Tao tao,Vec X,PetscReal *f, Vec G,void *ptr);

#endif
