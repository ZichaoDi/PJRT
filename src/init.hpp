#ifndef INIT_H
#define INIT_H
#include <petsctao.h>

extern int nid,np;
void initialize_pirt(int argc,char **args,detector_geometry& det,Mat *A, Vec *sample_x, Vec *data_b);
void read_synthetic_data(char filename[255],Mat *A,Vec *sample_x,
                         Vec *answer_x,PetscInt m[2],PetscInt *nTau,PetscInt numThetan);
void print_solution_to_file(Vec recon_sample,char filename[255]);
void finalize_pirt(Mat *A,Vec *recon_sample,Vec *data);
#endif
