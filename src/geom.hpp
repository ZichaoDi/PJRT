#ifndef GEOM_H
#define GEOM_H
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include "petsc.h"

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
typedef bg::model::polygon<point_t> polygon_t;
typedef bg::model::segment<point_t> segment_t;
typedef bg::model::linestring<point_t> linestring_t;

class detector_geometry;

std::vector<double> linspace(double a, double b, size_t N);
void setMatrixElements(int tau_idx, int theta_idx, detector_geometry det,Mat A_mat);
void intersectionSet(linestring_t line, detector_geometry det, double theta,
                     std::vector<int>& index,
                     std::vector<double>& Lvec);


// Detector / sample geometry class
class detector_geometry{

public:
  double omega[4],m[2],dz[2],Tol;
  std::vector<double> thetan,x,y; //Maybe make this point_t
  std::vector<point_t> DetKnot0,SourceKnot0;
  polygon_t rectangle;
  PetscBool synthetic;
  char solFile[255];
  int numThetan,nTau;
  detector_geometry(){
    omega[0] = -2;
    omega[1] =  2;
    omega[2] = -2;
    omega[3] = 2;
    Tol = 1e-2;
    numThetan = 60;
    m[0] = 30; m[1] = 30;
    nTau = 30;
  }

  detector_geometry(PetscReal omega_in[4], PetscReal tol_in,
                    PetscInt numThetan_in, PetscInt m_in[2], PetscInt nTau_in,
                    char solFile_in[255]);
};

#endif
