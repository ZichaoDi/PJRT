#include "geom.hpp"

std::vector<double> linspace(double a, double b, size_t N) {
  double h = (b - a) / static_cast<double>(N-1);
  std::vector<double> xs(N);
  typename std::vector<double>::iterator x;
  double val;
  for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
    *x = val;
  return xs;
}


detector_geometry::detector_geometry(double omega_in[4], double tol_in,
                                     int numThetan_in, int m_in[2], int nTau_in){
  double alpha,Tau,tol1;
  point_t detS0,detE0,SourceS0,SourceE0;
  std::vector<double> knot;

  // Copy in the passed values
  for (int i=0;i<4;i++){
    omega[i] = omega_in[i]*tol_in;
  }
  m[0] = m_in[0]; m[1] = m_in[1];
  Tol = tol_in;
  numThetan = numThetan_in;
  nTau = nTau_in;


  // angle of diagonal
  alpha = atan((omega[3] - omega[2])/(omega[1]-omega[0]));
  Tau = omega[1]-omega[0];
  dz[0] = (omega[1]-omega[0])/m[1];
  dz[1] = (omega[3]-omega[2])/m[0];
  tol1 = 0.5 * m[0];

  detS0 = point_t(Tau/2 *tan(alpha) + tol1*dz[0], -Tau/2-tol1*dz[0]);
  detE0 = point_t(Tau/2 *tan(alpha) + tol1*dz[0], Tau/2+tol1*dz[0]);

  knot = linspace(detS0.get<1>(),detE0.get<1>(),nTau);

  for(int i=0; i<knot.size(); i++){
    DetKnot0.push_back(point_t(detS0.get<0>(),knot[i]));
  }

  SourceS0 = point_t(-Tau/2 *tan(alpha) - tol1*dz[0], -Tau/2-tol1*dz[0]);
  SourceE0 = point_t(-Tau/2 *tan(alpha) - tol1*dz[0], Tau/2+tol1*dz[0]);

  knot = linspace(SourceS0.get<1>(),SourceE0.get<1>(),nTau);
  for(int i=0; i<knot.size(); i++){
    SourceKnot0.push_back(point_t(SourceS0.get<0>(),knot[i]));
  }
  thetan = linspace(1.,360.,numThetan);

  x = linspace(omega[0],omega[1],m[0]+1);
  y = linspace(omega[2],omega[3],m[1]+1);

  bg::append(rectangle, point_t(omega[0],omega[2]));
  bg::append(rectangle, point_t(omega[0],omega[3]));
  bg::append(rectangle, point_t(omega[1],omega[3]));
  bg::append(rectangle, point_t(omega[1],omega[2]));

}


