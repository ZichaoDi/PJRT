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


bool point_comparison(const point_t& p, const point_t& q){
  return p.get<0>()<q.get<0>();
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
  bg::append(rectangle, point_t(omega[0],omega[2]));

}


void setMatrixElements(int tau_idx, int theta_idx, detector_geometry det, Mat A_mat){
  PetscErrorCode ierr;
  linestring_t line;
  point_t pt_rot;
  double theta;
  std::vector<int> index;
  std::vector<double> Lvec;
  linestring_t A;
  std::vector<point_t> Q, QC;

  theta = det.thetan[theta_idx]/180. * M_PI;

  //Rotate point 1 of line
  pt_rot = point_t(cos(theta)*det.SourceKnot0[tau_idx].get<0>()
                  -sin(theta)*det.SourceKnot0[tau_idx].get<1>(),
                  sin(theta)*det.SourceKnot0[tau_idx].get<0>()
                  +cos(theta)*det.SourceKnot0[tau_idx].get<1>());
  bg::append(line, pt_rot);

  //Rotate point 2 of line
  pt_rot = point_t(cos(theta)*det.DetKnot0[tau_idx].get<0>()
                   -sin(theta)*det.DetKnot0[tau_idx].get<1>(),
                   sin(theta)*det.DetKnot0[tau_idx].get<0>()
                   +cos(theta)*det.DetKnot0[tau_idx].get<1>());
  bg::append(line, pt_rot);

  //Find intersection
  bg::intersection(line,det.rectangle,A);
  bg::unique(A);

  if (A.size()==0||A.size()==1){
    //Do nothing
  } else {
    if(theta==M_PI/2){
      //Only use first point in special cases because second point
      // is grabbed from direction
      //90 degrees - vertical line
      for (int i=0; i<det.y.size(); i++){
        Q.push_back(point_t(A[0].get<0>(),det.y[i]));
      }
    } else if(theta==0 || theta==2*M_PI){
      // 0, 360 - horizontal line
      for (int i=0; i<det.x.size(); i++){
        Q.push_back(point_t(det.x[i],A[0].get<1>()));
      }
    } else if (theta==M_PI){
      // 180 - source, detector switched location
      for (int i=(det.x.size()-1); i>=0; i--){
        Q.push_back(point_t(det.x[i],A[0].get<1>()));
      }
    } else if (theta==3*M_PI/2){
      // 270 - source, detector switched location
      for (int i=(det.y.size()-1); i>=0; i--){
        Q.push_back(point_t(A[0].get<0>(),det.y[i]));
      }
    } else {
      for (int i=0;i<det.x.size();i++){
        Q.push_back(point_t(det.x[i],
                            (A[1].get<1>() - A[0].get<1>())/(A[1].get<0>() - A[0].get<0>())*det.x[i]
                            + (A[0].get<1>()*A[1].get<0>()- A[1].get<1>()*A[0].get<0>())/
                            (A[1].get<0>() - A[0].get<0>())));
      }
      for (int i=0;i<det.y.size();i++){
        Q.push_back(point_t((det.y[i] - (A[0].get<1>()*A[1].get<0>() - A[1].get<1>()*A[0].get<0>())/
                             (A[1].get<0>() - A[0].get<0>()))/
                            ((A[1].get<1>() - A[0].get<1>())/(A[1].get<0>() - A[0].get<0>())),
                            det.y[i]));
      }
    }

    int j = 0,size = Q.size();
    for (int i=0;i<size;i++){
      if(!bg::covered_by(Q[j],det.rectangle)){
        Q.erase(Q.begin() + j);
      } else {
        j++;
      }
    }
    //Sort for ascending x
    std::sort(Q.begin(),Q.end(),point_comparison);

    //Get distance between neighbors
    for (int i=0; i<Q.size()-1; i++){
      point_t tmp_point = point_t( (Q[i].get<0>() + Q[i+1].get<0>())/2,
                                   (Q[i].get<1>() + Q[i+1].get<1>())/2);

      int x_ind = (tmp_point.get<0>() - det.omega[0])/det.dz[0];
      int y_ind = (tmp_point.get<1>() - det.omega[2])/det.dz[1];
      //m[1] and m[0] are (confusingly) correct
      if (x_ind>=0&&x_ind<det.m[1]&&y_ind>=0&&y_ind<det.m[0]){
        int j_mat = x_ind * det.m[0] + y_ind;
        int i_mat = tau_idx*det.numThetan + theta_idx;
        PetscScalar add_to_mat;
        add_to_mat = bg::distance(Q[i],Q[i+1]);
        MatSetValue(A_mat,i_mat,j_mat,add_to_mat,INSERT_VALUES);
      }
    }

    // Add Lvec unique here TODO
    // auto last = std::unique(index.begin(), index.end());
    // index.erase(last,index.end());
  }
}

void intersectionSet(linestring_t line, detector_geometry det, double theta,
                     std::vector<int>& index,
                     std::vector<double>& Lvec){
  linestring_t A;
  std::vector<point_t> Q, QC;

  bg::intersection(line,det.rectangle,A);
  bg::unique(A);

  if (A.size()==0||A.size()==1){
    //Return 'empty'
    return;
  } else {
    if(theta==M_PI/2){
      //Only use first point in special cases because second point
      // is grabbed from direction
      //90 degrees - vertical line
      for (int i=0; i<det.y.size(); i++){
        Q.push_back(point_t(A[0].get<0>(),det.y[i]));
      }
    } else if(theta==0 || theta==2*M_PI){
      // 0, 360 - horizontal line
      for (int i=0; i<det.x.size(); i++){
        Q.push_back(point_t(det.x[i],A[0].get<1>()));
      }
    } else if (theta==M_PI){
      // 180 - source, detector switched location
      for (int i=(det.x.size()-1); i>=0; i--){
        Q.push_back(point_t(det.x[i],A[0].get<1>()));
      }
    } else if (theta==3*M_PI/2){
      // 270 - source, detector switched location
      for (int i=(det.y.size()-1); i>=0; i--){
        Q.push_back(point_t(A[0].get<0>(),det.y[i]));
      }
    } else {
      for (int i=0;i<det.x.size();i++){
        Q.push_back(point_t(det.x[i],
                            (A[1].get<1>() - A[0].get<1>())/(A[1].get<0>() - A[0].get<0>())*det.x[i]
                            + (A[0].get<1>()*A[1].get<0>()- A[1].get<1>()*A[0].get<0>())/
                            (A[1].get<0>() - A[0].get<0>())));
      }
      for (int i=0;i<det.y.size();i++){
        Q.push_back(point_t((det.y[i] - (A[0].get<1>()*A[1].get<0>() - A[1].get<1>()*A[0].get<0>())/
                             (A[1].get<0>() - A[0].get<0>()))/
                            ((A[1].get<1>() - A[0].get<1>())/(A[1].get<0>() - A[0].get<0>())),
                            det.y[i]));
      }
    }

    int j = 0,size = Q.size();
    for (int i=0;i<size;i++){
      if(!bg::covered_by(Q[j],det.rectangle)){
        Q.erase(Q.begin() + j);
      } else {
        j++;
      }
    }
    //Sort for ascending x
    std::sort(Q.begin(),Q.end(),point_comparison);

    //Get distance between neighbors
    for (int i=0; i<Q.size()-1; i++){
      point_t tmp_point = point_t( (Q[i].get<0>() + Q[i+1].get<0>())/2,
                                   (Q[i].get<1>() + Q[i+1].get<1>())/2);

      int x_ind = (tmp_point.get<0>() - det.omega[0])/det.dz[0];
      int y_ind = (tmp_point.get<1>() - det.omega[2])/det.dz[1];
      //m[1] and m[0] are (confusingly) correct
      if (x_ind>=0&&x_ind<det.m[1]&&y_ind>=0&&y_ind<det.m[0]){
        index.push_back(y_ind * det.m[0] + x_ind);
        Lvec.push_back(bg::distance(Q[i],Q[i+1]));
      }
    }

    // Add Lvec unique here TODO
    // auto last = std::unique(index.begin(), index.end());
    // index.erase(last,index.end());

  }
}
