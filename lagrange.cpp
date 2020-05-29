#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;

class LagrangeInterpolation{
  MatrixXd W;
  MatrixXd a;

  MatrixXd Vandermonde(const MatrixXd &V, MatrixXd &W){
    // complete
    int n = V.rows();

    for(size_t i=0;i<n;i++){
      double x = V(i,0);

      for(size_t j=1;j<n;j++) W(i,j) = pow(x,j);
    }

    return W;
  }

public:
  LagrangeInterpolation(const MatrixXd &V1){
    W = MatrixXd::Ones(V1.rows(),V1.rows());
    W = Vandermonde(V1, W);
    // complete
    int n = V1.rows();
    a.resize(V1.rows(), 1);
    MatrixXd y(V1.rows(), 1);

    for(size_t i=0;i<n;i++) y(i,0) = V1(i,1);

    a = W.colPivHouseholderQr().solve(y);
  }

  // evaluate interpolated function at time t
  float eval_function(float t){
    // complete
    int n = a.rows();
    float ans = 0.0;

    for(size_t i=0;i<n;i++) ans += a(i,0)*pow(t,i);

    return ans;
  }
};
