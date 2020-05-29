#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;
using namespace std;

class CubicInterpolation{
  MatrixXd M;// coefficients
  MatrixXd a;// solution
  MatrixXd V;// points to interpolate
  MatrixXd y;// left size of linear system

  /*
  Initialize system constraints
  */
  void init_system(const MatrixXd &V, MatrixXd &M, MatrixXd &y){
    // complete
    int n = V.rows();
    int cont = 0;

    for(size_t i=0;i<=n-2;i++){
      y(cont, 0) = 0.0;
      y(cont+1, 0) = V(i,1);
      y(cont+2, 0) = V(i+1,1);
      y(cont+3, 0) = 0.0;
      cont += 4; 
    }

    M(0,0) = 0.0;
    M(0,1) = 1.0;
    M(0,2) = 2*V(0,0);
    M(0,3) = 3*V(0,0)*V(0,0);
    M(4*n-5, 4*n-8) = 0.0;
    M(4*n-5, 4*n-7) = 1.0;
    M(4*n-5, 4*n-6) = 2*V(n-1,0);
    M(4*n-5, 4*n-5) = 3*V(n-1,0)*V(n-1,0);
    M(4*n-6, 4*n-8) = 1.0;
    M(4*n-6, 4*n-7) = V(n-1,0);
    M(4*n-6, 4*n-6) = V(n-1,0)*V(n-1,0);
    M(4*n-6, 4*n-5) = V(n-1,0)*V(n-1,0)*V(n-1,0);
    M(4*n-7, 4*n-8) = 1.0;
    M(4*n-7, 4*n-7) = V(n-2,0);
    M(4*n-7, 4*n-6) = V(n-2,0)*V(n-2,0);
    M(4*n-7, 4*n-5) = V(n-2,0)*V(n-2,0)*V(n-2,0);

    int contRow = 1;
    int contCol = 0;

    for(size_t i=0;i<=n-3;i++){
      M(contRow,contCol) = 1.0;
      M(contRow,contCol+1) = V(i,0);
      M(contRow,contCol+2) = V(i,0)*V(i,0);
      M(contRow,contCol+3) = V(i,0)*V(i,0)*V(i,0);
      M(contRow+1,contCol) = 1.0;
      M(contRow+1,contCol+1) = V(i+1,0);
      M(contRow+1,contCol+2) = V(i+1,0)*V(i+1,0);
      M(contRow+1,contCol+3) = V(i+1,0)*V(i+1,0)*V(i+1,0);
      M(contRow+2,contCol) = 0.0;
      M(contRow+2,contCol+1) = 1.0;
      M(contRow+2,contCol+2) = 2*V(i+1,0);
      M(contRow+2,contCol+3) = 3*V(i+1,0)*V(i+1,0);
      M(contRow+2,contCol+4) = 0.0;
      M(contRow+2,contCol+5) = -1.0;
      M(contRow+2,contCol+6) = -2*V(i+1,0);
      M(contRow+2,contCol+7) = -3*V(i+1,0)*V(i+1,0);
      M(contRow+3,contCol) = 0.0;
      M(contRow+3,contCol+1) = 0.0;
      M(contRow+3,contCol+2) = 2.0;
      M(contRow+3,contCol+3) = 6*V(i+1,0);
      M(contRow+3,contCol+4) = 0.0;
      M(contRow+3,contCol+5) = 0.0;
      M(contRow+3,contCol+6) = -2.0;
      M(contRow+3,contCol+7) = -6*V(i+1,0);
      contRow += 4;
      contCol += 4;
    }
  }

public:
  CubicInterpolation(const MatrixXd &V1){
    M = MatrixXd::Zero(4*(V1.rows() - 1),4*(V1.rows() - 1));
    y = MatrixXd::Zero(4*(V1.rows() - 1),1);
    V = V1;
    init_system(V1, M, y );
    // complete
    a = M.colPivHouseholderQr().solve(y);
  }


  /*
  Evaluate tangent at step i
  */
  void eval_tangent(float i, MatrixXd &dX, float x){
    // complete
    float dx = 0.0001;
    dX(i,0) = x + dx;
    dX(i,1) = eval_function(dX(i,0));
  }

  /*
  Evaluate function at time t
  */
  float eval_function(float t){
    // complete
    int n= V.rows();
    int ind;

    for(size_t i=0;i<n-1;i++){
      if(t >= V(i,0) && t <= V(i+1,0)){
        ind = i;
        break;
      }
    }

    if(t < V(0,0)) ind = 0;
    if(t > V(n-1,0)) ind = n-2;

    float ans = 0.0;

    for(int i=0;i<4;i++) ans += a(4*ind+i,0)*pow(t,i);

    return ans;
  }
};
