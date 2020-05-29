#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
using namespace Eigen;
using namespace std;

class HermiteInterpolation{

  MatrixXd solutionx;// solution for x
  MatrixXd solutiony;// solution for y
  MatrixXd steps;// points to interpolate
  MatrixXd solvey;// left size of linear system
  MatrixXd solvex;// left size of linear system
  MatrixXd slopes;
  MatrixXd W;

void solve_x(){
  // complete
  solvex = MatrixXd::Zero(4*(steps.rows() - 1), 1);

  int cont = 0;

  for(size_t i=0;i<steps.rows() - 1;i++){
    solvex(cont,0) = steps(i,0);
    solvex(cont+1,0) = steps(i+1,0);
    solvex(cont+2,0) = slopes(i,0);
    solvex(cont+3,0) = slopes(i+1,0);
    cont += 4;
  }
}

void solve_y(){
  // complete
  solvey = MatrixXd::Zero(4*(steps.rows() - 1), 1);

  int cont = 0;

  for(size_t i=0;i<steps.rows() - 1;i++){
    solvey(cont,0) = steps(i,1);
    solvey(cont+1,0) = steps(i+1,1);
    solvey(cont+2,0) = slopes(i,1);
    solvey(cont+3,0) = slopes(i+1,1);
    cont += 4;
  }
}

public:
  HermiteInterpolation(const MatrixXd &V1){
    // complete
    steps = V1;
    W = MatrixXd::Zero(4*(steps.rows() - 1), 4*(steps.rows() - 1));

    for(size_t i=0;i<=W.rows() - 4;i+=4){
      W(i,i) = 1.0; W(i,i+1) = 0.0; W(i,i+2) = 0.0; W(i,i+3) = 0.0;
      W(i+1,i) = 1.0; W(i+1,i+1) = 1.0; W(i+1,i+2) = 1.0; W(i+1,i+3) = 1.0;
      W(i+2,i) = 0.0; W(i+2,i+1) = 1.0; W(i+2,i+2) = 0.0; W(i+2,i+3) = 0.0;
      W(i+3,i) = 0.0; W(i+3,i+1) = 1.0; W(i+3,i+2) = 2.0; W(i+3,i+3) = 3.0;  
    }

    slopes = MatrixXd::Zero(5,2);
    slopes(0,0) = 1; slopes(0,1) = 1;
    slopes(1,0) = 1; slopes(1,1) = 1;
    slopes(2,0) = 1.0; slopes(2,1) = 1;
    slopes(3,0) = -2.5; slopes(3,1) = 1;
    slopes(4,0) = 1.0; slopes(4,1) = 1;
    // complete
    solve_x();
    solve_y();
    solutionx = W.colPivHouseholderQr().solve(solvex);
    solutiony = W.colPivHouseholderQr().solve(solvey);
  }

  // complete linspace with corrdinates
  void eval_function(MatrixXd &linspace){
    // complete
    int n = steps.rows();
    int p = linspace.rows();
    int cont = 0;

    for(size_t i=0;i<n-1;i++){
      for(double par = 0.0; par <= 1.0; par += (1.0*(n-1))/(1.0*p)){
        if(cont < p){
          linspace(cont,0) = 0.0;
          linspace(cont,1) = 0.0;

          for(int j=0;j<4;j++){
            linspace(cont,0) += solutionx(4*i+j,0)*pow(par,j);
            linspace(cont,1) += solutiony(4*i+j,0)*pow(par,j);
          }

          cont++;
        }
      }
    }
  }

  MatrixXd get_slopes(){
    return slopes;
  }
};
