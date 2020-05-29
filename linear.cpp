#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/readPLY.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <ostream>
#include <bits/stdc++.h>

using namespace Eigen;
using namespace std;

typedef pair<float, float> ff;

const float eps = 1e-6f;

class LinearInterpolation{
  MatrixXd V;
  public:
    LinearInterpolation(const MatrixXd &V0){
      V = V0;
    }

    // evaluate function at time t
    float eval_function(float t){
      // complete
      vector<ff> points;
      int n = V.rows();

      for(size_t i=0;i<n;i++) points.push_back({V(i,0), V(i,1)});

      sort(points.begin(), points.end());

      int ind;
      for(size_t i=0;i<n-1;i++){
        if(t >= points[i].first && t <= points[i+1].first){
          ind = i;
          break;
        }
      }

      if(abs(points[ind+1].first - points[ind].first) < eps) return points[ind].second;
      else return ((points[ind+1].second - points[ind].second)/(points[ind+1].first - points[ind].first))*(t - points[ind].first) + points[ind].second;
    }
  };
