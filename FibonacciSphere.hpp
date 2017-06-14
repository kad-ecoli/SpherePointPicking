#include <vector>
#include <math.h>
#include "MathTools.hpp"

using namespace std;

/* return Fibonacci Sphere with 'num_point' on a sphere of radius 'r' using
 * spiral approximation */
void FibonacciSphere(int num_point, float r,
    vector<vector<float> >& cart_coor_array)
{
    float golden_angle = PI * (3 - sqrt(5));
    float radius,theta,z;
    for (int i=0;i<num_point;i++)
    {
        theta=golden_angle*i;
        z=(1-2.*i/num_point)*(1-1./num_point);
        radius=sqrt(1-z*z);

        cart_coor_array[i][0]=radius*cos(theta)*r;
        cart_coor_array[i][1]=radius*sin(theta)*r;
        cart_coor_array[i][2]=z*r;
    }
}
