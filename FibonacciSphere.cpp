const char* docstring=""
"FibonacciSphere num_point r point.pdb\n"
"    approximately uniformly sampling 'num_point' points on sphere of\n"
"    radius 'r' using spiral approximation.\n"
"    Write the coordinate in PDB format to 'point.pdb'\n"
;

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdlib.h>

#include "RandomSphereSampling.hpp"
#include "FibonacciSphere.hpp"

using namespace std;

int main(int argc, char **argv)
{
    /* parse commad line argument */
    if(argc<3)
    {
        cerr<<docstring;
        return 0;
    }
    
    int num_point=atoi(argv[1]); // number of points to sample
    float r=strtof(argv[2],NULL); // radius of sphere

    cout<<"REMARK "<<num_point<<" points on sphere of radius "<<r<<endl;

    vector<float> temp_float(3,0.);
    vector<vector<float> > cart_coor_array(num_point,temp_float);
    FibonacciSphere(num_point,r,cart_coor_array);

    string txt=cart_coor_to_pdb(cart_coor_array);

    if(argc>3)
    {
        ofstream fp(argv[3]);
        fp<<txt;
        fp.close();
    }
    else cout<<txt;

    return 0;
}
