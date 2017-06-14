const char* docstring=""
"UniformSphereSampling num_point r point.pdb\n"
"    uniformly and evenly sampling 'num_point' points on sphere of radius 'r'\n"
"    Write the coordinate in PDB format to 'point.pdb'\n"
;

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdlib.h>

#include "RandomSphereSampling.hpp"
#include "FibonacciSphere.hpp"
#include "UniformSphereSampling2.hpp"

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
    //RandomSphereSampling(num_point,r,cart_coor_array);

    int nstep=10000>10*num_point?10000-10*num_point:1;
    vector<vector<vector<float> > >mc_cart_coor_array(nstep,cart_coor_array);
    if (nstep>1)
        UniformSphereSampling_MonteCarlo(mc_cart_coor_array,r,nstep);

    //string txt=cart_coor_to_pdb(mc_cart_coor_array);
    string txt=cart_coor_to_pdb(mc_cart_coor_array[mc_cart_coor_array.size()-1]);

    if(argc>3)
    {
        ofstream fp(argv[3]);
        fp<<txt;
        fp.close();
    }
    else cout<<txt;

    return 0;
}
