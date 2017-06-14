const char* docstring=""
"icosphere num_point r point.pdb\n"
"    Create icosphere of radius 'r' with at least 'num_point' points\n"
"    Write the coordinate in PDB format to 'point.pdb'\n"
;

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <stdlib.h>

#include "icosphere.hpp"

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

    cout<<"REMARK at least "<<num_point<<" points on sphere"<<endl;

    vector<vector<float> > v; // vertices
    vector<vector<int> > f; // faces

    //if (num_point<=3) triangle(v,f);
    if (num_point<=4) tetrahedron(v,f);
    else if (num_point<=6) octahedron(v,f);
    else if (num_point<=8) skewed_cube(v,f);
    else icosahedron(v,f); // initizile a unit icosahedron

    subdivide_icosahedron(v,f,num_point);

    int i,j;
    for (i=0;i<v.size();i++)
        for (j=0;j<3;j++)
            v[i][j]*=r;

    cout<<"REMARK "<<v.size()<<" points on sphere of radius "<<r<<endl;
    string txt=cart_coor_to_pdb(v,f);

    if(argc>3)
    {
        ofstream fp(argv[3]);
        fp<<txt;
        fp.close();
    }
    else cout<<txt;

    return 0;
}
