#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <math.h>
#include <time.h>

#include "MathTools.hpp"

using namespace std;

#define RUNIF (float(rand())/RAND_MAX)
const float PI=3.14159265359;

/* 
 * randomly sample "num_point" points on a sphere of radius 'r'
 * by sampling points on cylinder and map it to sphere
 */
void RandomSphereSampling(int num_point,float r, 
    vector<vector<float> >& cart_coor_array)
{
    srand (time(NULL)); // random seed

    float z,theta;
    for (int i=0;i<num_point;i++)
    {
        z=2*RUNIF-1;
        theta=RUNIF*2*PI;
        cart_coor_array[i][0]=r*sqrt(1-z*z)*cos(theta); // x
        cart_coor_array[i][1]=r*sqrt(1-z*z)*sin(theta); // y
        cart_coor_array[i][2]=r*z;
    }
}

/* convert coordinate in cartesian space to PDB format */
/*
string cart_coor_to_pdb(vector<vector<float> >& cart_coor_array)
{
    stringstream buf;
    for (int i=0;i<cart_coor_array.size();i++)
    {
        buf<<setiosflags(ios::left)<<setw(6)<<"ATOM";
        buf<<resetiosflags(ios::left)<<setw(5)<<i+1<<"  CA  GLY S";
        buf<<setw(4)<<i+1<<"    "<<setiosflags(ios::fixed)<<setprecision(3);
        buf<<setw(8)<<cart_coor_array[i][0];
        buf<<setw(8)<<cart_coor_array[i][1];
        buf<<setw(8)<<cart_coor_array[i][2]<<'\n';
    }
    return buf.str();
}
*/

/* convert coordinate in cartesian space to PDB format */
string cart_coor_to_pdb(vector<vector<float> >& cart_coor_to_pdb)
{
    stringstream buf;
    int i,j;
    // list of insertion code
    string icode_list=" ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890abcdefghijklmnopqrstuvwxyz";
    for (i=0;i<cart_coor_to_pdb.size();i++)
    {
        buf<<setiosflags(ios::left)<<setw(6)<<"ATOM";
        buf<<resetiosflags(ios::left)<<setw(5)<<(i+1)%100000;
        buf<<"  CA  GLY "<<icode_list[1+int(i/100000)%icode_list.length()];
        buf<<setw(4)<<(i+1)%10000<<icode_list[int(i/10000)%icode_list.length()];
        buf<<"   "<<setiosflags(ios::fixed)<<setprecision(3);
        buf<<setw(8)<<cart_coor_to_pdb[i][0];
        buf<<setw(8)<<cart_coor_to_pdb[i][1];
        buf<<setw(8)<<cart_coor_to_pdb[i][2]<<'\n';
    }
    return buf.str();
}

/* convert coordinate in cartesian space to PDB format */
string cart_coor_to_pdb(vector<vector<vector<float> > >& mc_cart_coor_array)
{
    stringstream buf;
    int model_num=mc_cart_coor_array.size();
    int first_model_index=model_num-9999;
    if (first_model_index<0) first_model_index=0;

    int model_index=0;
    for (int i=first_model_index;i<model_num;i++)
    {
        model_index++;
        buf<<"MODEL  "<<setw(4)<<model_index<<'\n';
        buf<<cart_coor_to_pdb(mc_cart_coor_array[i])<<"ENDMDL\n";
    }
    buf<<"END\n";
    return buf.str();
}
