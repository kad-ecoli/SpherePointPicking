#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>

#include "MathTools.hpp"
using namespace std;


/* calculate midpoint between vertices 't1' & 't2' in vertices list 'v',
 * add it to 'v' and return index */
int getMidPoint(int t1,int t2, vector<vector<float> >& v)
{
    vector<float> p1=v[t1]; // coordinate of vertice t1
    vector<float> p2=v[t2]; // coordinate of vertice t1

    vector<float> pm(3,0);  // coordinate of midpoint
    for (int i=0;i<3;i++) pm[i]=(p1[i]+p2[i])/2;
    float norm=sqrt(pm[0]*pm[0]+pm[1]*pm[1]+pm[2]*pm[2]);
    for (int i=0;i<3;i++) pm[i]/=norm; // normalize onto unit sphere
    
    v.push_back(pm);
    return v.size()-1;
}


/* recursively subdivide triangle faces */
void subdivide_icosahedron(vector<vector<float> >& v,
    vector<vector<int> >& f, int num_point)
{
    vector<int> temp_int(3,0);
    int a,b,c;
    int i,j,m,n;
    
    while (v.size()<num_point)
    {
        vector<vector<int> > f_(f.size()*4,temp_int);

        for (i=0;i<f.size();i++)
        {
            // calculate mid points (add new points to v)
            a = getMidPoint(f[i][0],f[i][1],v);
            b = getMidPoint(f[i][1],f[i][2],v);
            c = getMidPoint(f[i][2],f[i][0],v);

            int nfc[4][3]={
                {f[i][0],a,c},
                {f[i][1],b,a},
                {f[i][2],c,b},
                {      a,b,c},
            };
            for (m=0;m<4;m++)
                for (n=0;n<3;n++)
                    f_[4*i+m][n]=nfc[m][n];
        }

        int face_already_exist=-1;
        for (i=0;i<f_.size();i++)
        {
            if (i<f.size()) f[i]=f_[i];
            else f.push_back(f_[i]);
        }
    }
}

/* return a triangle with 3 vertex & 1 face */
/*
void triangle(vector<vector<float> >& vertice_list,
    vector<vector<int> >& face_list)
{
    float v[3][3]={
        {         0,  1, 0}, // v0
        {-sqrt(3)/2,-.5, 0}, // v1
        { sqrt(3)/2,-.5, 0}, // v2
    };
    int f[1][3]={
        {0,1,2}, // f1
    };

    int i,j;
    vector<float> temp_float(3,0);
    for (i=0;i<3;i++) 
    {
        for (j=0;j<3;j++) temp_float[j]=v[i][j];
        vertice_list.push_back(temp_float);
    }
    face_list.push_back(vector<int>(f[0],f[0]+sizeof f[0]/sizeof f[0][0]));
}
*/

/* return a tetrahedron with 4 vertex & 4 face */
void tetrahedron(vector<vector<float> >& vertice_list,
    vector<vector<int> >& face_list)
{
    float t=1/sqrt(3);
    float v[4][3]={
        { t, t, t}, // v0
        { t,-t,-t}, // v1
        {-t, t,-t}, // v2
        {-t,-t, t}, // v3
    };
    int f[4][3]={
        {0,1,2}, // f0
        {0,1,3}, // f1
        {0,2,3}, // f2
        {1,2,3}, // f3
    };

    int i,j;
    vector<float> temp_float(3,0);
    for (i=0;i<4;i++) 
    {
        for (j=0;j<3;j++) temp_float[j]=v[i][j];
        vertice_list.push_back(temp_float);
    }
    for (i=0;i<4;i++) face_list.push_back(
        vector<int>(f[i],f[i]+sizeof f[i]/sizeof f[i][0]));
}


/* return an octahedron with 6 vertex & 8 face */
void octahedron(vector<vector<float> >& vertice_list,
    vector<vector<int> >& face_list)
{
    float v[6][3]={
        { 1, 0, 0}, // v0
        {-1, 0, 0}, // v1
        { 0, 1, 0}, // v2
        { 0,-1, 0}, // v3
        { 0, 0, 1}, // v4
        { 0, 0,-1}, // v5
    };
    int f[8][3]={
        {0,2,4}, // f0
        {0,3,4}, // f1
        {1,2,4}, // f2
        {1,3,4}, // f3
        {0,2,5}, // f4
        {0,3,5}, // f5
        {1,2,5}, // f6
        {1,3,5}, // f7
    };

    int i,j;
    vector<float> temp_float(3,0);
    for (i=0;i<6;i++) 
    {
        for (j=0;j<3;j++) temp_float[j]=v[i][j];
        vertice_list.push_back(temp_float);
    }
    for (i=0;i<8;i++) face_list.push_back(
        vector<int>(f[i],f[i]+sizeof f[i]/sizeof f[i][0]));
}

/* return a skewed cube with 8 vertex & 8 out of 10 face */
void skewed_cube(vector<vector<float> >& vertice_list,
    vector<vector<int> >& face_list)
{
    float y=sqrt((2*sqrt(2)-1)/7);
    float x=sqrt((4-sqrt(2))/7);
    float t=sqrt(2)*x;
    float v[8][3]={
        { x, x, y}, // v0
        { x,-x, y}, // v1
        {-x, x, y}, // v2
        {-x,-x, y}, // v3
        { t, 0,-y}, // v4
        {-t, 0,-y}, // v5
        { 0, t,-y}, // v6
        { 0,-t,-y}, // v7
    };
    int f[8][3]={
        {0,1,4}, // f0
        {0,2,6}, // f1
        {2,3,5}, // f2
        {1,3,7}, // f3
        {0,4,6}, // f4
        {2,5,6}, // f5
        {3,5,7}, // f6
        {1,4,7}, // f7
    };

    int i,j;
    vector<float> temp_float(3,0);
    for (i=0;i<8;i++) 
    {
        for (j=0;j<3;j++) temp_float[j]=v[i][j];
        vertice_list.push_back(temp_float);
    }
    for (i=0;i<8;i++) face_list.push_back(
        vector<int>(f[i],f[i]+sizeof f[i]/sizeof f[i][0]));
}

/* return a unit icosahedron with 12 vertex & 20 face */
void icosahedron(vector<vector<float> >& vertice_list,
    vector<vector<int> >& face_list)
{
    float t=(1+sqrt(5))/2;
    float r=sqrt(1+t*t);
    float v[12][3]={
        {-1, t, 0}, // v0
        { 1, t, 0}, // v1
        {-1,-t, 0}, // v2
        { 1,-t, 0}, // v3
        { 0,-1, t}, // v4
        { 0, 1, t}, // v5
        { 0,-1,-t}, // v6
        { 0, 1,-t}, // v7
        { t, 0,-1}, // v8
        { t, 0, 1}, // v9
        {-t, 0,-1}, // v10
        {-t, 0, 1}, // v11
    };
    int f[20][3]={
        { 0,11, 5}, // f0
        { 0, 5, 1}, // f1
        { 0, 1, 7}, // f2
        { 0, 7,10}, // f3
        { 0,10,11}, // f4
        { 1, 5, 9}, // f5
        { 5,11, 4}, // f6
        {11,10, 2}, // f7
        {10, 7, 6}, // f8
        { 7, 1, 8}, // f9
        { 3, 9, 4}, // f10
        { 3, 4, 2}, // f11
        { 3, 2, 6}, // f12
        { 3, 6, 8}, // f13
        { 3, 8, 9}, // f14
        { 4, 9, 5}, // f15
        { 2, 4,11}, // f16
        { 6, 2,10}, // f17
        { 8, 6, 7}, // f18
        { 9, 8, 1}, // f19
    };

    int i,j;
    vector<float> temp_float(3,0);
    for (i=0;i<12;i++) 
    {
        for (j=0;j<3;j++) temp_float[j]=v[i][j]/r;
        vertice_list.push_back(temp_float);
    }
    for (i=0;i<20;i++) face_list.push_back(
        vector<int>(f[i],f[i]+sizeof f[i]/sizeof f[i][0]));
}

/* convert coordinate in cartesian space to PDB format */
string cart_coor_to_pdb(vector<vector<float> >& v,vector<vector<int> >& f)
{
    stringstream buf;
    int i,j;
    // list of insertion code
    string icode_list=" ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890abcdefghijklmnopqrstuvwxyz";
    for (i=0;i<v.size();i++)
    {
        buf<<setiosflags(ios::left)<<setw(6)<<"ATOM";
        buf<<resetiosflags(ios::left)<<setw(5)<<(i+1)%100000;
        buf<<"  CA  GLY "<<icode_list[1+int(i/100000)%icode_list.length()];
        buf<<setw(4)<<(i+1)%10000<<icode_list[int(i/10000)%icode_list.length()];
        buf<<"   "<<setiosflags(ios::fixed)<<setprecision(3);
        buf<<setw(8)<<v[i][0];
        buf<<setw(8)<<v[i][1];
        buf<<setw(8)<<v[i][2]<<'\n';
    }

    vector<int> temp_int(1,0);
    vector<vector<int> > conect(v.size(),temp_int); // vertices to be connected
    for (i=0;i<v.size();i++) conect[i][0]=i+1;
    int a,b,c; // three vertices in a triangle
    for (int i=0;i<f.size();i++)
    {
        a=f[i][0]%100000;
        b=f[i][1]%100000;
        c=f[i][2]%100000;

        if(find(conect[a].begin(), conect[a].end(), b+1) == conect[a].end())
            conect[a].push_back(b+1);
        if(find(conect[a].begin(), conect[a].end(), c+1) == conect[a].end()) 
            conect[a].push_back(c+1);

        if(find(conect[b].begin(), conect[b].end(), a+1) == conect[b].end()) 
            conect[b].push_back(a+1);
        if(find(conect[b].begin(), conect[b].end(), c+1) == conect[b].end()) 
            conect[b].push_back(c+1);

        if(find(conect[c].begin(), conect[c].end(), a+1) == conect[c].end()) 
            conect[c].push_back(a+1);
        if(find(conect[c].begin(), conect[c].end(), b+1) == conect[c].end())
            conect[c].push_back(b+1);
    }

    for (i=0;i<conect.size();i++)
    {
        for (j=1;j<conect[i].size();j+=4)
        {
            buf<<setiosflags(ios::left)<<setw(6)<<"CONECT";
            buf<<resetiosflags(ios::left)<<setiosflags(ios::fixed);
            buf<<setw(5)<<conect[i][0]<<setw(5)<<conect[i][j];
            if (j+1<conect[i].size()) buf<<setw(5)<<conect[i][j+1];
            if (j+2<conect[i].size()) buf<<setw(5)<<conect[i][j+2];
            if (j+3<conect[i].size()) buf<<setw(5)<<conect[i][j+3];
            buf<<endl;
        }
    }
    
    buf<<"END\n";
    return buf.str();
}
