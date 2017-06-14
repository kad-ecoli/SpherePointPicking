/* UniformSphereSampling only move one point at a time */
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>

#include <math.h>
#include <time.h>
#include <float.h>

#include "MathTools.hpp"

using namespace std;

/* calculate sphere point distribution energy
 * E=1/N \sum_{i} max_{j!=i}(1/d_{ij})
 * dij is normalzied to sphere radius of 1
 */
float SphereDistrEnergy(const vector<vector<float> >&cart_coor_array,float r,
    vector<vector<float> >& reci_dist_array, vector<float>&point_energy_vec)
{
    int num_point=cart_coor_array.size();

    int i,j;
    float dx,dy,dz,d; // distance between i and j
    for (i=0;i<num_point-1;i++)
    {
        for (j=i+1;j<num_point;j++)
        {
            dx=cart_coor_array[i][0]-cart_coor_array[j][0];
            dy=cart_coor_array[i][1]-cart_coor_array[j][1];
            dz=cart_coor_array[i][2]-cart_coor_array[j][2];
            d=sqrt(dx*dx+dy*dy+dz*dz)/r;
            reci_dist_array[i][j]=reci_dist_array[j][i]=(d==0)?FLT_MAX:1/d;
        }
    }

    float energy=0;
    for (i=0;i<num_point;i++)
    {
        point_energy_vec[i]=0;
        for (j=0;j<num_point;j++)
        {
            if (i!=j && reci_dist_array[i][j]>point_energy_vec[i])
                point_energy_vec[i]=reci_dist_array[i][j];
        }
        energy+=point_energy_vec[i];
    }
    energy/=num_point;
    return energy;
}

/* alterntive implementation to calculate sphere point distribution energy
 * E=1/N \sum_{i} max_{j!=i}(1/d_{ij})
 * dij is normalzied to sphere radius of 1
 *
 * if i and j belongs to the sane or adjacent lattice, dij is euclidian 
 * distance. otherwise, it is set to inf
 */
float GridSphereDistrEnergy(const vector<vector<float> >&cart_coor_array,
    vector<vector<float> >& reci_dist_array, vector<float>&point_energy_vec,
    float r, vector<vector<int> >&grid_coor_array)
{
    int num_point=cart_coor_array.size();

    int i,j;
    float dx,dy,dz,d; // distance between i and j
    for (i=0;i<num_point-1;i++)
    {
        for (j=i+1;j<num_point;j++)
        {
            if (abs(grid_coor_array[i][0]-grid_coor_array[j][0])<=1 &&
                abs(grid_coor_array[i][1]-grid_coor_array[j][1])<=1 &&
                abs(grid_coor_array[i][2]-grid_coor_array[j][2])<=1)
            {
                dx=cart_coor_array[i][0]-cart_coor_array[j][0];
                dy=cart_coor_array[i][1]-cart_coor_array[j][1];
                dz=cart_coor_array[i][2]-cart_coor_array[j][2];
                d=sqrt(dx*dx+dy*dy+dz*dz)/r;
                reci_dist_array[i][j]=reci_dist_array[j][i]=
                    (d==0)?FLT_MAX:1/d;
            }
            else reci_dist_array[i][j]=reci_dist_array[j][i]=0;
        }
    }

    float energy=0;
    for (i=0;i<num_point;i++)
    {
        point_energy_vec[i]=0;
        for (j=0;j<num_point;j++)
        {
            if (i!=j && reci_dist_array[i][j]>point_energy_vec[i])
                point_energy_vec[i]=reci_dist_array[i][j];
        }
        energy+=point_energy_vec[i];
    }
    energy/=num_point;
    return energy;
}


/* put points in cart_coor_array onto a grid in cartesian space.
 * return the coordinate of points in the grid in grid_coor_array.
 * both arrays must be in the same dimension
 *
 * lattice_size is the minimum distance between a pair of adjacent grid nodes
 */
void PutPointsOnGrid(vector<vector<float> >&cart_coor_array,
    vector<vector<int> >&grid_coor_array,float r,float lattice_size=0)
{
    int num_point=cart_coor_array.size();
    if (lattice_size==0) lattice_size=4/sqrt(num_point);

    for (int i=0;i<num_point;i++)
    {
        grid_coor_array[i][0]=int(cart_coor_array[i][0]/lattice_size);
        grid_coor_array[i][1]=int(cart_coor_array[i][1]/lattice_size);
        grid_coor_array[i][2]=int(cart_coor_array[i][2]/lattice_size);
    }
}

/* make random movement for point on sphere
 * cart_coor_array is N*L array of initial conformation, N is number of point
 * move_size is the movement size in relative to radius of sphere
 */
float SphereMove(vector<vector<float> >&cart_coor_array,int &resi,
    float move_size,float r,vector<vector<float> >&reci_dist_array,
    vector<float>&point_energy_vec)
{
    /** perform movement **/
    int num_point=cart_coor_array.size();

    float x,y,z; // cartesian coordinate of moving point
    float rho,theta,phi; // spherical coordinate of moving point
    float new_x,new_y,new_z,new_theta; // new location after movememnt

    resi=int(RUNIF*num_point); // moving point index
    x=cart_coor_array[resi][0];
    y=cart_coor_array[resi][1];
    z=cart_coor_array[resi][2];
        
    // cartesian coordinate to spherical coordniate
    rho=sqrt(x*x+y*y+z*z);
    //if (x==0) theta=0;
    //else theta=atan(y/x)+(x<0)*PI+(x>0&&y<0)*2*PI;
    theta=atan2(y,x);
    phi=acos(z/rho);
    
    new_z    =z    +  r*(2*RUNIF-1)*move_size;
    //new_theta=theta+ PI*(2*RUNIF-1)*move_size;
    new_theta=theta+ PI*(2*RUNIF-1)*move_size+PI*(new_z<-r || new_z>r);

    if (new_z<-r)
        new_z=-2*r-new_z;
    else if (new_z>r) 
        new_z=2*r-new_z;

    new_x=sqrt(r*r-new_z*new_z)*cos(new_theta);
    new_y=sqrt(r*r-new_z*new_z)*sin(new_theta);

    cart_coor_array[resi][0]=new_x;
    cart_coor_array[resi][1]=new_y;
    cart_coor_array[resi][2]=new_z;

    /** calculate energy **/
    int i,j;
    float dx,dy,dz,d; // distance between i and resi
    float energy=0;
    for (j=0;j<num_point;j++)
    {
        if (j==resi) continue;

        dx=cart_coor_array[resi][0]-cart_coor_array[j][0];
        dy=cart_coor_array[resi][1]-cart_coor_array[j][1];
        dz=cart_coor_array[resi][2]-cart_coor_array[j][2];
        d=sqrt(dx*dx+dy*dy+dz*dz)/r;
        reci_dist_array[resi][j]=reci_dist_array[j][resi]=(d==0)?FLT_MAX:1/d;
    }

    for (i=0;i<num_point;i++)
    {
        point_energy_vec[i]=0; //reci_dist_array[i][resi];
        for (j=0;j<num_point;j++)
        {
            if (i!=j && reci_dist_array[i][j]>point_energy_vec[i])
                point_energy_vec[i]=reci_dist_array[i][j];
        }
        energy+=point_energy_vec[i];
    }
    energy/=num_point;
    return energy;
}

/* make random movement for point on sphere
 * cart_coor_array is N*L array of initial conformation, N is number of point
 * move_size is the movement size in relative to radius of sphere
 * grid_coor_array is N*L array of grid coordinate of initial conformation.
 * lattice_size is the size of grid lattice
 */
float GridSphereMove(vector<vector<float> >&cart_coor_array,int &resi,
    float move_size,float r,vector<vector<float> >&reci_dist_array,
    vector<float>&point_energy_vec, vector<vector<int> >&grid_coor_array,
    float lattice_size)
{
    /** perform movement **/
    int num_point=cart_coor_array.size();

    float x,y,z; // cartesian coordinate of moving point
    float rho,theta,phi; // spherical coordinate of moving point
    float new_x,new_y,new_z,new_theta; // new location after movememnt

    resi=int(RUNIF*num_point); // moving point index
    x=cart_coor_array[resi][0];
    y=cart_coor_array[resi][1];
    z=cart_coor_array[resi][2];
        
    // cartesian coordinate to spherical coordniate
    rho=sqrt(x*x+y*y+z*z);
    theta=atan2(y,x);
    phi=acos(z/rho);
    
    new_z    =z    +  r*(2*RUNIF-1)*move_size;
    new_theta=theta+ PI*(2*RUNIF-1)*move_size+PI*(new_z<-r || new_z>r);

    if (new_z<-r) new_z=-2*r-new_z;
    else if (new_z>r) new_z=2*r-new_z;

    new_x=sqrt(r*r-new_z*new_z)*cos(new_theta);
    new_y=sqrt(r*r-new_z*new_z)*sin(new_theta);

    cart_coor_array[resi][0]=new_x;
    cart_coor_array[resi][1]=new_y;
    cart_coor_array[resi][2]=new_z;
    // put point on the grid
    grid_coor_array[resi][0]=int(new_x/lattice_size);
    grid_coor_array[resi][1]=int(new_y/lattice_size);
    grid_coor_array[resi][2]=int(new_x/lattice_size);

    /** calculate energy **/
    int i,j;
    float dx,dy,dz,d; // distance between i and resi
    float energy=0;
    for (j=0;j<num_point;j++)
    {
        if (j==resi) continue;

        if (abs(grid_coor_array[resi][0]-grid_coor_array[j][0])<=1 &&
            abs(grid_coor_array[resi][1]-grid_coor_array[j][1])<=1 &&
            abs(grid_coor_array[resi][2]-grid_coor_array[j][2])<=1)
        {
            dx=cart_coor_array[resi][0]-cart_coor_array[j][0];
            dy=cart_coor_array[resi][1]-cart_coor_array[j][1];
            dz=cart_coor_array[resi][2]-cart_coor_array[j][2];
            d=sqrt(dx*dx+dy*dy+dz*dz)/r;
            reci_dist_array[resi][j]=reci_dist_array[j][resi]=
                (d==0)?FLT_MAX:1/d;
        }
        else reci_dist_array[resi][j]=reci_dist_array[j][resi]=0;
    }

    for (i=0;i<num_point;i++)
    {
        point_energy_vec[i]=0; //reci_dist_array[i][resi];
        for (j=0;j<num_point;j++)
        {
            if (i!=j && reci_dist_array[i][j]>point_energy_vec[i])
                point_energy_vec[i]=reci_dist_array[i][j];
        }
        energy+=point_energy_vec[i];
    }
    energy/=num_point;
    return energy;
}

/* 
 * uniformly sample "num_point" points on a sphere of radius 'r'
 * by sampling points on cylinder and map it to sphere using Monte Carlo
 */
void UniformSphereSampling_MonteCarlo(
    vector<vector<vector<float> > >& mc_cart_coor_array,float r, int nstep)
{
    int num_point=mc_cart_coor_array[0].size();
    float T=1./num_point;    // temperature
    float sim_ann=0.999;     //factor by which temprature decreases each step
    float move_size=1./num_point; //movement size in relative to sphere radius
    float move_shrink=0.999; //factor by which move_size decreases each step
    float max_reject=1e4;    //max num reject

    float min_move_size=1./num_point; // minimum move size
    if (min_move_size>1e-2) min_move_size=1e-2;

    // discretize cartesian space into grids
    float lattice_size=4*r/sqrt(num_point);
    vector<int> grid_coor(3,0);
    vector<vector<int> > grid_coor_array(num_point,grid_coor);
    vector<vector<int> > prev_grid_coor_array(num_point,grid_coor);
    PutPointsOnGrid(mc_cart_coor_array[0],grid_coor_array,r,lattice_size);

    // energy at each step
    vector <float> energy_vec(nstep,FLT_MAX);
    // energy at each point
    vector<float> point_energy_vec(num_point,0.);
    vector<float> prev_point_energy_vec(num_point,0.);

    vector<float> temp_float(num_point,0.); // reciprocol of distance
    vector<vector<float> > reci_dist_array(num_point,temp_float);
    vector<vector<float> > prev_reci_dist_array(num_point,temp_float);

    energy_vec[0]=SphereDistrEnergy(mc_cart_coor_array[0],r,
        reci_dist_array,prev_point_energy_vec);
    energy_vec[0]=GridSphereDistrEnergy(mc_cart_coor_array[0],
        reci_dist_array,prev_point_energy_vec,r,grid_coor_array);
    cout<<"0\t"<<setiosflags(ios::fixed)<<setprecision(6)<<energy_vec[0]<<endl;

    int step,reject_count,i;
    int resi=0; // moving point index
    float energy;
    for (step=1;step<nstep;step++)
    {
        reject_count=0;
        mc_cart_coor_array[step]=mc_cart_coor_array[step-1];
        prev_grid_coor_array[resi]=grid_coor_array[resi];
        for (i=0;i<num_point;i++) prev_reci_dist_array[i]=reci_dist_array[i];

        while(1)
        {
            mc_cart_coor_array[step][resi]=mc_cart_coor_array[step-1][resi];
            grid_coor_array[resi]=prev_grid_coor_array[resi];
            for (i=0;i<num_point;i++)
            {
                reci_dist_array[i][resi]=reci_dist_array[resi][i
                    ]=prev_reci_dist_array[resi][i];
                point_energy_vec[i]=prev_point_energy_vec[i];
            }
             
            energy=SphereMove(mc_cart_coor_array[step],resi,move_size,r,
                reci_dist_array,point_energy_vec);
            //energy=GridSphereMove(mc_cart_coor_array[step],resi,move_size,r,
                //reci_dist_array,point_energy_vec,grid_coor_array,lattice_size);
            if (RUNIF<exp(-(energy-energy_vec[step-1])/T)) break;
            reject_count++;
            if (reject_count>max_reject) break; // rejection rate > 0.9999
        }
        if (reject_count>max_reject) break; // trapped in minimum
        energy_vec[step]=energy;
        for (i=0;i<num_point;i++)
            prev_point_energy_vec[i]=point_energy_vec[i];
        
        T*=sim_ann;
        move_size*=move_shrink;
        if (move_size<min_move_size) move_size=min_move_size;

        cout<<step<<'\t'<<setiosflags(ios::fixed);
        cout<<setprecision(6)<<energy_vec[step]<<endl;

    }

    mc_cart_coor_array.resize(step);
}
