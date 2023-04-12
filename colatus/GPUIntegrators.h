//
//  GPUIntegrators.h
//  NBodySim
//
//  Created by José Ricardo da Silva  Júnior on 08/02/12.
//  Modified by Diego H. Stalder Diaz   18/09/12
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//
extern "C" {
#include "Direct.h"
__host__ __device__ double3 bodyBodyInteraction(double3,double3, double3,double);
__host__ __device__ double bodyBodypot(double,double3, double3);

void setSofteningValue(double epsilon);
void setL(double Lb);
void setNumBodies(int);
void setm(double mass);
void setmG(double massGrav);
void setisnap(int*,int );
void seteta(double);
void flip_pos(double3*,double3*,double3*,double3*);

void set_address(double3*,double3*,double3*,double3*,double3*,double*,double*,double*,double*);
void set_addressi(double3*,double3*,double3*,double3*,double3*,double3*,double3*,double3*,
		          double3*,double3*,double3*,double3*,double3*,double3*,double3*,double3*);

void g_EPot(int);
void g_EKin(int mNumBodies);
void g_compute_acc_andtimestep(double, double,double ,int);
void g_compute_potential(double,double,double,int);

void g_save_isnap(int ni,int);
void g_compute_acc(double,double,double,int);
void g_compute_stepsize(int);

void g_integrate_lf_1_2(int,double,double,double,int);
void g_integrate(int,double,double,double,int);
}
