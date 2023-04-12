//
//  Direct.h
//  NBodySim
//
//  Created by José Ricardo da Silva  Júnior on 06/02/12.
//  Modified by Diego H. Stalder Diaz   18/09/12
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#ifndef __DIRECT_H__

#define __DIRECT_H__

#include <string>

#include <vector_types.h>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <sstream>
#include <iostream>
//#include <boost/filesystem.hpp>
#include "GPUIntegrators.h"
#include <libxml/xmlwriter.h>
#include <libxml/xmlstring.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <sys/time.h>
//#include <tr1/random>

#include <helper_timer.h>
#include <cuda_runtime.h>

//#define TXT_OUTPUT 1
#define WORKSIZE 100000


double lenght(double3);
struct para{
float omegaM;
float omegaL;
float omegaR;
};

namespace NBody {
namespace Simulators {
struct SnapData{
	std::string filename;
	float elapsedTime;
};
extern "C"{
	double fun(double,void*);

//	double Grav_const(void);
}

/* Classe utilizada para simulacao de nbodies utilizando interacao p2p (complexidade n^2) */
class Direct {
public:
	Direct(char *filename);
	~Direct();
	void Go();
private:

	int Load_sim_parameters();
	int LoadParticleData();
	void write_snapshots_info();
	void DoCPUSimulation();
	void DoGPUSimulation();
	// Integracao de euler utilizando CPU
	void save_i_snap();
	void Bodyinteration(int,double);

	void compute_acc(double,double,double);
	void compute_acc_stepsize(double,double,double);
	void integrate_lf_1_2(double,double,double);
	void integrate(double,double,double);
	void setSaveSnapshotRedshift();
	void compute_energyGPU();
	void compute_energyCPU();
	void inicialize_leap_frog_CPU();
	void inicialize_leap_frog_GPU();

	void integrate_edelta();
	void save_snapshot();

	void save_particle_info();
	void save_energy();
	double calcular_tempo(struct timeval t1, struct timeval t0);
	void copy_gputocpuP();

	void copy_cputogpu();

	void allocar_memoriaCPU();
	void allocar_memoriaGPU();
	void liberar_memoria_gpu();
	void liberar_memoriaCPU();
	void copy_gputocpu_particules();
	double EKin();
	double EPot();
   double Compute_DTau();


private:
	std::string sim_cfgFilename;
	std::string mFileOut;
	std::string mICFilename;

	double3* Pos;
	double3* Vel;



	int* i_snap;
	int* i_snapg;

	double3** Pos_i;
	double3** Vel_i;

	double3* Pos_ig0;
	double3* Vel_ig0;
	double3* Pos_ig1;
	double3* Vel_ig1;
	double3* Pos_ig2;
	double3* Vel_ig2;
	double3* Pos_ig3;
	double3* Vel_ig3;
	double3* Pos_ig4;
	double3* Vel_ig4;
	double3* Pos_ig5;
	double3* Vel_ig5;
	double3* Pos_ig6;
	double3* Vel_ig6;
	double3* Pos_ig7;
	double3* Vel_ig7;

	double3* Pos_g_a;
	double3* Pos_g_b;
	double3* Vel_g_a;
	double3* Vel_g_b;
	double3* acc;
	double3* acc_g;
	FILE *f;
	double3 sum_mGrji_dij2;
	double* norm_acc_g;
	double* norm_vel_g;
	double* ekin;
	double* epot;
	double* etotal;
	double* edelta;
	double* pot_g;
	double* pot;
	double* d_g;
	int snap;
    int it;
    int it_old;
    double eps_sq;
    double p0;
    double d[2];
    double a_up;
    double a_dow;
    double p;
    double dp;
	double dp_dt;
	double Htau;
	double dtau;
	double* tau;
	double dp_min;
	double dtau_min;
	double ek;
	double A;
	double B;
	double* a;
	double da_dtau;
	double d2a_dtau2;
    bool flip;
	double m;
	int mNumBodies;
    int ni_snap;
	double G;
	int max_steps;
	int max_snapshots;
	double epsilon;
	double omegaM;
	double omegaBM;
	double omegaL;
	double omegaR;
	double zin;
	double zf;
	double a0;
	double af;
	double H0;
	double eta;
	double alpha;
	double L;
	double z;
	bool mBinOutput;
	bool sav_snap;
	bool use_gpu;
	bool track_on;

	bool barionic_matter;
	bool sav_ener;
	bool Benchmark;
	double m_Mpc;
	double By_s;
	double _10km_m;
	double msol;
	bool mBinInput;
	char* buf;
	double* z_s;
	double tcom;
	double tproc;
};

}
}
#endif
