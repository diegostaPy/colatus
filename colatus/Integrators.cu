#include <vector_functions.h>
#include <math_functions.h>
#include <stdio.h>
#include <helper_timer.h>
#include <cuda_runtime.h>
#include "Direct.h"
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif
#define TILE_WIDTH 128

//#define __USE_SHARED__
//#define __USE_NBODY_CALC__

// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at 
// once.
#define WRAP(x,mi) (((x)<mi)?(x):(x-mi))  // Mod without divide, works on values from 0 up to 2m


__constant__ double eps_sq;
__constant__ double L;
__constant__ int mNumBodies;

__constant__ double mG;
__constant__ double m;
__constant__ int* i_snap;
__constant__ int ni_snap;
__constant__ double eta;
__constant__ double3* Pos;
__constant__ double3* Vel;
__constant__ double3* nPos;
__constant__ double3* nVel;


__constant__ double3* Pos_i0;
__constant__ double3* Vel_i0;

__constant__ double3* Pos_i1;
__constant__ double3* Vel_i1;

__constant__ double3* Pos_i2;
__constant__ double3* Vel_i2;

__constant__ double3* Pos_i3;
__constant__ double3* Vel_i3;

__constant__ double3* Pos_i4;
__constant__ double3* Vel_i4;

__constant__ double3* Pos_i5;
__constant__ double3* Vel_i5;

__constant__ double3* Pos_i6;
__constant__ double3* Vel_i6;

__constant__ double3* Pos_i7;
__constant__ double3* Vel_i7;
__constant__ double3* Acc;
__constant__ double* Norm_acc;
__constant__ double* Norm_vel;
__constant__ double* pot;
__constant__ double* d_g;

extern "C"{
__host__ __device__ double3 bodyBodyInteraction(double3 sum_mGrij_dij2, double3 ri, double3 rj, double eps_sqc){
	double3 r;
	double L_2=L*0.5;
	   r.x = fmod(rj.x+L+L_2-ri.x,L) - L_2; 	   r.y = fmod(rj.y+L+L_2-ri.y,L)-  L_2;       r.z = fmod(rj.z+L+L_2-ri.z,L) - L_2;

	  double distsqrt = r.x*r.x + r.y*r.y + r.z*r.z+eps_sqc;
	  double _mGrij_dij3 = mG*rsqrtf(distsqrt)/distsqrt;
	   sum_mGrij_dij2.x += r.x * _mGrij_dij3;
	   sum_mGrij_dij2.y += r.y * _mGrij_dij3;
	   sum_mGrij_dij2.z += r.z * _mGrij_dij3;
	   return sum_mGrij_dij2;
}
__host__ __device__ double bodyBodypot(double sum_mG_rij, double3 ri, double3 rj){
	double3 r;
	//r.x = rj.x - ri.x;	r.y = rj.y - ri.y;		r.z = rj.z - ri.z;
	double L_2=L*0.5;
    r.x = fmod(rj.x+L_2-ri.x,L) - L_2; 	   r.y = fmod(rj.y+L_2-ri.y,L)-  L_2;       r.z = fmod(rj.z+L_2-ri.z,L) - L_2;

	double distsqrt = r.x*r.x + r.y*r.y + r.z*r.z;
	double mg_rij = mG*rsqrtf(distsqrt);
	//printf("\n invDist %f",invDist);
	sum_mG_rij += mg_rij;
	return sum_mG_rij;
}
}

extern "C"{
void setL(double Lb){cudaMemcpyToSymbol(L, &Lb, sizeof(double));};

void setSofteningValue(double epsilon){
	double epsilonSquared = epsilon * epsilon;
	cudaMemcpyToSymbol(eps_sq, &epsilonSquared, sizeof(double));
}
void setNumBodies(int NumBodies){
cudaMemcpyToSymbol(mNumBodies, &NumBodies, sizeof(int));};
void setm(double mass){cudaMemcpyToSymbol(m, &mass, sizeof(double));};
void setmG(double massGrav){cudaMemcpyToSymbol(mG, &massGrav, sizeof(double));};
void setisnap(int* isnap,int nisnap){
	cudaMemcpyToSymbol(i_snap, &isnap, sizeof(int*));
	cudaMemcpyToSymbol(ni_snap, &nisnap, sizeof(int));}

void seteta(double eta_){cudaMemcpyToSymbol(eta, &eta_, sizeof(double));};

void set_address(double3* Pos_g_a,double3* Vel_g_a,double3* Pos_g_b,double3* Vel_g_b,double3* acc_g,double* norm_acc_g,double* norm_vel_g,double* pot_g,double* d_){

cudaMemcpyToSymbol(Pos, &Pos_g_a, sizeof(double3*));
cudaMemcpyToSymbol(Vel, &Vel_g_a, sizeof(double3*));
cudaMemcpyToSymbol(nPos, &Pos_g_b, sizeof(double3*));
cudaMemcpyToSymbol(nVel, &Vel_g_b, sizeof(double3*));
cudaMemcpyToSymbol(Acc, &acc_g, sizeof(double3*));
cudaMemcpyToSymbol(Norm_acc, &norm_acc_g, sizeof(double*));
cudaMemcpyToSymbol(Norm_vel, &norm_vel_g, sizeof(double*));
cudaMemcpyToSymbol(pot, &pot_g, sizeof(double*));
cudaMemcpyToSymbol(d_g, &d_, sizeof(double*));
};

void set_addressi(double3* Pos_ig0,double3* Vel_ig0,double3* Pos_ig1,double3* Vel_ig1,double3* Pos_ig2,double3* Vel_ig2,double3* Pos_ig3,double3* Vel_ig3,
		double3* Pos_ig4,double3* Vel_ig4,double3* Pos_ig5,double3* Vel_ig5,double3* Pos_ig6,double3* Vel_ig6,double3* Pos_ig7,double3* Vel_ig7){
cudaMemcpyToSymbol(Pos_i0, &Pos_ig0, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i0, &Vel_ig0, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i1, &Pos_ig1, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i1, &Vel_ig1, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i2, &Pos_ig2, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i2, &Vel_ig2, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i3, &Pos_ig3, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i3, &Vel_ig3, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i4, &Pos_ig4, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i4, &Vel_ig4, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i5, &Pos_ig5, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i5, &Vel_ig5, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i6, &Pos_ig6, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i6, &Vel_ig6, sizeof(double3*));
cudaMemcpyToSymbol(Pos_i7, &Pos_ig7, sizeof(double3*));
cudaMemcpyToSymbol(Vel_i7, &Vel_ig7, sizeof(double3*));
};
void flip_pos(double3* Pos_,double3* Vel_,double3* nPos_,double3* nVel_){
cudaMemcpyToSymbol(Pos, &Pos_, sizeof(double3*));
cudaMemcpyToSymbol(Vel, &Vel_, sizeof(double3*));
cudaMemcpyToSymbol(nPos, &nPos_, sizeof(double3*));
cudaMemcpyToSymbol(nVel, &nVel_, sizeof(double3*));
};

}

namespace NBody {
namespace Simulators {

__global__   void EPot(){
	int bx     = blockIdx.x;
	int tx     = threadIdx.x;
	int dimX   = blockDim.x;
	int idx    = bx * dimX + tx;
	double sum_mG_rij = 0.0;
	if (idx < mNumBodies){
		double3 pos = Pos[idx];
		for (int j=0; j<mNumBodies; j++)
			if(idx!=j)
				sum_mG_rij = bodyBodypot(sum_mG_rij, Pos[j],pos);
	//	printf("\n summgG %f",sum_mG_rij);
		pot[idx] = sum_mG_rij;
	}
}
__global__   void EKin(){
	int bx     = blockIdx.x;
	int tx     = threadIdx.x;
	int dimX   = blockDim.x;
	int idx    = bx * dimX + tx;
	if (idx < mNumBodies){
		double3 vel = Vel[idx];
		for (int j=0; j<mNumBodies; j++)
			pot[idx] = (vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);
	}
}


__global__ void compute_acc(double C1, double C2,double eps_comov) {
	int bx     = blockIdx.x;	int tx     = threadIdx.x;
	int dimX   = blockDim.x;	int idx    = bx * dimX + tx;
	//printf("\n dim block %d blid %d txid %d i %d",dimX,bx,tx,idx);
	if (idx <mNumBodies){
		double n_acc=0.0,n_vel=0.0;
	    double3 sum_mGrji_dij2 =make_double3(0,0,0);
		double3 acc =make_double3(0,0,0);
		double3 pos = Pos[idx];		double3 vel = Vel[idx];
		for (int j=0; j<idx; j++)
			sum_mGrji_dij2= bodyBodyInteraction(sum_mGrji_dij2, pos, Pos[j], eps_comov);
		for (int j=idx+1; j<mNumBodies; j++)
			sum_mGrji_dij2= bodyBodyInteraction(sum_mGrji_dij2, pos, Pos[j], eps_comov);
	acc.x =C1*(sum_mGrji_dij2.x)+C2*vel.x;
	acc.y =C1*(sum_mGrji_dij2.y)+C2*vel.y;
	acc.z =C1*(sum_mGrji_dij2.z)+C2*vel.z;
	n_acc=(acc.x*acc.x+acc.y*acc.y+acc.z*acc.z);
	n_vel=(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z);

	Acc[idx]=acc;
	Norm_acc[idx]=n_acc;
	Norm_vel[idx]=n_vel;

}
}

__global__ void compute_stepsize() {
			double acc_max=0.0,vel_max=0.0,acc,vel;
			for (int j=0; j<mNumBodies; j++)
				{
				acc=Norm_acc[j];
				vel=Norm_vel[j];
				acc_max=acc>acc_max?acc:acc_max;
				vel_max=vel>vel_max?vel:vel_max;
				}
			d_g[0]=acc_max;
			d_g[1]=vel_max;
		//	printf("\n accmax %f velmax %f",acc_max,vel_max);

}


__global__ void integrate_lf_1_2(int it,double C1,double C2,double C3) {
	int bx     = blockIdx.x;	int tx     = threadIdx.x;
	int dimX   = blockDim.x;	int idx    = bx * dimX + tx;

	if (idx < mNumBodies){
		double3 pos = Pos[idx];	double3 vel = Vel[idx];double3 velo = Vel[idx];    double3 acc= Acc[idx];
		velo.x =velo.x*C1 ;
		velo.y =velo.y*C1 ;
		velo.z =velo.z*C1 ;
		nVel[idx] = velo;
		vel.x =vel.x*C1  +acc.x * C2;
		vel.y =vel.y*C1  +acc.y * C2;
		vel.z =vel.z*C1  +acc.z * C2;
		// Atualizar posicao
		pos.x = fmod (pos.x+vel.x*C3+L,L);
		pos.y = fmod (pos.y+vel.y*C3+L,L);
		pos.z = fmod (pos.z+vel.z*C3+L,L);
		// Atualizar posicao global
		nPos[idx] = pos;

}
}
__global__ void integrate(int it,double C1,double C2,double C3) {
	int bx     = blockIdx.x;	int tx     = threadIdx.x;
	int dimX   = blockDim.x;	int idx    = bx * dimX + tx;

	double3 acc =make_double3(0,0,0);

	if (idx < mNumBodies){
		double3 pos = Pos[idx]; 	double3 vel = Vel[idx];		double3 acc = Acc[idx];
		// Atualizar velocidade
		vel.x =vel.x*C1  +acc.x * C2;
		vel.y =vel.y*C1  +acc.y * C2;
		vel.z =vel.z*C1  +acc.z * C2;
		// Atualizar posicao
		pos.x = fmod (pos.x+vel.x*C3+L,L);
		pos.y = fmod (pos.y+vel.y*C3+L,L);
		pos.z = fmod (pos.z+vel.z*C3+L,L);

		// Atualizar posicao global
		nVel[idx] = vel;
		nPos[idx] = pos;

}
}
__global__ void save_isnap(int it) {
	int bx     = blockIdx.x;	int tx     = threadIdx.x;
	int dimX   = blockDim.x;	int idx    = bx * dimX + tx;

	double3 pos = nPos[i_snap[idx]];
	double3 vel = nVel[i_snap[idx]];

	if (idx == 0){
			Pos_i0[it]=pos;		Vel_i0[it]=vel;
	}
	if (idx == 1){
				Pos_i1[it]=pos;		Vel_i1[it]=vel;
		}
	if (idx == 2){
				Pos_i2[it]=pos;		Vel_i2[it]=vel;
		}
	if (idx == 3){
			Pos_i3[it]=pos;		Vel_i3[it]=vel;
	}
	if (idx == 4){
		Pos_i4[it]=pos;		Vel_i4[it]=vel;
    }
	if (idx == 5){
				Pos_i5[it]=pos;		Vel_i5[it]=vel;
		}
	if (idx == 6){
					Pos_i6[it]=pos;		Vel_i6[it]=vel;
			}
	if (idx == 7){
					Pos_i7[it]=pos;		Vel_i7[it]=vel;
			}

}

}
}


extern "C" {
void g_EPot(int mNumBodies) {
	dim3 dimBlock(TILE_WIDTH, 1, 1);
	dim3 dimGrid((mNumBodies/dimBlock.x)+1, 1);

	NBody::Simulators::EPot<<<dimGrid, dimBlock>>>();
}
void g_EKin(int mNumBodies) {
	dim3 dimBlock(TILE_WIDTH, 1, 1);
	dim3 dimGrid((mNumBodies/dimBlock.x)+1, 1);

	NBody::Simulators::EKin<<<dimGrid, dimBlock>>>();
}


void g_compute_acc(double C1, double C2,double eps_comov,int mNumBodies){
	dim3 dimBlock(TILE_WIDTH, 1, 1);
	dim3 dimGrid((mNumBodies/dimBlock.x)+1, 1);
	NBody::Simulators::compute_acc<<<dimGrid, dimBlock>>>(C1,C2,eps_comov);
}

void g_save_isnap(int ni_snap,int it)
{
	dim3 dimBlock(ni_snap, 1, 1);
	dim3 dimGrid(1, 1);
	NBody::Simulators::save_isnap<<<dimGrid, dimBlock>>>(it);
}
void g_compute_stepsize(int mNumBodies){
	NBody::Simulators::compute_stepsize<<<1, 1>>>();
}


void g_integrate_lf_1_2(int it,double C1,double C2,double C3,int mNumBodies){
	dim3 dimBlock(TILE_WIDTH, 1, 1);
	dim3 dimGrid((mNumBodies/dimBlock.x)+1, 1);
	NBody::Simulators::integrate_lf_1_2<<<dimGrid, dimBlock>>>(it,C1,C2,C3);

}
void g_integrate(int it,double C1,double C2,double C3,int mNumBodies){
	dim3 dimBlock(TILE_WIDTH, 1, 1);
	dim3 dimGrid((mNumBodies/dimBlock.x)+1, 1);
	NBody::Simulators::integrate<<<dimGrid, dimBlock>>>(it,C1,C2,C3);

}


}
