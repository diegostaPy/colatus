#include <vector_functions.h>
#include <math_functions.h>
#include <stdio.h>

#define TILE_WIDTH 64
#define SOFTENING_SQUARED 0.001

#define __USE_SHARED__
//#define __USE_NBODY_CALC__


// WRAP is used to force each block to start working on a different 
// chunk (and wrap around back to the beginning of the array) so that
// not all multiprocessors try to read the same memory locations at 
// once.
#define WRAP(x,m) (((x)<m)?(x):(x-m))  // Mod without divide, works on values from 0 up to 2m



namespace NBody {
    namespace Kernels {
        

        
        
        __global__ void show(float4* position,
                                        float4* pos) {

            int bx     = blockIdx.x;
            int tx     = threadIdx.x;
            int dimX   = blockDim.x;
            int idx    = bx * dimX + tx;
            float4 pos_v = pos[idx];

            position[idx]= pos_v;



        }
    }
}


   
extern "C" {
    void g_show(float4* position, float4* pos,int numbodies) {
        // Processamento do delta
        dim3 dimBlock(TILE_WIDTH, 1, 1);
        dim3 dimGrid((numbodies/dimBlock.x), 1);

        NBody::Kernels::show<<<dimGrid, dimBlock>>>(position, pos);


    }
}




