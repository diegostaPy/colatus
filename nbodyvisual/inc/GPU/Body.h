#ifndef BODY_H
#define BODY_H

#include <cuda.h>
#include <helper_cuda.h>
#include <vector_functions.h>
#include "nbody_kernels.cuh"
#include "jrfxgl.h"
#include "../GLCudaBuffer.h"

namespace NBody {
    /*
     * Classe de body para ser usada na GPU
     */
    class Body {
    private:
    	size_t num_bodies;
        cuda::GLCudaBufferFloat4 cudaGLBuffer;
		JRFXGL::Graphics::GLVertexBuffer positionVertexBuffer;
        bool initialized;

    public:
        float4* posi;


        Body();
        ~Body();

        static Body* load(const char *file_name,int numbodies,bool file_bin);

        void clear();

        int show(const char *file_name,int numbodies,bool file_bin);

        int getNumBodies(){return num_bodies;}

        JRFXGL::Graphics::GLVertexBuffer& getGLPositionBuffer(){ return positionVertexBuffer; }
        
        

    private:

        
        void createBuffers(int size);
    };
}
#endif // BODY_H
