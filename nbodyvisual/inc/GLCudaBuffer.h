#ifndef __GLCudaBuffer_h__
#define __GLCudaBuffer_h__

#include <jrfxgl.h>
#include <cuda_gl_interop.h>
#include <helper_cuda.h>
#include <helper_cuda.h>
#include <helper_cuda_drvapi.h>

#include <cuda_runtime.h>
namespace cuda
{
	template <class T>
	class GLCudaBuffer
	{
	public:
		GLCudaBuffer() :
			mBuffer(NULL),
			mMapped(false)
		{
		}

		~GLCudaBuffer()
		{
			if (mBuffer)
			{
				unmapBuffer();

				// Remover registro
				cudaGraphicsUnregisterResource(this->mBufferCudaMap);
				mBuffer = NULL;
			}
		}

		void unregisterBuffer()
		{
			if (mBuffer)
			{
				unmapBuffer();

				// Remover registro
				cudaGraphicsUnregisterResource(this->mBufferCudaMap);
				mBuffer = NULL;
			}
		}

		void setBuffer(JRFXGL::Graphics::GLVertexBuffer*  vlBuffer)
		{
			unregisterBuffer();

			mBuffer = vlBuffer;

			// Registrar buffer
			cudaGraphicsGLRegisterBuffer(&this->mBufferCudaMap,mBuffer->id(), cudaGraphicsMapFlagsNone);
		}

		void mapBuffer()
		{
			assert(mBuffer && "JRFXGL::Graphics::GLVertexBuffer is NULL");

			if (!mMapped)
			{

			cudaGraphicsMapResources(1, &this->mBufferCudaMap, 0);

				size_t num_bytes;
				cudaGraphicsResourceGetMappedPointer((void**)&mDevice_data,
						&num_bytes, this->mBufferCudaMap);

				mMapped = true;
			}

		}

		void unmapBuffer()
		{
			if (mMapped)
			{
				assert(mBuffer && "JRFXGL::Graphics::GLVertexBuffer is NULL");

				cudaGraphicsUnmapResources(1, &this->mBufferCudaMap, 0);

				mMapped = false;
			}
		}

		T getData(){ return mDevice_data; }

	private:
		JRFXGL::Graphics::GLVertexBuffer* mBuffer;
		struct cudaGraphicsResource*      mBufferCudaMap;
		T                                 mDevice_data;
		bool                              mMapped;
	};
	
	typedef GLCudaBuffer<float4*> GLCudaBufferFloat4;
}
#endif
