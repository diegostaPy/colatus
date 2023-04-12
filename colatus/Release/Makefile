CUDATOOLKIT_DIR = /usr/local/cuda
CUDASDK_DIR = $(HOME)/NVIDIA_GPU_Computing_SDK
TARGETDIR = $(HOME)/workspace/binarios
#NBodyCuda/Release

CUDATOOLKIT_LIBDIR = $(CUDATOOLKIT_DIR)/lib64/
CUDASDK_INC = $(CUDASDK_DIR)/C/common/inc/
CUDASDK_LIBDIR = $(CUDASDK_DIR)/C/lib/

BOOST = -I/usr/local/include/boost


CC = g++
NVCC =  $(CUDATOOLKIT_DIR)/bin/nvcc

MACHINE= $(uname -s)
ifeq ($(MACHINE),Darwin)
    OPENGL_INC=
    #OPENGL_LIB= -framework OpenGL
    SDL_INC= -I/opt/local/include/SDL -D_GNU_SOURCE=1 -D_THREAD_SAFE
    SDL_LIB= -L/opt/local/lib -lSDLmain -lSDL_image -lSDL -framework Cocoa
else
    OPENGL_INC= -I/usr/X11R6/include
    OPENGL_LIB= -I/usr/lib64 -lGL -lGLU
    SDL_INC= -I/usr/include/SDL -D_GNU_SOURCE=1 -D_REENTRANT
    SDL_LIB= -L/usr/lib -lSDL -lpthread 
endif

INCLUDES        = $(SDL_INC) $(BOOST) -I/usr/local/include -I/usr/include/freetype2 -I/usr/include/libxml2 -I$(CUDASDK_INC) -I$(CUDATOOLKIT_DIR)/include

LIBS    = $(SDL_LIB) -lGLEW -lfreeimage -lftgl -L$(CUDASDK_LIBDIR) -lcutil_x86_64  -lboost_system -lxml2 -lgsl -lgslcblas
# add command line parameters so we can target multiple architectures
GENCODE_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE_SM11 = -gencode=arch=compute_11,code=\"sm_11,compute_11\"
GENCODE_SM12 = -gencode=arch=compute_12,code=\"sm_12,compute_12\"
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"

GENCODE_ARCH    = $(GENCODE_SM20)

CUDAFLAGS       = $(GENCODE_SM20) -m64 -v -O2  --compiler-options -fno-strict-aliasing

CFLAGS  = -O2 -fno-strict-aliasing 

# Add source files here
EXECUTABLE      := $(TARGETDIR)/NBodySim

# CUDA source files (compiled with cudacc)
CUFILES         :=      ../Integrators.cu

# CUDA dependency files
CU_DEPS         :=

# C/C++ source files (compiled with gcc / g++)
CCFILES  := ../main.cpp \
                        ../Direct.cpp \
                    


OBJS +=  $(patsubst %.cpp,%.cpp.o,$(CCFILES))
OBJS +=  $(patsubst %.cu,%.cu.o,$(CUFILES))

# Rules and targets

all:$(EXECUTABLE)

# regras para gerar o executavel
$(EXECUTABLE):$(OBJS)
		$(NVCC) $(LIBS) $(INCLUDES) -o $@ $(OBJS)

# regras de compilação
%.cu.o: %.cu
		$(NVCC) $(CUDAFLAGS) $(INCLUDES) $< -c -o $@
        
%.cpp.o: %.cpp
		$(CC) $(CFLAGS)  $(INCLUDES) $< -c -o $@     

clean:
		rm -rf `find .. -name '*.o'` $(EXECUTABLE)
