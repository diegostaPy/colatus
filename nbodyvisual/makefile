CUDATOOLKIT_DIR = /usr/local/cuda-6.5
CUDASDK_DIR =/usr/local/cuda-6.5
TARGETDIR = $(HOME)/cuda-workspaceb/binarios
JRFXGL = $(HOME)/cuda-workspaceb/jrfxgl

CUDATOOLKIT_LIBDIR = $(CUDATOOLKIT_DIR)/lib64
CUDASDK_INC = $(CUDASDK_DIR)/samples/common/inc
CUDASDK_LIBDIR = $(CUDASDK_DIR)/lib64

JRFXGL_INC = -I$(JRFXGL)
JRFXGL_LIB = -L$(TARGETDIR) -lJRFXGL
BOOST = -I/usr/local/include


CC = g++ -fpermissive
NVCC = $(CUDATOOLKIT_DIR)/bin/nvcc

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

INCLUDES	= $(JRFXGL_INC) $(SDL_INC) ${BOOST} -I/usr/local/include  -I/usr/include/freetype2 -I$(CUDATOOLKIT_DIR)/include -I$(CUDASDK_INC) -I/usr/include/libxml2 
LIBS	= $(SDL_LIB) $(JRFXGL_LIB) -lGLEW -lGL -lfreeimage -lftgl -L$(CUDASDK_LIBDIR)  -lxml2
    
# add command line parameters so we can target multiple architectures
GENCODE_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\"
GENCODE_SM11 = -gencode=arch=compute_11,code=\"sm_11,compute_11\"
GENCODE_SM12 = -gencode=arch=compute_12,code=\"sm_12,compute_12\"
GENCODE_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\"


GENCODE_ARCH    = $(GENCODE_SM20)

CUDAFLAGS	= $(GENCODE_SM20) -m64 -v -O2 --compiler-options -fno-strict-aliasing
CFLAGS	= -O2 -fno-strict-aliasing

# Add source files here
EXECUTABLE	:= $(TARGETDIR)/NBodyvisual

# CUDA source files (compiled with cudacc)
CUFILES		:= 	src/GPU/nbody_kernels.cu

# CUDA dependency files
CU_DEPS		:= 

# C/C++ source files (compiled with gcc / g++)
CCFILES	 := main.cpp \
						src/AppMain.cpp \
						src/GPU/Body.cpp \
						src/GPU/NBodySystem.cpp \
						src/PointSpriteMaterial.cpp
	
OBJS +=  $(patsubst %.cpp,%.cpp.o,$(CCFILES))
OBJS +=  $(patsubst %.cu,%.cu.o,$(CUFILES))

# Rules and targets	

all: $(EXECUTABLE)

# regras para gerar o executavel
$(EXECUTABLE) : $(OBJS) 
	    $(NVCC) $(LIBS) $(INCLUDES)-o $@ $(OBJS)

# regras de compilação
%.cu.o: %.cu
	$(NVCC) $(CUDAFLAGS) $(INCLUDES)$< -c -o $@ 

%.cpp.o: %.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $< -c -o $@

clean:
	rm -rf `find ./ -name '*.o'` $(EXECUTABLE)
