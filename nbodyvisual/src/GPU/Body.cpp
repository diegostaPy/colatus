#include "../../inc/GPU/Body.h"
#include <vector>
#include <string>
#include <sstream>
#include <cstring>
#include <iostream>
#include <vector_types.h>
namespace NBody {



/*
     * Construtor
     */
    Body::Body(){
    	num_bodies = 0;
        initialized = false;
    }

    /*
     * Destrutor
     */
    Body::~Body(){
        clear();
    }

    /*
     * Limpar memoria
     */
    void
    Body::clear(){
     	 if (num_bodies > 0) {
    	            cudaFree(posi);
    	 }

         num_bodies = 0;
    }
    
    void 
    Body::createBuffers(int size) {
      cudaMalloc((void**)&posi, sizeof(float4) * size);
        // Setar posicao no vertex buffer
        positionVertexBuffer.genBuffer();
        initialized = true;
        num_bodies=size;
    }

    /*
     * Le arquivo com corpos no formato
     * massa posx posy posz velx vely velz
     */
    Body*
    Body::load(const char *file_name,int numbodies,bool file_bin)
    {
       double massa; int flag;
        double3 Pos,Vel;
        Body* b = new Body();
        FILE *fp;
       // Alocar memoria
        float4 *pos = (float4*)malloc(sizeof(float4) *numbodies);
        float4 *vel = (float4*)malloc(sizeof(float4)*numbodies);
        if (file_bin) fp = fopen( file_name, "rb");
        else fp = fopen(  file_name, "r");

        if (!file_bin)
        {
        	for (int i = 0; i < numbodies; i++){
        		flag=fscanf(fp, "%f %f %f %f %f %f %f\n",  &pos[i].w,&pos[i].x,&pos[i].y,&pos[i].z, &vel[i].x, &vel[i].y, &vel[i].z);

        	}
        }
        else
        {
        	// Primeiramente, armazenamos as posicoes
        	for (int i = 0; i < numbodies; i++) {
        		flag=fread(&massa, sizeof(double), 1, fp);
        		flag=fread(&Pos, sizeof(double3),1, fp);
        		flag=fread(&Vel, sizeof(double3),1, fp);
        		pos[i].x = Pos.x;
        		pos[i].y = Pos.y;
        		pos[i].z = Pos.z;
        		pos[i].w = massa;
        		vel[i].x = Vel.x;
        		vel[i].y = Vel.y;
        		vel[i].z = Vel.z;
        	}

        }
        fclose(fp);
        
        // Inicializar buffers
        b->createBuffers(numbodies);

        cudaMemcpy(b->posi, pos, sizeof(float4)*numbodies, cudaMemcpyHostToDevice);

        // Setar posicao no vertex buffer
        b->positionVertexBuffer.bindBuffer();
        b->positionVertexBuffer.setData(sizeof(float4)*numbodies, pos, GL_DYNAMIC_DRAW);
        b->cudaGLBuffer.setBuffer(&b->positionVertexBuffer);
        b->positionVertexBuffer.unbindBuffer();

        free(pos);
        free(vel);

        return b;
    }



    int   Body::show(const char *file_name,int numbodies,bool file_bin){
    	FILE *fp;
    	double massa;
    	int i,flag;

    	double3 Pos,Vel;
    	// Alocar memoria
    	float4 *pos = (float4*)malloc(sizeof(float4) *numbodies);
    	float4 *vel = (float4*)malloc(sizeof(float4)*numbodies);
    	if (file_bin) fp = fopen(file_name, "rb");
    	else fp = fopen(file_name, "r");

    	if (!file_bin)
    	{
    		for (int i = 0; i < numbodies; i++)
    			flag=fscanf(fp, "%f %f %f %f %f %f %f\n",  &pos[i].w,&pos[i].x,&pos[i].y,&pos[i].z, &vel[i].x, &vel[i].y, &vel[i].z);

    	}
    	else
    	{
    		// Primeiramente, armazenamos as posicoes
    		for (int i = 0; i < numbodies; i++) {
    			flag=fread(&massa, sizeof(double), 1, fp);
    			flag=fread(&Pos, sizeof(double3),1, fp);
    			flag=fread(&Vel, sizeof(double3),1, fp);
    			pos[i].x = Pos.x;
    			pos[i].y = Pos.y;
    			pos[i].z = Pos.z;
    			pos[i].w = massa;
    			vel[i].x = Vel.x;
    			vel[i].y = Vel.y;
    			vel[i].z = Vel.z;
    		}


    	}
    		fclose(fp);
    		cudaMemcpy(posi, pos , sizeof(float4)*numbodies, cudaMemcpyHostToDevice);
    		// Mapear posicao
    		cudaGLBuffer.mapBuffer();
        g_show(cudaGLBuffer.getData(), posi,numbodies);
        cudaGLBuffer.unmapBuffer();
        free(pos);
        free(vel);
        return 0;
    }

}
