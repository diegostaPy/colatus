#ifndef NBODYSYSTEM_H
#define NBODYSYSTEM_H

#include "Body.h"
#include <iostream>
#include "../GLCudaBuffer.h"
#include "jrfxgl.h"
#include "../PointSpriteMaterial.h"
#include <iostream>
#include <sstream>
#include <string>

namespace NBody{
	namespace GPU {

		/*
		 * Classe para simulacao de nbodies em GPU
		 */
		class NBodySystem:public JRFXGL::Graphics::Mesh{
		public:
			NBodySystem();
			~NBodySystem();
			void render(int);
			int load(const char *,int,bool);
			int act(const char *,int,int,bool);
			int reload();
			bool isLoaded();
			int getNumBodies(){ return m_bodies->getNumBodies(); }
			//int save(const string& file_name);


			void clear();


		 private:

			Body* m_bodies;
			std::stringstream m_file_name;
			bool m_loaded;
						bool m_collision;
						bool m_showVel;

						double m_t;
						bool m_wait;

			};

	}

}
#endif // BODYSYSTEM_H
