#include "../../inc/GPU/NBodySystem.h"

namespace NBody {
	namespace GPU {

		NBodySystem::NBodySystem() {
					m_bodies = new Body();
		}

		NBodySystem::~NBodySystem() {
					delete m_bodies;
		}

		int NBodySystem::load(const char *file_name,int numbod,bool file_bin) {
			m_loaded = false;
			 std::stringstream ss;
			ss <<file_name<<0<<".dat";
		    std::string filename;
		            		filename=ss.str();
		            		const char *c = (const char*)filename.c_str();
		        	m_bodies = Body::load(c,numbod,file_bin);
			m_loaded = true;
			return 0;
		}

		int NBodySystem::act(const char *file_name,int numbod,int snaps,bool file_bin){
            std::stringstream ss;
            if (getNumBodies() > 0){
            		ss <<file_name<<snaps<<".dat";
            	    std::string filename;
            		filename=ss.str();
            		const char *c = (const char*)filename.c_str();

            		m_bodies->show(c,numbod,file_bin);
            	   m_t = m_t + 0.001;
               }

			return 0;
		}

		void NBodySystem::render(int numbod){
			 if (getNumBodies() == 0)
			                return;

			glEnableClientState(GL_VERTEX_ARRAY);
			glEnable(GL_POINT_SPRITE);
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
			m_bodies->getGLPositionBuffer().bindBuffer();
			glVertexPointer(4, GL_FLOAT, 0, 0);
			mMaterial->apply();
			glColor3f(1, 0, 0);
			glDrawArrays(GL_POINTS, 0, numbod);
			m_bodies->getGLPositionBuffer().unbindBuffer();
			mMaterial->disable();
			glDisableClientState(GL_VERTEX_ARRAY);
			glDisable(GL_POINT_SPRITE);
		}


        

	}
}
