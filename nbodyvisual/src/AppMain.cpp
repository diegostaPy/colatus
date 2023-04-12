

#include "../inc/AppMain.h"

#include <helper_cuda.h>
#include <helper_cuda_drvapi.h>

#include <cuda_runtime.h>



	//---------------
	// CLASSES
	//---------------

	/**
	 *
	 */
	AppMain::AppMain (params P)
	{
Par=P;
	}
		

	/**
	 *
	 */
	AppMain::~AppMain ()
	{ 

	}

	bool AppMain::OnInit()
	{
		//printf("\n OnInit start");
		cudaGLSetGLDevice(0);

		camera.getCamera().getProjection().setAspectRatio(800, 600);
        camera.getCamera().getProjection().setPlanes(0.01f, 10000);
		camera.mouseSpeed = 60;
		camera.speed = 100;


		Render::Materials::PointSpriteMaterialPtr psMaterial(
				new Render::Materials::PointSpriteMaterial());
		psMaterial->setPointScale((float)(800.0f / tanf(90.0f*0.5f*3.1415f/180.0f)));
		psMaterial->setPointRadius(0.05f);
		nbody_sytem.setMaterial(psMaterial);

	//

		nbody_sytem.load(Par.filename.c_str(),Par.numbodies,Par.BinInput);
        
		pauseSimulation = true;
       ipsTime = 0;
        ipsCount = 0;
        updateTime = 0;
       timerCounter = 0;

       std::stringstream ss;
       ss << "GPU NBody";
       std::string ms;
       ms= ss.str();
        setCaption((char*)ms.c_str());
        JRFXGL::Input::HideMouseCursor();

		return true;
	}

void AppMain::OnGUI() {
    std::stringstream stats;
 if (ipsCount > Par.snapshots-1)
    {
        updateTime = (float) ipsTime / (float) ipsCount;
        ipsCount = 0;
        ipsTime = 0;
	    timerCounter = 0;
    }
    
    stats << "Bodies: " << Par.numbodies <<  " Snaps: " << ipsCount;
    JRFXGL::Graphics::GUIManager::getInstance().label(stats.str(), 10, 20);
}

	void AppMain::OnUpdate(float dt)
	{
		camera.update(dt);
	if (!pauseSimulation)
        {
        timerCounter += dt;
     mStopWatch.reset();
            {

    	 	nbody_sytem.act(Par.filename.c_str(),Par.numbodies,ipsCount,Par.BinInput);
	    	ipsCount++;
	    	   	for (int i=1;i<1e20;i++)
	    	        		int j=i*i;
        }
	    ipsTime += mStopWatch.getElapsedMilliseconds();
      }
	}

	void AppMain::OnRender()
	{
		glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
		camera.getCamera().applyTransform();
       // Recuperar a dimensao do box
        JRFXGL::Graphics::Geometry::DrawWireBox(
            mc::math::Vec3f(0,0,0), mc::math::Vec3f(Par.L,Par.L,Par.L), mc::math::Vec3f(1,1,1));
        nbody_sytem.render(Par.numbodies);
	}




