//---------------------------------------------------------------------------
#ifndef __APP_MAIN_H__
#define __APP_MAIN_H__
//---------------------------------------------------------------------------

#include "jrfxgl.h"
#include "GPU/NBodySystem.h"
#include "GPU/Body.h"
#include "PointSpriteMaterial.h"
#include <iostream>
#include <sstream>
#include <string>

#include <helper_cuda.h>
struct params {
	int steps;
	int snapshots;
	int numbodies;
	std::string filename;
    bool BinInput;
	float L;
};

	class AppMain : public JRFXGL::App::AppRunner
	{
	  public:
		virtual ~AppMain ();

	  public:
		AppMain (params P);


	  protected:
		virtual void OnRender       ();
        virtual void OnGUI          ();
		virtual bool OnInit         ();
		virtual void OnUpdate       (float dt);
		virtual void OnKeyDown(const SDL_keysym & keysym)
		{
			if (keysym.sym == SDLK_ESCAPE)
				OnQuit();


			if (keysym.sym == SDLK_p)
				pauseSimulation = !pauseSimulation;
		}


	  private:
		JRFXGL::Graphics::FPSCameraController camera;
		NBody::GPU::NBodySystem nbody_sytem;
        params Par;
		bool pauseSimulation;
        float ipsTime;
        int ipsCount;
       // int lastIpsCount;
        float updateTime;
        float timerCounter;
         JRFXGL::Util::Stopwatch mStopWatch;
	};


//---------------------------------------------------------------------------
#endif

