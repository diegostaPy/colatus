/*
 * PointSpriteMaterial.h
 *
 *  Created on: 25/07/2011
 *      Author: josericardo
 */

#ifndef POINTSPRITEMATERIAL_H_
#define POINTSPRITEMATERIAL_H_

#include "jrfxgl.h"
#include <boost/shared_ptr.hpp>

namespace Render
{
	namespace Materials
	{
		class PointSpriteMaterial : public JRFXGL::Graphics::Material
		{
		public:
			PointSpriteMaterial();

			void setPointRadius(float _pointradius){ pointRadius = _pointradius; }
			void setPointScale(float _pointScale){ pointScale = _pointScale; }

			void apply();
			void disable();

		private:
			JRFXGL::Graphics::ProgramPtr program;
			float pointRadius;
			float pointScale;
		};

		typedef boost::shared_ptr<PointSpriteMaterial> PointSpriteMaterialPtr;

	}
}


#endif /* POINTSPRITEMATERIAL_H_ */
