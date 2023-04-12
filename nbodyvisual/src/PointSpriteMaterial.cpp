/*
 * PointSpriteMaterial.cpp
 *
 *  Created on: 25/07/2011
 *      Author: josericardo
 */

#include "../inc/PointSpriteMaterial.h"


namespace Render
{
	namespace Materials
	{
		PointSpriteMaterial::PointSpriteMaterial()
		{
			program.reset(new JRFXGL::Graphics::Program(
					"resources/PointSpriteRendering.vert",
					"resources/PointSpriteRendering.frag"));

		}

		void PointSpriteMaterial::apply()
		{
			program->use();
			program->setUniform("pointRadius", pointRadius);
			program->setUniform("pointScale", pointScale);
		}

		void PointSpriteMaterial::disable()
		{
			program->disable();
		}
	}
}
