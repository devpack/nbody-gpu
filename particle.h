#ifndef __PARTICLE_H_
#define __PARTICLE_H_

#include <string>
#include <GL/glew.h>

#include "camera.h"

#include "glm/glm.hpp"

class Body
{
	public:
		glm::vec3 pos;
		glm::vec3 vel;
		glm::vec3 acc;

		float mass;

	public:
		Body(glm::vec3 pos, glm::vec3 vel, glm::vec3 acc, float mass) 
		{
			this->pos = pos;
			this->vel = vel;
			this->acc = acc;

			this->mass = mass;
		}

		virtual ~Body() {}
};

#endif
