#ifndef __RENDER_H_
#define __RENDER_H_

#define GLM_ENABLE_EXPERIMENTAL

#include <GL/glew.h>
#include <glm/glm.hpp>

#include "particle.h"

#include <string>
#include <vector>

/*---------------------------------------------------------------------------*/

class Render
{
	public:
		std::vector<Body> particles;

	public:
		Render() {}
		Render(std::vector<glm::vec3> vertices);
		virtual ~Render();

		void Bind(std::vector<glm::vec3> vertices);
		void Bind();

		void Draw();
	private:

		enum
		{
			VERTEX_VBO,
			INDEX_VBO,
			NB_BUFFERS
		};

		GLuint vao_id;
		GLuint vbo[NB_BUFFERS];
		unsigned int nb_vertices;
};

#endif
