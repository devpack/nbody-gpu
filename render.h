#ifndef __RENDER_H_
#define __RENDER_H_

#define GLM_ENABLE_EXPERIMENTAL

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <string>
#include <vector>

/*---------------------------------------------------------------------------*/

class Render
{
	public:
		Render(std::vector<glm::vec3> vertices);
		virtual ~Render();

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
