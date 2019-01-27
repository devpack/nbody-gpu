#ifndef __SHADER_H_
#define __SHADER_H_

#include <string>
#include <GL/glew.h>

#include "camera.h"

#include "glm/glm.hpp"

class Shader
{
	public:
		Shader(const std::string& vertexShaderFilename, const std::string& fragmentShaderFilename);
		virtual ~Shader();

		void Bind();

		// apply the MVP "transform" matrix to the vertices in our vertex shader
		void Update(glm::mat4 model, Camera *camera);

	private:
		static const unsigned int NUM_SHADERS = 2;

		std::string LoadShader(const std::string& fileName);
		void CheckShaderError(GLuint shader, GLuint flag, bool isProgram, const std::string& errorMessage);
		GLuint CreateShader(const std::string& text, unsigned int type);

		GLuint program;
		GLuint shaders[NUM_SHADERS];

		// our transformations
		enum
		{
			VIEW_PROJECTION,
			MODEL,
			NB_TRANSFORMATION
		};

		GLuint transformations[NB_TRANSFORMATION];
};

#endif
