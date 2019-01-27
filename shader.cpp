#include "shader.h"
#include <iostream>
#include <fstream>

/*---------------------------------------------------------------------------*/

Shader::Shader(const std::string& vertexShaderFilename, const std::string& fragmentShaderFilename)
{
    // assign our program handle a "name"
	program = glCreateProgram();

    // creates and compiles vertex/fragment shaders
	shaders[0] = CreateShader(LoadShader(vertexShaderFilename), GL_VERTEX_SHADER);
	shaders[1] = CreateShader(LoadShader(fragmentShaderFilename), GL_FRAGMENT_SHADER);

    // attach our shaders to our program
	for(unsigned int i = 0; i < NUM_SHADERS; i++)
		glAttachShader(program, shaders[i]);

    // bind attribute index 0 (coordinates) to "position"
    // attribute locations must be setup before calling glLinkProgram.
	glBindAttribLocation(program, 0, "position");

    // link: shader => binary code uploaded to the GPU, if there is no error
	glLinkProgram(program);
	CheckShaderError(program, GL_LINK_STATUS, true, "Error linking shader program");

    // checks to see whether the executables contained in program can execute given the current OpenGL state
	glValidateProgram(program);
	CheckShaderError(program, GL_LINK_STATUS, true, "Invalid shader program");

    transformations[VIEW_PROJECTION] = glGetUniformLocation(program, "view_projection");
    transformations[MODEL] = glGetUniformLocation(program, "model");
}

/*---------------------------------------------------------------------------*/

Shader::~Shader()
{
	for(unsigned int i = 0; i < NUM_SHADERS; i++)
    {
        glDetachShader(program, shaders[i]);
        glDeleteShader(shaders[i]);
    }

	glDeleteProgram(program);
}

/*---------------------------------------------------------------------------*/

GLuint Shader::CreateShader(const std::string& text, unsigned int type)
{
    // shader handle
    GLuint shader = glCreateShader(type);

    if(shader == 0)
		std::cerr << "Error compiling shader type " << type << std::endl;

    const GLchar* p[1];
    p[0] = text.c_str();
    GLint lengths[1];
    lengths[0] = text.length();

    glShaderSource(shader, 1, p, lengths);
    glCompileShader(shader);

    CheckShaderError(shader, GL_COMPILE_STATUS, false, "Error compiling shader!");

    return shader;
}

/*---------------------------------------------------------------------------*/

void Shader::Update(glm::mat4 model, Camera *camera)
{
    glm::mat4 mvp = camera->GetViewProjection();

    glUniformMatrix4fv(transformations[VIEW_PROJECTION], 1, GL_FALSE, &mvp[0][0]);
    glUniformMatrix4fv(transformations[MODEL], 1, GL_FALSE, &model[0][0]);
}

void Shader::Bind()
{
	glUseProgram(program);
}

/*---------------------------------------------------------------------------*/

std::string Shader::LoadShader(const std::string& fileName)
{
    std::ifstream file;
    file.open((fileName).c_str());

    std::string output;
    std::string line;

    if(file.is_open())
    {
        while(file.good())
        {
            getline(file, line);
			output.append(line + "\n");
        }
    }
    else
    {
		std::cerr << "Unable to load shader: " << fileName << std::endl;
    }

    return output;
}

/*---------------------------------------------------------------------------*/

void Shader::CheckShaderError(GLuint shader, GLuint flag, bool isProgram, const std::string& errorMessage)
{
    GLint success = 0;
    GLchar error[1024] = { 0 };

    if(isProgram)
        glGetProgramiv(shader, flag, &success);
    else
        glGetShaderiv(shader, flag, &success);

    if(success == GL_FALSE)
    {
        if(isProgram)
            glGetProgramInfoLog(shader, sizeof(error), NULL, error);
        else
            glGetShaderInfoLog(shader, sizeof(error), NULL, error);

        std::cerr << errorMessage << ": '" << error << "'" << std::endl;
    }
}
