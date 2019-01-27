#include "render.h"

/*---------------------------------------------------------------------------*/

Render::Render(std::vector<glm::vec3> vertices)
{
	this->nb_vertices = vertices.size();

	// used for INDEX_VBO
	//this->nb_vertices = nb_indices;

	// VAO: object which contains one or more VBO (buffer)

	// allocate and assign a VAO to vao_id
    glGenVertexArrays(1, &vao_id);

	// bind our VAO as the current used object: so any operation that would affect a VAO will affect this particular VAO
	glBindVertexArray(vao_id);

	// allocate and assign our VBO(s) (currently only one for vertex)
	glGenBuffers(NB_BUFFERS, vbo);

    // bind our first VBO as being the active buffer and storing vertex attributes (any operation that would affect a VBO will affect this particular VBO)
	glBindBuffer(GL_ARRAY_BUFFER, vbo[VERTEX_VBO]);

	// copy the vertex data to our VBO
	glBufferData(GL_ARRAY_BUFFER, this->nb_vertices * sizeof(vertices[0]), &vertices[0], GL_STATIC_DRAW);

	// used for INDEX_VBO
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[INDEX_VBO]);

	// used for INDEX_VBO
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, nb_indices * sizeof(indices[0]), indices, GL_STATIC_DRAW);

    // Enable attribute index 0 as being used (our vertex VBO)
    glEnableVertexAttribArray(0);

    // Specify that our coordinate data is going into attribute index 0, and contains 3 floats per vertex
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	// unbind our VAO as the current used object: so any operation that would affect a VAO will not affect this particular VAO anymore
	glBindVertexArray(0);
	//--glDisableVertexAttribArray(0);
}

/*---------------------------------------------------------------------------*/

Render::~Render()
{
	// VBO cleanup
	glDeleteBuffers(NB_BUFFERS, vbo);

	// VAO cleanup
	glDeleteVertexArrays(1, &vao_id);
}

/*---------------------------------------------------------------------------*/

void Render::Draw()
{
	// bind our VAO as the current used object: so any operation that would affect a VAO will affect this particular VAO
	glBindVertexArray(vao_id);
    //--glEnableVertexAttribArray(0);

	// used for INDEX_VBO
	//glDrawElements(GL_LINES, this->nb_vertices, GL_UNSIGNED_INT, 0);
	
	//glDrawArrays(GL_LINE_LOOP, 0, this->nb_vertices);
	
	glDrawArrays(GL_POINTS, 0, this->nb_vertices);

	// unbind our VAO as the current used object: so any operation that would affect a VAO will not affect this particular VAO anymore
	glBindVertexArray(0);
	//--glDisableVertexAttribArray(0);
}
