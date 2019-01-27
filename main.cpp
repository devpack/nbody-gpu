#include "display.h"
#include "camera.h"
#include "timer.h"

#include "shader.h"
#include "render.h"

#include <sstream>
#include <vector>

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

/*--------------------------------- LOOP ------------------------------------*/

void game_loop(Display *display, Render *render, Shader *shader, Camera *camera) {

	// FPS timer
    Timer fps_timer;
	fps_timer.start();
    int frame = 0;

    // title timer
    Timer title_timer;
    title_timer.start();

	// time between current frame and last frame (around 16ms for 60fps)
	// used for keyboard camera motion normatize
	float deltaTime = 0.0f;
	float lastFrameTime = 0.0f;

	// key camera motion
	bool forward = false;
	bool backward = false;
	bool left = false;
	bool right = false;
	bool up = false;
	bool down = false;

	bool stop_motion = false;

	// cube motion
	float motion_counter = 0.0f;

	// main loop
	bool loop = true;

	while (loop)
	{
		// frame duration in ms
        float currentFrame = fps_timer.get_ticks();
        deltaTime = currentFrame - lastFrameTime;
        lastFrameTime = currentFrame;
		
		// mouse
		float mouse_xoffset = 0;
		float mouse_yoffset = 0;

		SDL_Event event;

		while (SDL_PollEvent(&event))
		{

			if(event.type == SDL_QUIT) {
				loop = false;
				break;
			}

			// mouse
			if(event.type == SDL_MOUSEMOTION) {
				mouse_xoffset += (float)event.motion.xrel;
				mouse_yoffset += (float)event.motion.yrel;
			}

			// key up
			if(event.type == SDL_KEYDOWN) {

				if(event.key.keysym.sym == SDLK_ESCAPE) {
					loop = false;
					break;
				}

				if(event.key.keysym.sym == SDLK_SPACE) {
					stop_motion = !stop_motion;
				}

				if(event.key.keysym.sym == SDLK_UP)
					forward = true;
				if (event.key.keysym.sym == SDLK_DOWN)
					backward = true;

				if(event.key.keysym.sym == SDLK_LEFT)
					left = true;
				if(event.key.keysym.sym == SDLK_RIGHT)
					right = true;

				if(event.key.keysym.sym == SDLK_q)
					up = true;
				if(event.key.keysym.sym == SDLK_w)
					down = true;
			}

			// key down
			if(event.type == SDL_KEYUP) {

				if(event.key.keysym.sym == SDLK_UP)
					forward = false;
				if (event.key.keysym.sym == SDLK_DOWN)
					backward = false;

				if(event.key.keysym.sym == SDLK_LEFT)
					left = false;
				if(event.key.keysym.sym == SDLK_RIGHT)
					right = false;

				if(event.key.keysym.sym == SDLK_q)
					up = false;
				if(event.key.keysym.sym == SDLK_w)
					down = false;
			}

		} // end while poll event


		// clear
		display -> Clear(0.0f, 0.0f, 0.0f, 1.0f);

		// glUseProgram
		shader -> Bind();

		// compute the ViewProjection matrix (projection * lookAt)
		camera -> ProcessMouse(mouse_xoffset, mouse_yoffset, true);
		camera -> ProcessKeyboard(forward, backward, left, right, up, down, deltaTime);

		// compute a Model matrix (some motion for our cube)
		glm::vec3 pos = glm::vec3();
		pos.x = sinf(motion_counter);
		glm::mat4 tr_mx = glm::translate(pos);

		glm::vec3 rotx = glm::vec3();
		rotx.x = motion_counter;
		glm::mat4 rotx_mx = glm::rotate(rotx.x, glm::vec3(1.0f, 0.0f, 0.0f));

		glm::vec3 roty = glm::vec3();
		roty.y = -motion_counter;
		glm::mat4 rot_my = glm::rotate(roty.y, glm::vec3(0.0f, 1.0f, 0.0f));

		glm::vec3 rotz = glm::vec3();
		rotz.z = -motion_counter;
		glm::mat4 rot_mz = glm::rotate(rotz.z, glm::vec3(0.0f, 0.0f, 1.0f));

		// send our MVP matrix to the currently bound shader
		// camera = view_project matrix; model matrix = tr_mx * rotx_mx * rot_my * rot_mz 
		shader -> Update(tr_mx * rotx_mx * rot_my * rot_mz, camera);
		//shader -> Update(glm::mat4(1.0f), camera);

		// vao / vbo
		render -> Draw();

		// show back buffer
		display -> SwapBuffers();

		// for our cube motion
		if(!stop_motion) {
			motion_counter += 0.01f;
		}

		// FPS
        frame++;

        // update once per sec
        if( title_timer.get_ticks() > 1000 ) 
		{
            std::stringstream s;

            s << "FPS: " << frame / ( fps_timer.get_ticks() / 1000.f );

			SDL_SetWindowTitle(display->mainWindow, s.str().c_str());

            title_timer.start();
        }

	} // end while loop
} 

float xrand(float xl, float xh)
{
	return (xl + (xh - xl) *  drand48() ); 
}

/*---------------------------------------------------------------------------*/
/*--------------------------------- MAIN ------------------------------------*/
/*---------------------------------------------------------------------------*/

// screen globals
const int screen_width = 1280;
const int screen_height = 800;
bool fullscreen = false;
bool vsync = true;

// camera globals
float keyboard_sensitivity = 0.01f;
float mouse_sensitivity = 0.1f;
float znear = 0.01f;
float zfar = 1000.0f;
float fov = 70.0f;
glm::vec3 camera_pos = glm::vec3(0, 0, 5);

int main(int argc, char* argv[]) 
{     
	// SDL screen
    Display *display = new Display(screen_width, screen_height, fullscreen, vsync);

	srand(time(NULL));

	// cube points
	std::vector<glm::vec3> cube;

    for(unsigned int i = 0; i < 10000; i++)
		cube.push_back(glm::vec3(xrand(-1.0, 1.0), xrand(-1.0, 1.0), xrand(-1.0, 1.0)));

	// data vao/vbo
	Render *render = new Render(cube);

	// shaders
	Shader *shader = new Shader("../shader.vertex", "../shader.fragment");

	// camera
	Camera *camera = new Camera(camera_pos, fov, (float)display->screen_width/(float)display->screen_height, znear, zfar, mouse_sensitivity, keyboard_sensitivity);

	// screen size
	int actual_screen_width, actual_screen_height;
	SDL_GetWindowSize(display->mainWindow, &actual_screen_width, &actual_screen_height);

	std::cout << "Screen size: " << actual_screen_width << "x" << actual_screen_height << std::endl;

	// main loop
	game_loop(display, render, shader, camera);

    delete display;
    delete camera;
    delete render;

    return 0;
}
