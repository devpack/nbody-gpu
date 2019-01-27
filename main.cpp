#include "display.h"
#include "camera.h"
#include "timer.h"

#include "shader.h"
#include "render.h"

#include "particle.h"

#include <sstream>
#include <vector>

#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

/*---------------------------------------------------------------------------*/

float xrand(float xl, float xh)
{
	return (xl + (xh - xl) *  drand48() ); 
}

/*---------------------------------------------------------------------------*/

glm::vec3 pickball(float rad)
{
	glm::vec3 v;
	float rsq;

	do {
		rsq = 0.0;

		v.x = xrand(-1.0, 1.0);
		v.y = xrand(-1.0, 1.0);
		v.z = xrand(-1.0, 1.0);

		rsq = (v.x * v.x) + (v.y * v.y) + (v.z * v.z);

	} while (rsq > 1.0);

	v.x = v.x * rad;
	v.y = v.y * rad;
	v.z = v.z * rad;

	return(v);
}

/*---------------------------------------------------------------------------*/

void compute_particles_motion(const float eps, const float dt, const int nb_particles, Render *render) 
{
	float  eps2, dr2, drabs, phii, half_dt;

	glm::vec3 dr_v3;
	glm::vec3 acci_v3;
	glm::vec3 dvel_v3;
	glm::vec3 dpos_v3;

	eps2 = eps * eps;
	half_dt = 0.5 * dt;

	// Integration time leap-frog 1/2

	for(unsigned int i = 0; i < nb_particles; i++)
	{			
		// dvel = mulvs(Acc(p), 0.5 * dt);
		dvel_v3 = render->particles[i].acc * half_dt;

		// Vel(p) = addv(Vel(p), dvel);
		render->particles[i].vel = render->particles[i].vel + dvel_v3;

		// dpos = mulvs(Vel(p), dt);
		dpos_v3 = render->particles[i].vel * dt;

		// Pos(p) = addv(Pos(p), dpos);
		render->particles[i].pos = render->particles[i].pos + dpos_v3;
	}

	// Computer PP forces O(N^2)

	for(unsigned int pi = 0; pi < nb_particles; pi++) {
		for(unsigned int pj = 0; pj < nb_particles; pj++) {
			if(pi != pj)
			{
				// dr = subv(Pos(pj), Pos(pi));
				dr_v3 = render->particles[pj].pos - render->particles[pi].pos;

				// dr2 = dotvp(dr, dr);
				dr2 = glm::dot(dr_v3, dr_v3);
				dr2 += eps2;
				drabs = sqrt(dr2);

				phii = render->particles[pj].mass / (drabs*dr2);

				// acci = mulvs(dr, phii);
				acci_v3 = dr_v3 * phii;

				// Acc(pi) = addv(Acc(pi), acci);
				render->particles[pi].acc = render->particles[pi].acc + acci_v3;
			}
		}
	}

	// Integration time leap-frog 1/2

	for(unsigned int i = 0; i < nb_particles; i++)
	{
		// dvel = mulvs(Acc(p), 0.5 * dt);
		dvel_v3 = render->particles[i].acc * half_dt;

		// Vel(p) = addv(Vel(p), dvel);
		render->particles[i].vel = render->particles[i].vel + dvel_v3;

		// Acc(p) = clrv();
		render->particles[i].acc = glm::vec3(0.0);
	}
}

/*--------------------------------- LOOP ------------------------------------*/

void game_loop(const float eps, const float dt, const int nb_particles, Display *display, Render *render, Shader *shader, Camera *camera) {

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

	// particles
	std::vector<Body> particles;

	for(unsigned int i = 0; i < nb_particles; i++)
	{
		//glm::vec3 pos = glm::vec3(xrand(-1.0, 1.0), xrand(-1.0, 1.0), xrand(-1.0, 1.0));
		glm::vec3 pos = pickball(1.0);

		glm::vec3 vel = glm::vec3(0.0, 0.0, 0.0);
		glm::vec3 acc = glm::vec3(0.0, 0.0, 0.0);
		particles.push_back( Body(pos, vel, acc, 1.0) );
	}

	// vbo / vba
	render->particles = particles;

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

		// send our MVP matrix to the currently bound shader
		// camera = view_project matrix; model matrix = tr_mx * rotx_mx * rot_my * rot_mz 
		//shader -> Update(tr_mx * rotx_mx * rot_my * rot_mz, camera);
		shader -> Update(glm::mat4(1.0f), camera);

		// particules motion
		if(!stop_motion) {
			compute_particles_motion(eps, dt, nb_particles, render);
		}

		// bind particles data to VBO
		render->Bind();

		// vao / vbo => glDrawArrays()
		render -> Draw();

		// show back buffer
		display -> SwapBuffers();

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

// particles globals
const int nb_particles = 500;  // nb body
const float eps = 0.36;        // softening
const float dt = 1.0/128.0;     // time step

// main 
int main(int argc, char* argv[]) 
{     
	srand(time(NULL));

	// SDL screen
    Display *display = new Display(screen_width, screen_height, fullscreen, vsync);

	// data vao/vbo
	Render *render = new Render();

	// shaders
	Shader *shader = new Shader("../shader.vertex", "../shader.fragment");

	// camera
	Camera *camera = new Camera(camera_pos, fov, (float)display->screen_width/(float)display->screen_height, znear, zfar, mouse_sensitivity, keyboard_sensitivity);

	// screen size
	int actual_screen_width, actual_screen_height;
	SDL_GetWindowSize(display->mainWindow, &actual_screen_width, &actual_screen_height);

	std::cout << "Screen size: " << actual_screen_width << "x" << actual_screen_height << std::endl;

	// main loop
	game_loop(eps, dt, nb_particles, display, render, shader, camera);

    delete display;
    delete camera;
    delete render;

    return 0;
}
