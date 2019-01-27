#ifndef __DISPLAY_H_
#define __DISPLAY_H_

#include <string>
#include <iostream>

#define GL3_PROTOTYPES 1
#include <GL/glew.h>
#include <SDL2/SDL.h>

/*---------------------------------------------------------------------------*/

class Display
{
public:
       int screen_width;
       int screen_height;

        // Our SDL_Window
        SDL_Window *mainWindow;

        // Our opengl context handle
        SDL_GLContext mainContext;

public:
       Display(int screen_width, int screen_height, bool fullscreen, bool vsync);
       virtual ~Display();

       void Clear(float r, float g, float b, float a);
       void SwapBuffers();

};

#endif
