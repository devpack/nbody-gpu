cmake_minimum_required(VERSION 3.5)

project("nbody")

find_package(SDL2 REQUIRED)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(SRC render.cpp shader.cpp timer.cpp display.cpp main.cpp)
  
add_executable(nbody ${SRC} )

target_link_libraries(nbody GLEW GL SDL2)
#target_include_directories(nbody BEFORE PUBLIC /usr/include)




