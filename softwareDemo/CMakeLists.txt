cmake_minimum_required(VERSION 3.10)
project(softwareDemo)

set(CMAKE_CXX_STANDARD 11)
set(TARGET_LIB
        "/usr/lib/x86_64-linux-gnu/libGL.so"
        "/usr/lib/x86_64-linux-gnu/libGLU.so"
        "/usr/lib/x86_64-linux-gnu/libglut.so"
        )

add_executable(softwareDemo main.cpp math/myMath.cpp math/myMath.h pipline/myPiple.cpp pipline/myPiple.h)

target_link_libraries(softwareDemo ${TARGET_LIB})
