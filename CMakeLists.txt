cmake_minimum_required(VERSION 3.10)
project(fix_perspective)

set(CMAKE_CXX_STANDARD 14)

find_package( OpenCV REQUIRED )

add_executable(fix_perspective main.cpp)
target_link_libraries( fix_perspective ${OpenCV_LIBS} )
