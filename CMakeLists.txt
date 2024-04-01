cmake_minimum_required(VERSION 3.1.0)
project(correlationTest VERSION 1.0.0 LANGUAGES CXX)

set(EIGEN3_INCLUDE_DIRS /usr/include/eigen3)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -march=native -fopenmp -Ofast -Wno-class-memaccess -fvisibility=hidden") #-Wno-class-memaccess to avoid warning in Eigen 

include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)

find_library(FFTW3_LIBRARIES fftw3)
if(FFTW3_LIBRARIES)
	cmake_path(GET FFTW3_LIBRARIES PARENT_PATH FFTW3_LIBRARY_DIRS)
	set(FFTW3_INCLUDE_DIRS ${FFTW3_LIBRARY_DIRS}/../include)
	message(STATUS "FFTW3_LIBRARIES:" ${FFTW3_LIBRARIES})
	message(STATUS "FFTW3_LIBRARY_DIRS:" ${FFTW3_LIBRARY_DIRS})
	message(STATUS "FFTW3_INCLUDE_DIRS:" ${FFTW3_INCLUDE_DIRS})
	include_directories(${FFTW3_INCLUDE_DIRS})
else()
  message(WARNING "library ${FFTW3_LIBRARIES} not found, noise generator is disabled")
endif()

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} ${FFTW3_LIBRARIES})

