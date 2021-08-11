cmake_minimum_required(VERSION 3.8)

project(RRT3D)

set(CMAKE_CXX_STANDARD 17)

find_package (OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

find_package(autodiff)
find_package(OpenMP)

add_library(Particle src/particle.hpp src/particle.cpp)

add_library(PostProcess src/postprocess.hpp src/postprocess.cpp)
target_link_libraries(PostProcess ${OpenCV_LIBS})

add_library(Physics src/physics.hpp src/physics.cpp)
target_link_libraries(Physics Particle)

add_library(KSMetric src/metric_Kerr_Schild.hpp src/metric_Kerr_Schild.cpp src/metric.hpp)
target_link_libraries(KSMetric autodiff::autodiff)

add_library(KSMetricSchwarzschild src/metric_KS_Schwarzschild.hpp src/metric_KS_Schwarzschild.cpp src/metric.hpp)
target_link_libraries(KSMetricSchwarzschild autodiff::autodiff)

add_library(Draw src/draw.hpp src/draw.cpp)
target_link_libraries(Draw ${OpenCV_LIBS} Particle)

add_library(Raytracing src/raytracing.hpp src/raytracing.cpp)
target_link_libraries(Raytracing ${OpenCV_LIBS} Particle KSMetric KSMetricSchwarzschild Physics)
target_compile_options(Raytracing PRIVATE "${OpenMP_CXX_FLAGS}")

add_executable(example_schwarzschild src/example_schwarzschild.cpp)
target_link_libraries(example_schwarzschild ${OpenCV_LIBS} "${OpenMP_CXX_FLAGS}" Raytracing KSMetricSchwarzschild Particle Draw Physics PostProcess)
target_compile_options(example_schwarzschild PRIVATE "${OpenMP_CXX_FLAGS}")

add_executable(example_kerr src/example_kerr.cpp)
target_link_libraries(example_kerr ${OpenCV_LIBS} "${OpenMP_CXX_FLAGS}" Raytracing KSMetric Particle Draw Physics PostProcess)
target_compile_options(example_kerr PRIVATE "${OpenMP_CXX_FLAGS}")

add_executable(example_kerr_dust src/example_kerr_dust.cpp)
target_link_libraries(example_kerr_dust ${OpenCV_LIBS} "${OpenMP_CXX_FLAGS}" Raytracing KSMetric Particle Draw Physics PostProcess)
target_compile_options(example_kerr_dust PRIVATE "${OpenMP_CXX_FLAGS}")



