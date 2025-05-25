# CMakeList.txt: StructureLightAlgorithms 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required(VERSION 3.8)

project(StructureLightAlgorithm)

#查找OpenMP包，它是一个内置于编译器的并行编程标准
find_package(OpenMP REQUIRED)

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(include)

message("status: ${CMAKE_CURRENT_SOURCE_DIR}")


if(CMAKE_COMPILER_IS_GNUCXX)
    # 使用 -O0 选项禁用优化
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -fopenmp")
elseif(MSVC)
    # 使用 /Od 选项禁用优化
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Od /openmp")
endif()



#因为 GLOB 的变化不会触发 CMake 重新生成构建文件。也就是说：新增文件后，可能需要你手动重新 cmake 一次。
file(GLOB_RECURSE SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_executable(StructureLightAlgorithm)

#使用C++17标准
target_compile_features(StructureLightAlgorithm PRIVATE cxx_std_17)


target_sources(StructureLightAlgorithm PRIVATE 
                ${SOURCE} 
                "include/CloudPointsGeneration.h" 
                "CloudPointsGeneration.cpp" 
                "include/MulitExposure.h" 
                "MulitExposure.cpp" "monocularStructureLight.cpp" "include/monocularStructureLight.h" "include/ImagepreProcessing.h" "ImagepreProcessing.cpp" "include/ImageReprocess.h" "ImageReprocess.cpp")



target_link_libraries(${PROJECT_NAME} 
                      ${OpenCV_LIBRARIES}
                      OpenMP::OpenMP_CXX)
