# CMakeList.txt: GrayCodeGenerate 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})


# 将源代码添加到此项目的可执行文件。
add_executable (GrayCodeGenerate "GrayCodeGenerate.cpp" "GrayCodeGenerate.h" "main.cpp")

target_link_libraries(
	GrayCodeGenerate
	${OpenCV_LIBRARIES}
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET GrayCodeGenerate PROPERTY CXX_STANDARD 20)
endif()

# TODO: 如有需要，请添加测试并安装目标。
