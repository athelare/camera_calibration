cmake_minimum_required(VERSION 3.15)
project(camera_calibration)
MESSAGE(STATUS "operation system is ${CMAKE_SYSTEM}")
set(CMAKE_CXX_STANDARD 17)

IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
    add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
    add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")
    set(OpenCV_DIR D:/opencv/build)
ENDIF (CMAKE_SYSTEM_NAME MATCHES "Windows")

FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

add_executable(camera_calibration main.cpp DvpFrameCapture.h include/DVPCamera.h)

TARGET_LINK_LIBRARIES(camera_calibration  ${OpenCV_LIBS} ${CMAKE_CURRENT_SOURCE_DIR}/lib/DVPCamera64.lib)