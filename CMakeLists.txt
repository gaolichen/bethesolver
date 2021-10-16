cmake_minimum_required (VERSION 3.1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#SET(MYDEV D:/dev)
SET(MYDEV /home/gaolichen/mydev)
SET(BOOST_ROOT ${MYDEV}/boost_1_70_0)

SET(EIGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../eigen/)
find_package(Boost 1.70)

project (bethesolver)
SET (CMAKE_BUILD_TYPE Release)
# to use gdb enable the following two flags.
#set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

add_subdirectory (src)
