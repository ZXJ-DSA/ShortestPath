cmake_minimum_required(VERSION 3.14)
project(ShortestPath)

set(CMAKE_CXX_STANDARD 11)

# boost
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()
# openMP
find_package(OpenMP)
if (OPENMP_FOUND)
    message(STATUS "OpenMP find.")
    message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    message(STATUS "OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP_EXE_LINKER_FLAGS: ${OpenMP_EXE_LINKER_FLAGS}")
    message(STATUS "OPENMP_INCLUDE_DIRS: ${OPENMP_INCLUDE_DIRS}")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


add_definitions(-w) # omit warning information

add_executable(SP main.cpp Founda.hpp IndexFree.hpp)
add_executable(Preprocess mainPre.cpp Preprocess.hpp)

if (Boost_FOUND)
    if(UNIX AND NOT APPLE)
        target_link_libraries(SP pthread boost_system boost_thread)# For linux
        target_link_libraries(Preprocess pthread boost_system boost_thread)# For linux
    endif()
    if(APPLE)
        target_link_libraries(SP boost_system boost_thread-mt)# For mac
        target_link_libraries(Preprocess boost_system boost_thread-mt)# For mac
    endif()
endif ()
