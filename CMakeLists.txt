cmake_minimum_required(VERSION 3.14)
project(ShortestPath)

set(CMAKE_CXX_STANDARD 17)

# boost
find_package(Boost 1.71.0 COMPONENTS system filesystem REQUIRED)
if (Boost_FOUND)
    message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Boost_VERSION: ${Boost_VERSION}")
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif ()

add_definitions(-w) # omit warning information

add_executable(SP main.cpp Founda.hpp IndexFree.hpp IndexBased.hpp)
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
