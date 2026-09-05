cmake_minimum_required(VERSION 4.1)
project(Gmap_fruits)

set(CMAKE_CXX_STANDARD 20)

find_package(CGAL REQUIRED COMPONENTS Qt6)

add_executable(Gmap_fruits main.cpp
        tomato/grammar_1.h)

target_link_libraries(Gmap_fruits PRIVATE CGAL::CGAL CGAL::CGAL_Basic_viewer Qt6::Widgets)