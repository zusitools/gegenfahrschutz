cmake_minimum_required (VERSION 3.1)
project (gegenfahrschutz)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(parser)

find_package(Boost COMPONENTS filesystem program_options REQUIRED)

add_executable(gegenfahrschutz gegenfahrschutz.cpp nowide/src/iostream.cpp)
# TODO: target_link_libraries(Boost::filesystem) funktioniert derzeit unter MXE nicht
target_include_directories(gegenfahrschutz PRIVATE ${Boost_INCLUDE_DIRS} nowide/include)
target_link_libraries(gegenfahrschutz PRIVATE zusi_parser ${Boost_LIBRARIES})

install(TARGETS gegenfahrschutz RUNTIME DESTINATION bin)
