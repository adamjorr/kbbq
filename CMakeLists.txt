cmake_minimum_required(VERSION 3.10)

# set the project name
project(cbbq)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# library stuff
add_subdirectory(src/minionrng)
#add_library(hts SHARED IMPORTED)
add_subdirectory(src/yak)
add_subdirectory(src/kbbq)

# add the executable
add_executable(cbbq src/kbbq/cbbq.cc)

# add libraries for linking
target_link_libraries(cbbq kbbq)

target_include_directories(cbbq PUBLIC
                          "${PROJECT_SOURCE_DIR}/include/kbbq"
                          )
