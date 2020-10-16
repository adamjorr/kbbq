cmake_minimum_required(VERSION 3.10)

# set the project name
project(kbbq)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# library stuff
add_subdirectory(src/minionrng)
#add_library(hts SHARED IMPORTED)
add_subdirectory(src/kbbq)

# add the executable
add_executable(kbbq-bin src/kbbq/kbbq.cc)

# force native build to enable vectorization
add_compile_options(-march=native)

# add libraries for linking
target_link_libraries(kbbq-bin kbbq)

target_include_directories(kbbq-bin PUBLIC
                          "${PROJECT_SOURCE_DIR}/include/kbbq"
                          )

set_target_properties(kbbq-bin PROPERTIES OUTPUT_NAME kbbq)

install(TARGETS kbbq-bin)
