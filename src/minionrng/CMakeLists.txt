add_library(minionrng SHARED
    minion.cc
)

target_include_directories(minionrng PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/minionrng"
)
