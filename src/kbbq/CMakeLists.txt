add_library(kbbq STATIC
    bloom.cc
    covariateutils.cc
    recalibrateutils.cc
    readutils.cc
    htsiter.cc
)

target_include_directories(kbbq PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/kbbq"
)

target_link_libraries(kbbq PUBLIC
    hts
    minionrng
)
