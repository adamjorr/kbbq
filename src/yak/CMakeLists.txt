add_library(yak SHARED
    kthread.c
    bbf.c
    htab.c
    bseq.c
    misc.c
    sys.c
    6gjdn.c
    count.c
    qv.c
    triobin.c
    trioeval.c
    inspect.c
)

target_include_directories(yak PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/../../include/yak"
)

#add_library(m SHARED IMPORTED)
#add_library(z SHARED IMPORTED)
#add_library(pthread SHARED IMPORTED)

target_link_libraries(yak PUBLIC
    m
    z
    pthread
)
