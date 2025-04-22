if(ADVANCED_OUTPUT)
    add_compile_options(-DADVANCED_OUTPUT)
endif()

add_executable(get_espe gauss.c ../ran_pcg.c get_espe.c)

target_include_directories(get_espe PRIVATE
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/pcg/include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>)

target_link_libraries(get_espe
        PRIVATE jibal
        PRIVATE pcg
        "$<$<BOOL:${UNIX}>:m>")
