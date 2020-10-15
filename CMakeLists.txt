cmake_minimum_required(VERSION 3.13)
if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
  message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

project(mcerd
        VERSION 1.0.3
        DESCRIPTION "Monte Carlo ERD"
        LANGUAGES C
        )

set(CMAKE_REQUIRED_LIBRARIES "m")

set(CMAKE_C_STANDARD 99)

find_package(Jibal 0.2.6 REQUIRED)
if(Jibal_FOUND)
    message(STATUS "JIBAL ${Jibal_VERSION} headers found at ${Jibal_INCLUDE_DIR}")
endif()

option(DEVELOPER_MODE_ENABLE "Enable developer mode" ON)
if(DEVELOPER_MODE_ENABLE)
    set(MCERD_DATADIR "${PROJECT_SOURCE_DIR}/data")
else()
    set(MCERD_DATADIR "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}/")
endif()

configure_file(mcerd_config.h.in mcerd_config.h @ONLY)

add_executable(mcerd 
    main.c cross_section.c read_input.c scattering_angle.c
    init_params.c misc.c potential.c init_simu.c random.c ran_pcg.c
    read_target.c read_detector.c init_detector.c print_data.c
    ion_simu.c finalize.c finish_ion.c output.c erd_scattering.c
    rotate.c pre_simulation.c erd_detector.c elsto.c
    debug.c virtual_detector.c rough_surface.c ion_stack.c)

add_subdirectory(pcg)

target_include_directories(mcerd PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/pcg/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>) #Because of mcerd_config.h

if(APPLE)
    set_target_properties(mcerd PROPERTIES LINK_FLAGS -Wl,-stack_size,0x100000) #Hexadecimal 8 MB
elseif(MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608") #8 MB
endif()

target_link_libraries(mcerd PRIVATE Jibal PRIVATE pcg
    "$<$<BOOL:${UNIX}>:m>"
    )

add_subdirectory(get_espe)

INSTALL(TARGETS mcerd get_espe
        RUNTIME DESTINATION bin)
