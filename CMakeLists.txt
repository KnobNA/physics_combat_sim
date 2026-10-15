cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0091 NEW)

project(RobotCombatSim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
foreach(_config IN ITEMS Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER ${_config} _config_upper)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_config_upper} ${CMAKE_BINARY_DIR}/bin/${_config})
endforeach()

# ---------------------------------------------------------------------------
# Eigen (required by Project Chrono; header-only)
# ---------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_GetProperties(eigen)
if(NOT eigen_POPULATED)
    FetchContent_Populate(eigen)
endif()
set(EIGEN3_INCLUDE_DIR "${eigen_SOURCE_DIR}" CACHE PATH "Eigen headers" FORCE)

# ---------------------------------------------------------------------------
# Intel oneMKL (Pardiso is the FEA linear solver; MUMPS+AOCL is a later option)
# ---------------------------------------------------------------------------
function(rcs_locate_mkl_dir out_var)
    if(MKL_DIR AND EXISTS "${MKL_DIR}/MKLConfig.cmake")
        set(${out_var} "${MKL_DIR}" PARENT_SCOPE)
        return()
    endif()
    set(_hints)
    foreach(_env IN ITEMS MKL_DIR MKLROOT MKL_ROOT ONEAPI_ROOT)
        if(DEFINED ENV{${_env}} AND NOT "$ENV{${_env}}" STREQUAL "")
            list(APPEND _hints
                "$ENV{${_env}}"
                "$ENV{${_env}}/lib/cmake/mkl"
                "$ENV{${_env}}/mkl/latest/lib/cmake/mkl")
        endif()
    endforeach()
    foreach(_root IN ITEMS
            "C:/Program Files (x86)/Intel/oneAPI"
            "C:/Program Files/Intel/oneAPI"
            "/opt/intel/oneapi")
        list(APPEND _hints "${_root}/mkl/latest/lib/cmake/mkl")
        file(GLOB _ver_cfgs "${_root}/mkl/*/lib/cmake/mkl/MKLConfig.cmake")
        foreach(_cfg IN LISTS _ver_cfgs)
            get_filename_component(_dir "${_cfg}" DIRECTORY)
            list(APPEND _hints "${_dir}")
        endforeach()
    endforeach()
    foreach(_hint IN LISTS _hints)
        if(EXISTS "${_hint}/MKLConfig.cmake")
            set(${out_var} "${_hint}" PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${out_var} "" PARENT_SCOPE)
endfunction()

rcs_locate_mkl_dir(_rcs_mkl_dir)
if(NOT _rcs_mkl_dir)
    message(FATAL_ERROR
        "Intel oneMKL was not found. This project builds Chrono FEA against Pardiso from MKL.\n"
        "Install it, then reconfigure:\n"
        "  winget install Intel.oneMKL\n"
        "Or pass -DMKL_DIR=<oneAPI>/mkl/latest/lib/cmake/mkl")
endif()
set(MKL_DIR "${_rcs_mkl_dir}" CACHE PATH "Intel oneMKL CMake config directory" FORCE)
get_filename_component(MKL_ROOT "${MKL_DIR}/../../.." ABSOLUTE)
set(MKL_ROOT "${MKL_ROOT}" CACHE PATH "Intel oneMKL root" FORCE)
set(MKL_ARCH "intel64" CACHE STRING "MKL architecture")
set(MKL_INTERFACE "lp64" CACHE STRING "MKL interface")
set(MKL_LINK "dynamic" CACHE STRING "MKL linking mode")
# MSVC Chrono FEA uses vcomp OpenMP. Intel OpenMP from MKL's intel_thread layer
# cannot share that process. oneTBB is the parallel MKL path that coexists.
if(MSVC)
    set(MKL_THREADING "tbb_thread" CACHE STRING "MKL threading lib" FORCE)
else()
    set(MKL_THREADING "intel_thread" CACHE STRING "MKL threading lib")
endif()

# oneMKL's tbb_thread layer requires find_package(TBB CONFIG). The TBB that
# ships with standalone oneMKL often has libs but no TBBConfig.cmake.
set(_tbb_root_hints
    "${MKL_ROOT}/../../tbb/latest"
    "C:/Program Files (x86)/Intel/oneAPI/tbb/latest"
    "C:/Program Files/Intel/oneAPI/tbb/latest"
)
if(DEFINED ENV{TBBROOT})
    list(PREPEND _tbb_root_hints "$ENV{TBBROOT}")
endif()
set(_tbb_root "")
foreach(_hint IN LISTS _tbb_root_hints)
    if(EXISTS "${_hint}/lib/tbb12.lib" AND EXISTS "${_hint}/bin/tbb12.dll")
        get_filename_component(_tbb_root "${_hint}" ABSOLUTE)
        break()
    endif()
endforeach()
if(_tbb_root AND MKL_THREADING STREQUAL "tbb_thread")
    set(_tbb_cfg "${CMAKE_BINARY_DIR}/cmake/tbb")
    file(MAKE_DIRECTORY "${_tbb_cfg}")
    file(WRITE "${_tbb_cfg}/TBBConfig.cmake"
"if(NOT TARGET TBB::tbb)
  add_library(TBB::tbb SHARED IMPORTED)
  set_target_properties(TBB::tbb PROPERTIES
    IMPORTED_LOCATION \"${_tbb_root}/bin/tbb12.dll\"
    IMPORTED_LOCATION_RELEASE \"${_tbb_root}/bin/tbb12.dll\"
    IMPORTED_LOCATION_DEBUG \"${_tbb_root}/bin/tbb12_debug.dll\"
    IMPORTED_IMPLIB \"${_tbb_root}/lib/tbb12.lib\"
    IMPORTED_IMPLIB_RELEASE \"${_tbb_root}/lib/tbb12.lib\"
    IMPORTED_IMPLIB_DEBUG \"${_tbb_root}/lib/tbb12_debug.lib\")
endif()
set(TBB_FOUND TRUE)
")
    set(TBB_DIR "${_tbb_cfg}" CACHE PATH "oneTBB CMake config directory" FORCE)
    message(STATUS "oneTBB: ${_tbb_root}")
elseif(MKL_THREADING STREQUAL "tbb_thread")
    message(WARNING "oneTBB was not found; using sequential MKL so Pardiso does not load Intel OpenMP next to MSVC vcomp.")
    set(MKL_THREADING "sequential" CACHE STRING "MKL threading lib" FORCE)
endif()
find_package(MKL CONFIG REQUIRED)
if(NOT TARGET MKL::MKL)
    message(FATAL_ERROR "MKL::MKL target was not created by find_package(MKL).")
endif()
message(STATUS "Intel oneMKL config: ${MKL_DIR}")
message(STATUS "Intel oneMKL root:   ${MKL_ROOT}")

# ---------------------------------------------------------------------------
# Project Chrono (configured as its own CMake project, then find_package)
# ---------------------------------------------------------------------------
set(CHRONO_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/chrono")
set(CHRONO_BINARY_DIR "${CMAKE_BINARY_DIR}/chrono" CACHE PATH "Chrono build tree")

if(NOT Chrono_DIR)
    if(EXISTS "${CHRONO_BINARY_DIR}/cmake/ChronoConfig.cmake")
        set(Chrono_DIR "${CHRONO_BINARY_DIR}/cmake")
    endif()
endif()

message(STATUS "Configuring Project Chrono in ${CHRONO_BINARY_DIR}")
set(_chrono_cmake_args
    -S "${CHRONO_SOURCE_DIR}"
    -B "${CHRONO_BINARY_DIR}"
    -G "${CMAKE_GENERATOR}"
    -DCH_ENABLE_MODULE_FEA=ON
    -DCH_ENABLE_MODULE_PARDISO_MKL=ON
    -DCH_ENABLE_MODULE_MUMPS=OFF
    -DCH_ENABLE_MODULE_ROS=OFF
    -DCHRONO_GPU_VENDOR=NONE
    -DBUILD_DEMOS=OFF
    -DBUILD_TESTING=OFF
    -DBUILD_BENCHMARKING=OFF
    -DEIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR}
    -DEigen3_ROOT=${EIGEN3_INCLUDE_DIR}
    -DMKL_DIR=${MKL_DIR}
    -DMKL_ROOT=${MKL_ROOT}
    -DMKL_ARCH=intel64
    -DMKL_INTERFACE=lp64
    -DMKL_LINK=dynamic
    -DMKL_THREADING=${MKL_THREADING}
    -DCMAKE_PREFIX_PATH=${MKL_ROOT}
    -DCMAKE_CXX_FLAGS=/Zc:preprocessor
)
if(TBB_DIR)
    list(APPEND _chrono_cmake_args -DTBB_DIR=${TBB_DIR})
endif()
if(CMAKE_GENERATOR_PLATFORM)
    list(APPEND _chrono_cmake_args -A "${CMAKE_GENERATOR_PLATFORM}")
endif()
if(CMAKE_GENERATOR_TOOLSET)
    list(APPEND _chrono_cmake_args -T "${CMAKE_GENERATOR_TOOLSET}")
endif()
if(CMAKE_BUILD_TYPE)
    list(APPEND _chrono_cmake_args -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
endif()
if(CMAKE_CXX_COMPILER)
    list(APPEND _chrono_cmake_args -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER})
endif()
if(CMAKE_C_COMPILER)
    list(APPEND _chrono_cmake_args -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER})
endif()

execute_process(
    COMMAND ${CMAKE_COMMAND} ${_chrono_cmake_args}
    RESULT_VARIABLE _chrono_cfg_result
)
if(NOT _chrono_cfg_result EQUAL 0)
    message(FATAL_ERROR "Failed to configure Project Chrono. See output above.")
endif()
set(Chrono_DIR "${CHRONO_BINARY_DIR}/cmake")

# ChronoConfig.cmake requires the libraries to already exist on disk.
file(GLOB _chrono_built_libs
    "${CHRONO_BINARY_DIR}/lib/Release/*Chrono_pardisomkl*"
    "${CHRONO_BINARY_DIR}/bin/Release/*Chrono_pardisomkl*"
    "${CHRONO_BINARY_DIR}/Release/*Chrono_pardisomkl*"
    "${CHRONO_BINARY_DIR}/lib/*Chrono_pardisomkl*"
    "${CHRONO_BINARY_DIR}/bin/*Chrono_pardisomkl*"
)
if(NOT _chrono_built_libs)
    message(STATUS "Building Project Chrono with PardisoMKL (Release). This can take a while.")
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build "${CHRONO_BINARY_DIR}" --config Release
        RESULT_VARIABLE _chrono_build_result
    )
    if(NOT _chrono_build_result EQUAL 0)
        message(FATAL_ERROR "Failed to build Project Chrono. See output above.")
    endif()
endif()

# ChronoConfig.cmake emits an unquoted prefix path, which breaks on directories
# that contain spaces (this repo lives under "Robot Sim").
file(READ "${CHRONO_BINARY_DIR}/cmake/ChronoConfig.cmake" _chrono_config_text)
string(REGEX REPLACE
    "set\\(CHRONO_PACKAGE_PREFIX_DIR ([^\n\"]+)\\)"
    "set(CHRONO_PACKAGE_PREFIX_DIR \"\\1\")"
    _chrono_config_text "${_chrono_config_text}")
file(WRITE "${CHRONO_BINARY_DIR}/cmake/ChronoConfig.cmake" "${_chrono_config_text}")

find_package(Chrono REQUIRED COMPONENTS PardisoMKL CONFIG)

add_custom_target(chrono_build
    COMMAND ${CMAKE_COMMAND} --build "${CHRONO_BINARY_DIR}" --config $<CONFIG>
    COMMENT "Building Project Chrono"
    VERBATIM
)

# ---------------------------------------------------------------------------
# SDL2 (needed by Magnum Sdl2Application; fetched if not already available)
# ---------------------------------------------------------------------------
find_package(SDL2 CONFIG QUIET)
if(NOT TARGET SDL2::SDL2 AND NOT TARGET SDL2 AND NOT TARGET SDL2-static)
    include(FetchContent)
    FetchContent_Declare(SDL2
        GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
        GIT_TAG release-2.32.8
    )
    set(SDL_SHARED ON CACHE BOOL "" FORCE)
    set(SDL_STATIC OFF CACHE BOOL "" FORCE)
    set(SDL_TEST OFF CACHE BOOL "" FORCE)
    set(SDL_TESTS OFF CACHE BOOL "" FORCE)
    set(SDL2_DISABLE_INSTALL ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(SDL2)
endif()

# ---------------------------------------------------------------------------
# Magnum (designed for add_subdirectory)
# ---------------------------------------------------------------------------
set(CORRADE_BUILD_STATIC ON CACHE BOOL "" FORCE)
set(MAGNUM_BUILD_STATIC ON CACHE BOOL "" FORCE)
set(MAGNUM_BUILD_STATIC_PIC ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_SDL2APPLICATION ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_GL ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_SCENEGRAPH ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_PRIMITIVES ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_SHADERS ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_MESHTOOLS ON CACHE BOOL "" FORCE)
set(MAGNUM_WITH_TRADE ON CACHE BOOL "" FORCE)

add_subdirectory(external/corrade EXCLUDE_FROM_ALL)
add_subdirectory(external/magnum EXCLUDE_FROM_ALL)

set(IMGUI_DIR "${CMAKE_SOURCE_DIR}/external/imgui" CACHE PATH "Dear ImGui source" FORCE)
set(MAGNUM_WITH_IMGUIINTEGRATION ON CACHE BOOL "" FORCE)
set(MAGNUM_BUILD_STATIC ON CACHE BOOL "" FORCE)
add_subdirectory(external/magnum-integration EXCLUDE_FROM_ALL)

# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------
add_library(tetgen STATIC
    external/tetgen/tetgen.cxx
    external/tetgen/predicates.cxx)
target_include_directories(tetgen PUBLIC external/tetgen)
target_compile_definitions(tetgen PUBLIC TETLIBRARY)
if(MSVC)
    target_compile_options(tetgen PRIVATE /w /bigobj)
endif()

set(COMBAT_SIM_PHYSICS_SOURCES
    src/physics_core/part.cpp
    src/physics_core/heat.cpp
    src/physics_core/reissner_j2.cpp
    src/physics_core/thermal_tet.cpp
    src/physics_core/contact_load.cpp
    src/physics_core/mesh_split.cpp
    src/physics_core/assembly.cpp
    src/physics_core/cad.cpp
    src/physics_core/simulator.cpp
    src/physics_core/physics_units.cpp
)

add_executable(combat_sim
    src/main.cpp
    ${COMBAT_SIM_PHYSICS_SOURCES}
    src/render_core/orbit_camera.cpp
    src/render_core/renderer.cpp
    src/render_core/gizmos.cpp
    src/editor/ui.cpp
    src/editor/faces.cpp
)

target_include_directories(combat_sim PRIVATE src)

target_link_libraries(combat_sim PRIVATE
    tetgen
    ${CHRONO_TARGETS}
    MKL::MKL
    Magnum::Magnum
    Magnum::Sdl2Application
    Magnum::GL
    Magnum::Primitives
    Magnum::Shaders
    Magnum::MeshTools
    Magnum::Trade
    MagnumIntegration::ImGui
)

add_dependencies(combat_sim chrono_build)

add_executable(combat_sim_physics_smoke EXCLUDE_FROM_ALL
    src/physics_core/physics_smoke.cpp
    ${COMBAT_SIM_PHYSICS_SOURCES}
)
target_include_directories(combat_sim_physics_smoke PRIVATE src)
target_link_libraries(combat_sim_physics_smoke PRIVATE tetgen ${CHRONO_TARGETS} MKL::MKL)
add_dependencies(combat_sim_physics_smoke chrono_build)
set_target_properties(combat_sim_physics_smoke PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}
)

# ---------------------------------------------------------------------------
# GoogleTest unit suite (known I/O per physics calculation; no Simulator::step)
# ---------------------------------------------------------------------------
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)
enable_testing()
include(GoogleTest)

add_executable(combat_sim_unit_tests
    tests/unit/test_loads.cpp
    tests/unit/test_materials.cpp
    tests/unit/test_constitutive.cpp
    tests/unit/test_fracture.cpp
    tests/unit/test_welds.cpp
    tests/unit/test_mesh_geom.cpp
    tests/unit/test_contact_gains.cpp
    tests/unit/test_recording.cpp
    tests/unit/test_joints.cpp
    tests/unit/test_idealization_policy.cpp
    tests/unit/test_constitutive_compat.cpp
    ${COMBAT_SIM_PHYSICS_SOURCES}
)
target_include_directories(combat_sim_unit_tests PRIVATE src)
target_link_libraries(combat_sim_unit_tests PRIVATE
    GTest::gtest_main
    tetgen
    ${CHRONO_TARGETS}
    MKL::MKL
)
add_dependencies(combat_sim_unit_tests chrono_build)
set_target_properties(combat_sim_unit_tests PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE}
    RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG}
)
gtest_discover_tests(combat_sim_unit_tests
    DISCOVERY_MODE PRE_TEST
    WORKING_DIRECTORY $<TARGET_FILE_DIR:combat_sim_unit_tests>
)

if(DEFINED CHRONO_DATA_DIR)
    target_compile_definitions(combat_sim PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    target_compile_definitions(combat_sim_physics_smoke PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
    target_compile_definitions(combat_sim_unit_tests PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"")
endif()

if(MSVC)
    target_compile_options(combat_sim PRIVATE /bigobj /utf-8)
    target_compile_options(combat_sim_physics_smoke PRIVATE /bigobj /utf-8)
    target_compile_options(combat_sim_unit_tests PRIVATE /bigobj /utf-8)
    if(DEFINED CHRONO_MSVC_RUNTIME_LIBRARY)
        set_target_properties(combat_sim PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
        set_target_properties(combat_sim_physics_smoke PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
        set_target_properties(combat_sim_unit_tests PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
    endif()
    set_property(TARGET combat_sim PROPERTY VS_DEBUGGER_ENVIRONMENT
        "MKL_INTERFACE_LAYER=LP64\nMKL_THREADING_LAYER=TBB")
    set_property(TARGET combat_sim_physics_smoke PROPERTY VS_DEBUGGER_ENVIRONMENT
        "MKL_INTERFACE_LAYER=LP64\nMKL_THREADING_LAYER=TBB")
    set_property(TARGET combat_sim_unit_tests PROPERTY VS_DEBUGGER_ENVIRONMENT
        "MKL_INTERFACE_LAYER=LP64\nMKL_THREADING_LAYER=TBB")
endif()

if(COMMAND add_CHRONO_DLLS_copy_command)
    add_CHRONO_DLLS_copy_command()
endif()

# Always copy Chrono (and SDL2 if shared) next to the executable.
foreach(_exe IN ITEMS combat_sim combat_sim_unit_tests)
    foreach(_chrono_tgt IN ITEMS Chrono::Chrono_core Chrono::Chrono_pardisomkl Chrono::pardisomkl)
        if(TARGET ${_chrono_tgt})
            add_custom_command(TARGET ${_exe} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    $<TARGET_FILE:${_chrono_tgt}>
                    $<TARGET_FILE_DIR:${_exe}>
                COMMENT "Copying ${_chrono_tgt} next to ${_exe}"
            )
        endif()
    endforeach()
endforeach()

# Pardiso loads MKL kernels (and Intel OpenMP) by DLL name at runtime.
set(_mkl_dll_dirs
    "${MKL_ROOT}/bin"
    "${MKL_ROOT}/redist/intel64"
    "${MKL_ROOT}/../../compiler/latest/bin"
    "${MKL_ROOT}/../../compiler/latest/windows/redist/intel64_win/compiler"
    "${MKL_ROOT}/../../tbb/latest/bin"
    "${MKL_ROOT}/../../tbb/latest/bin/vc_mt"
    "${MKL_ROOT}/../compiler/latest/bin"
)
set(_mkl_runtime_dlls)
foreach(_dir IN LISTS _mkl_dll_dirs)
    get_filename_component(_dir "${_dir}" ABSOLUTE)
    if(EXISTS "${_dir}")
        file(GLOB _found
            "${_dir}/mkl_core*.dll"
            "${_dir}/mkl_intel_thread*.dll"
            "${_dir}/mkl_tbb_thread*.dll"
            "${_dir}/mkl_sequential*.dll"
            "${_dir}/mkl_rt*.dll"
            "${_dir}/mkl_avx*.dll"
            "${_dir}/mkl_def*.dll"
            "${_dir}/mkl_mc*.dll"
            "${_dir}/mkl_vml_*.dll"
            "${_dir}/libiomp5md.dll"
            "${_dir}/libimalloc.dll"
            "${_dir}/tbb12.dll"
            "${_dir}/tbbmalloc.dll"
            "${_dir}/tbbbind.dll"
            "${_dir}/tbbbind_2_0.dll"
            "${_dir}/tbbbind_2_5.dll"
        )
        list(APPEND _mkl_runtime_dlls ${_found})
    endif()
endforeach()
if(_mkl_runtime_dlls)
    list(REMOVE_DUPLICATES _mkl_runtime_dlls)
    foreach(_exe IN ITEMS combat_sim combat_sim_unit_tests)
        foreach(_dll IN LISTS _mkl_runtime_dlls)
            add_custom_command(TARGET ${_exe} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "${_dll}"
                    $<TARGET_FILE_DIR:${_exe}>
                COMMENT "Copying MKL runtime DLL next to ${_exe}"
            )
        endforeach()
    endforeach()
endif()

if(EXISTS "${CHRONO_BINARY_DIR}/bin/Release/yaml-cpp.dll")
    foreach(_exe IN ITEMS combat_sim combat_sim_unit_tests)
        add_custom_command(TARGET ${_exe} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CHRONO_BINARY_DIR}/bin/$<CONFIG>/yaml-cpp.dll"
                $<TARGET_FILE_DIR:${_exe}>
            COMMENT "Copying yaml-cpp next to ${_exe}"
        )
    endforeach()
endif()

if(TARGET SDL2::SDL2)
    get_target_property(_sdl2_type SDL2::SDL2 TYPE)
    if(_sdl2_type STREQUAL "SHARED_LIBRARY")
        add_custom_command(TARGET combat_sim POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:SDL2::SDL2>
                $<TARGET_FILE_DIR:combat_sim>
            COMMENT "Copying SDL2 next to combat_sim"
        )
    endif()
endif()
