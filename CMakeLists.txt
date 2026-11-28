cmake_minimum_required(VERSION 3.28)
project(metriko)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_OPENMP "Enable OpenMP support" ON)
option(METRIKO_BUILD_DIST "Build the distributable targets (metriko_core, metriko-cli, python extension)" OFF)
option(METRIKO_BUILD_CLI "Build metriko-cli (needs METRIKO_BUILD_DIST)" ON)
option(METRIKO_BUILD_PYTHON "Build the nanobind extension metriko._core (needs METRIKO_BUILD_DIST)" ON)
option(METRIKO_REQUIRE_SUITESPARSE "Fail when CHOLMOD is not found instead of falling back to Eigen's solvers" ${METRIKO_BUILD_DIST})

#========= libigl (also provides Eigen) ================================================================================
if (NOT TARGET igl::core)
    include(FetchContent)
    FetchContent_Declare(
            libigl
            GIT_REPOSITORY https://github.com/libigl/libigl.git
            GIT_TAG 7100764c2a2833284a8bdfa9b948d2b3bf124624 # v2.6.4
    )
    FetchContent_MakeAvailable(libigl)
endif()

# Work with non-standard homebrew installations
if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    find_program(HOMEBREW_EXECUTABLE brew)
    mark_as_advanced(FORCE HOMEBREW_EXECUTABLE)
    if (HOMEBREW_EXECUTABLE)
        # Detected a Homebrew install, query for its install prefix.
        execute_process(COMMAND ${HOMEBREW_EXECUTABLE} --prefix
                OUTPUT_VARIABLE HOMEBREW_INSTALL_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Detected Homebrew with install prefix: "
                "${HOMEBREW_INSTALL_PREFIX}, adding to CMake search paths.")
        list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_INSTALL_PREFIX}")
    endif()
endif()

#========= SuiteSparse (CHOLMOD, via the CMake config shipped with SuiteSparse >= 7) ===================================
set(SUITESPARSE_HELP "Enable SuiteSparse (CHOLMOD). Falls back to Eigen's solvers when off or not found.")
option(SUITESPARSE "${SUITESPARSE_HELP}" ON)
set(METRIKO_HAVE_SUITESPARSE FALSE)
if (SUITESPARSE)
    if (METRIKO_REQUIRE_SUITESPARSE)
        find_package(CHOLMOD CONFIG REQUIRED)
    else()
        find_package(CHOLMOD CONFIG QUIET)
    endif()
    if (CHOLMOD_FOUND)
        message(STATUS "Found SuiteSparse (CHOLMOD ${CHOLMOD_VERSION}), building with SuiteSparse.")
        set(METRIKO_HAVE_SUITESPARSE TRUE)
        if (TARGET SuiteSparse::CHOLMOD_static)
            set(METRIKO_CHOLMOD_TARGET SuiteSparse::CHOLMOD_static)
        else()
            set(METRIKO_CHOLMOD_TARGET SuiteSparse::CHOLMOD)
        endif()
        if (APPLE)
            set(BLA_VENDOR Apple)
        endif()
        find_package(OpenBLAS CONFIG QUIET)
        if (TARGET OpenBLAS::OpenBLAS)
            set(METRIKO_BLAS_TARGETS OpenBLAS::OpenBLAS)
        else()
            find_package(BLAS REQUIRED)
            find_package(LAPACK REQUIRED)
            set(METRIKO_BLAS_TARGETS LAPACK::LAPACK BLAS::BLAS)
        endif()
        message(STATUS "CHOLMOD target: ${METRIKO_CHOLMOD_TARGET}, BLAS/LAPACK: ${METRIKO_BLAS_TARGETS}")
    else()
        message(STATUS "SuiteSparse (CHOLMOD) not found, falling back to Eigen's solvers.")
        set(SUITESPARSE OFF CACHE BOOL "${SUITESPARSE_HELP}" FORCE)
    endif()
else()
    message(STATUS "Building without SuiteSparse, using Eigen's solvers.")
endif()

#========= Try finding OpenMp ==========================================================================================

set(METRIKO_HAVE_OPENMP FALSE)
if (ENABLE_OPENMP OR METRIKO_BUILD_DIST)
    if (METRIKO_BUILD_DIST)
        find_package(OpenMP REQUIRED COMPONENTS CXX)
    else()
        find_package(OpenMP)
    endif()
    if (OpenMP_CXX_FOUND)
        set(METRIKO_HAVE_OPENMP TRUE)
        message(STATUS "OpenMP found. CXX flags: ${OpenMP_CXX_FLAGS}")
        # Apple clang: FindOpenMP finds the flags (-Xclang -fopenmp) but not the runtime library, and linking
        # without it leaves __kmpc_* undefined in the extension module and the CLI. The caller may pass the
        # absolute path (conan's llvm-openmp, e.g. -DMETRIKO_OPENMP_LIBRARY=.../lib/libomp.a); otherwise look for
        # it on the usual search paths. find_library alone is not enough: a conan toolchain overwrites
        # CMAKE_PREFIX_PATH, so the package folder is not on the search path.
        if (APPLE AND NOT METRIKO_OPENMP_LIBRARY)
            find_library(METRIKO_OPENMP_LIBRARY omp)
        endif()
        if (APPLE AND METRIKO_OPENMP_LIBRARY)
            message(STATUS "OpenMP runtime: ${METRIKO_OPENMP_LIBRARY}")
        elseif (APPLE)
            message(WARNING "OpenMP flags found but no libomp: the link will fail on __kmpc_* symbols")
        endif()
    else()
        message(WARNING "OpenMP not found, building without OpenMP support.")
    endif()
endif()
#========== Done finding OpenMp ========================================================================================

add_library(${PROJECT_NAME} INTERFACE)

target_link_libraries(${PROJECT_NAME} INTERFACE igl::core)
target_include_directories(${PROJECT_NAME} INTERFACE include)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    # Eigen parallelizes some sparse products with its own OpenMP temporaries, whose accumulation order (and so
    # the last bits) depends on the scheduling: a run-to-run difference there cascades into a different
    # extraction. the explicit OpenMP loops stay parallel; Eigen's internal threading is turned off.
    target_compile_definitions(${PROJECT_NAME} INTERFACE NOMINMAX _USE_MATH_DEFINES EIGEN_DONT_PARALLELIZE)
    target_compile_options(${PROJECT_NAME} INTERFACE /utf-8 /bigobj /EHsc /Zc:preprocessor)
endif()

if(METRIKO_HAVE_SUITESPARSE)
    target_compile_definitions(${PROJECT_NAME} INTERFACE METRIKO_HAVE_SUITESPARSE GC_HAVE_SUITESPARSE)
    target_link_libraries(${PROJECT_NAME} INTERFACE ${METRIKO_CHOLMOD_TARGET} ${METRIKO_BLAS_TARGETS})
endif()

if(METRIKO_HAVE_OPENMP)
    target_link_libraries(${PROJECT_NAME} INTERFACE OpenMP::OpenMP_CXX)
    if (APPLE AND METRIKO_OPENMP_LIBRARY)
        # absolute path: propagates to _core and metriko-cli without -L/-l juggling
        target_link_libraries(${PROJECT_NAME} INTERFACE ${METRIKO_OPENMP_LIBRARY})
    endif()
endif()

#========= Distributable targets (static core shared by the CLI and the python extension) ==============================
if(METRIKO_BUILD_DIST)
    #--------- pmp-library (vendored in deps/), static, 64-bit scalars -----------------------------------------------
    set(PMP_BUILD_EXAMPLES OFF)
    set(PMP_BUILD_TESTS    OFF)
    set(PMP_BUILD_DOCS     OFF)
    set(PMP_BUILD_VIS      OFF)
    set(PMP_INSTALL        OFF)
    set(BUILD_SHARED_LIBS  OFF)
    set(PMP_SCALAR_TYPE    64)
    if(NOT EIGEN3_INCLUDE_DIR)
        FetchContent_GetProperties(eigen SOURCE_DIR EIGEN3_INCLUDE_DIR)
        if(NOT EIGEN3_INCLUDE_DIR AND TARGET Eigen3::Eigen)
            get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
        endif()
    endif()
    add_subdirectory(deps/pmp-library)
    target_include_directories(pmp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/deps/pmp-library/src)

    add_library(metriko_core STATIC src/metriko_api.cpp src/decimate.cpp)
    target_compile_definitions(metriko_core PRIVATE PMP_SCALAR_TYPE_64)
    set_target_properties(metriko_core PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_include_directories(metriko_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(metriko_core PUBLIC ${PROJECT_NAME} PRIVATE pmp)

    if(METRIKO_BUILD_CLI)
        find_package(CLI11 CONFIG REQUIRED)
        add_executable(metriko-cli src/cli/main.cpp)
        target_link_libraries(metriko-cli PRIVATE metriko_core CLI11::CLI11)
        install(TARGETS metriko-cli RUNTIME DESTINATION bin COMPONENT cli)

        enable_testing()
        add_test(NAME cli.fertility
                 COMMAND metriko-cli -i ${CMAKE_CURRENT_SOURCE_DIR}/example/models/fertility.obj
                                     -o ${CMAKE_CURRENT_BINARY_DIR}/fertility_quad.obj --scale 0.01 --quiet)
    endif()

    if(METRIKO_BUILD_PYTHON)
        find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module Development.SABIModule)
        find_package(nanobind CONFIG REQUIRED)
        nanobind_add_module(_core STABLE_ABI NB_STATIC python/bindings.cpp)
        target_link_libraries(_core PRIVATE metriko_core)
        install(TARGETS _core LIBRARY DESTINATION metriko COMPONENT python)
    endif()
endif()
