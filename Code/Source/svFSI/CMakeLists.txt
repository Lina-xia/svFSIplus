set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${SV_SOURCE_DIR}/ThirdParty/parmetis_svfsi/simvascular_parmetis_svfsi/ParMETISLib)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tetgen/simvascular_tetgen)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tinyxml/simvascular_tinyxml)
include_directories(${MPI_C_INCLUDE_PATH})

# Find Trilinos package if requested
if(SV_USE_TRILINOS)

  # source: https://github.com/trilinos/Trilinos_tutorial/wiki/CMakeFindPackageTrilinosExample

  # It may be required to add a path for Trilinos
  # This can be accomplished by adding the Trilinos
  # -DCMAKE_PREFIX_PATH=~/sfw/Trilinos_build/lib/cmake/Trilinos
  FIND_PACKAGE(Trilinos)

  # If FIND_PACKAGE successfully found your Trilinos install, it will
  # set the Boolean flag Trilinos_FOUND.  The following IF statement
  # fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
  # found, it prints out the values of some Trilinos configuration
  # details.  You may find them useful for building your application
  # that uses Trilinos.
  IF(Trilinos_FOUND)
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
    MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
    MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
    MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
    MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
    MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
    MESSAGE("End of Trilinos details\n")

    if(NOT Trilinos_LIBRARIES)
      MESSAGE(WARNING "Trilinos_FOUND is true but Trilinos_LIBRARIES is empty. This suggests that FIND_PACKAGE(Trilinos) failed quietly. Proceed with caution.")
    endif()

    set(USE_TRILINOS 1)
  ELSE()
    MESSAGE(WARNING "Could not find Trilinos. Compiling svFSI without Trilinos.")
  ENDIF()

endif()

# add trilinos flags and defines
if(USE_TRILINOS)
  ADD_DEFINITIONS(-DWITH_TRILINOS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Build with the PETSc linear algebra package.
#
if(NOT "${SV_PETSC_DIR}" STREQUAL "")

   if(NOT EXISTS ${SV_PETSC_DIR})
      message(ERROR "The PETSc directory ${SV_PETSC_DIR} could not be found.")
   endif()

   if(NOT EXISTS "${SV_PETSC_DIR}/include")
      message(ERROR "The PETSc include directory ${SV_PETSC_DIR}/include could not be found.")
   endif()

   if(NOT EXISTS "${SV_PETSC_DIR}/lib")
      message(ERROR "The PETSc library directory ${SV_PETSC_DIR}/lib could not be found.")
   endif()

  set(PETSC_INCLUDE_DIRS "${SV_PETSC_DIR}/include;${SV_PETSC_DIR}/../include")
  set(PETSC_LIBRARY_DIRS "-L${SV_PETSC_DIR}/lib -lpetsc")
  message("\nBuilding with the PETSc package; the include and libraries directories are: ")
  message("  PETSC_LIBRARY_DIRS = ${PETSC_LIBRARY_DIRS}")
  message("  PETSC_INCLUDE_DIRS = ${PETSC_INCLUDE_DIRS}")

  # Set PETSc include directory.
  include_directories(${PETSC_INCLUDE_DIRS})

  # Set C++ directive to use PETSc in svFSIplus code.
  ADD_DEFINITIONS(-DWITH_PETSC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic -std=c99")

  set(USE_PETSC 1)

else()

  set(USE_PETSC 0)

endif()

if(ENABLE_ARRAY_INDEX_CHECKING)
  ADD_DEFINITIONS(-DENABLE_ARRAY_INDEX_CHECKING)
endif()

# svFSI requires LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Include VTK either from a local build using SV_LOCAL_VTK_PATH
# or from a default installed version.
#
if("${SV_LOCAL_VTK_PATH}" STREQUAL "")
  find_package(VTK REQUIRED ) 
else()
  find_package(VTK PATHS ${SV_LOCAL_VTK_PATH} NO_DEFAULT_PATH REQUIRED) 
endif()

include(${VTK_USE_FILE})

set(lib ${SV_LIB_SVFSI_NAME})

set(CSRCS 
  Array3.h Array3.cpp 
  Array.h Array.cpp
  Tensor4.h Tensor4.cpp
  Vector.h Vector.cpp 

  lapack_defs.h

  DebugMsg.h 
  Parameters.h Parameters.cpp
  Simulation.h Simulation.cpp
  SimulationLogger.h
  VtkData.h VtkData.cpp

  all_fun.h all_fun.cpp
  baf_ini.h baf_ini.cpp
  bf.h bf.cpp
  cep.h cep.cpp
  cep_ion.h cep_ion.cpp
  cmm.h cmm.cpp
  consts.h consts.cpp
  contact.h contact.cpp
  distribute.h distribute.cpp
  eq_assem.h eq_assem.cpp
  fluid.h fluid.cpp
  fsi.h fsi.cpp
  fs.h fs.cpp
  fft.h fft.cpp
  heatf.h heatf.cpp
  heats.h heats.cpp
  initialize.h initialize.cpp
  l_elas.h l_elas.cpp
  lhsa.h lhsa.cpp
  ls.h ls.cpp
  main.cpp
  mat_fun.h mat_fun.cpp
  mat_fun_carray.h mat_fun_carray.cpp
  mat_models.h mat_models.cpp
  mesh.h mesh.cpp
  nn.h nn.cpp
  output.h output.cpp
  load_msh.h load_msh.cpp
  pic.h pic.cpp
  post.h post.cpp
  read_files.h read_files.cpp
  read_msh.h read_msh.cpp
  remesh.h remesh.cpp
  remeshTet.cpp
  set_bc.h set_bc.cpp
  shells.h shells.cpp
  stokes.h stokes.cpp
  sv_struct.h sv_struct.cpp
  txt.h txt.cpp
  utils.h utils.cpp
  ustruct.h ustruct.cpp
  vtk_xml.h vtk_xml.cpp
  vtk_xml_parser.h vtk_xml_parser.cpp

  CepMod.h CepMod.cpp
  CepModAp.h CepModAp.cpp
  CepModBo.h CepModBo.cpp
  CepModFn.h CepModFn.cpp
  CepModTtp.h CepModTtp.cpp
  ChnlMod.h ChnlMod.cpp
  CmMod.h CmMod.cpp
  ComMod.h ComMod.cpp
  Timer.h

  SPLIT.c
  
  cpl1DType.h cpl1DType.cpp
)

  # Set PETSc interace code.
  if(USE_PETSC)
     #set(CSRCS ${CSRCS} petsc_linear_solver.c)
  endif()

if(USE_TRILINOS)
  set(CSRCS ${CSRCS} trilinos_linear_solver.cpp)

  # trilinos directories and libraries
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  link_directories(${Trilinos_TPL_LIBRARY_DIRS})
endif()


# xia
  FILE(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/../../cvOneD/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/*.cxx")
  FILE(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/*.h")

  set(SRC_C ${SRC_C} CACHE INTERNAL "List of cvOneD source files")
  set(SRC_H ${SRC_H} CACHE INTERNAL "List of cvOneD head files")

  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD)
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD)
  # message("CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

IF(sparseSolverType STREQUAL "skyline")
  simvascular_add_executable(${SV_SVFSI_EXE}
      SRCS ${CSRCS} ${SRC_C} ${SRC_H}
      DEV_SCRIPT_NAME "mysvfsi"
      INSTALL_SCRIPT_NAME "svfsi"
      INSTALL_COMP CoreExecutables
      INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})

  target_compile_definitions(${SV_SVFSI_EXE} PRIVATE USE_SKYLINE)

ELSEIF(sparseSolverType STREQUAL "superlu")

  # FIND BLAS
  FIND_PACKAGE(BLAS REQUIRED)
  IF(BLAS_FOUND)
    INCLUDE_DIRECTORIES(${BLAS_INCLUDE_DIRS})
  ENDIF()

  # CHECK IF SUPERLU_DIR IS STILL UNDEFINED
  IF(${SUPERLU_DIR} STREQUAL " ")
    MESSAGE(FATAL_ERROR "The SUPERLU_DIR folder is undefined.\
                         Please enter the path to the superLU_MT library.")
  ENDIF()

  # CHECK IF SUPERLU_DIR EXISTS
  IF(EXISTS ${SUPERLU_DIR})

    SET(SUPERLU_INCLUDE_DIRS "${SUPERLU_DIR}")
    SET(SUPERLU_LIBRARY "${SUPERLU_DIR}/lib/libsuperlu_mt_PTHREAD.a")

    MESSAGE(STATUS "SuperLU Includes: ${SUPERLU_INCLUDE_DIRS}")
    MESSAGE(STATUS "SuperLU Library: ${SUPERLU_LIBRARY}")

  ELSE()
    MESSAGE(FATAL_ERROR "The SUPERLU_DIR folder is invalid.\
                         Please enter the path to the superLU_MT library.")
  ENDIF()

  # ADD FILES IN ADDITIONAL FOLDERS
  FILE(GLOB SUPERLU_SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.cxx"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/superlu/*.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/superlu/*.cxx")
  FILE(GLOB SUPERLU_SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.h"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/superlu/*.h")

  set(SUPERLU_SRC_C ${SUPERLU_SRC_C} CACHE INTERNAL "List of cvOneD supeLU source files")
  set(SUPERLU_SRC_H ${SUPERLU_SRC_H} CACHE INTERNAL "List of cvOneD superLU head files")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD
                      ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse 
                      ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/superlu
                      ${SUPERLU_INCLUDE_DIRS})

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/
                   ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse 
                   ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/superlu)

  simvascular_add_executable(${SV_SVFSI_EXE}
            SRCS ${CSRCS} ${SRC_C} ${SRC_H} ${SUPERLU_SRC_C} ${SUPERLU_SRC_H}
            DEV_SCRIPT_NAME "mysvfsi"
            INSTALL_SCRIPT_NAME "svfsi"
            INSTALL_COMP CoreExecutables
            INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})
  target_compile_definitions(${SV_SVFSI_EXE} PRIVATE USE_SUPERLU)


ELSEIF(sparseSolverType STREQUAL "csparse")

  # ADD FILES IN ADDITIONAL FOLDERS
  FILE(GLOB CSPARSE_SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.cxx"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/csparse/*.c"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/csparse/*.cxx")
  FILE(GLOB CSPARSE_SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/*.h"
                          "${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/csparse/*.h")

  set(CSPARSE_SRC_C ${CSPARSE_SRC_C} CACHE INTERNAL "List of cvOneD supeLU source files")
  set(CSPARSE_SRC_H ${CSPARSE_SRC_H} CACHE INTERNAL "List of cvOneD superLU head files")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD
                      ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse 
                      ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/csparse)

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD
                   ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse 
                   ${CMAKE_CURRENT_SOURCE_DIR}/../cvOneD/sparse/csparse)

  simvascular_add_executable(${SV_SVFSI_EXE}
        SRCS ${CSRCS} ${SRC_C} ${SRC_H} ${CSPARSE_SRC_C} ${CSPARSE_SRC_H}
        DEV_SCRIPT_NAME "mysvfsi"
        INSTALL_SCRIPT_NAME "svfsi"
        INSTALL_COMP CoreExecutables
        INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})
  target_compile_definitions(${SV_SVFSI_EXE} PRIVATE USE_CSPARSE)

ENDIF()

target_link_libraries(${SV_SVFSI_EXE}
  ${GLOBAL_LIBRARIES}
  ${INTELRUNTIME_LIBRARIES}
  ${ZLIB_LIBRARY}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${METIS_SVFSI_LIBRARY_NAME}
  ${PARMETIS_SVFSI_LIBRARY_NAME}
  ${TETGEN_LIBRARY_NAME}
  ${TINYXML_LIBRARY_NAME}
  ${SV_LIB_SVFSILS_NAME}
  ${SV_MPI_NAME_EXT}
  ${VTK_LIBRARIES}
  )

# extra MPI libraries only if there are not set to NOT_FOUND or other null
if(SV_MPI_EXTRA_LIBRARY)
  target_link_libraries(${SV_SVFSI_EXE} ${SV_MPI_EXTRA_LIBRARY})
endif()

if(USE_TRILINOS)
  target_link_libraries(${SV_SVFSI_EXE} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()

if(USE_PETSC)
  target_link_libraries(${SV_SVFSI_EXE} ${PETSC_LIBRARY_DIRS})
endif()

# coverage
if(ENABLE_COVERAGE)
  # set compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -coverage")

  # find required tools
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)

  # add coverage target
  add_custom_target(coverage
    # gather data
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info
    # exclude externals
    COMMAND ${LCOV} --remove coverage.info -o coverage.info '/usr/*' '/opt/*' '${PROJECT_SOURCE_DIR}/ThirdParty/*' '/Library/*' 'v1/*'
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# unit tests and Google Test
if(ENABLE_UNIT_TEST)

  # link pthread on ubuntu20
  find_package(Threads REQUIRED)

  # install Google Test
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  include(GoogleTest)

  # add test.cpp for unit test

  # remove the main.cpp and add test.cpp
  list(APPEND CSRCS "../../../tests/unitTests/test.cpp")
  list(REMOVE_ITEM CSRCS "main.cpp")
  
  # include source files (same as what svFSI does except for main.cpp)
  add_executable(run_all_unit_tests ${CSRCS})
  
  # libraries
  target_link_libraries(run_all_unit_tests
    ${GLOBAL_LIBRARIES}
    ${INTELRUNTIME_LIBRARIES}
    ${ZLIB_LIBRARY}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${METIS_SVFSI_LIBRARY_NAME}
    ${PARMETIS_SVFSI_LIBRARY_NAME}
    ${TETGEN_LIBRARY_NAME}
    ${TINYXML_LIBRARY_NAME}
    ${SV_LIB_SVFSILS_NAME}
    ${SV_MPI_NAME_EXT}
    ${VTK_LIBRARIES}
  )

  # link Google Test
  target_link_libraries(
    run_all_unit_tests
    gtest
    GTest::gtest_main 
    pthread   # link pthread on ubuntu20
  )

  # gtest_discover_tests(runUnitTest)
  add_test(NAME all_unit_tests COMMAND run_all_unit_tests)

endif()

