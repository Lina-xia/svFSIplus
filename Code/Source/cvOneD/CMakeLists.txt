set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${SV_SOURCE_DIR}/ThirdParty/parmetis_svfsi/simvascular_parmetis_svfsi/ParMETISLib)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tetgen/simvascular_tetgen)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tinyxml/simvascular_tinyxml)
include_directories(${MPI_C_INCLUDE_PATH})

# COMPILE CODE IN SPARSE SUBDIRECTORY
IF(sparseSolverType STREQUAL "skyline")
  FILE(GLOB SRC_C "${CMAKE_CURRENT_SOURCE_DIR}/Source/cvOneD/*.c"
                  "${CMAKE_CURRENT_SOURCE_DIR}/Source/cvOneD/*.cxx")
  FILE(GLOB SRC_H "${CMAKE_CURRENT_SOURCE_DIR}/Source/cvOneD/*.h")
  # # JUST ADD THE SOURCE FILES
  # ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H})

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
  FILE(GLOB SUPERLU_SRC_C "${CMAKE_SOURCE_DIR}/sparse/*.c"
                          "${CMAKE_SOURCE_DIR}/sparse/*.cpp"
                          "${CMAKE_SOURCE_DIR}/sparse/superlu/*.c"
                          "${CMAKE_SOURCE_DIR}/sparse/superlu/*.cpp")
  FILE(GLOB SUPERLU_SRC_H "${CMAKE_SOURCE_DIR}/sparse/*.h"
                          "${CMAKE_SOURCE_DIR}/sparse/superlu/*.h"
                          "${SUPERLU_INCLUDE_DIRS}/slu_mt_ddefs.h")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}
                      ${CMAKE_SOURCE_DIR}/sparse 
                      ${CMAKE_SOURCE_DIR}/sparse/superlu
                      ${SUPERLU_INCLUDE_DIRS})

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}
                   ${CMAKE_SOURCE_DIR}/sparse 
                   ${CMAKE_SOURCE_DIR}/sparse/superlu)

  # ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H} ${SUPERLU_SRC_C} ${SUPERLU_SRC_H})
  # TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${BLAS_LIBRARIES} ${SUPERLU_LIBRARY} pthread)  

ELSEIF(sparseSolverType STREQUAL "csparse")

  # ADD FILES IN ADDITIONAL FOLDERS
  FILE(GLOB CSPARSE_SRC_C "${CMAKE_SOURCE_DIR}/sparse/*.c"
                          "${CMAKE_SOURCE_DIR}/sparse/*.cpp"
                          "${CMAKE_SOURCE_DIR}/sparse/csparse/*.c"
                          "${CMAKE_SOURCE_DIR}/sparse/csparse/*.cpp")
  FILE(GLOB CSPARSE_SRC_H "${CMAKE_SOURCE_DIR}/sparse/*.h"
                          "${CMAKE_SOURCE_DIR}/sparse/csparse/*.h")

  # INCLUDE DIRECTORIES
  INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}
                      ${CMAKE_SOURCE_DIR}/sparse 
                      ${CMAKE_SOURCE_DIR}/sparse/csparse)

  # LINK DIRECTORIES
  LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}
                   ${CMAKE_SOURCE_DIR}/sparse 
                   ${CMAKE_SOURCE_DIR}/sparse/csparse)

  # # ADD SOURCE ON SPARSE FOLDER
  # ADD_EXECUTABLE(${PROJECT_NAME} ${SRC_C} ${SRC_H} ${CSPARSE_SRC_C} ${CSPARSE_SRC_H})

ENDIF()


