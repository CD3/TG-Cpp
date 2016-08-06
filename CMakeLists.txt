cmake_minimum_required(VERSION 2.8)

project(TGModel)


# FIND DEPENDENCIES

find_package( OpenMP )
if( OPENMP_FOUND )
  set( CMAKE_C_FLAGS   ${CMAKE_C_FLAGS}   ${OpenMP_C_FLAGS}   )
  set( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} )
endif()


# VERSION NUMBER

# this will extract version information from the git tags in the repository
set( GIT_COMMIT_DESC "X.X.X-UNKNOWN" )
find_package( Git )
if( GIT_FOUND )
  execute_process( COMMAND ${GIT_EXECUTABLE} status
                   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                   RESULT_VARIABLE IsGitRepo
                   OUTPUT_VARIABLE Trash
                   ERROR_VARIABLE  Trash)

  if( ${IsGitRepo} EQUAL 0 )
    execute_process( COMMAND ${GIT_EXECUTABLE} describe --tags --dirty HEAD
                     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                     OUTPUT_VARIABLE GIT_COMMIT_DESC
                     ERROR_VARIABLE Trash
                     # OUTPUT_STRIP_TRAILING_WHITESPACE
                     )
  endif()
endif()

# BUILD

include_directories(src)

set(SOURCES src/TGModel.cpp src/TGlib.cpp)

add_executable(TGModel ${SOURCES})
