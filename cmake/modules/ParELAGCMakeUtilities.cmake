# Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at the
# Lawrence Livermore National Laboratory. LLNL-CODE-745557. All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the ParElag library. For more information and source code
# availability see http://github.com/LLNL/parelag.
#
# ParElag is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

# A handy function to add the current source directory to a local
# filename. To be used for creating a list of sources.
function(convert_filenames_to_full_paths NAMES)
  unset(tmp_names)
  foreach(name ${${NAMES}})
    list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
  endforeach()
  set(${NAMES} ${tmp_names} PARENT_SCOPE)
endfunction()

# A function to add each executable in the list to the build with the
# correct flags, includes, and linkage
function(add_parelag_executables EXE_SRCS)
  # Add one executable per cpp file
  foreach(SRC_FILE IN LISTS ${EXE_SRCS})
    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME)

    string(REPLACE ".cpp" ".exe" EXE_NAME ${SRC_FILENAME})
    add_parelag_executable(${EXE_NAME}
      MAIN ${SRC_FILE}
      LIBRARIES ${PROJECT_NAME})
  endforeach(SRC_FILE)
endfunction()

# A slightly more versitile function for adding executables to ParELAG
function(add_parelag_executable PARELAG_EXE_NAME)
  # Parse the input arguments looking for the things we need
  set(POSSIBLE_ARGS "MAIN" "EXTRA_SOURCES" "EXTRA_HEADERS" "EXTRA_OPTIONS" "EXTRA_DEFINES" "LIBRARIES")
  set(CURRENT_ARG)
  foreach(arg ${ARGN})
    list(FIND POSSIBLE_ARGS ${arg} is_arg_name)
    if (${is_arg_name} GREATER -1)
      set(CURRENT_ARG ${arg})
      set(${CURRENT_ARG}_LIST)
    else()
      list(APPEND ${CURRENT_ARG}_LIST ${arg})
    endif()
  endforeach()

  # Actually add the exe
  add_executable(${PARELAG_EXE_NAME} ${MAIN_LIST}
    ${EXTRA_SOURCES_LIST} ${EXTRA_HEADERS_LIST})

  # Append the additional libraries and options
  target_link_libraries(${PARELAG_EXE_NAME} PRIVATE ${LIBRARIES_LIST})
  target_compile_options(${PARELAG_EXE_NAME} PRIVATE ${EXTRA_OPTIONS_LIST})
  target_compile_definitions(${PARELAG_EXE_NAME} PRIVATE ${EXTRA_DEFINES_LIST})

  # Handle the MPI separately
  target_link_libraries(${PARELAG_EXE_NAME} PRIVATE ${MPI_CXX_LIBRARIES})

  target_include_directories(${PARELAG_EXE_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH})
  if (MPI_CXX_COMPILE_FLAGS)
    target_compile_options(${PARELAG_EXE_NAME} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
  endif()

  if (MPI_CXX_LINK_FLAGS)
    set_target_properties(${PARELAG_EXE_NAME} PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif()

  # FIXME: This should be generalized or scrapped entirely
  if(VALGRIND_LIB_DIR)
    target_link_libraries(${EXE_NAME} debug "${VALGRIND_LIB_DIR}/libmpiwrap-amd64-darwin.so")
  endif()

endfunction()

# Function to add unit tests to ctest. This will just use the cpp
# name (without the .cpp) as the test name registered with CMake.
function(add_parelag_unittests EXE_SRCS)
  foreach(SRC_FILE IN LISTS ${EXE_SRCS})
    get_filename_component(SRC_FILENAME ${SRC_FILE} NAME)

    string(REPLACE ".cpp" ".exe" EXE_NAME ${SRC_FILENAME})
    string(REPLACE ".cpp" "" TEST_NAME ${SRC_FILENAME})

    add_test(unit_${TEST_NAME} ${EXE_NAME})
    set_tests_properties(unit_${TEST_NAME}
      PROPERTIES
      PASS_REGULAR_EXPRESSION
      "Unit test passed.")
  endforeach()
endfunction()

# Function that uses "dumb" logic to try to figure out if a library
# file is a shared or static library. This won't work on Windows; it
# will just return "unknown" for everything.
function(parelag_determine_library_type lib_name output_var)

  # Test if ends in ".a"
  string(REGEX MATCH "\.a$" _static_match ${lib_name})
  if (_static_match)
    set(${output_var} STATIC PARENT_SCOPE)
    return()
  endif (_static_match)

  # Test if ends in ".so(.version.id.whatever)"
  string(REGEX MATCH "\.so($|\..*$)" _shared_match ${lib_name})
  if (_shared_match)
    set(${output_var} SHARED PARENT_SCOPE)
    return()
  endif (_shared_match)

  # Test if ends in ".dylib(.version.id.whatever)"
  string(REGEX MATCH "\.dylib($|\..*$)" _mac_shared_match ${lib_name})
  if (_mac_shared_match)
    set(${output_var} SHARED PARENT_SCOPE)
    return()
  endif (_mac_shared_match)

  set(${output_var} "UNKNOWN" PARENT_SCOPE)
endfunction(parelag_determine_library_type lib_name output)
