include(cmake/modules/common-compiler-flags.cmake)
tfel_enable_cxx_compiler_flag(OPTIMISATION_FLAGS "fno-fast-stdlib")
tfel_enable_cxx_compiler_flag(OPTIMISATION_FLAGS "march=auto")

set(OPTIMISATION_FLAGS "-DNO_RUNTIME_CHECK_BOUNDS ${OPTIMISATION_FLAGS}")

if((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))
  set(OPTIMISATION_FLAGS "-O2 -DNDEBUG -OPT:Olimit=0 ${OPTIMISATION_FLAGS}")
endif((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE STREQUAL "Release"))

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_definitions("-g")
endif(CMAKE_BUILD_TYPE STREQUAL "Debug")

if(HAVE_FORTRAN)
  if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PathScale")
    set(PATHSCALE_FORTRAN_COMPILER ON)
    if(MSYS OR APPLE OR ("${CMAKE_SYSTEM_NAME}" STREQUAL "FreeBSD"))
      message(FATAL_ERROR "pathscal compiler is not supported under MSYS")
    else(MSYS OR APPLE OR ("${CMAKE_SYSTEM_NAME}" STREQUAL "FreeBSD"))
      add_definitions("-D'F77_FUNC(X,Y)=X\\#\\#_' -D'F77_FUNC_(X,Y)=X\\#\\#_'")
    endif(MSYS OR APPLE OR ("${CMAKE_SYSTEM_NAME}" STREQUAL "FreeBSD"))
  else(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PathScale")
    message(FATAL_ERROR "unsupported fortran compiler ${CMAKE_Fortran_COMPILER_NAME}")
  endif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "PathScale")
endif(HAVE_FORTRAN)

