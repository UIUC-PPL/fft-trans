#     This file is part of LLCMc++ (a.k.a. llcmcpp)
#     Copyright (C) 2011 INRIA Julien Bigot <julien.bigot@inria.fr>
#     Copyright (C) 2011 INRIA Christian Pérez <christian.perez@inria.fr>
# 
#     LLCMc++ is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     LLCMc++ is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public License
#     along with LLCMc++.  If not, see <http://www.gnu.org/licenses/>.

# - Find file for Charm++
#
# It set the following variables:
#  Charm_FOUND              - When false, don't try to use Charm++
#  Charm_SHARED             - is Charm++ build with shared library support
#  Charm_INCLUDE_DIR        - Charm++ include directory
#  Charm_INCLUDE_DIRS       - All directories to include to use Charm++ (including dependancies)
#  Charm_LIBRARY            - The Charm++ library
#  Charm_LIBRARIES          - All libraries to link against to use Charm++ (including dependancies)
#  Charm_COMPILE_FLAGS      - Compile flags required to build a Charm++ based object file 
#  Charm_CHARMXI_EXECUTABLE - The charm interface compiler

#### Find paths
#
# Charm_DIR can be used to make it simpler to find the various include
# directories and compiled libraries when omniORB was not installed in the
# usual/well-known directories (e.g. because you made an in tree-source
# compilation or because you installed it in an "unusual" directory).
# Just set Charm_DIR to point to your specific installation directory.
#
# the Charm_DIR environment variable (as opposed to the Charm_DIR cmake variable)
# is also used to search for an omniORB installation.
#
# WARNING: The order of precedence is the following
#    1/ when set Charm_DIR (the cmake variable)
#    2/ when set Charm_DIR (the environment variable)
#    3/ the default system path
# This precedence order goes against the usual pratice and default behavior
# of cmake's FIND_* macros. For more on this debate see e.g.
#     http://www.mail-archive.com/kde-buildsystem@kde.org/msg00589.html

cmake_minimum_required(VERSION 2.8)

include(FindPackageHandleStandardArgs)

set(Charm_DIR "$ENV{Charm_DIR}" CACHE PATH "Charm++ directory")

get_filename_component(_CHARM_CMAKE_DIR  "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_CHARM_CMAKE_DIR  "${_CHARM_CMAKE_DIR}" ABSOLUTE)

# dependancy on pthread
unset(_FIND_PACKAGE_ARGS)
if ( DEFINED PACKAGE_FIND_VERSION )
	list(APPEND _FIND_PACKAGE_ARGS "${PACKAGE_FIND_VERSION}")
endif()
if ( Charm_FIND_QUIETLY )
	list(APPEND _FIND_PACKAGE_ARGS QUIET)
endif()
if ( Charm_FIND_REQUIRED )
	list(APPEND _FIND_PACKAGE_ARGS REQUIRED)
endif()
find_package(Threads ${_FIND_PACKAGE_ARGS})
unset(_FIND_PACKAGE_ARGS)

# Include dir
find_path(Charm_INCLUDE_DIR
	NAMES "charm++.h"
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/include" "$ENV{Charm_DIR}/include" NO_DEFAULT_PATH
	DOC "Charm++ include directory"
)
find_path(Charm_INCLUDE_DIR
	NAMES "charm++.h"
	DOC "Charm++ include directory"
)
list(APPEND Charm_INCLUDE_DIRS "${Charm_INCLUDE_DIR}")

# Include dir
find_file(Charm_MEMORY_DEFAULT
	NAMES "libmemory-default.o"
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Charm++ default memory handler"
)
find_file(Charm_THREADS_DEFAULT
	NAMES "libthreads-default.o"
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Charm++ default threads handler"
)
find_file(Charm_LDB_RAND
	NAMES "libldb-rand.o"
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Charm++ random load balancer"
)

# Library
find_library(Charm_LIBRARY
	NAMES ck
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Charm++ library"
)
find_library(Charm_LIBRARY
	NAMES ck
	DOC "Charm++ library"
)
list(APPEND Charm_LIBRARIES "${Charm_LIBRARY}")

find_library(Charm_CONVCPLUSY_LIBRARY
	NAMES conv-cplus-y
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Converse C++ library"
)
find_library(Charm_CONVCPLUSY_LIBRARY
	NAMES conv-cplus-y
	DOC "Converse C++ library"
)
list(APPEND Charm_LIBRARIES "${Charm_CONVCPLUSY_LIBRARY}")

find_library(Charm_CONVCORE_LIBRARY
	NAMES conv-core
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Converse core library"
)
find_library(Charm_CONVCORE_LIBRARY
	NAMES conv-core
	DOC "Converse core library"
)
list(APPEND Charm_LIBRARIES "${Charm_CONVCORE_LIBRARY}")

find_library(Charm_CONVUTIL_LIBRARY
	NAMES conv-util
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Converse utility library"
)
find_library(Charm_CONVUTIL_LIBRARY
	NAMES conv-util
	DOC "Converse utility library"
)
list(APPEND Charm_LIBRARIES "${Charm_CONVUTIL_LIBRARY}")

find_library(Charm_CKQT_LIBRARY
	NAMES conv-util
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/lib_so" "$ENV{Charm_DIR}/lib_so" "${Charm_DIR}/lib" "$ENV{Charm_DIR}/lib" NO_DEFAULT_PATH
	DOC "Charm++ ckqt library"
)
find_library(Charm_CKQT_LIBRARY
	NAMES conv-util
	DOC "Charm++ ckqt library"
)
list(APPEND Charm_LIBRARIES "${Charm_CKQT_LIBRARY}") 

find_library(Charm_MATH_LIBRARY
	NAMES m
	DOC "standard math library"
)
list(APPEND Charm_LIBRARIES ${CMAKE_THREAD_LIBS_INIT} ${CMAKE_DL_LIBS} ${Charm_MATH_LIBRARY})

# Charmxi
find_program(Charm_CHARMXI
	NAMES charmxi
	PATHS "${Charm_DIR}" "$ENV{Charm_DIR}" "${Charm_DIR}/bin" "$ENV{Charm_DIR}/bin" NO_DEFAULT_PATH
	DOC "path to charmxi (Charm++ interface compiler)"
)
find_program(Charm_CHARMXI
	NAMES charmxi
	DOC "path to charmxi (Charm++ interface compiler)"
)

if ( EXISTS ${Charm_INCLUDE_DIR} )
	execute_process(COMMAND ${_CHARM_CMAKE_DIR}/ReadCharmConfig.sh ${Charm_INCLUDE_DIR} CMK_NO_BUILD_SHARED OUTPUT_VARIABLE _CMK_NO_BUILD_SHARED)
	if ( ${_CMK_NO_BUILD_SHARED} STREQUAL "true" )
		set(Charm_SHARED FALSE)
	else()
		set(Charm_SHARED TRUE)
	endif()
	unset(_CMK_NO_BUILD_SHARED)
endif()

#TODO: find about version

mark_as_advanced(CLEAR Charm_DIR)

#Default arguments handling
find_package_handle_standard_args(Charm REQUIRED_VARS Charm_INCLUDE_DIR Charm_LIBRARY Charm_MEMORY_DEFAULT Charm_THREADS_DEFAULT Charm_SHARED)
set(Charm_FOUND ${CHARM_FOUND})

if ( ${Charm_FOUND} )
	include(${_CHARM_CMAKE_DIR}/CharmMacros.cmake)
	set(Charm_COMPILE_FLAGS "-D__CHARMC__=1 -D_REENTRANT")
endif()

unset(_CHARM_CMAKE_DIR)

#opts: -cs -language charm++ -cc-option -fPIC -shared ${_CHARM_CPP_ARGS} -o ${CMAKE_CURRENT_BINARY_DIR}/${_OBJ} -c ${CMAKE_CURRENT_SOURCE_DIR}/${S}
#CMK_CXX (CMK_PRODUCTION CMK_CXX_PRODUCTION)?

#opts: -cs -language charm++ -cc-option -fPIC -shared {_CHARM_CPP_ARGS} -o {CMAKE_CURRENT_BINARY_DIR}/libmodule{NAME}.so {_OBJECTS}
#CMK_LDXX (CMK_PRODUCTION CMK_C_PRODUCTION)? CMK_LD_SHARED CMK_LD_SHARED_LIBS
