cmake_minimum_required(VERSION 2.8)

set(CHARMDIR "" CACHE PATH "Charm++ installation directory")

set(CHARMC ${CHARMDIR}/bin/charmc CACHE FILEPATH "Charm++ compiler" FORCE)
set(CHARMXI ${CHARMDIR}/bin/charmxi CACHE FILEPATH "Charm++ ci compiler" FORCE)
set(Charm_INCLUDE_PATH ${CHARMDIR}/include CACHE FILEPATH "Charm++ include path" FORCE)
set(Charm_COMPILE_FLAGS "-D__CHARMC__=1 -D_REENTRANT" CACHE FILEPATH "Charm++ include path" FORCE)
set(Charm_LIBRARY_DIRS ${CHARMDIR}/lib_so CACHE FILEPATH "Charm++ library dir path" FORCE)
set(Charm_LINK_FLAGS "-lck -lconv-cplus-y -lconv-core -lconv-util -lpthread -lckqt -ldl -lm" CACHE STRING "Charm++ linking flags" FORCE)
mark_as_advanced(CHARMC CHARMXI Charm_INCLUDE_PATH Charm_COMPILE_FLAGS Charm_LIBRARY_DIRS Charm_LINK_FLAGS)

get_filename_component(_CHARM_CMAKE_DIR  "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_CHARM_CMAKE_DIR  "${_CHARM_CMAKE_DIR}" ABSOLUTE)
include(${_CHARM_CMAKE_DIR}/Charm++Macros.cmake)

#reqs:
# CMK_NO_BUILD_SHARED==false
execute_process(COMMAND ${_CHARM_CMAKE_DIR}/ReadCharmConfig.sh ${CHARMDIR}/include CMK_NO_BUILD_SHARED OUTPUT_VARIABLE Charm_NO_BUILD_SHARED)
if (${Charm_NO_BUILD_SHARED} STREQUAL "true")
	message(SEND_ERROR "Charm++ build without shared lib support")
else()
	message(STATUS "Charm++ Found")
endif(${Charm_NO_BUILD_SHARED} STREQUAL "true")

#opts: -cs -language charm++ -cc-option -fPIC -shared ${_CHARM_CPP_ARGS} -o ${CMAKE_CURRENT_BINARY_DIR}/${_OBJ} -c ${CMAKE_CURRENT_SOURCE_DIR}/${S}
#CMK_CXX (CMK_PRODUCTION CMK_CXX_PRODUCTION)?

#opts: -cs -language charm++ -cc-option -fPIC -shared {_CHARM_CPP_ARGS} -o {CMAKE_CURRENT_BINARY_DIR}/libmodule{NAME}.so {_OBJECTS}
#CMK_LDXX (CMK_PRODUCTION CMK_C_PRODUCTION)? CMK_LD_SHARED CMK_LD_SHARED_LIBS

unset(_CHARM_CMAKE_DIR)