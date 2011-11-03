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

cmake_minimum_required(VERSION 2.8)

macro(ci_files PUBLIC_HEADERS)
	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	foreach(_IT ${ARGN})
		if (IS_ABSOLUTE ${_IT})
			file(RELATIVE_PATH ${_IT} ${CMAKE_CURRENT_SOURCE_DIR} ${_IT})
		endif(IS_ABSOLUTE ${_IT})
		get_filename_component(_CP ${_IT} NAME_WE)
		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_CP}.decl.h ${CMAKE_CURRENT_BINARY_DIR}/${_CP}.def.h
			COMMAND ${Charm_CHARMXI} ${CMAKE_CURRENT_SOURCE_DIR}/${_IT}
			MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${_IT}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			VERBATIM
		)
		list(APPEND ${PUBLIC_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/${_CP}.decl.h)
		unset(_CP)
	endforeach(_IT)
	unset(_IT)
endmacro(ci_files)

macro(add_charm_executable NAME)
	include_directories(${Charm_INCLUDE_PATH})
	add_executable(${NAME} ${Charm_MEMORY_DEFAULT} ${Charm_THREADS_DEFAULT} ${Charm_LDB_RAND} ${ARGN})
	target_link_libraries(${NAME} ${Charm_LIBRARIES})
	set_property(TARGET ${NAME} PROPERTY INSTALL_RPATH_USE_LINK_PATH true)
endmacro(add_charm_executable)