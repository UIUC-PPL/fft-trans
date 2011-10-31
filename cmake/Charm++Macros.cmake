#     This file is part of Gluon.
#     Gluon is a component model build on top of Charm++
#     Copyright © 2011  Julien Bigot <julien.bigot@inria.fr>
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Lesser General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#      ( at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Lesser General Public License for more details.
# 
#     You should have received a copy of the GNU Lesser General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

macro(ci_files PUBLIC_HEADERS)
	include_directories(${CMAKE_CURRENT_BINARY_DIR})
	foreach(_IT ${ARGN})
		if (IS_ABSOLUTE ${_IT})
			file(RELATIVE_PATH ${_IT} ${CMAKE_CURRENT_SOURCE_DIR} ${_IT})
		endif(IS_ABSOLUTE ${_IT})
		get_filename_component(_CP ${_IT} NAME_WE)
		add_custom_command(
			OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_CP}.decl.h ${CMAKE_CURRENT_BINARY_DIR}/${_CP}.def.h
			COMMAND ${CHARMXI} ${CMAKE_CURRENT_SOURCE_DIR}/${_IT}
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
	link_directories(${Charm_LIBRARY_DIRS})
	
	add_executable(${NAME} ${CHARMDIR}/lib/libmemory-default.o ${CHARMDIR}/lib/libthreads-default.o ${CHARMDIR}/lib/libldb-rand.o ${ARGN})
endmacro(add_charm_executable)