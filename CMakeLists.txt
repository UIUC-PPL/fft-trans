project(fft_trans)
cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Charm)

include_directories(include ${Charm_INCLUDE_DIRS})
link_directories(${CMAKE_CURRENT_BINARY_DIR} ${Charm_LIBRARY_DIRS})



ci_files(Fft_HEADERS src/Fft/Fft.ci)
add_library(Fft SHARED src/Fft/fft.cpp ${Fft_HEADERS})
set_property(TARGET Fft APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET Fft APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS Fft LIBRARY DESTINATION lib)
install(FILES ${Fft_HEADERS} DESTINATION include/fft_trans)

ci_files(transpose_HEADERS src/transpose/transpose.ci)
add_library(transpose SHARED src/transpose/transpose.cpp ${transpose_HEADERS})
set_property(TARGET transpose APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET transpose APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS transpose LIBRARY DESTINATION lib)
install(FILES ${transpose_HEADERS} DESTINATION include/fft_trans)

ci_files(Validation_HEADERS src/Validation/Validation.ci)
add_library(Validation SHARED src/Validation/validation.cpp ${Validation_HEADERS})
set_property(TARGET Validation APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET Validation APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS Validation LIBRARY DESTINATION lib)
install(FILES ${Validation_HEADERS} DESTINATION include/fft_trans)



ci_files(FftComputer_HEADERS src/FftComputer/FftComputer.ci)
add_library(FftComputer MODULE src/FftComputer/fftcomputer.cpp ${FftComputer_HEADERS})
target_link_libraries(FftComputer Fft transpose -lfftw3)
set_property(TARGET FftComputer APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET FftComputer APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS FftComputer LIBRARY DESTINATION lib)

ci_files(FftValidator_HEADERS src/FftValidator/FftValidator.ci)
add_library(FftValidator MODULE src/FftValidator/fftvalidator.cpp ${FftValidator_HEADERS})
target_link_libraries(FftValidator Fft Validation -lfftw3)
set_property(TARGET FftValidator APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET FftValidator APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS FftValidator LIBRARY DESTINATION lib)

ci_files(p2p_transpose_HEADERS src/p2p_transpose/p2p_transpose.ci)
add_library(p2p_transpose MODULE src/p2p_transpose/p2p_transpose.cpp ${p2p_transpose_HEADERS})
target_link_libraries(p2p_transpose transpose)
set_property(TARGET p2p_transpose APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET p2p_transpose APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS p2p_transpose LIBRARY DESTINATION lib)

ci_files(ValidatorMain_HEADERS src/ValidatorMain/ValidatorMain.ci)
add_library(ValidatorMain MODULE src/ValidatorMain/validatormain.cpp ${ValidatorMain_HEADERS})
target_link_libraries(ValidatorMain Validation)
set_property(TARGET ValidatorMain APPEND PROPERTY COMPILE_FLAGS ${Charm_COMPILE_FLAGS})
set_property(TARGET ValidatorMain APPEND PROPERTY LINK_FLAGS ${Charm_LINK_FLAGS})
install(TARGETS ValidatorMain LIBRARY DESTINATION lib)

install(DIRECTORY include/ DESTINATION include/fft_trans FILES_MATCHING PATTERN "*.h")
install(DIRECTORY gad/ DESTINATION lib/gad/fft_trans FILES_MATCHING PATTERN "*.gad")