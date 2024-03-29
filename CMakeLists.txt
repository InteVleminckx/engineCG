cmake_minimum_required(VERSION 3.5)

project(engine)

############################################################
# Set compiler flags
############################################################
set(OWN_GXX_FLAGS "-std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS} ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(engine_sources
	easy_image.cc
	engine.cc
	ini_configuration.cc
		draw2DLine.cpp
		draw2DLine.h
		hulpFile.cpp
		hulpFile.h
		l_parser.cc
		l_parser.h
		_2DLsystems.cpp
		_2DLsystems.h
        drieDdrawings.cpp
        drieDdrawings.h
		vector3d.cc
		vector3d.h
		DrieDLsystems.cpp
		DrieDLsystems.h
        Z_buffer_driehoeken.cpp Z_buffer_driehoeken.h)

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
