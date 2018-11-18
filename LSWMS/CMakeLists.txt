##############################
PROJECT( lineSegment )
##############################
cmake_minimum_required (VERSION 2.6)

# Packages
find_package( OpenCV REQUIRED )


# Includes
include_directories(	
	${OPENCV_INCLUDE_DIR}
	${PROJECT_BINARY_DIR}
)

set(SOURCE_FILES 
	main.cpp
	LSWMS.cpp
	LSWMS.h
)


# Add executable and target link libraries
ADD_EXECUTABLE( lineSegment ${SOURCE_FILES})
TARGET_LINK_LIBRARIES( lineSegment ${OpenCV_LIBS})


