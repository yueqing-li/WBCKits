cmake_minimum_required(VERSION 3.0)


set ( WBCKits_VERSION_MAJOR 0 )
set ( WBCKits_VERSION_MINOR 1 )
set ( WBCKits_VERSION_PATCH 2 )

set ( WBCKits_VERSION
  ${WBCKits_VERSION_MAJOR}.${WBCKits_VERSION_MINOR}.${WBCKits_VERSION_PATCH}
)
set ( WBCKits_SO_VERSION
  ${WBCKits_VERSION_MAJOR}.${WBCKits_VERSION_MINOR}
)
set (PROJECT_VERSION ${WBCKits_VERSION})

# Setup the project.
project(WBCKITS VERSION ${PROJECT_VERSION})
set(CMAKE_BUILD_TYPE "Release") # Debug Release

# findXXX.cmake
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)


#find packages
find_package(Eigen3 REQUIRED)

IF (NOT EIGEN3_FOUND)
  message(WARNING "Could not find Eigen3 on your system.")
ENDIF (NOT EIGEN3_FOUND)

IF (EIGEN3_FOUND)
  INCLUDE_DIRECTORIES (SYSTEM ${EIGEN3_INCLUDE_DIR})
  message(STATUS "Find eigen3 directory : " ${EIGEN3_INCLUDE_DIR})
ENDIF (EIGEN3_FOUND)

# Options
option(WBCKITS_TEST "If ON, build wbc_test project" ON)
option(WBCKITS_BUILD_STATIC "If ON, build static library instead of shared" OFF)

# Source files for WBCKITS
set ( WBCKITS_SOURCES
  src/wbc_base.cpp
  src/task.cpp
  src/constraint.cpp
  src/options.cpp
  src/wqp.cpp
  src/rhp.cpp
  src/hqp.cpp
  src/ghc.cpp
 )

######## build library #############

IF (WBCKITS_BUILD_STATIC)
  include_directories(${PROJECT_SOURCE_DIR}/include)
  add_library( wbckits-static STATIC ${WBCKITS_SOURCES})
  set_target_properties(wbckits-static PROPERTIES 
      OUTPUT_NAME "wbckits")
  get_target_property(OUTPUT_VALUE wbckits-static OUTPUT_NAME)
  message(STATUS "The wbckits-static OUTPUT_NAME is:" ${OUTPUT_VALUE})
ELSE (WBCKITS_BUILD_STATIC)
  add_library( wbckits SHARED ${WBCKITS_SOURCES})
  target_include_directories( wbckits 
    PUBLIC   
      ${PROJECT_SOURCE_DIR}/include)
  target_link_libraries( wbckits
      qpOASES)
  set_target_properties(wbckits PROPERTIES 
      VERSION ${WBCKits_VERSION}
      SOVERSION ${WBCKits_SO_VERSION})
ENDIF (WBCKITS_BUILD_STATIC)


######## INSTALL ###########

IF (WBCKITS_BUILD_STATIC)
  install(TARGETS wbckits-static
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
ELSE (WBCKITS_BUILD_STATIC)
  install(TARGETS wbckits
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
ENDIF (WBCKITS_BUILD_STATIC)

#headers
install(FILES include/WBCKits.hpp
  DESTINATION include)
install(DIRECTORY include/WBCKits
  DESTINATION include
  FILES_MATCHING PATTERN "*.h")
  install(DIRECTORY include/WBCKits
  DESTINATION include
  FILES_MATCHING PATTERN "*.hpp")


#######  BUILD TEST ###########

IF (WBCKITS_TEST)
    add_subdirectory(test)
ELSE (WBCKITS_TEST)

ENDIF (WBCKITS_TEST)


