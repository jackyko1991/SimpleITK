set(CTEST_SITE "circleci")

set(CTEST_DASHBOARD_ROOT "$ENV{SIMPLEITK_BUILD_DIR}")

string(SUBSTRING $ENV{CIRCLE_SHA1} 0 7 commit)

set(what "#$ENV{CIRCLE_PR_NUMBER}")
if("$ENV{CIRCLE_PR_NUMBER}" STREQUAL "")
  set(what "$ENV{CIRCLE_BRANCH}")
endif()
set(CTEST_BUILD_NAME "_${what}_${commit}")
set(CTEST_CONFIGURATION_TYPE Release)
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(CTEST_SOURCE_DIRECTORY "$ENV{SIMPLEITK_SRC_DIR}")
set(dashboard_no_update 1)

if( DEFINED ENV{MAKE_J} )
  set(CTEST_BUILD_FLAGS "-j$ENV{MAKE_J}")
endif()


#set(CTEST_TEST_ARGS PARALLEL_LEVEL 8)

set(dashboard_model Experimental)
set(dashboard_track Circle-CI)

set(dashboard_cache "${dashboard_cache}
  ENABLE_SHARED:BOOL=ON
  WRAP_LUA:BOOL=OFF
  WRAP_PYTHON:BOOL=OFF
  WRAP_JAVA:BOOL=OFF
  WRAP_CSHARP:BOOL=OFF
  WRAP_TCL:BOOL=OFF
  WRAP_R:BOOL=OFF
  WRAP_RUBY:BOOL=OFF
  ITK_REPOSITORY:PATH=$ENV{ITK_REPOSITORY}
  USE_SYSTEM_SWIG:BOOL=ON
  USE_SYSTEM_LUA:BOOL=ON
  CMAKE_BUILD_TYPE=Release
")

# Include driver script
include(${CTEST_SCRIPT_DIRECTORY}/simpleitk_common.cmake)
