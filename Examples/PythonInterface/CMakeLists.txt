cmake_minimum_required(VERSION 2.8)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(SimpleITK REQUIRED)
include(${SimpleITK_USE_FILE})

find_package( PythonLibs REQUIRED)

include_directories( ${PYTHON_INCLUDE_DIR} )

add_library( SimpleITKPythonInterface MODULE PythonInterface.cxx )
target_link_libraries( SimpleITKPythonInterface
  ${ITK_LIBRARIES}
  ${SimpleITK_LIBRARIES}
  ${PYTHON_LIBRARIES}
  /dev/shm/SimpleITK/SimpleITK-build/lib/_SimpleITK.so  )
set_target_properties( SimpleITKPythonInterface  PROPERTIES PREFIX "")
