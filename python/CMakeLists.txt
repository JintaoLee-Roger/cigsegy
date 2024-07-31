pybind11_add_module(cigsegy MODULE
  ${CMAKE_SOURCE_DIR}/src/segy.cpp
  PySegy.cpp 
)

target_include_directories(cigsegy PRIVATE ${CMAKE_SOURCE_DIR}/src/include)
target_include_directories(cigsegy PRIVATE ${PYBIND11_INCLUDE_DIR})
# target_link_libraries(cigsegy PRIVATE segy)
target_compile_definitions(cigsegy PRIVATE USE_PYBIND11)

install(TARGETS cigsegy 
  LIBRARY DESTINATION python/cigsegy
)

FILE(GLOB pyfiles ${CMAKE_SOURCE_DIR}/python/*.py)
message(${pyfiles})

set(COPYRIGHT 
  ${CMAKE_SOURCE_DIR}/LICENSE
  ${CMAKE_SOURCE_DIR}/README.rst
  ${CMAKE_SOURCE_DIR}/Attributions.md
  ${CMAKE_SOURCE_DIR}/python/setup.py)
set(PYFILE
  ${CMAKE_SOURCE_DIR}/python/cigsegy.pyi
  ${pyfiles})
install(FILES ${COPYRIGHT} DESTINATION python)
install(FILES ${PYFILE} DESTINATION python/cigsegy)