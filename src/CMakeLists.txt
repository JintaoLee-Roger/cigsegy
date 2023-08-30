set(PROJECT_INCLUDE_DIR "include")

set(SOURCE_FILES
  segy.cpp
)

add_library(segy STATIC ${SOURCE_FILES})

target_include_directories(segy PUBLIC ${PROJECT_INCLUDE_DIR})

if (OpenMP_CXX_FOUND)
  target_link_libraries(segy PRIVATE OpenMP::OpenMP_CXX)
endif()

target_link_libraries(segy PRIVATE fmt::fmt)

install(TARGETS segy 
  LIBRARY DESTINATION lib
)
