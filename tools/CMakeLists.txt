add_executable(SEGYRead SEGYRead.cpp cxxopts.hpp)
target_link_libraries(SEGYRead PUBLIC segy fmt::fmt)

if(NOT MSVC)
  set(compilecmd -Wall -Wextra -pedantic -Werror -Wno-unused-parameter -fms-extensions)
else()
  set(compilecmd "")
endif()

if(NOT MSVC)
target_compile_options(SEGYRead PRIVATE -Werror=return-type)
endif()
set_target_properties(SEGYRead PROPERTIES FOLDER tools)
target_compile_options(SEGYRead PRIVATE ${compilecmd})

install(TARGETS SEGYRead 
  RUNTIME DESTINATION bin
)
if(APPLE)
  set_target_properties(SEGYRead PROPERTIES INSTALL_RPATH "@executable_path/../lib${LIBSUFFIX}")
else()
  set_target_properties(SEGYRead PROPERTIES INSTALL_RPATH "$ORIGIN/../lib${LIBSUFFIX}")
endif()



add_executable(SEGYCreate SEGYCreate.cpp cxxopts.hpp)
target_link_libraries(SEGYCreate PUBLIC segy fmt::fmt)

if(NOT MSVC)
  target_compile_options(SEGYCreate PRIVATE -Werror=return-type)
endif()
set_target_properties(SEGYCreate PROPERTIES FOLDER tools)
target_compile_options(SEGYCreate PRIVATE ${compilecmd})

install(TARGETS SEGYCreate
  RUNTIME DESTINATION bin
)
if(APPLE)
  set_target_properties(SEGYCreate PROPERTIES INSTALL_RPATH "@executable_path/../lib")
else()
  set_target_properties(SEGYCreate PROPERTIES INSTALL_RPATH "$ORIGIN/../lib")
endif()
