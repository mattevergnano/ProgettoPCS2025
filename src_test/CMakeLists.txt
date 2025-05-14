message(${CMAKE_CURRENT_BINARY_DIR}) 
message(${CMAKE_BINARY_DIR}) 
message(${CMAKE_SOURCE_DIR}) 
message(${CMAKE_CURRENT_SOURCE_DIR}) 
message(${PROJECT_NAME}) 
message(${CMAKE_PROJECT_NAME}) 

list(APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/test_platonic_solids.hpp") 
list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR}) 

set(src_headers ${src_headers} PARENT_SCOPE)
set(src_sources ${src_sources} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)