cmake_minimum_required(VERSION 3.3)

project(LP)
# include(CTest)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(lp
        problem_1.cpp)

target_link_libraries(lp ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples lp )
endif()

# enable_testing()
