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

add_executable(lp_problem_0 "problem_0.cpp")
add_executable(lp_problem_0_2 "problem_0_2.cpp")
add_executable(lp_problem_1 "problem_1.cpp")
add_executable(lp_problem_2 "problem_2.cpp")
add_executable(lp_problem_3 "problem_3.cpp")

target_link_libraries(lp_problem_0 ${SCIP_LIBRARIES})
target_link_libraries(lp_problem_0_2 ${SCIP_LIBRARIES})
target_link_libraries(lp_problem_1 ${SCIP_LIBRARIES})
target_link_libraries(lp_problem_2 ${SCIP_LIBRARIES})
target_link_libraries(lp_problem_3 ${SCIP_LIBRARIES})

if( TARGET examples )
    add_dependencies( examples lp_problem_0 )
endif()

# enable_testing()
