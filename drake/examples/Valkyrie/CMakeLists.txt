# drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieBalancing REQUIRES gurobi COMMAND runValkyrieBalancing)  # FIXME: see #2839

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug" AND UNIX) # FIXME: see #3147.
  drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieBalancingPerturb REQUIRES gurobi lcm libbot yaml-cpp COMMAND runValkyrieBalancingPerturb)
endif()

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug" AND UNIX) # FIXME: see #3147.
  drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieWalking REQUIRES gurobi lcm libbot yaml-cpp COMMAND runValkyrieWalking SIZE large)
endif()

drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieVisualize OPTIONAL bullet gurobi snopt COMMAND runValkyrieVisualize)

add_library(robot_state_translator robot_state_translator.cc)
pods_install_libraries(robot_state_translator)
drake_install_headers(robot_state_translator.h)
target_link_libraries(robot_state_translator
    drakeLCMSystem2
    drakeLCMUtil
    drakeRBM)

add_executable(valkyrie_simulation valkyrie_simulation.cc)
add_dependencies(valkyrie_simulation drake_lcmtypes_hpp)
target_link_libraries(valkyrie_simulation
    drake_rbp
    drakeLCMSystem2
    robot_state_translator
    drakeSystemAnalysis
    drakeSystemFramework)

add_subdirectory(test)
