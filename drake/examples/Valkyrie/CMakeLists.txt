# drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieBalancing REQUIRES gurobi COMMAND runValkyrieBalancing)  # FIXME: see #2839

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug" AND UNIX) # FIXME: see #3147.
  drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieBalancingPerturb REQUIRES gurobi lcm libbot yaml-cpp COMMAND runValkyrieBalancingPerturb)
endif()

if (NOT ${CMAKE_BUILD_TYPE} STREQUAL "Debug" AND UNIX) # FIXME: see #3147.
  drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieWalking REQUIRES gurobi lcm libbot yaml-cpp COMMAND runValkyrieWalking SIZE large)
endif()

drake_add_matlab_test(NAME examples/Valkyrie/runValkyrieVisualize OPTIONAL bullet gurobi snopt COMMAND runValkyrieVisualize)

#add_library_with_exports(LIB_NAME drakeRobotCommandToDesiredEffortConverter SOURCE_FILES
#    simple_robot_command_to_desired_effort_converter.cc)
#target_link_libraries(drakeRobotCommandToDesiredEffortConverter
#    drakeLCMSystem2)
#pods_install_libraries(drakeRobotCommandToDesiredEffortConverter)
#drake_install_headers(simple_robot_command_to_desired_effort_converter.h)

add_library_with_exports(LIB_NAME drakeRobotStateEncoder SOURCE_FILES
    robot_state_encoder.cc)
add_dependencies(drakeRobotStateEncoder drake_lcmtypes_hpp)
target_link_libraries(drakeRobotStateEncoder
    drakeLCMSystem2
    drakeLCMUtil
    drakeRigidBodyPlant
    drakeSide
    )
pods_install_libraries(drakeRobotStateEncoder)
drake_install_headers(robot_state_encoder.h)

add_executable(
    valkyrie_simulation valkyrie_simulation.cc
    simple_robot_command_to_desired_effort_converter.cc
    simple_joint_sensors.cc simple_joint_sensors.h)
add_dependencies(valkyrie_simulation drake_lcmtypes_hpp)
target_link_libraries(valkyrie_simulation
    drakeRigidBodyPlant
    drakeLCMSystem2
    drakeRobotStateEncoder
    drakeSystemAnalysis
    drakeSystemFramework)
drake_install_headers(simple_robot_command_to_desired_effort_converter.h)

add_subdirectory(test)
