include_directories(${CMAKE_INSTALL_PREFIX}/include/lcmtypes) #TODO(#3351)


add_library_with_exports(LIB_NAME drakeRobotStateDecoder SOURCE_FILES
    robot_state_decoder.cc)
add_dependencies(drakeRobotStateDecoder drake_lcmtypes_hpp)
target_link_libraries(drakeRobotStateDecoder
    drakeLCMSystem2
    drakeLCMUtil
    drakeRigidBodyPlant
    drakeRobotStateLcmTypeUtil
    )
pods_install_libraries(drakeRobotStateDecoder)
drake_install_headers(robot_state_decoder.h)
