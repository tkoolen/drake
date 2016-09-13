include(${LCM_USE_FILE})

if(PYTHONINTERP_FOUND)
  set(python_args PYTHON_SOURCES python_install_sources)
endif()
if(JAVA_FOUND AND TARGET lcm-java)
  set(java_args JAVA_SOURCES java_sources)
endif()

lcm_wrap_types(
    C_SOURCES sources
    C_HEADERS c_install_headers
    CPP_HEADERS cpp_install_headers
    ${python_args}
    ${java_args}
    valkyrie_effort_command_t.lcm
)

add_custom_target(drake_valkyrie_lcmtypes_hpp DEPENDS ${cpp_install_headers})

add_library(drake_valkyrie_lcmtypes STATIC ${sources})
target_link_libraries(drake_valkyrie_lcmtypes lcm)

install(FILES ${c_install_headers} DESTINATION include/lcmtypes)
install(FILES ${cpp_install_headers} DESTINATION include/lcmtypes/drake)

install(TARGETS drake_valkyrie_lcmtypes
    ARCHIVE DESTINATION lib${LIB_SUFFIX}
    INCLUDES DESTINATION include
    )
