project(tests)

list(APPEND Tests test_matrix)

foreach(TEST IN LISTS Tests)
    add_executable(${TEST} ${TEST}.cpp)
    target_link_libraries(${TEST} CMatrix)
    target_include_directories(${TEST} PRIVATE ../inc)
    set_target_properties(${TEST} PROPERTIES
        CXX_STANDARD 17
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/bin")
    add_test(NAME ${TEST} COMMAND ${TEST})
endforeach()