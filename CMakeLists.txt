INCLUDE_DIRECTORIES(
        ${PROJECT_SOURCE_DIR}
        ${PROJECT_BINARY_DIR}
)

LINK_DIRECTORIES(
        ${PROJECT_BINARY_DIR}/lemon
)


ADD_EXECUTABLE(vrp main.cc vrp.cc)
TARGET_LINK_LIBRARIES(vrp lemon)

#ADD_EXECUTABLE(in_gen input_generator.cc)
#TARGET_LINK_LIBRARIES(in_gen lemon)