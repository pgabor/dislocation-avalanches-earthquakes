# Sources
add_subdirectory("fields")

set(SRC
    simulation_data.cpp
    adiabatic_stress.cpp
    argreader.cpp
    precision_handler.cpp
    simulation.cpp
)

add_library(sdddstEQCore ${SRC})

target_link_libraries(sdddstEQCore sdddstEQFields)
target_link_libraries(sdddstEQCore ${UMFPACK_LIBRARIES})

add_executable(sdddstEQ main.cpp)
target_link_libraries(sdddstEQ sdddstEQCore)
