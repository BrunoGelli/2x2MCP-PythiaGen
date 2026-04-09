cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(mesonGenProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT
find_package(ROOT CONFIG REQUIRED COMPONENTS Core RIO Tree Physics)

# Threads
find_package(Threads REQUIRED)

# Pythia8
# Prefer the environment variable you already exported:
#   export PYTHIA8=/path/to/pythia-install
find_path(PYTHIA8_INCLUDE_DIR
    NAMES Pythia8/Pythia.h
    HINTS
        $ENV{PYTHIA8}/include
)

find_library(PYTHIA8_LIBRARY
    NAMES pythia8
    HINTS
        $ENV{PYTHIA8}/lib
        $ENV{PYTHIA8}/lib64
)

if(NOT PYTHIA8_INCLUDE_DIR OR NOT PYTHIA8_LIBRARY)
    message(FATAL_ERROR
        "Could not find Pythia8. Make sure PYTHIA8 is set, e.g.\n"
        "  export PYTHIA8=$SCRATCH/pythia/pythia-install\n"
        "and that include/ and lib/ exist under it."
    )
endif()

add_library(Pythia8::Pythia8 UNKNOWN IMPORTED)
set_target_properties(Pythia8::Pythia8 PROPERTIES
    IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR}"
)

# Main executable
add_executable(allInOne_multimeson allInOne_multimeson.cc)

target_link_libraries(allInOne_multimeson
    PRIVATE
        ROOT::Core
        ROOT::RIO
        ROOT::Tree
        ROOT::Physics
        Threads::Threads
        Pythia8::Pythia8
)

# Optional: make config files available from the build directory too
configure_file(beam.config beam.config COPYONLY)
configure_file(momentum.config momentum.config COPYONLY)

message(STATUS "Using PYTHIA8 include dir: ${PYTHIA8_INCLUDE_DIR}")
message(STATUS "Using PYTHIA8 library    : ${PYTHIA8_LIBRARY}")
