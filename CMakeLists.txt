cmake_minimum_required(VERSION 3.16)
project(da VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# --- Options ---
option(WITH_SYMBOLIC "Build symbolic (SDA) support" ON)
set(DA_CHECK_ENV 1 CACHE STRING "Enable env guard (1=on, 0=off)")

add_compile_definitions(DA_CHECK_ENV=${DA_CHECK_ENV})

# --- Symbolic support ---
if(WITH_SYMBOLIC)
    find_package(SymEngine REQUIRED CONFIG)
    add_compile_definitions(DA_WITH_SYMBOLIC)
    set(DA_SYMENGINE_INCLUDES ${SYMENGINE_INCLUDE_DIRS})
    set(DA_SYMENGINE_LIBS ${SYMENGINE_LIBRARIES})
else()
    set(DA_SYMENGINE_INCLUDES "")
    set(DA_SYMENGINE_LIBS "")
endif()

# --- Library sources ---
set(DA_SOURCES
    src/layout.cpp
    src/env.cpp
    src/engine.cpp
    src/davector.cpp
    src/functions.cpp
    src/base.cpp
    src/symbolic.cpp
)

# --- Shared library ---
add_library(daShared SHARED ${DA_SOURCES})
target_include_directories(daShared
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${DA_SYMENGINE_INCLUDES}
)
if(WITH_SYMBOLIC)
    target_link_libraries(daShared PUBLIC ${DA_SYMENGINE_LIBS})
endif()
set_target_properties(daShared PROPERTIES
    VERSION   ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# --- Static library ---
add_library(daStatic STATIC ${DA_SOURCES})
target_include_directories(daStatic
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        ${DA_SYMENGINE_INCLUDES}
)
if(WITH_SYMBOLIC)
    target_link_libraries(daStatic PUBLIC ${DA_SYMENGINE_LIBS})
endif()

# --- Install rules ---
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS daShared daStatic
    EXPORT daTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/da/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/da
)

install(EXPORT daTargets
    FILE      daTargets.cmake
    NAMESPACE da::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/daConfigVersion.cmake"
    VERSION       ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_SOURCE_DIR}/cmake/daConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/daConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/daConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/daConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/da
)

# --- Tests ---
enable_testing()
add_subdirectory(test)

# --- Examples ---
add_subdirectory(examples)
