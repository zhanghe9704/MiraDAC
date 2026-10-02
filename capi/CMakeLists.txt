# capi/CMakeLists.txt — C API library miradac_c (DA_BUILD_CAPI).

enable_language(C)

file(GLOB MDAC_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Symbol-list generator inputs: nm, preferred as llvm-nm (format-agnostic, so
# it reads Mach-O objects when cross-building for Apple platforms).
find_program(MDAC_NM NAMES llvm-nm nm)
if(NOT MDAC_NM)
    message(FATAL_ERROR "capi: nm not found (needed to generate the Apple export list)")
endif()

function(mdac_library name)
    # The sources are compiled as an OBJECT library: on Apple the exported-
    # symbol list is generated from the compiled objects (see below), which
    # requires the .o files as a unit.
    add_library(${name}_objs OBJECT ${MDAC_SOURCES})
    target_include_directories(${name}_objs PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )
    target_link_libraries(${name}_objs PRIVATE daStatic)
    target_compile_features(${name}_objs PRIVATE cxx_std_17)
    target_compile_definitions(${name}_objs PRIVATE MDAC_VERSION="${PROJECT_VERSION}")
    set_target_properties(${name}_objs PROPERTIES
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
    )

    add_library(${name} SHARED $<TARGET_OBJECTS:${name}_objs>)
    target_include_directories(${name} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(${name} PRIVATE daStatic)
    set_target_properties(${name} PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH ON
    )
    # Export only the MDAC_API (mdac_*) symbols. The version script also hides
    # std:: template instantiations, which keep default visibility under
    # -fvisibility=hidden; --exclude-libs covers daStatic.
    if(APPLE)
        # Apple ld has no --version-script and no --exclude-libs. The
        # equivalent is -exported_symbols_list with every mdac_* symbol,
        # generated from the compiled objects: everything else — daStatic
        # internals, std:: template instantiations, and in the JLL build the
        # statically linked SymEngine — stays hidden, as `local: *` does on
        # GNU ld. The POST_BUILD command below uses $<TARGET_OBJECTS:...>,
        # which CMake supports in TARGET build events only from 3.21.
        if(CMAKE_VERSION VERSION_LESS 3.21)
            message(FATAL_ERROR
                "capi: building miradac_c for APPLE requires CMake >= 3.21 "
                "(the export list is generated from the object files with "
                "a target-dependent generator expression)")
        endif()
        set(_exports "${CMAKE_CURRENT_BINARY_DIR}/exports.${name}.darwin")
        add_custom_command(TARGET ${name}_objs POST_BUILD
            BYPRODUCTS "${_exports}"
            COMMAND ${CMAKE_COMMAND}
                    -DNM=${MDAC_NM}
                    "-DOBJS=$<TARGET_OBJECTS:${name}_objs>"
                    -DOUT=${_exports}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/gen_export_list.cmake
            COMMENT "mdac: generating ${_exports}"
            VERBATIM)
        target_link_options(${name} PRIVATE
            -Wl,-exported_symbols_list,${_exports})
        set_property(TARGET ${name} APPEND PROPERTY LINK_DEPENDS ${_exports})
    else()
        target_link_options(${name} PRIVATE -Wl,--exclude-libs,ALL
            -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/src/exports.map)
        set_property(TARGET ${name} APPEND PROPERTY LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/exports.map)
    endif()
endfunction()

mdac_library(miradac_c)

install(TARGETS miradac_c
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES include/miradac.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(DA_BUILD_TESTS)
    add_executable(test_c_header test/test_c_header.c)
    set_target_properties(test_c_header PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
    target_compile_definitions(test_c_header PRIVATE MDAC_EXPECTED_VERSION="${PROJECT_VERSION}")
    target_link_libraries(test_c_header PRIVATE miradac_c)
    add_test(NAME capi_c_header COMMAND test_c_header)

    # C API against direct C++ results. The test links its own daStatic,
    # separate from the hidden copy inside the C API library.
    mdac_library(miradac_c_asan)
    target_compile_options(miradac_c_asan_objs PRIVATE -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(miradac_c_asan PRIVATE -fsanitize=address)
    foreach(variant IN ITEMS "" "_asan")
        set(t test_capi${variant})
        add_executable(${t} ${CMAKE_SOURCE_DIR}/test/catch_main.cc test/test_capi.cc)
        target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/test)
        target_compile_options(${t} PRIVATE -Wall -Wextra)
        target_link_libraries(${t} PRIVATE miradac_c${variant} daStatic)
        add_test(NAME capi_core${variant} COMMAND ${t})
    endforeach()
    target_compile_options(test_capi_asan PRIVATE -fsanitize=address -fno-omit-frame-pointer)
    target_link_options(test_capi_asan PRIVATE -fsanitize=address)
endif()
