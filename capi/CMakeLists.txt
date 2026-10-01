# capi/CMakeLists.txt — C API library miradac_c (DA_BUILD_CAPI).

enable_language(C)

file(GLOB MDAC_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

function(mdac_library name)
    add_library(${name} SHARED ${MDAC_SOURCES})
    target_include_directories(${name} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
    target_link_libraries(${name} PRIVATE daStatic)
    target_compile_features(${name} PRIVATE cxx_std_17)
    target_compile_definitions(${name} PRIVATE MDAC_VERSION="${PROJECT_VERSION}")
    # Export only the MDAC_API (mdac_*) symbols. The version script also hides
    # std:: template instantiations, which keep default visibility under
    # -fvisibility=hidden; --exclude-libs covers daStatic.
    set_target_properties(${name} PROPERTIES
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
        INSTALL_RPATH_USE_LINK_PATH ON
    )
    target_link_options(${name} PRIVATE -Wl,--exclude-libs,ALL
        -Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/src/exports.map)
    set_property(TARGET ${name} APPEND PROPERTY LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/exports.map)
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
    target_compile_options(miradac_c_asan PRIVATE -fsanitize=address -fno-omit-frame-pointer)
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
