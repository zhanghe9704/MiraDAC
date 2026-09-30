# python/CMakeLists.txt — nanobind module miradac._core (added by the root
# CMakeLists.txt when built through scikit-build-core).

find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE
    OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core NB_STATIC NOMINSIZE
    src/module.cpp
    src/bind_env.cpp
    src/bind_nda.cpp
    src/bind_cnda.cpp
)
target_link_libraries(_core PRIVATE daStatic)

if(WITH_SYMBOLIC)
    set_target_properties(_core PROPERTIES INSTALL_RPATH "${DA_SYMENGINE_PREFIX}/lib")
    target_sources(_core PRIVATE src/bind_expr.cpp src/bind_sda.cpp src/bind_csda.cpp)
    # Constants for the interop checks (plan A.8): the libsymengine we link,
    # its sha256 from the pin stamp, its soname, and the pinned symengine.py.
    get_target_property(_se_lib symengine LOCATION)
    get_filename_component(_se_lib "${_se_lib}" REALPATH)
    string(REGEX MATCH "^[0-9]+\\.[0-9]+" _se_major_minor "${DA_SYMENGINE_VERSION}")
    target_compile_definitions(_core PRIVATE
        DA_SYMENGINE_LIB="${_se_lib}"
        DA_SYMENGINE_LIB_SHA256="${DA_SYMENGINE_LIB_SHA256}"
        DA_SYMENGINE_SONAME="libsymengine.so.${_se_major_minor}"
        DA_SYMENGINE_PY_VERSION="${_da_pin_SYMENGINE_PY_VERSION}")
endif()

install(TARGETS _core LIBRARY DESTINATION miradac COMPONENT python)
