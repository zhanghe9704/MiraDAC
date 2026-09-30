# CheckSymEnginePin.cmake — find SymEngine and verify it is the pinned build.
#
# The pin is cmake/symengine_pin.txt. scripts/setup_symengine.sh builds the
# pinned commit and writes the stamp <prefix>/share/symengine/miradac-pin.txt;
# the pinned commit reports the same version as the v0.14.0 tag, so the stamp is
# what tells them apart. This file checks that the stamp's commit and build
# options equal the pin.
#
# With DA_IGNORE_SYMENGINE_PIN=ON any SymEngine version is accepted and a
# mismatch only warns. That exists solely to benchmark against another
# SymEngine release; never use it for a build that talks to symengine.py.
#
# Sets: DA_SYMENGINE_VERSION (pinned version), DA_SYMENGINE_PREFIX,
#       DA_SYMENGINE_LIB_SHA256 (from the stamp, empty when there is none),
#       DA_SYMENGINE_FIND_VERSION (version arguments for daConfig.cmake's
#       find_dependency: "<version> EXACT", or empty when the pin is ignored).

set(_da_pin_file "${CMAKE_CURRENT_LIST_DIR}/symengine_pin.txt")
file(STRINGS "${_da_pin_file}" _da_pin REGEX "^[A-Za-z0-9_]+=")
foreach(_line IN LISTS _da_pin)
    string(REGEX REPLACE "=.*" "" _key "${_line}")
    string(REGEX REPLACE "^[^=]*=" "" _val "${_line}")
    set(_da_pin_${_key} "${_val}")
endforeach()
set(DA_SYMENGINE_VERSION "${_da_pin_SYMENGINE_VERSION}")

set(_da_pin_fix "Build the pinned SymEngine and point CMake at it:\n  eval \"$(scripts/setup_symengine.sh --print-env)\"\nthen reconfigure (with a fresh build directory or -DSymEngine_DIR=\$SymEngine_DIR).")

if(DA_IGNORE_SYMENGINE_PIN)
    set(_da_pin_severity WARNING)
    set(DA_SYMENGINE_FIND_VERSION "")
    find_package(SymEngine REQUIRED CONFIG)
else()
    set(_da_pin_severity FATAL_ERROR)
    set(DA_SYMENGINE_FIND_VERSION "${DA_SYMENGINE_VERSION} EXACT")
    find_package(SymEngine ${DA_SYMENGINE_VERSION} EXACT CONFIG)
    if(NOT SymEngine_FOUND)
        set(_da_found "")
        foreach(_cfg IN LISTS SymEngine_CONSIDERED_CONFIGS)
            list(FIND SymEngine_CONSIDERED_CONFIGS "${_cfg}" _idx)
            list(GET SymEngine_CONSIDERED_VERSIONS ${_idx} _ver)
            get_filename_component(_cfg "${_cfg}" DIRECTORY)
            list(APPEND _da_found "  ${_cfg} (version ${_ver})\n")
        endforeach()
        list(REMOVE_DUPLICATES _da_found)
        string(REPLACE ";" "" _da_found "${_da_found}")
        message(FATAL_ERROR
            "SymEngine ${DA_SYMENGINE_VERSION} (commit ${_da_pin_SYMENGINE_COMMIT}) not found. "
            "Found instead:\n${_da_found}\n${_da_pin_fix}")
    endif()
endif()

get_filename_component(DA_SYMENGINE_PREFIX "${SymEngine_DIR}/../../.." REALPATH)
set(_da_stamp "${DA_SYMENGINE_PREFIX}/share/symengine/miradac-pin.txt")
set(DA_SYMENGINE_LIB_SHA256 "")

set(_da_pin_problems "")
if(NOT EXISTS "${_da_stamp}")
    set(_da_pin_problems "no stamp file ${_da_stamp}")
else()
    file(STRINGS "${_da_stamp}" _da_stamp_lines REGEX "^[A-Za-z0-9_]+=")
    foreach(_line IN LISTS _da_pin)
        string(REGEX REPLACE "=.*" "" _key "${_line}")
        # Download hashes and the symengine.py version are not build properties.
        if(_key MATCHES "^SYMENGINE_(TARBALL_SHA256|PY_VERSION|PY_SDIST_SHA256)$")
            continue()
        endif()
        list(FIND _da_stamp_lines "${_line}" _idx)
        if(_idx EQUAL -1)
            string(APPEND _da_pin_problems "stamp does not have ${_line}\n")
        endif()
    endforeach()
    foreach(_line IN LISTS _da_stamp_lines)
        if(_line MATCHES "^LIB_SHA256=(.*)")
            set(DA_SYMENGINE_LIB_SHA256 "${CMAKE_MATCH_1}")
        endif()
    endforeach()
endif()

if(_da_pin_problems)
    message(${_da_pin_severity}
        "SymEngine at ${SymEngine_DIR} is not the pinned build "
        "(commit ${_da_pin_SYMENGINE_COMMIT}):\n${_da_pin_problems}\n${_da_pin_fix}")
else()
    message(STATUS "SymEngine pin OK: ${_da_pin_SYMENGINE_COMMIT} in ${DA_SYMENGINE_PREFIX}")
endif()
