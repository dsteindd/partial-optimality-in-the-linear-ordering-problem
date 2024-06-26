set(GUROBI_ROOT_DIR "" CACHE PATH "GUROBI root directory.")

STRING(REGEX MATCH "^[0-9]+" GUROBI_VERSION "$ENV{GUROBI_HOME}")

find_path(GUROBI_INCLUDE_DIR gurobi_c++.h HINTS "$ENV{GUROBI_HOME}/include")
find_library(GUROBI_LIBRARY libgurobi110.so HINTS $ENV{GUROBI_HOME}/lib)
find_library(GUROBI_CPP_LIBRARY libgurobi_c++.a HINTS $ENV{GUROBI_HOME}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_CPP_LIBRARY GUROBI_INCLUDE_DIR)

if(GUROBI_FOUND)
    set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
    set(GUROBI_LIBRARIES ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(GUROBI_LIBRARIES "${GUROBI_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(GUROBI_FOUND)

mark_as_advanced(GUROBI_LIBRARY GUROBI_CPP_LIBRARY GUROBI_INCLUDE_DIR)
