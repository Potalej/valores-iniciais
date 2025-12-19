# ============================================================================
# Arquivo CMake para o gerador de valores iniciais
# ============================================================================
cmake_minimum_required(VERSION 3.15)
project(ncorpos_vi
        VERSION 0.2.3
        DESCRIPTION "Geracao de valores iniciais"
        LANGUAGES Fortran)
enable_language(Fortran)
set(NCORPOS_VI_PROJECT_VERSION ${PROJECT_VERSION} CACHE INTERNAL "")

# ============================================================================
# Biblioteca
# ----------------------------------------------------------------------------
if (NOT EXISTS "${VI_UTILIDADES_LOCAL}")
        set(VI_UTILIDADES_LOCAL ON)
endif()
set(VI_CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
include(libvi.cmake)