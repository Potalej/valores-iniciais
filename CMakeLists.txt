# ============================================================================
# Arquivo CMake para o gerador de valores iniciais
# ============================================================================
cmake_minimum_required(VERSION 3.15)
project(ncorpos_vi
        VERSION 0.2.3
        DESCRIPTION "Geracao de valores iniciais"
        LANGUAGES Fortran)
enable_language(Fortran)

# ============================================================================
# Biblioteca
# ----------------------------------------------------------------------------
set(VI_UTILIDADES_LOCAL ON)
set(VI_CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
include(libvi.cmake)