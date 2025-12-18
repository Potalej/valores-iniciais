# ============================================================================
# Arquivo CMake para o gerador de valores iniciais
# ============================================================================
cmake_minimum_required(VERSION 3.15)
project(ncorpos_vi 
        VERSION 0.2.3
        DESCRIPTION "Geracao de valores iniciais"
        LANGUAGES Fortran)
enable_language(Fortran)

# Para baixar as coisas
include(FetchContent)
include(repos.cmake)

# ============================================================================
# Biblioteca
# ----------------------------------------------------------------------------
set(VI_DIR "")
set(VI_CURRENT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
include(libvi.cmake)