# ============================================================================
# Arquivo CMake para o gerador de valores iniciais
# ============================================================================
cmake_minimum_required(VERSION 3.15)
project(valores_iniciais 
        VERSION 0.2.2
        DESCRIPTION "Geracao de valores iniciais"
        LANGUAGES Fortran)
enable_language(Fortran)

# Para baixar as coisas
include(FetchContent)
include(repos.cmake)

# ============================================================================
# Biblioteca
# ----------------------------------------------------------------------------
add_library(valores_iniciais STATIC)

set_target_properties(valores_iniciais PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_sources(valores_iniciais
  PRIVATE
    src/aleatorio.f90
    src/condicionamento.f90
    src/condicoes_iniciais.f90
    src/sorteio.f90
    include/tipos.F90
)

target_include_directories(valores_iniciais
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/mod
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ============================================================================
# Dependencia: utilidades
# ----------------------------------------------------------------------------
if (NOT TARGET utilidades)
  message(STATUS "Biblioteca 'utilidades' nao encontrada! baixando do repositorio...")
  FetchContent_Declare(
    utilidades
    GIT_REPOSITORY ${UTILIDADES_REPO}
    GIT_TAG ${UTILIDADES_TAG}
  )
  FetchContent_MakeAvailable(utilidades)
else()
  message(STATUS "Usando a 'utilidades' do programa.")
endif()

target_link_libraries(valores_iniciais
  PUBLIC utilidades
)