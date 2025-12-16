cmake_minimum_required(VERSION 3.15)
project(valores_iniciais LANGUAGES Fortran)

include(FetchContent)

# -----------------------------------------
# Dependência: utilidades
# -----------------------------------------
if (NOT TARGET utilidades)

  message(STATUS "utilidades not found — fetching")

  FetchContent_Declare(
    utilidades
    GIT_REPOSITORY https://github.com/potalej/utilidades.git
    GIT_TAG main   # ou tag / commit
  )

  FetchContent_MakeAvailable(utilidades)

else()
  message(STATUS "Using utilidades provided by parent project")
endif()



add_library(valores_iniciais STATIC)

option(USE_EXTERNAL_TIPOS
       "Use module tipos provided by the application"
       OFF)

set_target_properties(valores_iniciais PROPERTIES
  Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/mod
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

target_link_libraries(valores_iniciais
  PUBLIC utilidades
)

if (USE_EXTERNAL_TIPOS)
  target_compile_definitions(valores_iniciais PRIVATE HAVE_TIPOS)
endif()