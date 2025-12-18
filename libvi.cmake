add_library(valores_iniciais STATIC)

set_target_properties(valores_iniciais PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)

target_sources(valores_iniciais
  PRIVATE
    ${VI_DIR}include/tipos.F90
    ${VI_DIR}src/aleatorio.f90
    ${VI_DIR}src/condicionamento.f90
    ${VI_DIR}src/condicoes_iniciais.f90
    ${VI_DIR}src/sorteio.f90
)

target_include_directories(valores_iniciais
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/mod
    ${VI_CURRENT_DIR}/include
)

# ============================================================================
# Dependencia: utilidades
# ============================================================================
include(repos.cmake)

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