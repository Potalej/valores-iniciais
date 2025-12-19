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
add_subdirectory(libs/ncorpos_utilidades)

target_link_libraries(valores_iniciais
  PUBLIC utilidades
)