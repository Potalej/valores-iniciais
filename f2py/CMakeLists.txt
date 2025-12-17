# ============================================================================
# Projeto
# ============================================================================
cmake_minimum_required(VERSION 3.15)

project(valores_iniciais
  VERSION 0.1.0
  DESCRIPTION "Geracao de valores iniciais"
  LANGUAGES Fortran
)

enable_language(Fortran)

include(FetchContent)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

# ============================================================================
# Biblioteca Fortran (core)
# ============================================================================
add_library(valores_iniciais SHARED)

set_target_properties(valores_iniciais PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
  POSITION_INDEPENDENT_CODE ON
)

target_sources(valores_iniciais
  PRIVATE
    ../include/tipos.F90
    ../src/aleatorio.f90
    ../src/condicionamento.f90
    ../src/condicoes_iniciais.f90
    ../src/sorteio.f90
)

target_include_directories(valores_iniciais
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/mod
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ============================================================================
# Dependencia: utilidades
# ============================================================================
if (NOT TARGET utilidades)
  message(STATUS "Biblioteca 'utilidades' nao encontrada! baixando...")
  FetchContent_Declare(
    utilidades
    GIT_REPOSITORY https://github.com/potalej/utilidades.git
    GIT_TAG main
  )
  FetchContent_MakeAvailable(utilidades)
else()
  message(STATUS "Usando 'utilidades' do projeto principal")
endif()
# add_subdirectory(utilidades/)
target_link_libraries(valores_iniciais
  PUBLIC utilidades
)

# ============================================================================
# Wrapper f2py
# ============================================================================
set(F2PY_WRAPPER
  ${CMAKE_CURRENT_SOURCE_DIR}/valores_iniciais.f90
)

set(F2PY_MODULE_NAME _core)

# Descobrir sufixo do módulo Python (.so)
execute_process(
  COMMAND ${Python3_EXECUTABLE}
          -c
          "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
  OUTPUT_VARIABLE PY_EXT_SUFFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(F2PY_OUTPUT
  ${CMAKE_CURRENT_BINARY_DIR}/${F2PY_MODULE_NAME}${PY_EXT_SUFFIX}
)

# add_custom_command(
#   OUTPUT ${F2PY_OUTPUT}
#   COMMAND ${Python3_EXECUTABLE} -m numpy.f2py
#           -c
#           -m ${F2PY_MODULE_NAME}
#           ${F2PY_WRAPPER}
#           -L$<TARGET_FILE_DIR:valores_iniciais>
#           -lvalores_iniciais
#           -L$<TARGET_FILE_DIR:utilidades>
#           -lutilidades
#           -I${CMAKE_CURRENT_BINARY_DIR}/mod
#           --quiet
#   DEPENDS
#     valores_iniciais
#     ${F2PY_WRAPPER}
#   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
#   COMMENT "Gerando modulo Python via f2py"
# )

add_custom_command(
  OUTPUT ${F2PY_OUTPUT}
  COMMAND
    ${CMAKE_COMMAND} -E env
      LDFLAGS=-Wl,-rpath,'$$ORIGIN'
    ${Python3_EXECUTABLE} -m numpy.f2py
      -c
      -m ${F2PY_MODULE_NAME}
      ${F2PY_WRAPPER}
      -L$<TARGET_FILE_DIR:valores_iniciais>
      -lvalores_iniciais
      -L$<TARGET_FILE_DIR:utilidades>
      -lutilidades
      -I${CMAKE_CURRENT_BINARY_DIR}/mod
      --quiet
  DEPENDS
    valores_iniciais
    utilidades
    ${F2PY_WRAPPER}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Gerando modulo Python via f2py com RPATH=$ORIGIN"
)

add_custom_target(python_module ALL
  DEPENDS ${F2PY_OUTPUT}
)

# ============================================================================
# Mensagem final
# ============================================================================
message(STATUS "Modulo Python sera gerado em:")
message(STATUS "  ${F2PY_OUTPUT}")

# ============================================================================
# Pasta do pacote Python
# ============================================================================
set(PYTHON_PACKAGE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/../python/valores_iniciais
)

# ============================================================================
# Arquivos a serem copiados
# ============================================================================
set(PYTHON_CORE_SO
  ${PYTHON_PACKAGE_DIR}/_core.so
)

set(PYTHON_VALORES_SO
  ${PYTHON_PACKAGE_DIR}/libvalores_iniciais.so
)

set(PYTHON_UTILIDADES_SO
  ${PYTHON_PACKAGE_DIR}/libutilidades.so
)

# ============================================================================
# Comando de cópia
# ============================================================================
add_custom_command(
  OUTPUT
    ${PYTHON_CORE_SO}
    ${PYTHON_VALORES_SO}
    ${PYTHON_UTILIDADES_SO}

  COMMAND ${CMAKE_COMMAND} -E make_directory ${PYTHON_PACKAGE_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy
          ${F2PY_OUTPUT}
          ${PYTHON_CORE_SO}

  COMMAND ${CMAKE_COMMAND} -E copy
          $<TARGET_FILE:valores_iniciais>
          ${PYTHON_VALORES_SO}

  COMMAND ${CMAKE_COMMAND} -E copy
          $<TARGET_FILE:utilidades>
          ${PYTHON_UTILIDADES_SO}

  DEPENDS
    ${F2PY_OUTPUT}
    valores_iniciais
    utilidades

  COMMENT "Copiando _core.so, libvalores_iniciais.so e libutilidades.so para o pacote Python"
)

# ============================================================================
# Target agregador
# ============================================================================
add_custom_target(install_python_module ALL
  DEPENDS
    ${PYTHON_CORE_SO}
    ${PYTHON_VALORES_SO}
    ${PYTHON_UTILIDADES_SO}
)