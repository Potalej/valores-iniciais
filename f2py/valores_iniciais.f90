module valores_iniciais
  use iso_c_binding
  use sorteio_mod, only: sorteio
  USE tipos
  use condicoes_iniciais
  IMPLICIT NONE
  
  PRIVATE
  type(sorteio), TARGET, SAVE :: sorteio_global

  PUBLIC parametros, parametros_massas, parametros_posicoes, parametros_momentos, gerar, &
         massas_c, posicoes_c, momentos_c
  
  REAL(pf), ALLOCATABLE :: massas(:), posicoes(:,:), momentos(:,:)
  REAL(c_double), ALLOCATABLE :: massas_c(:), posicoes_c(:,:), momentos_c(:,:)

contains

  SUBROUTINE parametros (N, G, eps, modo, ed, jd, pd)
    INTEGER, INTENT(IN) :: N
    REAL(8), INTENT(IN) :: G, eps, ed
    CHARACTER(LEN=*), INTENT(IN) :: modo
    REAL(8), INTENT(IN), OPTIONAL :: jd(3), pd(3)
    REAL(8) :: angular(3), linear(3)

    angular = 0.0
    linear = 0.0
    IF (PRESENT(jd)) angular = jd
    IF (PRESENT(pd)) linear = pd

    sorteio_global % N = N
    sorteio_global % G = G
    sorteio_global % amortecedor = eps
    sorteio_global % modo = modo
    sorteio_global % ed = ed
    sorteio_global % jd = angular
    sorteio_global % pd = linear

    IF (ALLOCATED(massas_c)) DEALLOCATE(massas_c)
    IF (ALLOCATED(posicoes_c)) DEALLOCATE(posicoes_c)
    IF (ALLOCATED(momentos_c)) DEALLOCATE(momentos_c)

    ALLOCATE(massas_c(N))
    ALLOCATE(posicoes_c(N,3))
    ALLOCATE(momentos_c(N,3))

    massas_c = 0.0
    posicoes_c = 0.0
    momentos_c = 0.0
  END SUBROUTINE

  SUBROUTINE parametros_massas (distribuicao, regiao, intervalo, normalizado)
    CHARACTER(LEN=*), INTENT(IN) :: distribuicao, regiao
    REAL(8), INTENT(IN) :: intervalo(2)
    LOGICAL, INTENT(IN), OPTIONAL :: normalizado

    sorteio_global % massas % distribuicao = distribuicao
    sorteio_global % massas % regiao = regiao
    sorteio_global % massas % intervalo = intervalo
    sorteio_global % massas % normalizado = .FALSE.
    IF (PRESENT(normalizado)) sorteio_global % massas % normalizado = normalizado
  END SUBROUTINE

  SUBROUTINE parametros_posicoes (distribuicao, regiao, intervalo)
    CHARACTER(LEN=*), INTENT(IN) :: distribuicao, regiao
    REAL(8), INTENT(IN) :: intervalo(2)
    
    sorteio_global % posicoes % distribuicao = distribuicao
    sorteio_global % posicoes % regiao = regiao
    sorteio_global % posicoes % intervalo = intervalo
  END SUBROUTINE

  SUBROUTINE parametros_momentos (distribuicao, regiao, intervalo)
    CHARACTER(LEN=*), INTENT(IN) :: distribuicao, regiao
    REAL(8), INTENT(IN) :: intervalo(2)
    
    sorteio_global % momentos % distribuicao = distribuicao
    sorteio_global % momentos % regiao = regiao
    sorteio_global % momentos % intervalo = intervalo
  END SUBROUTINE

  SUBROUTINE gerar ()
    type(sorteio), pointer :: sorteio_local

    sorteio_local => sorteio_global

    ! Primeiro, aloca os valores
    IF (ALLOCATED(massas)) DEALLOCATE(massas)
    IF (ALLOCATED(posicoes)) DEALLOCATE(posicoes)
    IF (ALLOCATED(momentos)) DEALLOCATE(momentos)
    
    CALL condicionar(sorteio_local, massas, posicoes, momentos)

    massas_c = massas
    posicoes_c = posicoes
    momentos_c = momentos
  END SUBROUTINE

  FUNCTION get_massas ()
    REAL(c_double) :: get_massas(SIZE(massas))
    get_massas = massas
  END FUNCTION

  SUBROUTINE get_posicoes (posicoes)
    REAL(pf), INTENT(OUT), ALLOCATABLE:: posicoes(:,:)
  END SUBROUTINE

  SUBROUTINE get_momentos (momentos)
    REAL(pf), INTENT(OUT), ALLOCATABLE :: momentos(:,:)
  END SUBROUTINE

  ! subroutine sorteio_set_G(G) bind(C)
  !   real(c_double), value :: G
  !   sorteio_global%G = G
  ! end subroutine

  ! subroutine sorteio_run(massas, pos, mom)
  !   USE ico_c_binding
  !   real(c_double), allocatable, intent(out) :: massas(:)
  !   real(c_double), allocatable, intent(out) :: pos(:,:)
  !   real(c_double), allocatable, intent(out) :: mom(:,:)

  !   call condicionar(sorteio_global, massas, pos, mom)
  ! end subroutine

end module
