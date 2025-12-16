! *****************************************************************
!! METODOS DE CONDICIONAMENTO DE VALORES INICIAIS
!
! Objetivos:
!   Funcoes para o condicionamento de valores iniciais.
! 
! Modificado:
!   16 de dezembro de 2025
! 
! Autoria:
!   oap
! 
MODULE condicoes_iniciais
  USE tipos
  USE utilidades
  USE condicionamento ! Metodos gerais de restricao
  USE sorteio_mod
  IMPLICIT NONE
CONTAINS

! ************************************************************
!! Condicionamento geral
!
! Objetivos:
!   Subrotina geral para condicionamento. Todas as outras sao
!   chamadas a partir desta.
!
! Modificado:
!   16 de dezembro de 2025
!
! Autoria:
!   oap
!
SUBROUTINE condicionar (sorteio_infos, massas, posicoes, momentos)
  TYPE(sorteio), POINTER, INTENT(INOUT) :: sorteio_infos
  REAL(pf), INTENT(INOUT), ALLOCATABLE  :: massas(:), posicoes(:,:), momentos(:,:) ! out
  INTEGER  :: N
  REAL(pf) :: G, eps, ed
  REAL(pf), DIMENSION(:), ALLOCATABLE :: pd, jd

  N = sorteio_infos % N
  G = sorteio_infos % G
  eps = sorteio_infos % amortecedor

  WRITE(*,'(a)') "GERACAO DAS VALORES INICIAIS ("//sorteio_infos % modo//")"

  ! Sorteio dos valores na regiao e com distribuicao desejadas
  WRITE(*,'(a)') "  > gerando valores..."
  ALLOCATE(massas(N))
  ALLOCATE(posicoes(N,3))
  ALLOCATE(momentos(N,3))
  CALL gerar_valores(N, sorteio_infos, massas, posicoes, momentos)

  ! Captura as integrais primeiras
  ed = sorteio_infos % ed
  pd = sorteio_infos % jd
  jd = sorteio_infos % pd

  ! Agora faz o condicionamento de acordo com o metodo desejado
  WRITE(*,*) ""
  WRITE(*,'(a)') "  > condicionando..."
  
  SELECT CASE (TRIM(sorteio_infos % modo))
    ! Sem nenhum condicionamento
    CASE("sorteio")
      ! Apenas nao faz nada

    ! Condicionamento por integrais primeiras iterativamente
    CASE("sorteio_ip_iterativo")
      CALL condicionar_ip_iterativo(G, massas, posicoes, momentos, eps, ed, jd, pd, 50)

    ! Condicionamento por integrais primeiras diretamente
    CASE("sorteio_ip_direto")
      CALL condicionar_ip_direto(G, massas, posicoes, momentos, eps, ed, jd, pd)

    ! Condicionamento de Aarseth
    CASE("sorteio_aarseth")
      CALL condicionar_aarseth(G, massas, posicoes, momentos, eps)

    ! Condicionamento de Aarseth Modificado
    CASE("sorteio_aarseth_modificado")
      CALL condicionar_aarseth_modificado(G, massas, posicoes, momentos, eps, jd)

    ! Outro caso: erro
    CASE DEFAULT
      ERROR STOP "!! Modo de condicionamento desconhecido. !!"
  END SELECT

  ! Outputs com informacoes sobre os dados iniciais sorteados
  CALL condicionamento_outputs(G, massas, posicoes, momentos, eps)
END SUBROUTINE

SUBROUTINE condicionamento_outputs (G, massas, posicoes, momentos, eps)
  REAL(pf), INTENT(IN) :: G, massas(:), posicoes(:,:), momentos(:,:), eps
  
  ! informacoes para exibir
  REAL(pf) :: potencial, cinetica
  REAL(pf) :: lintot(3), angtot(3)
  REAL(pf) :: inercia, dilatacao
  REAL(pf) :: anitenine ! anisotropia do tensor de inercia
  REAL(pf) :: f_prod_q, virial
  INTEGER  :: a
  REAL(pf) :: dist, mais_distante

  cinetica  = energia_cinetica(massas, momentos)
  
  virial = cinetica + cinetica
  IF (eps == 0) THEN
    potencial = energia_potencial(G, massas, posicoes, eps)
    virial = virial + potencial
  ELSE
    f_prod_q = virial_potencial_amortecido(G, massas, posicoes, eps, potencial)
    virial = virial + f_prod_q
  ENDIF

  lintot    = momento_linear_total(momentos)
  angtot    = momento_angular_total(posicoes, momentos)
  
  inercia   = momento_inercia(massas, posicoes)
  dilatacao = momento_dilatacao(posicoes, momentos)

  anitenine = anisotropia_tensor_inercia(massas, posicoes)

  mais_distante = 0
  DO a = 1, SIZE(massas)
    dist = NORM2(posicoes(a,:))
    IF (dist > mais_distante) mais_distante = dist
  END DO

  WRITE(*,*) ""
  WRITE(*,'(a)') "  > valores iniciais condicionados!"
  WRITE(*,*) "     * V   = ", potencial
  WRITE(*,*) "     * T   = ", cinetica
  WRITE(*,*) "     * E   = ", potencial + cinetica
  WRITE(*,*) "     * Vir = ", virial
  WRITE(*,*) "     * J   = ", angtot
  WRITE(*,*) "     * P   = ", lintot
  WRITE(*,*) "     * I   = ", inercia
  WRITE(*,*) "     * D   = ", dilatacao
  WRITE(*,*) "     * A_I = ", anitenine
  WRITE(*,*) "     * MD  = ", mais_distante

END SUBROUTINE

END MODULE condicoes_iniciais