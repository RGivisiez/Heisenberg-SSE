Module Variables

  Integer*4 :: lx, ly, lz, d, L
  Integer*4 :: Nb, N, NH, NH_max = 0
  Integer*4 :: Nbins, mcsteps, term_steps

  Real*8 :: beta, ener, n_opH = 0.0d0

  Integer*4, Allocatable, Dimension (:) :: spin
  Integer*4, Allocatable, Dimension (:) :: opstring
  Integer*4, Allocatable, Dimension (:,:) :: bound
  Integer*4, Allocatable, Dimension (:) :: first_vertex_visitted
  Integer*4, Allocatable, Dimension (:) :: last_vertex_visitted
  Integer*4, Allocatable, Dimension (:) :: vertex_link

End Module Variables

Program Main

  Use Variables
  Implicit None

  Integer*4 seed(33), i, j

  seed(:33) = 42

  Call random_seed(put=seed)

  beta = 16.0d0  ! Inverse temperature.

  lx = 16 ! Number of spins in x.
  ly = 16  ! Number of spins in y.
  lz = 1  ! Number of spins in z.

  NH = 0  ! Number of operator different of the Identity one. 

  N = lx * ly ! Total number of spins.

  L = max(4, N / 4)  ! Truncated system size.

  Nbins = 20

  term_steps = 1e4

  mcsteps = 1e4

  Allocate(opstring(0:L-1))
  opstring = 0

  Allocate(first_vertex_visitted(N), last_vertex_visitted(N))
  Allocate(vertex_link(0:4 * L - 1))

  Call lattice
  Call init

  Do i = 1, term_steps

    Call diagonalupdate
    Call linkvertices
    Call loopupdate
    Call adjustcutoff(i)

  end do

  ! Close(30)

  Do j = 1, Nbins

    Do i = 1, mcsteps

      Call diagonalupdate
      Call linkvertices
      Call loopupdate
      Call measure

    end do

    Call write_results

  end do 

  Close(10)
  Close(20)

  Call results

  Call free_memory

End Program Main

Subroutine lattice

  Use Variables
  Implicit None

  Integer*4 :: i, j, bound_idx, spin_idx

  ! mod(x, L) returns 0 if x == Lx or returns x+1 if x != Lx

  if ( lx > 1 .and. ly > 1) then

    print*, 'Sistema fechado 2D'

    d = 2                     ! System's dimension
    Nb = d * lx*ly*lz         ! Number o bonds between spins.
    Allocate(bound(2, Nb))

    Do i = 0, lx - 1
      Do j = 0, ly - 1
        spin_idx = 1 + i + j * lx

        bound_idx = 1 + i + j * lx

        bound(1, bound_idx) = spin_idx
        bound(2, bound_idx) = 1 + mod(i + 1, lx) + j * lx

        bound(1, bound_idx + N) = spin_idx
        bound(2, bound_idx + N) = 1 + i + mod(j + 1, ly) * lx
      end do
    end do

  elseif ( (lx > 1) .and. (ly == 1 .and. lz == 1) ) then

    print*, 'Sistema fechado 1D'

    d = 1                     ! System's dimension
    Nb = d * lx*ly*lz         ! Number o bonds between spins.
    Allocate(bound(2, Nb))

    Do i = 0, lx - 1
       spin_idx = i + 1

       bound_idx = i + 1

       bound(1, bound_idx) = spin_idx
       bound(2, bound_idx) = 1 + mod(i + 1, lx)
    end do

  end if

End Subroutine lattice

Subroutine init

  Use Variables
  Implicit None

  Integer*4 :: i

  Allocate(spin(N))

  ! Spins with random values of -1 or 1.

  Do i = 1, N
    spin(i) = (-1)**(mod(i - 1, lx) + (i - 1) / lx)
  End do

  Open(10, file='results_bin.dat')
  Open(20, file='results_raw.dat')

End Subroutine init

Subroutine diagonalupdate()

  Use Variables
  Implicit None

  Integer*4 :: p, bound_idx, op
  Real*8 :: ran

  ! opstring = 2 * b(p) + a(p) - 1
  ! Como 2 * b(p) sempre é um número inteiro par,
  ! é o (a(p) - 1) que define qual tipo de operador temos.
  ! Quando a(p) = 1, opstring vai ser necessariamente
  ! par. Quando a(p) = 2, opstring é ímpar.
  ! Os valores de a(p) e b(p) pode ser obtidos a partir
  ! do opstring usando as equações,
  ! a(p)=MOD(opstring[p],2)+1, e
  ! b(p)=opstring[p]/2.

  Do p = 0, L - 1

    op = opstring(p)

    if (op == 0) then ! 0 == Operador Identidade

      ! Adiciona um operador diagonal,
      ! se os spins da ligação não estão em
      ! paralelo. 

      Call random_number(ran)

      bound_idx = int(ran * Nb) + 1

      if (spin(bound(1, bound_idx)) /= spin(bound(2, bound_idx))) then

        Call random_number(ran)

        if (ran <= (0.5d0 * beta * Nb) / (L - NH)) then

          opstring(p) = 2 * bound_idx
          NH = NH + 1
        end if

      end if

    elseif (mod(op, 2) == 0) then ! Par == operador diagonal

      ! Remove um operador diagonal,
      
      Call random_number(ran)

      if (ran <= (L - NH + 1) / (0.5d0 * beta * Nb)) then

        opstring(p) = 0
        NH = NH - 1
      end if

    else !Ímpar == Operador fora da diagonal

      ! Aqui não mudamos o número de operadores,
      ! só os spins são flipados.
      ! Adicionar ou remover operadores fora da diagonal
      ! é mais complicado, sendo necessário uma atenção
      ! maior. 

      ! Esse passo é permitido pq se existe um operador
      ! fora diagonal em p=1, vai existir um outro operador
      ! fora diagonal em p>1. Então é garantido que vamos flipar
      ! dois operadores fora da diagonal. Se isso não acontecer
      ! é pq estamos em uma configuração invalida para o sistema.

      bound_idx = op / 2

      spin(bound(1, bound_idx)) = -spin(bound(1, bound_idx))
      spin(bound(2, bound_idx)) = -spin(bound(2, bound_idx))

    end if

  end do

End Subroutine diagonalupdate

Subroutine linkvertices

  Use Variables
  Implicit None

  Integer*4 :: v0, p, op, bound_idx, first_vertex, last_vertex
  Integer*4 :: last_vertex1, last_vertex2
  Integer*4 :: spin_idx, spin_idx1, spin_idx2

  !Trocar boun por spin_idx_of_bound

  first_vertex_visitted(:) = -1
  last_vertex_visitted(:) = -1

  Do v0 = 0, 4 * L - 1, 4
    
    p = v0 / 4

    op = opstring(p)
    
    if (op /= 0) then
       
       bound_idx = op / 2
       
       spin_idx1 = bound(1, bound_idx)
       spin_idx2 = bound(2, bound_idx)
       
       last_vertex1 = last_vertex_visitted(spin_idx1)
       last_vertex2 = last_vertex_visitted(spin_idx2)
       
       if (last_vertex1 /= -1) then

          vertex_link(last_vertex1) = v0
          vertex_link(v0) = last_vertex1

       else

          first_vertex_visitted(spin_idx1) = v0

       endif

       if (last_vertex2 /= -1) then

          vertex_link(last_vertex2) = v0 + 1
          vertex_link(v0 + 1) = last_vertex2

       else

          first_vertex_visitted(spin_idx2) = v0 + 1

       endif

       last_vertex_visitted(spin_idx1) = v0 + 2
       last_vertex_visitted(spin_idx2) = v0 + 3

    else

       vertex_link(v0:v0 + 3) = -1

    endif
  enddo

  Do spin_idx = 1, N
    
    first_vertex = first_vertex_visitted(spin_idx)
    
    if (first_vertex /= -1) then

        last_vertex = last_vertex_visitted(spin_idx)

        vertex_link(last_vertex) = first_vertex
        vertex_link(first_vertex) = last_vertex

    endif
  enddo
  
End Subroutine linkvertices

Subroutine loopupdate

  Use Variables
  Implicit None

  Real*8 :: ran
  Integer*4 :: spin_idx, v0, vertex_in, vertex_out

  ! -1 para vertex visitados porém não flipados.
  ! -2 para vertex visitados e flipados.

  Do v0 = 0, 4 * L - 1, 2

    ! Se o vertex já foi visitado, pula.
    if (vertex_link(v0) < 0) cycle

    vertex_in = v0

    Call random_number(ran)

    if (ran < 0.5d0) then

      do 

        ! flip
        opstring(vertex_in / 4) = ieor(opstring(vertex_in / 4), 1)
        
        ! marked as visitted
        vertex_link(vertex_in) = -2
        
        ! next vertex
        vertex_out = ieor(vertex_in, 1)

        vertex_in = vertex_link(vertex_out)

        vertex_link(vertex_out) = -2
        
        ! Se voltar onde começou o loop está completo.
        if (vertex_in == v0) exit

      end do

    else

      do 

        vertex_link(vertex_in) = -1

        vertex_out = ieor(vertex_in, 1)

        vertex_in = vertex_link(vertex_out)

        vertex_link(vertex_out) = -1

        if (vertex_in == v0) exit

      end do

    end if

  end do

  ! flip das linhas de spins sem operadores.

  Do spin_idx = 1, N

    if (first_vertex_visitted(spin_idx) /= -1) then

      if (vertex_link(first_vertex_visitted(spin_idx)) == -2) then
        spin(spin_idx) = -spin(spin_idx)
      end if

    else
      
      Call random_number(ran)

      if (ran < 0.5d0) spin(spin_idx) = -spin(spin_idx)

    end if

  end do 

End Subroutine loopupdate

Subroutine adjustcutoff(step)

  Use Variables
  Implicit None

  Integer*4, Allocatable, Dimension (:) :: copy_opstring
  Integer*4 :: L_new, step

  L_new = NH + NH / 3

  ! Open(30, file='adjustcutoff.dat', position='append')

  ! if ( NH > NH_max) NH_max = NH

  ! if (L_new <= L) then

  !   write(30, *) step, L, NH, NH_max
  !   return

  ! else 

  !   write(30, *) step, L_new, NH, NH_max

  ! end if

  if (L_new <= L) return

  Allocate(copy_opstring(0:L-1))
  copy_opstring(:)=opstring(:)

  Deallocate(opstring)
  Allocate(opstring(0:L_new-1))
  opstring(0:L-1)=copy_opstring(:)
  opstring(L:L_new-1)=0
  Deallocate(copy_opstring)

  L=L_new
  Deallocate (vertex_link)
  Allocate(vertex_link(0:4*L-1))

End Subroutine adjustcutoff

Subroutine measure

  Use Variables
  Implicit None

  n_opH = n_opH + dble(NH)

  write(20, *) NH

End Subroutine measure

Subroutine write_results

  Use Variables
  Implicit None

  n_opH = n_opH / dble(mcsteps)

  ener = - ( n_opH / (beta * N) - 0.25d0 * dble(Nb) / dble(N))

  write(10, *) ener, n_opH

  n_opH = 0.0d0

End Subroutine write_results

Subroutine results

  Use Variables
  Implicit None

  Real*8 :: r1(Nbins)
  Integer*4 :: i

  Open(10, file='ener.dat')

  Do i = 1, Nbins

    read(10, *) r1(i)

  end do 

  Close(10)

  Open(20, file='results.dat', position='append')

  write(20, *) beta, sum(r1) / dble(Nbins)

  Close(20)

End Subroutine results

Subroutine free_memory

  Use Variables
  Implicit None

  Deallocate(spin)
  Deallocate(opstring)
  Deallocate(bound)
  Deallocate(first_vertex_visitted)
  Deallocate(last_vertex_visitted)
  Deallocate(vertex_link)

End Subroutine free_memory