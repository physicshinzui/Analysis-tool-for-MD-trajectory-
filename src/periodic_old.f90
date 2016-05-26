!***************************************************************
!This program is needed to obtain correct protein conformations
!when periodic boundary condition is used in simulation.
!
!Algorithm:
!If the distance between adjacent atom is larger than a box side, 
!the atom goes beyond the box side.Thus the atom is needed to return an appropiate position.
!
!README: 
!tot_atms is the number of total atoms of a system.
!x, y and z are coordinate of a system.
!xcell, ycell and zcell are box size of a system.
!narray is array size of these
! Memo:
! rdiff / xcell < 1 absolutely
!***************************************************************
module periodic 
  implicit none
contains

  subroutine DetectResidueNumber(total_atom, residue_number, total_residue,first_atom,last_atom) 
    integer :: i  
    integer, intent(in)  :: total_atom
    integer, intent(in)  :: residue_number(total_atom)
    integer, intent(in)  :: total_residue
    integer, intent(out) :: first_atom(:)
    integer, intent(out) :: last_atom(:)

    do i = 1, total_atom
      if( i == 1 ) then
        first_atom(1) = 1
      elseif( i /= 1 .and. residue_number(i) /= residue_number(i-1)) then
        first_atom(residue_number(i)) = i
        last_atom(residue_number(i-1)) = i - 1
      endif
    enddo
    last_atom(residue_number(i-1)) = i - 1

  end subroutine


  subroutine ReturnAtom(narray,nchain,first_atom,last_atom,tot_res,codx,cody,codz,cell)
    integer :: i, ii, j, jj, k, kk, tot_res, imove
    integer :: iatmst, iatmen
    integer,intent(in) :: narray
    integer,intent(in) :: first_atom(tot_res), last_atom(tot_res)
    integer :: nchain, iland(1), ichst(1), ichen(1) 
    real(4), intent(inout) :: codx(narray), cody(narray), codz(narray)
    real(4) :: rdiff
    double precision, intent(in) :: cell(3) 
!    print*,"cell",cell
  
    !print*,"first_atom,last_atom",lbound(first_atom),ubound(first_atom),lbound(last_atom),ubound(last_atom)
    !print*,"first_atom,last_atom",first_atom(1), last_atom(1)
    !print*,"codx",ubound(codx)
    !print*,"cody",ubound(cody)
    !print*,"codz",ubound(codz)
    !print*,"chain",nchain
    !print*,"codx",codx(1)
  
    !print*,"aa",narray
    !print*,"nchanin",nchain
    !print*,"fit atm",first_atom
    !print*,last_atom
    !print*,tot_res
    !print*,codx,cody,codz
    !print*,cell


  !*****************************************************
  !Move atoms at each residues if atoms go beyond a box.
  !*****************************************************
    do ii = 1, tot_res
      iatmst = first_atom(ii)                    !First atom number of a residues 
      iatmen = last_atom(ii)                     !Last atom number of a residues
      do jj = iatmst + 1, iatmen              !From next atom to last atom in a residue.
        rdiff = codx(iatmst) - codx(jj)       !calc. diff. betw. coordinate of first atom and jjth atom.
        imove = idnint(rdiff / cell(1))         !0 or 1. calc. how rdiff is large, comparing with a box side.
        codx(jj) = codx(jj) + (imove * cell(1)) !Update a coordinate of atom
      enddo
  
      do jj = iatmst + 1, iatmen        !From next atom to last atom in a residue.
        rdiff = cody(iatmst) - cody(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove = idnint(rdiff / cell(2))
        cody(jj) = cody(jj) + (imove * cell(2)) 
      enddo
  
      do jj = iatmst + 1, iatmen        !From next atom to last atom in a residue.
        rdiff = codz(iatmst) - codz(jj) !calc. diff. betw. coordinate of first atom and jjth atom.
        imove = idnint(rdiff / cell(3))
        codz(jj) = codz(jj) + (imove * cell(3)) 
      enddo
  
    enddo
  
  !*****************************************************
  !Move residue if a entire residue goes beyond a box.
  !*****************************************************
  
  !@@@@@@@@
    iland(1) = 9  !I still haven't understood why 9 is ok in iland. I think iland means intermediate of a chain.
    ichst(1) = 1 
    ichen(1) = 19
  
    do i = 1, nchain
      do ii = iland(i) - 1, ichst(i), -1
        iatmst = first_atom(ii)
        iatmen = last_atom(ii)
  
        !***All of atoms in a res. are moved***
        rdiff  = codx(first_atom(ii+1)) - codx(first_atom(ii)) !If neiboring res. is far to ii+1 res., rdiff is large. 
        imove  = idnint(rdiff / cell(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove * cell(1))
        enddo
  
        rdiff  = cody(first_atom(ii + 1)) - cody(first_atom(ii))
        imove  = idnint(rdiff / cell(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove * cell(2))
        enddo
        
        rdiff  = codz(first_atom(ii + 1)) - codz(first_atom(ii))
        imove  = idnint(rdiff / cell(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove * cell(3))
        enddo
  
      enddo
    enddo
  
    do i = 1, nchain
      do ii = iland(i)+1, ichen(i)
        iatmst = first_atom(ii)
        iatmen = last_atom(ii)
  
        rdiff = codx(first_atom(ii-1)) - codx(first_atom(ii))
        imove = idnint(rdiff / cell(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove * cell(1))
        enddo
  
        rdiff = cody(first_atom(ii-1)) - cody(first_atom(ii))
        imove = idnint(rdiff / cell(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove * cell(2))
        enddo
  
        rdiff = codz(first_atom(ii-1)) - codz(first_atom(ii))
        imove = idnint(rdiff / cell(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove * cell(3))
        enddo
  
      enddo
    enddo
    
    write(*,'(a)') "Note) Atoms have been returned (module periodic was uesd)."

  end subroutine
end module
