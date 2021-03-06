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
    integer :: i, iresi, j, jj, k, kk, tot_res, imove
    integer :: iatmst, iatmen
    integer,intent(in) :: narray
    integer,intent(in) :: first_atom(tot_res), last_atom(tot_res)
    integer :: nchain
    integer :: iland(nchain), ichst(nchain), ichen(nchain) 
    real(4), intent(inout) :: codx(narray), cody(narray), codz(narray)
    real(4) :: rdiff
    double precision, intent(in) :: cell(3) 

  !*****************************************************
  !Move atoms at each residues if atoms go beyond a box.
  !*****************************************************
    do iresi = 1, tot_res
      iatmst = first_atom(iresi)                    !First atom number of a residues 
      iatmen = last_atom(iresi)                     !Last atom number of a residues
      do jj = iatmst + 1, iatmen              !From next atom to last atom in a residue.
        rdiff = codx(iatmst) - codx(jj)       !calc. diff. betw. coordinate of first atom and jjth atom.
        imove = idnint(rdiff / cell(1))       !imove is 0 or 1. if 1, move back atom. if 0, remain position of atom. 
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
    !**For ET1
    iland(1) = 5
    iland(2) = 22
    ichst(1) = 1 
    ichst(2) = 21 
    ichen(1) = 20
    ichen(2) = 40
    
    !For p53
    !iland(1) = 5
    !ichst(1) = 1 
    !ichen(1) = 18

  
    do i = 1, nchain
      do iresi = iland(i)-1, ichst(i), -1
        iatmst = first_atom(iresi)
        iatmen = last_atom(iresi)
  
        !*************************************************************************************************
        !All atoms in a residue will be moved, 
        !if the distance btw 1st atom of residue(i+1) and residue(i) is over half of the length of x-cell.
        !*************************************************************************************************
        rdiff  = codx(first_atom(iresi+1)) - codx(first_atom(iresi)) !If neiboring res. is far to iresi+1 res., rdiff is large. 
        imove  = idnint(rdiff / cell(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove * cell(1))
        enddo
  
        rdiff  = cody(first_atom(iresi + 1)) - cody(first_atom(iresi))
        imove  = idnint(rdiff / cell(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove * cell(2))
        enddo
        
        rdiff  = codz(first_atom(iresi + 1)) - codz(first_atom(iresi))
        imove  = idnint(rdiff / cell(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove * cell(3))
        enddo
  
      enddo
    enddo

    do i = 1, nchain
      do iresi = iland(i)+1, ichen(i)
        iatmst = first_atom(iresi)
        iatmen = last_atom(iresi)
  
        rdiff = codx(first_atom(iresi-1)) - codx(first_atom(iresi))
        imove = idnint(rdiff / cell(1))
        do jj = iatmst, iatmen
          codx(jj) = codx(jj) + (imove * cell(1))
        enddo
  
        rdiff = cody(first_atom(iresi-1)) - cody(first_atom(iresi))
        imove = idnint(rdiff / cell(2))
        do jj = iatmst, iatmen
          cody(jj) = cody(jj) + (imove * cell(2))
        enddo
  
        rdiff = codz(first_atom(iresi-1)) - codz(first_atom(iresi))
        imove = idnint(rdiff / cell(3))
        do jj = iatmst, iatmen
          codz(jj) = codz(jj) + (imove * cell(3))
        enddo
      enddo

    enddo
    
    write(*,'(a)') "Note) Atoms have been returned (module periodic was uesd)."
  end subroutine
end module
