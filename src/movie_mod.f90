module movie
  implicit none
contains
  
  subroutine make_movie(TotAtoms, Ref_AtmNum, Ref_AtmName, Ref_ResName, Ref_ResNum, x_trj, y_trj, z_trj)
  !This routine constructs trajectory movie.
  !Explanation of variables:
  !  1st(INPUT): Total atoms (integer)
  !  2nd(INPUT): atom number (1D array) of reference
  !  3rd(INPUT): atom name (1D array) of reference (ex. CA, H, O,...)
  !  4th(INPUT): residue name (1D array) of reference 
  !  5th(INPUT): residue number (1D array) of reference
  !  6th(INPUT): x-coordinate (1D array) of trajectory
  !  7th(INPUT): y-coordinate (1D array) of trajectory
  !  8th(INPUT): z-coordinate (1D array) of trajectory
    integer :: j
    integer :: TotAtoms
    integer, allocatable :: AtomNum(:)
    character(len=4), allocatable :: AtomName(:)
    character(len=3), allocatable :: ResName(:)
    character(len=1), allocatable :: ChainId(:)
    integer, allocatable :: ResNum(:)
    double precision, allocatable :: x(:),y(:),z(:)


      write(1122,'("MODEL ", i7)') i
      write(1122,'("REMARK ENERGY [kcal/mol] ", f10.3)') potential(i)
      do j = 1, TotAtoms
        write(1122,'(3f10.3)') &



              trj_x(j), trj_y(j), trj_z(j)
      enddo
      write(1122,"(a6)") "ENDMDL"
  
  end subroutine


end module
