module calculation 
  implicit none
contains
  !------------------------------
  subroutine Calc_CenterMass(n,x,y,z,CenterOfMass)
    !  1st(INPUT) : The N of protein atoms 
    integer :: i
    integer,intent(in) :: n !Num of atom
    double precision,intent(in)  :: x(n),y(n),z(n)
    double precision,intent(out) :: CenterOfMass(3)
    double precision :: SumX,SumY,SumZ 
    real(4) :: xtmp(n),ytmp(n),ztmp(n)
    SumX = sum(x) 
    SumY = sum(y) 
    SumZ = sum(z) 
    CenterOfMass(1) = SumX/n
    CenterOfMass(2) = SumY/n
    CenterOfMass(3) = SumZ/n
    !CenterOfMass(1) = SumX/dble(n)
    !CenterOfMass(2) = SumY/dble(n)
    !CenterOfMass(3) = SumZ/dble(n)
  end subroutine
  !------------------------------
  !------------------------------
  subroutine Move_COM_to_origin(n,x,y,z,CenterOfMass)
    !  1st(INPUT): The N of atoms 
    integer,intent(in) :: n 
    double precision,intent(inout) :: CenterOfMass(3)
    double precision,intent(inout) :: x(n), y(n), z(n)
    character(len=6) :: type_vari
    
    x(:) = x(:) - CenterOfMass(1)
    y(:) = y(:) - CenterOfMass(2)
    z(:) = z(:) - CenterOfMass(3)

    call Calc_CenterMass(n,x,y,z,CenterOfMass)
     
  end subroutine
  !------------------------------

  !------------------------------
  subroutine Calc_RadiusOfGyration(n,x,y,z,COM,Rg)
  !  1st(INPUT) : N of atoms 
  !  5th(INPUT) : center of mass
  !  6th(OUTPUT): Radius of gyration
  integer, intent(in) :: n
  double precision, intent(in) :: x(n), y(n), z(n)
  double precision, intent(in) :: COM(3) 
  double precision, intent(out) :: Rg
  double precision :: xtmp(n), ytmp(n), ztmp(n)
  integer :: i, j
  double precision :: ti, tf

  call cpu_time(ti)
  !Initialize
  Rg      = 0.0d0
  xtmp(:) = 0.0d0
  ytmp(:) = 0.0d0
  ztmp(:) = 0.0d0

  !(1) COM is subtructed from each coordinate (r - r_av) 
  xtmp(:) = x(:) - COM(1)
  ytmp(:) = y(:) - COM(2)
  ztmp(:) = z(:) - COM(3)

  !(2) ∑(ri - r_av)^2
  do i = 1, n  
    Rg = Rg + xtmp(i)**2 + ytmp(i)**2 + ztmp(i)**2
  enddo

  !(3) 1/N ∑(ri - r_av)^2
  Rg = Rg / n

  call cpu_time(tf)

  end subroutine
!------------------------------

!------------------------------
  subroutine Calc_EndToEnd(NumberOfAtomPair,x,y,z,TotalAtoms,ETE)
    integer,          intent(in)  :: NumberOfAtomPair(2) !Specify pair of atom number 
    integer,          intent(in)  :: TotalAtoms          !total atoms you read 
    double precision, intent(in)  :: x(:),y(:),z(:)
    double precision, intent(out) :: ETE
    double precision :: dx, dy, dz
    

    !(1) Calculate difference bet. ri- rj
    dx = x(NumberOfAtomPair(1)) - x(NumberOfAtomPair(2)) 
    dy = y(NumberOfAtomPair(1)) - y(NumberOfAtomPair(2)) 
    dz = z(NumberOfAtomPair(1)) - z(NumberOfAtomPair(2)) 
    ETE = sqrt(dx**2 + dy**2 + dz**2) 
  end subroutine


end module
