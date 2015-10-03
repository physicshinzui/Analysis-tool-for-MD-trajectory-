module rmsd_module
  implicit none
contains
  
  subroutine Jacobi(a,e,v,n,nn)
  integer :: i, j, aaa
  integer,intent(in) :: n, nn              !dimension of array (a, e, and v)
  integer :: p, q, kaisuu                  !Subscripts of a 
  integer, parameter :: kmax = 30000       !Max val for calculation of Jacobi loop
  real(8), intent(inout)  :: a(nn,nn)      !Symmetric matrix
  real(8), intent(out) :: e(nn), v(n,n)    !eigen val and vector
  real(8) :: eps = 0.000001, bunbo,R,T,S,C
  real(8) :: apq, apqmax
  real(8) :: aip, aiq, apj, aqj, vip, viq

!***Judge weather symmetric matrix or not
  do i = 1, n
    do j = 1, n
      if (a(j,i) /= a(i,j)) then 
       !print*,i,j, a(j,i), a(i,j)
        stop "???Non-symmetric matrix was inputed.???"
      endif
    enddo
  enddo

!***Construct Unit matrix
  do i = 1, N
    do j = 1, N
      v(i,j) = 0.0d0
    enddo
    v(i,i) = 1.0d0
  enddo

!***Iteration of Jacobi
  do kaisuu = 1, kmax
    do p = 1, N -1
      do q = p+1, N
!***generate rotation angle
        bunbo = a(p,p) - a(q,q)
        if (bunbo /= 0.0d0) then
          R = 2.0d0 * a(p,q)/bunbo
          T = 0.5d0 * atan(R)
        else
          T = 0.78539818d0  !0.5 * (0.5*pi) <= lmit of arctan.   
        endif
        s = sin(T)
        c = cos(T)
!***Multiplication from left
        do j = 1, n
          apj = a(p,j)
          aqj = a(q,j)
          a(p,j) = apj*c + aqj*s
          a(q,j) = -apj*s + aqj *c
        enddo
!***Multiplication from right
        do i = 1, n
          aip = a(i,p)
          aiq = a(i,q)
          a(i,p) =  aip*c + aiq*s
          a(i,q) = -aip*s + aiq*c
          vip = v(i,p)
          viq = v(i,q)
          v(i,p) =  vip*c + viq*s
          v(i,q) = -vip*s + viq*c
        enddo
      enddo
    enddo
!    print*,kaisuu,v
!***Judgement of convergence
    apqmax = 0.0d0
    do p = 1, n-1
      do q = p+1, n
        apq = abs(a(p,q))
        if(apq > apqmax) apqmax = apq
      enddo
    enddo
    if(apqmax <= eps) goto 100 
  enddo
100 continue
!write(*,'("N of loops(Jacobi): ",i5)') kaisuu
!***trace eigenvalue
  forall (i = 1:n) e(i) = a(i,i)

  end subroutine
!-----------------------------------------

!-----------------------------------------
  subroutine MakeSymMat(n,rA,rB,S)
  integer :: i, j, k
  integer, intent(in) ::  n               !Array size
  real(8), intent(in) :: rA(3,n), rB(3,n)
  real(8), intent(out) :: S(4,4)          !Symmetric matrix
  real(8) :: a(3,n), b(3,n)

  !***
  do i = 1, n
    a(1,i) = rA(1,i) + rB(1,i)
    a(2,i) = rA(2,i) + rB(2,i)
    a(3,i) = rA(3,i) + rB(3,i)

    b(1,i) = rB(1,i) - rA(1,i)
    b(2,i) = rB(2,i) - rA(2,i)
    b(3,i) = rB(3,i) - rA(3,i)
  enddo

  !***Construct symmetric matrix
  S(:,:) = 0.0d0
  do j = 1, n
    S(1,1) = S(1,1) + b(1,j)**2 + b(2,j)**2 + b(3,j)**2
    S(2,1) = S(2,1) + a(3,j)*b(2,j) - a(2,j)*b(3,j)
    S(3,1) = S(3,1) - a(3,j)*b(1,j) + a(1,j)*b(3,j)
    S(4,1) = S(4,1) + a(2,j)*b(1,j) - a(1,j)*b(2,j)

    S(1,2) = S(1,2) + a(3,j)*b(2,j) - a(2,j)*b(3,j)
    S(2,2) = S(2,2) + b(1,j)**2 + a(2,j)**2 + a(3,j)**2
    S(3,2) = S(3,2) + b(1,j)*b(2,j) - a(1,j)*a(2,j)
    S(4,2) = S(4,2) + b(1,j)*b(3,j) - a(1,j)*a(3,j)

    S(1,3) = S(1,3) - a(3,j)*b(1,j) + a(1,j)*b(3,j)
    S(2,3) = S(2,3) + b(1,j)*b(2,j) - a(1,j)*a(2,j)
    S(3,3) = S(3,3) + a(1,j)**2 + b(2,j)**2 + a(3,j)**2 
    S(4,3) = S(4,3) + b(2,j)*b(3,j) - a(2,j)*a(3,j)

    S(1,4) = S(1,4) + a(2,j)*b(1,j) - a(1,j)*b(2,j)
    S(2,4) = S(2,4) + b(1,j)*b(3,j) - a(1,j)*a(3,j)
    S(3,4) = S(3,4) + b(2,j)*b(3,j) - a(2,j)*a(3,j)
    S(4,4) = S(4,4) + a(1,j)**2 + a(2,j)**2 + b(3,j)**2
  enddo
  S(:,:) = S(:,:)/n

!***Judge weather symmetric matrix or not
  do i = 1, 4
    do j = 1, 4
      if(S(j,i) == S(i,j)) then
!        write(*,'("   Symmetry! ","S(j,i): ", 2i4, f10.3, ", S(i,j): ", 2i4, f10.3)') i, j, S(j,i), j, i, S(i,j) 
      else
        write(*,'("No Symmetry! ","S(j,i): ", 2i4, f10.3, ", S(i,j): ", 2i4, f10.3)') i, j, S(j,i), j, i, S(i,j) 
        stop
      endif
    enddo
  enddo
  end subroutine
!--------------------------------------------

!--------------------------------------------
  subroutine MakeRotationMat(q, R)
    integer :: i
    real(8), intent(in)  :: q(0:3)     !quartanion
    real(8), intent(out) :: R(1:3,1:3) !Rotation matrix
    R(1,1) = 2*q(0)**2   + 2*q(1)**2 - 1
    R(2,1) = 2*q(1)*q(2) + 2*q(0)*q(3)
    R(3,1) = 2*q(1)*q(3) - 2*q(0)*q(2)
    R(1,2) = 2*q(1)*q(2) - 2*q(0)*q(3)
    R(2,2) = 2*q(0)**2   + 2*q(2)**2 - 1
    R(3,2) = 2*q(2)*q(3) + 2*q(0)*q(1)
    R(1,3) = 2*q(1)*q(3) + 2*q(0)*q(2)
    R(2,3) = 2*q(2)*q(3) - 2*q(0)*q(1)
    R(3,3) = 2*q(0)**2   + 2*q(3)**2 - 1 
  end subroutine
!--------------------------------------------

!--------------------------------------------
  subroutine CalcRMSD(n, R, rA, rB, RmsdVal)
    integer :: n                            !N of atoms of Ca, back bone, or all atoms  
    real(8), intent(in) :: R(1:3,1:3)       !Rotation matrix which was optimized
    real(8), intent(in) :: rA(3,n), rB(3,n)     !Coordinates of conformation A(<=target) and B(<=trajectry)
    real(8), intent(out) :: RmsdVal         !RMSD value 
    real(8) :: RrA(3,n) 
    real(8) :: xdiff, ydiff, zdiff, rdiff
    integer :: i, j
    
    xdiff = 0.0d0
    ydiff = 0.0d0
    zdiff = 0.0d0

    RrA  = matmul(R,rA)
    do i = 1, n
      xdiff = xdiff + (rB(1,i) - RrA(1,i))**2
      ydiff = ydiff + (rB(2,i) - RrA(2,i))**2
      zdiff = zdiff + (rB(3,i) - RrA(3,i))**2
    enddo
    rdiff = xdiff + ydiff + zdiff

    RmsdVal = sqrt( rdiff / dble(n) )
    !print*,"RMSD= ",RmsdVal
  end subroutine
!--------------------------------------------

end module

