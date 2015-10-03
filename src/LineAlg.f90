module line_alg
  implicit none
contains

!  subroutine Jacobi(a,b,x,IterationMax,error0)

!  end subroutine

  subroutine gauss_seidel(a,b,x,n,IterationMax,error0)
    integer, intent(in) :: n, IterationMax
    real(8), intent(in) :: a(n,n), b(n), error0
    real(8), intent(out) :: x(n)
    real(8) :: s, er, rd(n),r(n)
    integer :: i, itr
    do i = 1, n
      if (a(i,i) == 0.0d0) stop "a(i,i) == 0.0d0"
      rd(i) = 1.0d0/a(i,i)
    enddo
    x(1:n) = 0.0d0
    print*,"x0=",x
    do itr = 1, IterationMax
      do i = 1, n
        s = dot_product(a(i,1:i-1), x(1:i-1))
        s = s + dot_product(a(i,i+1:n), x(i+1:n))
        x(i) = rd(i) * (b(i) -s)
      enddo
      r(1:n) = b(1:n) - matmul(a,x) 
      er= dot_product(r,r)
      write(*,*) "itr = ", itr, "error = ", er
      if(er <= error0) then
        write(*,*) "#converged#"
        exit
      endif
    enddo
  end subroutine

end module


program main
  use line_alg
  implicit none
  integer,parameter :: n = 2
  integer :: IterationMax=100
  real(8) :: a(n,n), b(n), error0 = 0.0000000001 
  real(8) :: x(n)
  real(8) :: s, er, rd(n),r(n)
  integer :: i, itr
  a(1,1) = 3.0d0
  a(2,1) = 2.0d0
  a(1,2) = 1.0d0
  a(2,2) = 1.0d0
  b(1:n) = (/9.0d0, 7.0d0/)

  call gauss_seidel(a,b,x,n,IterationMax,error0)
  print*,"x=",x


end program
