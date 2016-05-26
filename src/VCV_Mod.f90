module vcv_mod
  implicit none
contains

!------------------------------------------
  subroutine SaveDist(n,Npairs,AtmNums,x,y,z,DistsOfPair)
    integer,              intent(in)  :: n
    integer,              intent(in)  :: Npairs
    integer,              intent(in)  :: AtmNums(n) 
    real(8),              intent(in)  :: x(n),y(n),z(n)
    real(8),              intent(out) :: DistsOfPair(:) 

    integer :: i, j, icount
    real(8) distX,distY,distZ,dist

    icount = 0
    DistsOfPair(:) = 0.0d0
    do i = 1, n - 1
      do j = i + 1, n 
        icount = icount + 1
        distX = x(i) - x(j) 
        distY = y(i) - y(j)
        distZ = z(i) - z(j)
        dist  = sqrt(distX**2 + distY**2 + distZ**2)  
        DistsOfPair(icount) = dist
      enddo
    enddo
  end subroutine

!------------------------------------------
!------------------------------------------
!  subroutine SaveCoord

!  end subroutine


!------------------------------------------
!------------------------------------------
  subroutine MakeTermsOfVCV(prob,Npairs, DistsOfPair,q,qq)
  !Memo:
  !    q and qq are updated iteratively.
    real(8), intent(in)    :: prob
    integer, intent(in)    :: Npairs
    real(8), intent(in)    :: DistsOfPair(Npairs) 
    real(8), intent(inout) :: q(Npairs),qq(Npairs,Npairs) 
    integer :: i, j, k
    integer :: icount

    !<q_i>
    do i = 1, Npairs
        !q(i) = q(i) + DistsOfPair(i) 
        q(i) = q(i) + prob * DistsOfPair(i) 
    enddo

    !<qi*qj>
    do i = 1, Npairs
      do j = 1, Npairs
        !qq(j,i) = qq(j,i) + DistsOfPair(j) * DistsOfPair(i)
        qq(j,i) = qq(j,i) + prob * DistsOfPair(j) * DistsOfPair(i)
      enddo
    enddo
  end subroutine

!------------------------------------------
!------------------------------------------

  subroutine MakeVCV(Npairs,q,qq,vcv)
    integer, intent(in)    :: Npairs
    real(8), intent(in)    :: q(Npairs),qq(Npairs,Npairs)
    real(8), intent(out)   :: vcv(Npairs,Npairs)
    integer :: i, j, k
    do i = 1, Npairs
      do j = 1, Npairs
        vcv(j,i) = qq(j,i) - q(j) * q(i)
      enddo
    enddo
  end subroutine
!------------------------------------------
!------------------------------------------


end module
