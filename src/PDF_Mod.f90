module Read_Assign_PDF
  implicit none
contains

  subroutine CountPDF(FilNamPDF, ncount)
    integer :: i, ios
    integer :: UnitPDF = 1022
    integer           , intent(out) :: ncount
    character(len=120), intent(in)  :: FilNamPDF
    real(8)  :: energy, pdf

    open(UnitPDF,file = FilNamPDF )
    print*,FilNamPDF
    ncount = 0
    ios = 0
    do
      read(UnitPDF,*,iostat=ios) energy, pdf
      if (ios /= 0)  exit
      ncount = ncount + 1
    enddo
    close(UnitPDF)
  end subroutine
  
  subroutine ReadPDF(FilNamPDF,ncount,energy,pdf)
    integer :: i, j, k
    integer, intent(in) ::  ncount
    integer :: UnitPDF =1001, UnitOut = 123
    integer :: ios
    character(len=120), intent(in)  :: FilNamPDF
    real(8)           , intent(out) :: energy(:), pdf(:)
    real(8) :: MaxPDF
    real(8), parameter :: LowLimit = 0.0000000001

    print*," File Name(PDF):", FilNamPDF
    open(UnitPDF,file = FilNamPDF )
    !open(UnitOut,file = "NonLogPdf.out" )

    !print*,"size(ene,pdf)=",size(energy,1), size(pdf,1)
!***Read PDF
    do i = 1, ncount
      read(UnitPDF,*) energy(i), pdf(i) 
      !write(1111,*)energy(i), pdf(i)
    enddo
    !print*,"n of pdf=", ncount

    MaxPDF = maxval(pdf)
    !print*,"Maximum PDF = ", MaxPDF

!***Normalization of PDF based on max value of PDF
    do i = 1, ncount
      pdf(i) = exp(pdf(i) - MaxPDF)
      if (pdf(i) <= LowLimit) pdf(i) = 0.0d0
      !write(1111, *) energy(i), pdf(i)
    enddo

    close(UnitPDF)
    !close(UnitOut)

  end subroutine


  subroutine AssignPDF(NdataPDF,enePDF,pdf,potential,prob)
    integer              :: i ,j ,k
    integer, intent(in)  :: NdataPDF
    real(8), intent(in)  :: enePDF(NdataPDF),pdf(NdataPDF)
    real(4), intent(in)  :: potential
    real(4)              :: MultProb
    real(8), intent(out) :: prob
!    print*,"Ndata,enePDF,pdf",NdataPDF, Nconf!enePDF(1564),pdf(1564)

    !do i = 1,nbin - 1
    do i = 1, NdataPDF - 1
      if(potential >= enePDF(i) .and. potential < enePDF(i+1)) then 
        MultProb = pdf(i) + pdf(i+1)
        MultProb = MultProb * 0.5d0
        prob     = MultProb
        write(*,'("Potential                          :", f15.3)') potential
        write(*,'("The range where a conf. is detected:", 2f15.3)') enePDF(i), enePDF(i+1)
        write(*,'("Neighbor Prob (i and i+1)          :", 2f12.8)') pdf(i), pdf(i+1) 
        write(*,'("Probability (Non-log scale)        :", f12.8)') prob
        exit
      endif
    enddo


  end subroutine



end module
