!Fri, Mar 13 2015
!

module reading
  use rmsd_module
  implicit none
contains
!------------------------------------------------------
  subroutine count_pdb(FilName,NumAtom)
    !This routine counts N of atoms in PDB.
    !Explanation of variables:
    !  1st(INPUT) : File name of PDB
    !  2nd(OUTPUT): N of atoms in PDB

    character(len=120) :: FilName
    character(len=6) :: atom
    integer, parameter :: FilUnit = 11
    integer :: i, ios, ter 
    integer :: NotAtomLine, NumAtom 
    integer :: NumAtomProtein

    open(FilUnit, file=FilName)

    write(*,'("  File Name: ",a80)') FilName

    NotAtomLine = 0
    do
      read(FilUnit,*) atom
      !---
      if(atom /= "ATOM" .and. atom /= "HETATM") then
        NotAtomLine = NotAtomLine + 1  
        cycle
      !---
      elseif(atom == "ATOM" .or. atom == "HETATM") then
        backspace(FilUnit)

        ios = 0
        ter = 0
        NumAtom = 0
!***Reading PDB
        do
          read(FilUnit,*, iostat=ios) atom

          select case(ios)
            case(0) !---ios equals zero(i.e. not finished reading)
              if(atom /= "ATOM" .and. atom /= "HETATM" .and. atom /= "TER") then
                write(*,'("    *In PDB reading, Non-ATOM line was detected, so stopped reading.")')
                goto 123
              !---I moniter TER just in case.
              elseif(atom == "TER") then
                ter = ter + 1
              !---When atom equals ATOM or HETATM, I count the number of atoms 
              else
                NumAtom = NumAtom + 1
              endif

            case default !---ios does not equal zero(i.e. finished reading)
              write(*,'("    *IO status      :",i5 )') ios
              goto 123
          end select

        enddo
      endif
    enddo
123 write(*,'("    *The N. of atoms: ",i8)') NumAtom
    write(*,'("    *TER lines      : ", i8)') ter

    close(FilUnit)
  end subroutine
!------------------------------------------------------

!------------------------------------------------------
  subroutine read_pdb(FilName,TotalAtoms,NumProteinAtom,&
                      AtomNum,AtomName,ResName,ChainId,ResNum,x,y,z)
    integer, intent(in) :: TotalAtoms
    character(len=120), intent(in) :: FilName

    integer,          intent(out) :: AtomNum(:)
    character(len=4), intent(out) :: AtomName(:)
    character(len=3), intent(out) :: ResName(:)
    character(len=1), intent(out) :: ChainId(:)
    integer         , intent(out) :: ResNum(:)
    double precision, intent(out) :: x(:),y(:),z(:)
    
    !I commented out at 7/18 Sat.
!    integer         , allocatable :: AtomNum(:)
!    character(len=4), allocatable :: AtomName(:)
!    character(len=3), allocatable :: ResName(:)
!    character(len=1), allocatable :: ChainId(:)
!    integer         , allocatable :: ResNum(:)
!    double precision, allocatable :: x(:),y(:),z(:)
    
    character(len=6) :: atom
    integer, parameter :: FilUnit = 11
    integer, parameter :: OutPDBUnit  = 12
    integer :: i, ios, ter 
    integer :: NotAtomLine,NumAtom 
    integer,intent(out) :: NumProteinAtom
    double precision :: ti, tf

    call cpu_time(ti)

    !I commented out at 7/18 Sat.
    !allocate(AtomNum(TotalAtoms),AtomName(TotalAtoms),ResName(TotalAtoms),&
    !         ChainId(TotalAtoms),ResNum(TotalAtoms),x(TotalAtoms),y(TotalAtoms),z(TotalAtoms))

    open(FilUnit, file=FilName)
    write(*,'("  File Name: ",a80)') FilName

!***Detect first atom line
    do
      read(FilUnit,*) atom
      if(atom /= "ATOM" .and. atom /= "HETATM") then
        NotAtomLine = NotAtomLine + 1  
        cycle
      elseif(atom == "ATOM" .or. atom == "HETATM") then
        exit
      endif
    enddo

    backspace(FilUnit)

!***If detected, start reading
    NumProteinAtom = 0
    NumAtom = 0
    !NumAtom = 1
    ios = 0
    ter = 0
    !do
    do i = 1, TotalAtoms 
      read(FilUnit,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)", iostat=ios) &
           atom       , &
           AtomNum(i) , &
           AtomName(i), &
           ResName(i) , &
           ChainId(i) , &
           ResNum(i)  , &
           x(i),        &
           y(i),        &
           z(i)  
      select case(ios)
       
        case(0) !---ios equals zero(i.e. not finished reading)
          if(atom /= "ATOM" .and. atom /= "HETATM" .and. atom /= "TER") then
            goto 123  !Finish this reading
          elseif(atom == "TER") then
            ter = ter + 1
          else
!***Count protein atoms
            if(ResName(i) /= "WAT" .and. &
               ResName(i) /= "NIO" .and. &
               ResName(i) /= "PIO" .and. &
               ResName(i) /= "CIM" .and. &
               ResName(i) /= "CIP" ) then
              NumProteinAtom = NumProteinAtom + 1
            endif
            NumAtom = NumAtom + 1 !Count the N of atoms 
          endif

        case default
          write(*,'("    *IO status: ",i8)') ios
          goto 123  !Finish this reading
      end select
    enddo
123 NumAtom = NumAtom 
    if (NumAtom /= TotalAtoms) then 
      print*,""
      print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      print*,"!  NumAtom /= TotalAtoms    !"
      print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      print*,"NumAtom =",NumAtom," TotalAtoms=", TotalAtoms
      stop
    endif
    print*,"************************************************************"
    write(*,'("    *The N. of atoms           : ",i8)') NumAtom
    write(*,'("    *The N. of atoms of protein: ",i8)') NumProteinAtom
    write(*,'("    *TER lines                 : ",i8)') ter

    call cpu_time(tf)
    write(*,'("    *CPU TIME(Reading PDB)     : ", f10.5," sec")') tf - ti

    close(FilUnit)
  end subroutine
!------------------------------------------------------
!------------------------------------------------------
  subroutine CountRmsdAtm(RmsdType,RegionCalcRMSD,N_TotalAtoms, AtomName,ResName,ResNum,N_RmsdAtoms)
    character(len=2)             , intent(in)  :: RmsdType
    integer                      , intent(in)  :: RegionCalcRMSD(2)
    integer                      , intent(in)  :: N_TotalAtoms 
    character(len=4)             , intent(in)  :: AtomName(N_TotalAtoms) 
    character(len=3)             , intent(in)  :: ResName(N_TotalAtoms) 
    integer                      , intent(in)  :: ResNum(N_TotalAtoms) 
!    character(len=4)             , intent(in)  :: AtomName(:) 
!    character(len=3)             , intent(in)  :: ResName(:) 
!    integer                      , intent(in)  :: ResNum(:) 
    integer                      , intent(out) :: N_RmsdAtoms
    integer :: i, icount

    N_RmsdAtoms = 0
    icount      = 0

    !write(*,'("  ***List of atoms used in RMSD calculation***")')
    do i = 1, N_TotalAtoms
      select case(RmsdType)
!***Count CA atoms 
        case("CA") 
          if(        (AtomName(i) == " CA ") &
            .and. (ResName(i) /= "ACE" .and. ResName(i) /= "NME" .and. ResName(i) /= "NH2") & 
            .and. (ResNum(i) >= RegionCalcRMSD(1) .and. ResNum(i) <= RegionCalcRMSD(2)) & 
            )then
              icount = icount + 1
!             print*, ResNum(i),AtomName(i), ResName(i)
          endif
!***Count backbone atoms 
        case("BB")
          if( (AtomName(i) == " CA " .or. AtomName(i) == " N " .or. AtomName(i) == " C ") .and. &
              (ResName(i) /= "ACE" .and. ResName(i) /= "NME" .and. ResName(i) /= "NH2")   .and. & 
              (ResNum(i) >= RegionCalcRMSD(1) .and. ResNum(i) <= RegionCalcRMSD(2)) ) then
              icount = icount + 1
!             print*, ResNum(i),AtomName(i), ResName(i)
          endif
!***Don't count 
        case("NO")  
          if(i == 1) then 
            print*," ??????????????????????????????????????????????????????????????????"
            print*," ??? Atoms used in  RMSD calc were not counted (parameter==NO). ???"
            print*," ??????????????????????????????????????????????????????????????????"
          endif
      end select  
    enddo
    N_RmsdAtoms = icount
    !write(*,'("    *Number of atoms used in RMSD calc.: ",i5)'), N_RmsdAtoms
  end subroutine
!------------------------------------------------------
!------------------------------------------------------
  subroutine SaveAtomNumForRmsd(RmsdType,RegionCalcRMSD,N_RmsdAtoms,N_TotalAtoms,AtomName, AtomNum,ResName,ResNum, AtomNumForRmsd)
    character(len=2)             , intent(in)  :: RmsdType
    integer                      , intent(in)  :: RegionCalcRMSD(2)
    integer                      , intent(in)  :: N_RmsdAtoms 
    integer                      , intent(in)  :: N_TotalAtoms
    character(len=4)             , intent(in)  :: AtomName(N_TotalAtoms)
    integer                      , intent(in)  :: AtomNum(N_TotalAtoms)
    character(len=3)             , intent(in)  :: ResName(N_TotalAtoms) 
    integer                      , intent(in)  :: ResNum(N_TotalAtoms) 
    integer                      , intent(out) :: AtomNumForRmsd(N_RmsdAtoms)
    integer :: i , icount

    icount = 1
    AtomNumForRmsd(:) = 0
    do i = 1, N_TotalAtoms
      select case(RmsdType)
!***Count CA atoms 
        case("CA") 
          if(        (AtomName(i) == " CA ") &
            .and. (ResName(i) /= "ACE" .and. ResName(i) /= "NME" .and. ResName(i) /= "NH2") & 
            .and. (ResNum(i) >= RegionCalcRMSD(1) .and. ResNum(i) <= RegionCalcRMSD(2)) & 
            )then
              AtomNumForRmsd(icount) = AtomNum(i)
              icount = icount + 1
          endif
!***Count backbone atoms 
        case("BB")
          if( (AtomName(i) == " CA " .or. AtomName(i) == " N " .or. AtomName(i) == " C ") .and. &
              (ResName(i) /= "ACE" .and. ResName(i) /= "NME" .and. ResName(i) /= "NH2")   .and. & 
              (ResNum(i) >= RegionCalcRMSD(1) .and. ResNum(i) <= RegionCalcRMSD(2)) ) then
              AtomNumForRmsd(icount) = AtomNum(i)
              icount = icount + 1
          endif
!***Don't count 
        case("NO")  
          if(i == 1) then 
            print*," ??????????????????????????????????????????????????????????????????"
            print*," ??? Atoms used in  RMSD calc were not saved (parameter==NO). ???"
            print*," ??????????????????????????????????????????????????????????????????"
          endif
      end select  
    enddo
!    print*,""
!    write(*,'("  ***List of atoms used in RMSD calculation***")')
!    do i = 1, N_RmsdAtoms 
!      print*,AtomNumForRmsd(i), AtomName(AtomNumForRmsd(i)),& 
!             ResNum(AtomNumForRmsd(i)),ResName(AtomNumForRmsd(i))
!    enddo
  end subroutine

!------------------------------------------------------
!------------------------------------------------------
  subroutine SaveCoordForRmsd(RmsdType, NumAtmsForRmsd, TotalAtoms,AtomNumForRmsd, Xcod, Ycod, Zcod, CoordForRmsd)
    integer :: i
    character(len=2)    , intent(in)  :: RmsdType
    integer             , intent(in)  :: NumAtmsForRmsd
    integer             , intent(in)  :: TotalAtoms
    integer             , intent(in)  :: AtomNumForRmsd(NumAtmsForRmsd) 
    real(8)             , intent(in)  :: Xcod(TotalAtoms), Ycod(TotalAtoms), Zcod(TotalAtoms) 
    real(8)             , intent(out) :: CoordForRmsd(:,:)
    !real(8), allocatable, intent(out) :: CoordForRmsd(:,:)
    !allocate(CoordForRmsd(3,NumAtmsForRmsd))
    do i = 1, NumAtmsForRmsd
      CoordForRmsd(1,i) = Xcod(AtomNumForRmsd(i))
      CoordForRmsd(2,i) = Ycod(AtomNumForRmsd(i))
      CoordForRmsd(3,i) = Zcod(AtomNumForRmsd(i))
    enddo
  end subroutine


!------------------------------------------------------
!------------------------------------------------------
  subroutine count_trj(FileName,tot_confs)
    !This routine counts trajectory file (binary), 
    !
    !Explanation of variables:
    !  *The 1st Variable: binary file name(input) 
    !  *The 2nd         : N of total atoms calculated by the subroutine(count_pdb) (input) 
    
    integer :: iconf                          !Loop variable 
    integer, intent(out) :: tot_confs         ! N of total confs
    integer :: ios, UnitTrjFile = 13    
    character(len=120),intent(in) :: FileName !Input trajectory file name 

    !***variables for reading trajectory which is binary format
    integer(4) :: istp,iyn15v,iyn15h
    !integer*4 :: istp,iyn15v,iyn15h
    real(4) :: sitime, sec, et, kinetic, temperature, rmsf, rmsd
    real(4), allocatable :: potential(:)
    write(*,*)
    write(*,'("#####Read trajectory#####")')

    open(UnitTrjFile, file = FileName, status="old",access="sequential",form="unformatted")
    write(*,'("  File(binary) Name: ",a80)') FileName

    !***Count conformations in trajectory file
    iconf = 0
    do 
      read(UnitTrjFile,iostat=ios)
      read(UnitTrjFile,iostat=ios)
      if (ios /= 0) then 
        print*," ios(gfortran=5001,ifort=-1):",ios 
        exit
      endif
      iconf = iconf + 1
    enddo
    tot_confs = iconf
    write(*,'("  N of conformations: ", i8 )') tot_confs 
    close(UnitTrjFile)
  end subroutine
!------------------------------------------------------



!------------------------------------------------------
  subroutine read_and_analy_trj(FileName,judge_output_pdb,tot_confs,iconfs,TotAtoms,&
                                total_residues,first_atom,last_atom,cell,ProteinAtoms,& 
                                AtomNum,AtomName,ResName,ChainId,ResNum,&
                                AtomNumForRmsd, Ref_AtomNumForRmsd,Ref_AtomNumForVCV, NumAtmsForRmsd,NumAtmsForVCV, AtomPair,&
                                trj_x,trj_y,trj_z,CoordForRmsd, &
                                TotalProb, NdataPDF, enepdf, pdf, Npairs,q,qq,vcv, &
                                atm_no_ete, targt_rsidu_atom, npair_for_targt_rsidu_atom)
    use calculation 
    use Read_Assign_PDF
    use vcv_mod
    use periodic
    integer :: i, j, ii, iii, ios
    character(len=126) :: blnk=""
    integer,           intent(in) :: TotAtoms, ProteinAtoms, tot_confs
    integer, intent(inout) :: iconfs

!@@@Periodic info.
    integer,           intent(in) :: total_residues
    integer,           intent(in) :: first_atom(total_residues)
    integer,           intent(in) :: last_atom(total_residues)
    real(8),intent(in) :: cell(3)
!@@@

!@@@Ref info.
    integer         , intent(in) :: AtomNum(TotAtoms)
    character(len=4), intent(in) :: AtomName(TotAtoms)
    character(len=3), intent(in) :: ResName(TotAtoms)
    character(len=1), intent(in) :: ChainId(TotAtoms)
    integer         , intent(in) :: ResNum(TotAtoms)
!@@@

    character(len=120),intent(in) :: FileName
    integer :: UnitTrjFile = 14, UnitOut = 1122

!@@@Trj. info.
    integer(4) :: istp,iyn15v,iyn15h    !No problem if you don't know
    real(4)   :: sitime, sec, et, kinetic, temperature, rmsf, rmsd
    real(4)   :: potential
    real(4), intent(out)   :: trj_x(:), trj_y(:), trj_z(:)
!@@@

    integer            , intent(in)    :: NumAtmsForRmsd
    integer            , intent(in)    :: NumAtmsForVCV
    integer            , intent(in)    :: AtomNumForRmsd(NumAtmsForRmsd) 
    integer            , intent(in)    :: Ref_AtomNumForRmsd(NumAtmsForRmsd) 
    integer            , intent(in)    :: Ref_AtomNumForVCV(NumAtmsForVCV) 

!For PDF 
    integer :: NdataPDF
    double precision :: enepdf(NdataPDF), pdf(NdataPDF)
    real(8) :: prob
    real(8), intent(inout) :: TotalProb

    double precision, allocatable :: tmp_x(:), tmp_y(:), tmp_z(:)
    double precision :: TrjCoordForRmsd(3,NumAtmsForRmsd)
!***Added Sep 28 2015
    double precision :: TrjDistForVCV(3,NumAtmsForVCV)

    double precision, intent(in) :: CoordForRmsd(NumAtmsForRmsd)
    double precision :: COM_trj(3)
    double precision :: COM_RMSD(3)
    double precision :: Rg_trj 
    double precision :: ETE_trj
    double precision :: SymMat(4,4)
    double precision :: EiganVal(4), EiganVec(4,4)
    double precision :: MinEigenVal, MinEigenVec(4)
    double precision :: Quartanion(0:3)
    double precision :: RotMat(1:3,1:3) 
    double precision :: RmsdVal 

    integer,intent(in) :: AtomPair(2)

!***For VCV
    integer, intent(in)    :: Npairs 
    real(8), allocatable   :: DistsOfPair(:) 
    real(8), intent(inout) :: q(Npairs),qq(Npairs,Npairs)
    real(8), intent(inout) :: vcv(Npairs,Npairs) 

!***
    integer :: targt_rsidu_atom
    integer :: npair_for_targt_rsidu_atom
    integer :: atm_no_ete(:)
    integer :: atom_pair_ete(2)
    integer :: unit_out

    integer :: No_of_chains

    !***
    character(len=130) :: fill_outpdbs
    character(len=3) :: judge_output_pdb

    !@@@For Rg which is a measure of hydrophobic core 
    real(8) :: HB_coord(19,3)
    !real(8) :: HB_coord(5,3)

    open(UnitTrjFile, file = FileName, status="old",access="sequential",form="unformatted")
    open(UnitOut,     file = "output.dat")
    open(74,          file = "each_dist.dat")
    open(63,          file = "q_dist.dat")
    open(62,          file = "vcv_mat.dat")

    allocate(tmp_x(TotAtoms),tmp_y(TotAtoms),tmp_z(TotAtoms) )
    tmp_x(:) = 0.0d0
    tmp_y(:) = 0.0d0
    tmp_z(:) = 0.0d0

    print*," Total Atoms:", TotAtoms

!***Read data in trajectry file
    ios = 0
    do i = 1, tot_confs
      read(UnitTrjFile) istp,sitime,sec,et,kinetic,temperature,&
                        potential,rmsf,iyn15v,iyn15h,rmsd
      read(UnitTrjFile, iostat=ios) (trj_x(j), trj_y(j), trj_z(j), j = 1, TotAtoms)
      if (ios > 0) then 
          print*,"NOTE) Error in reading trajectry file occurs. Iostat is positive:",ios
          stop
      endif

      iconfs = iconfs + 1
      print*,""
      write(*,'("#Number of conformation ",i8)') iconfs

      No_of_chains = 1   !For p53 
      !No_of_chains = 2   !For ET1 system
      call ReturnAtom(TotAtoms,No_of_chains,first_atom,last_atom,total_residues,trj_x,trj_y,trj_z,cell)
      print*, "aaa"

      !***Aim to making a structure possess a probability***
      call AssignPDF(NdataPDF,enepdf,pdf,potential,prob) 
      !*****************************************************
      TotalProb = TotalProb + prob 

!***Replace to convert real(4) to real(8) 
      do j = 1, TotAtoms
        tmp_x(j) = trj_x(j)
        tmp_y(j) = trj_y(j)
        tmp_z(j) = trj_z(j)
      enddo

      !call Calc_CenterMass(ProteinAtoms,tmp_x,tmp_y,tmp_z,COM_trj)
   !   print*,"Function COM =", COM(tmp_x,tmp_y,tmp_z)
    !  write(*,'("Center Of Mass (protein) [A]       : ",3f8.3)') COM_trj

    
!      call Calc_RadiusOfGyration(ProteinAtoms,tmp_x,tmp_y,tmp_z,COM_trj,Rg_trj)
!      write(*,'("Radius of gyration (protein)       : ", f10.3)') Rg_trj

      !@@@@
      HB_coord(:,:) = 0.0
      !***not used FOr Ac H+
      !HB_coord(1,:) = (/ tmp_x(139), tmp_y(139), tmp_z(139) /)
      !HB_coord(2,:) = (/ tmp_x(179), tmp_y(179), tmp_z(179) /)
      !HB_coord(3,:) = (/ tmp_x(205), tmp_y(205), tmp_z(205) /)
      !HB_coord(4,:) = (/ tmp_x(224), tmp_y(224), tmp_z(224) /)
      !HB_coord(5,:) = (/ tmp_x(241), tmp_y(241), tmp_z(241) /)

      !***Common atoms
      HB_coord(1,:) = (/ tmp_x(11), tmp_y(11), tmp_z(11) /)
      HB_coord(2,:) = (/ tmp_x(14), tmp_y(14), tmp_z(14) /)
      HB_coord(3,:) = (/ tmp_x(17), tmp_y(17), tmp_z(17) /)
      HB_coord(4,:) = (/ tmp_x(33), tmp_y(33), tmp_z(33) /)
      HB_coord(5,:) = (/ tmp_x(36), tmp_y(36), tmp_z(36) /)
      HB_coord(6,:) = (/ tmp_x(39), tmp_y(39), tmp_z(39) /)
      HB_coord(7,:) = (/ tmp_x(139), tmp_y(139), tmp_z(139) /)

      !***NonAc(H+) & Ac(H+)
      !HB_coord(8,:) = (/ tmp_x(157), tmp_y(157), tmp_z(157) /)
      !HB_coord(9,:) = (/ tmp_x(160), tmp_y(160), tmp_z(160) /)
      !HB_coord(10,:) = (/ tmp_x(163), tmp_y(163), tmp_z(163) /)
      !HB_coord(11,:) = (/ tmp_x(179), tmp_y(179), tmp_z(179) /)
      !HB_coord(12,:) = (/ tmp_x(182), tmp_y(182), tmp_z(182) /)
      !HB_coord(13,:) = (/ tmp_x(185), tmp_y(185), tmp_z(185) /)

      !***NonAc(H+)
      !HB_coord(14,:) = (/ tmp_x(201), tmp_y(201), tmp_z(201) /)
      !HB_coord(15,:) = (/ tmp_x(220), tmp_y(220), tmp_z(220) /)
      !HB_coord(16,:) = (/ tmp_x(237), tmp_y(237), tmp_z(237) /)
      !HB_coord(17,:) = (/ tmp_x(257), tmp_y(257), tmp_z(257) /)
      !HB_coord(18,:) = (/ tmp_x(260), tmp_y(260), tmp_z(260) /)
      !HB_coord(19,:) = (/ tmp_x(263), tmp_y(263), tmp_z(263) /)

      !***Ac(H+)
      !HB_coord(14,:) = (/ tmp_x(205), tmp_y(205), tmp_z(205) /)
      !HB_coord(15,:) = (/ tmp_x(224), tmp_y(224), tmp_z(224) /)
      !HB_coord(16,:) = (/ tmp_x(241), tmp_y(241), tmp_z(241) /)
      !HB_coord(17,:) = (/ tmp_x(261), tmp_y(261), tmp_z(261) /)
      !HB_coord(18,:) = (/ tmp_x(264), tmp_y(264), tmp_z(264) /)
      !HB_coord(19,:) = (/ tmp_x(267), tmp_y(267), tmp_z(267) /)

      !***NonAc & Ac
      HB_coord(8,:) = (/ tmp_x(156), tmp_y(156), tmp_z(156) /)
      HB_coord(9,:) = (/ tmp_x(159), tmp_y(159), tmp_z(159) /)
      HB_coord(10,:) = (/ tmp_x(162), tmp_y(162), tmp_z(162) /)
      HB_coord(11,:) = (/ tmp_x(178), tmp_y(178), tmp_z(178) /)
      HB_coord(12,:) = (/ tmp_x(181), tmp_y(181), tmp_z(181) /)
      HB_coord(13,:) = (/ tmp_x(184), tmp_y(184), tmp_z(184) /)

      !***NonAc
      !HB_coord(14,:) = (/ tmp_x(200), tmp_y(200), tmp_z(200) /)
      !HB_coord(15,:) = (/ tmp_x(219), tmp_y(219), tmp_z(219) /)
      !HB_coord(16,:) = (/ tmp_x(236), tmp_y(236), tmp_z(236) /)
      !HB_coord(17,:) = (/ tmp_x(256), tmp_y(256), tmp_z(256) /)
      !HB_coord(18,:) = (/ tmp_x(259), tmp_y(259), tmp_z(259) /)
      !HB_coord(19,:) = (/ tmp_x(262), tmp_y(262), tmp_z(262) /)

      !***Ac
      HB_coord(14,:) = (/ tmp_x(204), tmp_y(204), tmp_z(204) /)
      HB_coord(15,:) = (/ tmp_x(223), tmp_y(223), tmp_z(223) /)
      HB_coord(16,:) = (/ tmp_x(240), tmp_y(240), tmp_z(240) /)
      HB_coord(17,:) = (/ tmp_x(260), tmp_y(260), tmp_z(260) /)
      HB_coord(18,:) = (/ tmp_x(263), tmp_y(263), tmp_z(263) /)
      HB_coord(19,:) = (/ tmp_x(266), tmp_y(266), tmp_z(266) /)

      !do j = 1, 5! 19
      !  write(321,'(3f10.3)') HB_coord(j,:)
      !enddo
      call Calc_CenterMass(19,HB_coord(:,1),HB_coord(:,2),HB_coord(:,3),COM_trj)
      call Calc_RadiusOfGyration(19,HB_coord(:,1),HB_coord(:,2),HB_coord(:,3),COM_trj,Rg_trj)
      !call Calc_CenterMass(5,HB_coord(:,1),HB_coord(:,2),HB_coord(:,3),COM_trj)
      !call Calc_RadiusOfGyration(5,HB_coord(:,1),HB_coord(:,2),HB_coord(:,3),COM_trj,Rg_trj)
      write(321, *) Rg_trj, prob  
      !stop
      !@@@@

      !@@@ New line. Aug 13
      unit_out = 1001
      atom_pair_ete(:) = 0
      atom_pair_ete(1) = targt_rsidu_atom 
      do j = 1, npair_for_targt_rsidu_atom 
        atom_pair_ete(2) = atm_no_ete(j)
        call Calc_EndToEnd(atom_pair_ete,tmp_x,tmp_y,tmp_z,TotAtoms,ETE_trj)
        write(unit_out,*) iconfs, prob, ETE_trj
        unit_out = unit_out + 1 
      enddo
      !@@@

      call Calc_EndToEnd(AtomPair,tmp_x,tmp_y,tmp_z,TotAtoms,ETE_trj)
      write(*,'("End-to-end                         : ", f10.3," [A]")') ETE_trj

      !call Move_COM_to_origin(ProteinAtoms,tmp_x,tmp_y,tmp_z,COM_trj) !tmp cods are  exactly moved
      !write(*,'("COM moved to origin [A]      : ",3f8.3)') COM_trj

!***RMSD part

      !***Aim at saving coordinates of CA/BackBone atoms which are used in RMSD calc.***
      do ii = 1, NumAtmsForRmsd 
        TrjCoordForRmsd(1,ii) = tmp_x(Ref_AtomNumForRmsd(ii))
        TrjCoordForRmsd(2,ii) = tmp_y(Ref_AtomNumForRmsd(ii))
        TrjCoordForRmsd(3,ii) = tmp_z(Ref_AtomNumForRmsd(ii))
!        print*,"specified atom no. in RMSD: ", Ref_AtomNumForRmsd(ii)
      enddo
      !*********************************************************************************

!***RMSD part
!(1) Calculate center of mass of a conformation 
      call Calc_CenterMass(NumAtmsForRmsd, & 
                           TrjCoordForRmsd(1,:), & 
                           TrjCoordForRmsd(2,:), & 
                           TrjCoordForRmsd(3,:), & 
                           COM_RMSD)
      !write(*,'("Center Of Mass (RMSD) [A]    : ",3f8.3)') COM_RMSD


!(2) Move the center to origin in order to prepare for fitting a conformation into the target structure 
      call Move_COM_to_origin(NumAtmsForRmsd, & 
                              TrjCoordForRmsd(1,:), & 
                              TrjCoordForRmsd(2,:), & 
                              TrjCoordForRmsd(3,:), & 
                              COM_RMSD) 
      !write(*,'("COM moved to origin   [A]    : ",3f8.3)') COM_RMSD

!(3) Make symmetric matrix composed of coordinates of a conformation and a target, respectively
      call MakeSymMat(NumAtmsForRmsd,TrjCoordForRmsd,CoordForRmsd,SymMat)

!(4) Solve eigen value problem for the symmetric matrix 
      call Jacobi(SymMat,EiganVal,EiganVec,4,4)
      MinEigenVal = 10000000
      do iii = 1, 4
      !  write(*,'(i5,", Eigen Value: ",f10.6)') iii, EiganVal(iii)
      !  write(*,'(5x,"  Eigen vector",4f15.3)') (EiganVec(j,iii), j =1, 4)

!(5) Detect minimum eigen value and eigen vector
        if(EiganVal(iii) < MinEigenVal) then
          MinEigenVal = EiganVal(iii)
          MinEigenVec = EiganVec(:,iii)
        endif
      enddo
      !write(*,'("Minimum eigen vale: ", f10.6)') MinEigenVal
      !write(*,'("Minimum eigen vector: ")')
      !write(*,'(4f10.6)') (MinEigenVec(iii), iii = 1, 4)

!(6) Make matrix of quartanion 
      Quartanion(0:3) = MinEigenVec(1:4)
      call MakeRotationMat(Quartanion, RotMat)

!(7) Calculate RMSD

      call CalcRMSD(NumAtmsForRmsd, RotMat, TrjCoordForRmsd, CoordForRmsd, RmsdVal)
      write(*,'("RMSD                               : ",f10.3," [A]")') RmsdVal


!---Making variance-covariance matrix from here---

!***Save the coordinates used in the calc. of variance-covariance matrix
!      if (prob >= 0.2 .and. prob <= 1.0) then !***Pick only conformations at a certain T 

            do ii = 1, NumAtmsForVCV
              TrjDistForVCV(1,ii) = tmp_x(Ref_AtomNumForVCV(ii))
              TrjDistForVCV(2,ii) = tmp_y(Ref_AtomNumForVCV(ii))
              TrjDistForVCV(3,ii) = tmp_z(Ref_AtomNumForVCV(ii))
            enddo
      
      !(1) Save distances betw. atom pairs
            if (.not. allocated(DistsOfPair)) allocate(DistsOfPair(Npairs))
      
            call SaveDist(NumAtmsForVCV, & 
                          Npairs,Ref_AtomNumForVCV, &
                          TrjDistForVCV(1,:), &
                          TrjDistForVCV(2,:), & 
                          TrjDistForVCV(3,:), & 
                          DistsOfPair)
      
            write(74,*) iconfs, prob
            write(74,*) (DistsOfPair(j), j = 1, Npairs)
      
      !(2) Make terms of Variance-Covariance Matrix 
            call MakeTermsOfVCV(prob,Npairs,DistsOfPair,q,qq)
      
!      else
!        print*,"Not picked. conf.NO= ", iconfs, "Probability= ", prob
!      endif
      
        !***Aim at generating PDBs
        ios = 0
        select case(judge_output_pdb)
          case("YES")
            !The lower write line are temporary. this should be rewritten.2015/7/24
            write(fill_outpdbs, & 
            '("/Volumes/data_shinji/Endth_dat/TestPDBs/",i6.6,".pdb")') iconfs 
  
            print*," Output PDB name: ",fill_outpdbs
            open(10, file = fill_outpdbs, iostat=ios)
            if (ios > 0) stop "error"
            write(10,"('MODEL', i7)") iconfs
            write(10,"('# RMSD value: ',f10.3)") RmsdVal
            do j = 1, ProteinAtoms 
              write(10,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3,a26)") &
                "ATOM  ", &
                AtomNum(j), &
                AtomName(j),&
                ResName(j),&
                ChainId(j), &
                ResNum(j), &
                trj_x(j), trj_y(j), trj_z(j),&
                blnk
            enddo
            write(10,"(a6)") "ENDMDL"
            close(10)

  
          case("NO")
            !PDBs are not outputted 
        end select 
        !*****************************************************

        write(UnitOut,*) iconfs, RmsdVal, potential, ETE_trj, Rg_trj, prob
    enddo!Finish reading
    write(*,'("#####END#####")') 

    deallocate(tmp_x,tmp_y,tmp_z)
    close(UnitTrjFile)

  end subroutine
!------------------------------------------------------


end module
