program main
  use calculation
  use reading 
  use Read_Assign_PDF
  use vcv_mod
  use periodic
  implicit none
  integer :: i,ii, j
  integer :: iconfs
  integer :: iensmbl
  integer :: itrj
  integer :: count_for_error
  integer, parameter :: Threshold = 5

  character(len=3) judge_output_pdb

  !***For PDF
  character(len=120), allocatable :: FilNameOfPdf(:)

  !***FOR PERIODIC BOUNDARY ANALYSIS
  real(8), allocatable :: cell_ens(:,:) 

  !***VARIABLES OF REFERENCE PDB DATA
  integer :: total_residue
  integer, allocatable :: first_atom_all(:,:) !used in periodic module
  integer, allocatable :: last_atom_all(:,:)  !used in periodic module
  integer, allocatable :: total_residue_all(:)

  integer, allocatable          :: AtomNum_ref_all(:,:)
  character(len=4), allocatable :: AtomName_ref_all(:,:)
  character(len=3), allocatable :: ResName_ref_all(:,:)
  character(len=1), allocatable :: ChainID_ref_all(:,:)
  integer, allocatable          :: ResNum_ref_all(:,:) 
  double precision, allocatable :: x_ref_all(:,:), y_ref_all(:,:), z_ref_all(:,:)    
  integer :: MaxNumOfAtoms  !Left side array size of the above ref. parameters
  integer :: max_value_of_res_num


  !***Variables for target PDB data
  integer, allocatable :: AtomNumOfTarget(:)
  character(len=4), allocatable :: AtomNameOfTarget(:)
  character(len=3), allocatable :: ResNameOfTarget(:)
  character(len=1), allocatable :: ChainIdOfTarget(:)
  integer, allocatable          :: ResNumOfTarget(:)
  double precision, allocatable :: XcoordOfTarget(:),YcoordOftarget(:),ZcoordOftarget(:)
  double precision, allocatable :: CoordForRmsd(:,:) 
  integer :: max_number

  !---variables of the PDBs plotted on PC space
  integer,          allocatable :: AtomNumOfPlt(:)
  character(len=4), allocatable :: AtomNameOfPlt(:)
  character(len=3), allocatable :: ResNameOfPlt(:)
  character(len=1), allocatable :: ChainIdOfPlt(:)
  integer, allocatable          :: ResNumOfPlt(:)
  double precision, allocatable :: XcoordOfPlt(:),YcoordOfPlt(:),ZcoordOfPlt(:)
  double precision, allocatable :: CoordForPlt(:,:) 


  !***Variables for binary file
  real(4), allocatable :: potential(:)
  real(4), allocatable :: trj_x(:), trj_y(:), trj_z(:)
  integer :: more_max_Natoms

  character(len=120), allocatable :: FilNameOfRef(:) !Input reference PDB name
  character(len=120) :: FilNameOfTag !Input target PDB name

  character(len=120), allocatable :: FilNameOfTrj(:,:) !Input binary name

  character(len=120), allocatable :: FilNameOfPltOnPCAspace(:)
  integer, allocatable :: NoOfTrj(:)
  integer :: max_no_trj
  integer :: NoOfPlt

  character(len=2) :: Rmsd_Type, PCA_Type
  integer          :: N_RmsdAtoms 
  integer, allocatable :: N_RmsdAtomsRef(:)
  integer :: N_AtomsPlt

  integer, allocatable  :: N_VCVAtomsRef(:)
  integer, allocatable :: Ref_AtomNumForVCV(:,:)

  !***ATOM NNUMBERS ARE SAVED IN:
  integer, allocatable :: Tag_AtomNumForRmsd(:)
  integer, allocatable :: Ref_AtomNumForRmsd(:,:)
  integer, allocatable :: AtomNums_for_plt(:) 

  integer, allocatable :: NumOfAtomsOfRef(:)          
  integer :: NumOfAtomsOfTag            
  integer, allocatable :: NumOfAromsOfProteinOfRef(:)   
  integer :: NumOfAtomsOfProteinOfTarget  
  integer, allocatable :: NumOfAtomsOfProteinOfPlt(:)     
  integer :: NumOfConfs                !The Number of conformations in binary
  integer, allocatable :: NumOfAtomsOfPlt(:)
 
  integer, allocatable :: RegionCalc_RMSD_Ref(:,:)
  integer, allocatable :: RegionCalc_VCV_Ref(:,:)
  integer :: RegionCalc_Tag(2)
  integer :: RegionCalc_proj(2)

  double precision :: ComOfTaget(3) 
  double precision :: ComOfAtomsOfRmsd(3) 

  double precision :: RgOfTarget

  integer :: NumberOfAtomPair(2)
  integer :: targt_rsidu_atom
  integer :: npair_for_targt_rsidu_atom
  integer, allocatable :: atm_no_ete(:)
  double precision :: EteOfTrj

!**RMSD variables
  real(8) :: R(1:3,1:3) !ROTATION MATRIX

!For PDF
  real(8) :: TotalProb
  integer :: No_of_ensembles
  integer, allocatable :: NdataPDF(:)
  integer :: max_NdataPDF
  double precision, allocatable :: enepdf_ens(:,:), pdf_ens(:,:)
  real(8) :: prob

!For VCV
  real(8), allocatable :: q(:),qq(:,:),vcv(:,:)
  integer :: Npairs, Npair_plt
  real(8), allocatable :: DistsOfPair_plt(:)
  real(8), allocatable :: x_proj(:), y_proj(:),z_proj(:)

  double precision :: ti, tf

  call cpu_time(ti)

!***INPUT PARAMETERS***

  read(*,*) judge_output_pdb
  
  !***READ NO. OF ENSEMBLES***
  read(*,*) No_of_ensembles
  allocate(FilNameOfRef(No_of_ensembles))
  do i = 1, No_of_ensembles 
    read(*,*) FilNameOfRef(i)
  enddo
  !***************************

  read(*,*) FilNameOfTag

  read(*,*) NoOfPlt

  !***READ PDB PROJECTED ON PC SPACE***
  select case(NoOfPlt)
    case(0)
    case(1)
      allocate(FilNameOfPltOnPCAspace(NoOfPlt), &
               NumOfAtomsOfPlt(NoOfPlt),NumOfAtomsOfProteinOfPlt(NoOfPlt) ) 
      NumOfAtomsOfProteinOfPlt(:) = 0
      do i = 1, NoOfPlt 
        read(*,*) FilNameOfPltOnPCAspace(i)
      enddo
    case default
      stop "?????The no.of files have to be one only. ?????"
  end select
  !************************************

  !***READ CELL SIZE***
  allocate(cell_ens(3,No_of_ensembles))
  do iEnsmbl = 1, No_of_ensembles
    read(*,*) cell_ens(1:3,iEnsmbl)
  enddo
  !********************

  !@@@ loop of ref of each ensemble is going to put here!   Fri, 7/17 2015
  read(*,*) targt_rsidu_atom, npair_for_targt_rsidu_atom
  allocate(atm_no_ete(npair_for_targt_rsidu_atom))
  read(*,*) (atm_no_ete(i), i = 1, npair_for_targt_rsidu_atom)
  !@@@
  read(*,*) NumberOfAtomPair(1), NumberOfAtomPair(2)

  read(*,'(a2)') Rmsd_Type
  read(*,'(a2)') PCA_Type

  allocate(RegionCalc_RMSD_Ref(2,No_of_ensembles))
  do iensmbl = 1, No_of_ensembles 
    read(*,*) RegionCalc_RMSD_Ref(1,iensmbl), RegionCalc_RMSD_Ref(2,iensmbl)
  enddo

  allocate(RegionCalc_VCV_Ref(2,No_of_ensembles))
  do iensmbl = 1, No_of_ensembles
    read(*,*) RegionCalc_VCV_Ref(1,iensmbl), RegionCalc_VCV_Ref(2,iensmbl)
  enddo

  read(*,*) RegionCalc_Tag(1), RegionCalc_Tag(2)
  read(*,*) RegionCalc_proj(1), RegionCalc_proj(2)

  allocate( FilNameOfPdf(No_of_ensembles) )
  do i = 1, No_of_ensembles
    read(*,*) FilNameOfPdf(i) 
  enddo

  allocate(NoOfTrj(No_of_ensembles))
  read(*,*) (NoOfTrj(i), i = 1, No_of_ensembles)
  max_no_trj = maxval(NoOfTrj)

  allocate(FilNameOfTrj(max_no_trj, No_of_ensembles))
  do j = 1, No_of_ensembles
    do i = 1, NoOfTrj(j)
      read(*,*) FilNameOfTrj(i,j)
    enddo
  enddo
!**********************

  
  write(*,'("#####Input parameters#####")')
  write(*,'("1.  Do you output pdb files?                      : ",a3)') judge_output_pdb 
  write(*,'("2.  The No. of ensembles                          : ", i5)') No_of_ensembles 
  write(*,'("3.  Reference PDBs are : ")') 
  do i = 1, No_of_ensembles
    write(*,'("     * Ensemble No.", i1, " :", a80)') i, FilNameOfRef(i)
  enddo
  write(*,'("4.  Target PDB                                    : ", a80)') FilNameOfTag
  write(*,'("5.  The No. of PDBs to plot it on PC space        : ", i2)') NoOfPlt
  select case(NoOfPlt)
  case(0)
    write(*,'("    *NOTE) The PDBs for plotting on PC space are not inputed.")')
  case(1:)
    do i = 1, NoOfPlt
      write(*,'("     *", a80)') FilNameOfPltOnPCAspace(i)
    enddo
  end select
  
  write(*,'("6.  Cell size                                     : ", 3f10.3)') cell_ens(1:3,:)

  write(*,'("7.  Atom No. of target residue (end-to-end)       : ",2i5)') targt_rsidu_atom, npair_for_targt_rsidu_atom 
  write(*,'("    *** List of atom No. of partner residue *** ")') 
  do i = 1, npair_for_targt_rsidu_atom
  write(*,'("    *                 ",i5,"                   *")') atm_no_ete(i) 
  enddo
  write(*,'("    *******************************************")')

  write(*,'("7.  No. of atom pair to calc. end-to-end          : ",2i5)') NumberOfAtomPair
  write(*,'("8.  Atom type for RMSD calc.                      : ", a2)') Rmsd_Type
  write(*,'("9.  Atom type for VCV calc.                       : ", a2)') PCA_Type

  write(*,'("10. The range of residues Ref. for RMSD calc.     : ")') 
  do iensmbl = 1, No_of_ensembles
    write(*,'("      * Ref. No.",i1,": ", 1x,2i5)') iensmbl, RegionCalc_RMSD_Ref(1:2,iensmbl) 
  enddo

  write(*,'("11. The range of residues Ref. for VCV calc.      : ")')
  do iensmbl = 1, No_of_ensembles
    write(*,'("      * Ref. No.",i1,": ", 1x,2i5)') iensmbl, RegionCalc_VCV_Ref(1:2,iensmbl) 
  enddo
  write(*,'("12. The range of Tag. resi. for RMSD              : ",2i5)') RegionCalc_Tag(1:2) 
  write(*,'("13. The range of residues of Ref. projection      : ", 2i5)') RegionCalc_proj(1:2) 

  write(*,'("14. Input PDFs are : " )')
  do i = 1, No_of_ensembles
    write(*,'("     *"," Ensemble No.",i1, " :",1x, a80)' ) i, FilNameOfPdf(i) 
  enddo

  write(*,'("15. The No. of files of trajectory are : ")') 
  do i = 1, No_of_ensembles
    write(*,'("     *"," Ensemble No.",i1, " :", 1x, i5)') i, NoOfTrj(i) 
  enddo

  write(*,'("16. Input trajectories are : ")') 
  do j = 1, No_of_ensembles
  write(*,'("    --Ensemble No.",i1,"--")') j 
    do i = 1, NoOfTrj(j)
      write(*,'("     *", a80)') FilNameOfTrj(i,j) 
    enddo
  enddo

  print*,""
  write(*,'("###Count atoms in reference PDB###")')

  allocate(NumOfAtomsOfRef(No_of_ensembles))
  do i = 1, No_of_ensembles
    call count_pdb(FilNameOfRef(i),NumOfAtomsOfRef(i))
  enddo
  MaxNumOfAtoms = maxval(NumOfAtomsOfRef)
  
  !***Aim for error detection***
  if ( NumberOfAtomPair(2) > minval(NumOfAtomsOfRef)) then
    print*,""
    print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print*,"!  Last atom no. for ETE distance calc. is bigger !" 
    print*,"!  than total atoms. The no. should be smaller.   !"
    print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print*,"                   NumOfAtoms =",minval(NumOfAtomsOfRef)
    print*,"AtomPair(2) should be smaller)=",NumberOfAtomPair(2)
    stop
  endif
  !*****************************
  print*,""

  write(*,'("###Count atoms in target PDB###")')
  call count_pdb(FilNameOfTag,NumOfAtomsOfTag)
  print*,""

  write(*,'("###Count atoms in the PDBs plotted on PC space###")')
  do i = 1, NoOfPlt 
    call count_pdb(FilNameOfPltOnPCAspace(i), NumOfAtomsOfPlt(i))
  enddo

!----------For reference PDB----------
  print*,""
  write(*,'("###Store reference PDB###")')

!**Aim for allocation of data of reference PDBs*** 
  allocate(NumOfAromsOfProteinOfRef(No_of_ensembles)      ,&
           AtomNum_ref_all(MaxNumOfAtoms,No_of_ensembles) ,&
           AtomName_ref_all(MaxNumOfAtoms,No_of_ensembles),&
           ResName_ref_all(MaxNumOfAtoms,No_of_ensembles) ,&
           ChainID_ref_all(MaxNumOfAtoms,No_of_ensembles) ,&
           ResNum_ref_all(MaxNumOfAtoms,No_of_ensembles)  ,&
           x_ref_all(MaxNumOfAtoms,No_of_ensembles)       ,&
           y_ref_all(MaxNumOfAtoms,No_of_ensembles)       ,&
           z_ref_all(MaxNumOfAtoms,No_of_ensembles)  )
           
  NumOfAromsOfProteinOfRef(:) = 0
  ResNum_ref_all(:,:) = 0
  x_ref_all(:,:) = 0.0d0
  y_ref_all(:,:) = 0.0d0
  z_ref_all(:,:) = 0.0d0
!************************************************ 

  do i = 1, No_of_ensembles
    call read_pdb(FilNameOfRef(i)            , &
                  NumOfAtomsOfRef(i)         , &
                  NumOfAromsOfProteinOfRef(i), &
                  AtomNum_ref_all(:,i)       , &
                  AtomName_ref_all(:,i)      , &
                  ResName_ref_all(:,i)       , &
                  ChainID_ref_all(:,i)       , &
                  ResNum_ref_all(:,i)        , &
                  x_ref_all(:,i), y_ref_all(:,i), z_ref_all(:,i)) 
  enddo

!***Aim at making sure that ref.s were read correctly by output. *** 
  open(2233, file = "in_pdb.log")
  do i = 1, No_of_ensembles 
    write(2233,*) "#Reference NO.",i
    do j = 1, NumOfAtomsOfRef(i)
       write(2233,"(a6,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)") &
         "ATOM  ", &
         AtomNum_ref_all(j,i), &
         AtomName_ref_all(j,i),&
         ResName_ref_all(j,i),&
         ChainID_ref_all(j,i), &
         ResNum_ref_all(j,i), &
         x_ref_all(j,i), y_ref_all(j,i), z_ref_all(j,i)
     enddo
   enddo
!***************************************************************** 

!@@@Preparatoin for modification of periodic boundary effect
  allocate(total_residue_all(No_of_ensembles))

!***Aim for allocation of first_atom_all & last_atom_all using more maximum value 
!***of total residues
  do iEnsmbl = 1, No_of_ensembles
    total_residue_all(iEnsmbl) = maxval(ResNum_ref_all(:,iEnsmbl))
  enddo
  max_value_of_res_num = maxval(total_residue_all)

  !print*,"total_residue_all=",total_residue_all !max_value_of_res_num

  allocate(first_atom_all(max_value_of_res_num, No_of_ensembles), &
            last_atom_all(max_value_of_res_num, No_of_ensembles) ) 
!******************************************************************************** 

  do i  = 1, No_of_ensembles
    call DetectResidueNumber(NumOfAromsOfProteinOfRef(i) ,ResNum_ref_all(:,i),total_residue_all(i),&
                             first_atom_all(:,i), last_atom_all(:,i)) 
  enddo

  do i = 1, No_of_ensembles  
    print*,""
    print*,"     Ensemble No.",i, "Total res.: ",total_residue_all(i)
    do j = 1, total_residue_all(i)
      if (first_atom_all(j,i) == 0 .and. last_atom_all(j,i) ==0 ) then 
        print*,"       ??? Exit this loop cuz first atom and last atom are 0 ???"
        exit
      endif
      print*,"       * Res. No., First, Last :", j, first_atom_all(j,i), last_atom_all(j,i)
    enddo
  enddo
!@@@

!***Specify atom numbers to calc. RMSD***
  allocate(N_RmsdAtomsRef(No_of_ensembles))
  do i = 1, No_of_ensembles 
    call CountRmsdAtm(Rmsd_Type,RegionCalc_RMSD_Ref(:,i),NumOfAtomsOfRef(i),& 
                      AtomName_ref_all(:,i), ResName_ref_all(:,i), ResNum_ref_all(:,i), N_RmsdAtomsRef(i) ) 
    write(*,'("    *Number of atoms of ref used in RMSD calc.: ",i5)') N_RmsdAtomsRef(i)
  enddo
  print*,""
  if (No_of_ensembles >= 2) then
    if (N_RmsdAtomsRef(1) /= N_RmsdAtomsRef(2) ) then 
      print*,"??? The No. of RMSD atoms differs. ???"
      print*,"N_RmsdAtomsRef(1) /= N_RmsdAtomsRef(2): ",N_RmsdAtomsRef(1), N_RmsdAtomsRef(2)
      stop
    endif
  endif
!****************************************

!***Saving atom numbers of each references which will be used in RMSD calc.***
  allocate(Ref_AtomNumForRmsd(N_RmsdAtomsRef(1),No_of_ensembles))
  do i = 1, No_of_ensembles
    call SaveAtomNumForRmsd(Rmsd_Type               , & 
                            RegionCalc_RMSD_Ref(:,i), &
                            N_RmsdAtomsRef(i)       , &
                            NumOfAtomsOfRef(i)      , &
                            AtomName_ref_all(:,i)   , &
                            AtomNum_ref_all(:,i)    , &
                            ResName_ref_all(:,i)    , &
                            ResNum_ref_all(:,i)     , &
                            Ref_AtomNumForRmsd(:,i) )
  enddo
!*****************************************************************************

  write(*,'("      ***LIST OF ATOMS USED IN RMSD CALCULATION***")')
  do iEnsmbl = 1, No_of_ensembles
    print*,"     *Ensemble No.",iEnsmbl,"                 *"
    do j = 1, N_RmsdAtomsRef(iEnsmbl) 
      print*,"     * ",&
             Ref_AtomNumForRmsd(j,iEnsmbl)                            , &
             AtomName_ref_all(Ref_AtomNumForRmsd(j,iEnsmbl), iEnsmbl ), &
             ResNum_ref_all(Ref_AtomNumForRmsd(j,iEnsmbl), iEnsmbl)   , &
             ResName_ref_all(Ref_AtomNumForRmsd(j,iEnsmbl), iEnsmbl)  , &
             "        *"
    enddo
    write(*,'("      ********************************************")')
  enddo

  print*,""
  write(*,'("    *****PCA alalysis info.*****")')
!***SPECIFY ATOM NUMBERS TO CALC. VCV MATRIX***
  allocate(N_VCVAtomsRef(No_of_ensembles)) 
  do i = 1, No_of_ensembles
    call CountRmsdAtm(PCA_Type               , &
                      RegionCalc_VCV_Ref(:,i), &
                      NumOfAtomsOfRef(i)     , & 
                      AtomName_ref_all(:,i)  , &
                      ResName_ref_all(:,i)   , &
                      ResNum_ref_all(:,i)    , & 
                      N_VCVAtomsRef(i) ) 
    write(*,'("    *The number of atoms of ref. used in VCV calc.: ",i5)') N_VCVAtomsRef(i) 
  enddo
  print*,""
!**********************************************

!***Saving atom numbers of each references which will be used in VCV calc.***
  allocate(Ref_AtomNumForVCV(N_VCVAtomsRef(No_of_ensembles), No_of_ensembles))
  do i = 1, No_of_ensembles
    call SaveAtomNumForRmsd(PCA_Type,&
                            RegionCalc_VCV_Ref(:,i), &
                            N_VCVAtomsRef(i), &
                            NumOfAtomsOfRef(i), &
                            AtomName_ref_all(:,i), &
                            AtomNum_ref_all(:,i),&
                            ResName_ref_all(:,i), &
                            ResNum_ref_all(:,i), &
                            Ref_AtomNumForVCV(:,i) )
  enddo
!****************************************************************************

  write(*,'("      ***LIST OF ATOMS USED IN VCV CALCULATION***")')
  do iEnsmbl = 1, No_of_ensembles
    print*,"     *Ensemble No.",iEnsmbl, "                *"
    do j = 1, N_VCVAtomsRef(iEnsmbl)
      print*,"     *",&
             Ref_AtomNumForVCV(j,iEnsmbl),&
             AtomName_ref_all(Ref_AtomNumForVCV(j,iEnsmbl), iEnsmbl), &
             ResNum_ref_all(Ref_AtomNumForVCV(j,iEnsmbl), iEnsmbl)  , &
             ResName_ref_all(Ref_AtomNumForVCV(j,iEnsmbl), iEnsmbl) , &
             "        *"
    enddo
    print*,"     *******************************************"
  enddo

!***Calculation of the number of atom pairs that is used in VCV calc.*** 
  Npairs = ( N_VCVAtomsRef(1) * (N_VCVAtomsRef(1) - 1) ) * 0.5 
  print*,"   *The no. of atom pairs of Ref. : ", Npairs
!*********************************************************************** 


!----------For target PDB----------
  print*,""
  write(*,'("###Store target PDB###")')
  allocate(AtomNumOfTarget(NumOfAtomsOfTag) , &
           AtomNameOfTarget(NumOfAtomsOfTag), &
           ResNameOfTarget(NumOfAtomsOfTag) , &
           ChainIdOfTarget(NumOfAtomsOfTag) , &
           ResNumOfTarget(NumOfAtomsOfTag)  , &
           XcoordOfTarget(NumOfAtomsOfTag)  , &
           YcoordOftarget(NumOfAtomsOfTag)  , &
           ZcoordOftarget(NumOfAtomsOfTag) )

  call read_pdb(FilNameOfTag,NumOfAtomsOfTag,NumOfAtomsOfProteinOfTarget,AtomNumOfTarget,AtomNameOfTarget,&
                ResNameOfTarget,ChainIdOfTarget,ResNumOfTarget,&
                XcoordOfTarget,YcoordOftarget,ZcoordOftarget) 


  if ( AtomNumOfTarget(1) /= 1) then 
    print*,""
    print*,"?????????????????????????????????????????????????????"
    print*,"??? Atom number must be listed in acending order. ???"
    print*,"?????????????????????????????????????????????????????"
    print*,"1st Atom number is ", AtomNumOfTarget(1)
    stop
  endif

  call Calc_CenterMass(NumOfAtomsOfProteinOfTarget,XcoordOfTarget,YcoordOftarget,ZcoordOftarget,ComOfTaget)
  print*,""
  write(*,'("    *Center Of Mass (protein) [A]  : ",3f10.3)') ComOfTaget(1:3)

  call Calc_RadiusOfGyration(NumOfAtomsOfProteinOfTarget,XcoordOfTarget,YcoordOftarget,ZcoordOftarget,ComOfTaget,RgOfTarget)
  write(*,'("    *Radius of gyration (protein)  : ", f10.3)') sqrt(RgOfTarget)


!----------For RMSD part of target----------
  call CountRmsdAtm(Rmsd_Type,RegionCalc_Tag,NumOfAtomsOfTag,& 
                    AtomNameOfTarget,ResNameOfTarget,ResNumOfTarget,N_RmsdAtoms) 
  write(*,'("    *Number of atoms used in RMSD calc.: ",i5)') N_RmsdAtoms

  allocate(Tag_AtomNumForRmsd(N_RmsdAtoms))

  call SaveAtomNumForRmsd(Rmsd_Type,RegionCalc_Tag,N_RmsdAtoms,NumOfAtomsOfTag,AtomNameOfTarget, AtomNumOfTarget,&
                          ResNameOfTarget,ResNumOfTarget,Tag_AtomNumForRmsd)

  allocate(CoordForRmsd(3,N_RmsdAtoms))
  call SaveCoordForRmsd(Rmsd_Type,N_RmsdAtoms,NumOfAtomsOfTag,Tag_AtomNumForRmsd, &
                         XcoordOfTarget,YcoordOftarget,ZcoordOftarget, CoordForRmsd) 

  !???????????????@@@@@
  write(*,'("      ***LIST OF ATOMS of TAG USED IN RMSD CALCULATION***")')
  do j = 1, N_RmsdAtoms 
    print*,"     * ",&
           Tag_AtomNumForRmsd(j)                  , &
           AtomNameOfTarget(Tag_AtomNumForRmsd(j)), &
           ResNumOfTarget(Tag_AtomNumForRmsd(j))  , &
           ResNameOfTarget(Tag_AtomNumForRmsd(j)) , &
           !XcoordOfTarget(Tag_AtomNumForRmsd(j))   , &
           "        *"
  enddo
  write(*,'("      ********************************************")')
  !???????????????@@@@@


  call Calc_CenterMass(N_RmsdAtoms,CoordForRmsd(1,:),CoordForRmsd(2,:),CoordForRmsd(3,:),ComOfAtomsOfRmsd)
  write(*,'("    *Center Of Mass (",a2,")","      [A]  : " ,3f10.3)') Rmsd_Type,ComOfAtomsOfRmsd(1:3)
  call Move_COM_to_origin(N_RmsdAtoms,CoordForRmsd(1,:),CoordForRmsd(2,:),CoordForRmsd(3,:),ComOfAtomsOfRmsd)
  write(*,'("    *COM moved to origin      [A]  : ",3f10.3)') ComOfAtomsOfRmsd(1:3)
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!---for PDBs plotted on PC space
  select case(NoOfPlt)
    case(0)

    case(1)
      print*,""
      write(*,'("###Store the PDBs plotted on PC space###")')
      allocate(AtomNumOfPlt(NumOfAtomsOfPlt(NoOfPlt)), &
               AtomNameOfPlt(NumOfAtomsOfPlt(NoOfPlt)), &
               ResNameOfPlt(NumOfAtomsOfPlt(NoOfPlt)) , &
               ChainIdOfPlt(NumOfAtomsOfPlt(NoOfPlt)), &
               ResNumOfPlt(NumOfAtomsOfPlt(NoOfPlt)), & 
               XcoordOfPlt(NumOfAtomsOfPlt(NoOfPlt)), &
               YcoordOfPlt(NumOfAtomsOfPlt(NoOfPlt)) , &
               ZcoordOfPlt(NumOfAtomsOfPlt(NoOfPlt)))
      do i = 1, NoOfPlt
        call read_pdb(FilNameOfPltOnPCAspace(i),NumOfAtomsOfPlt(i),NumOfAtomsOfProteinOfPlt(i),&
                      AtomNumOfPlt,AtomNameOfPlt ,ResNameOfPlt ,ChainIdOfPlt ,ResNumOfPlt ,&
                      XcoordOfPlt, YcoordOfPlt ,ZcoordOfPlt) 

!***WARNING DETECTION***
        if ( AtomNumOfTarget(1) /= 1) then 
          print*,""
          print*,"?????????????????????????????????????????????????????"
          print*,"??? Atom number must be listed in acending order. ???"
          print*,"?????????????????????????????????????????????????????"
          print*,"1st Atom number is ", AtomNumOfTarget(1)
          stop
        endif

        call CountRmsdAtm(PCA_Type,RegionCalc_proj,NumOfAtomsOfPlt(i),& 
                          AtomNameOfPlt,ResNameOfPlt,ResNumOfPlt,N_AtomsPlt) 
        write(*,'("    *No. of atoms of a region of the PDB plotted on PC space : ",i5)') N_AtomsPlt

!***WARNING DETECTION
        if (N_AtomsPlt /= N_VCVAtomsRef(1) ) then
          print*,""
          print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print*,"! The no. of atoms which is      !"
          print*,"! projected on PC space does not !"
          print*,"! correspond to reference one.   !"
          print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          print*,"No. (Ref) :",N_VCVAtomsRef
          print*,"No. (Proj):",N_AtomsPlt
          stop
        endif


        allocate(AtomNums_for_plt(N_AtomsPlt), x_proj(N_AtomsPlt), y_proj(N_AtomsPlt), z_proj(N_AtomsPlt))
        AtomNums_for_plt(:) = 0

        call SaveAtomNumForRmsd(PCA_Type,RegionCalc_proj,N_AtomsPlt ,NumOfAtomsOfPlt(i),AtomNameOfPlt, AtomNumOfPlt,&
                                ResNameOfPlt,ResNumOfPlt,AtomNums_for_plt)

!***SAVING COODRDINATES OF ATOMS WHICH ARE USED IN VCV CALC.***
        do ii = 1, N_AtomsPlt 
          x_proj(ii) = XcoordOfPlt(AtomNums_for_plt(ii))
          y_proj(ii) = YcoordOfPlt(AtomNums_for_plt(ii))
          z_proj(ii) = ZcoordOfPlt(AtomNums_for_plt(ii))
        enddo
!**************************************************************
                            
        Npair_plt = N_AtomsPlt * (N_AtomsPlt - 1) *0.5d0  
        print*,"   *Atom pairs of a PDB projected on PC space :",Npair_plt

!***SAVING DISTANCES BETW. THE ATOM PAIRS***
        allocate(DistsOfPair_plt(Npair_plt))
        call SaveDist(N_AtomsPlt, Npair_plt, AtomNums_for_plt, &
                      x_proj,y_proj,z_proj,DistsOfPair_plt)
!*******************************************

!***WRITE THE DISTANCES***
        open(123, file = "proj.dist")
        write(123,*) 0
        write(123,*) (DistsOfPair_plt(j), j = 1, Npair_plt)
!*************************
      enddo
   end select

!***CHECKING WHETHER RESIDUE NAME CORRESPONDS OR NOT***
  count_for_error = 0
  do i = 1, N_VCVAtomsRef(1) 
    if ((ResNameOfPlt(AtomNums_for_plt(i))) /= ResName_ref_all(Ref_AtomNumForVCV(i,1),1)) then 
      print*,""
      print*,""
      print*,"!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!"
      print*,"!  Both residue names are mismatch. !"
      print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      print*, "Res. No.          : ", i
      print*, "Res. Name in Ref  : ", ResName_ref_all(Ref_AtomNumForVCV(i,1),1)
      print*, "Res. name for proj: ", ResNameOfPlt(AtomNums_for_plt(i))
      print*,"-------------------------------"
      print*,""
      print*,""
      count_for_error = count_for_error + 1
      if (count_for_error > Threshold) then
        print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*,"!! Mismatch residues are too much !!"
        print*,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        stop
      endif
    endif
  enddo
!******************************************************



!-----For trajectory-----
  allocate(q(Npairs),qq(Npairs,Npairs),vcv(Npairs,Npairs))

  TotalProb = 0.0d0
  q(:)      = 0.0d0
  qq(:,:)   = 0.0d0
  vcv(:,:)  = 0.0d0

  print*,""
  print*,"###Read probability distribution info.###"

!***Aim to allocate enepdf_ens & pdf_ens using the max_NdataPDF****
  allocate(NdataPDF(No_of_ensembles))
  do iEnsmbl = 1, No_of_ensembles
    call CountPDF(FilNameOfPdf(iEnsmbl), NdataPDF(iEnsmbl))
  enddo
  write(*,*) NdataPDF(:)
  max_NdataPDF = maxval(NdataPDF)
  !print*,"More max n data of pdf for use of allocation: ",max_NdataPDF
  !print*,""

  allocate(enepdf_ens(max_NdataPDF,No_of_ensembles), pdf_ens(max_NdataPDF,No_of_ensembles))
!******************************************************************

  enepdf_ens(:,:) = 0.0d0
  pdf_ens(:,:)    = 0.0d0

!***Aim to save PDF data***
  do i = 1, No_of_ensembles 
    call ReadPDF(FilNameOfPdf(i), NdataPDF(i), enepdf_ens(:,i), pdf_ens(:,i))
    print*,""
  enddo
!**************************

  !***Aim for allocation of coordinates of trajectory 
  more_max_Natoms = maxval(NumOfAtomsOfRef)
  !************************************************** 

  allocate(trj_x(more_max_Natoms),trj_y(more_max_Natoms),trj_z(more_max_Natoms))

  iconfs = 0

  do iEnsmbl = 1, No_of_ensembles 
    do itrj = 1, NoOfTrj(iEnsmbl)
      call count_trj(FilNameOfTrj(itrj,iEnsmbl), NumOfConfs)
      call read_and_analy_trj(FilNameOfTrj(itrj,iEnsmbl),&
                              judge_output_pdb,&
                              NumOfConfs,&
                              iconfs    ,&
                              NumOfAtomsOfRef(iEnsmbl),&
                              total_residue_all(iEnsmbl),&
                              first_atom_all(:, iEnsmbl),&
                              last_atom_all(:,iEnsmbl),&
                              cell_ens(:,iEnsmbl),&
                              NumOfAromsOfProteinOfRef(iEnsmbl) ,&
                              AtomNum_ref_all(:,iEnsmbl),&
                              AtomName_ref_all(:,iEnsmbl),&
                              ResName_ref_all(:,iEnsmbl),&
                              ChainID_ref_all(:,iEnsmbl),&
                              ResNum_ref_all(:,iEnsmbl),&
                              Tag_AtomNumForRmsd,&
                              Ref_AtomNumForRmsd(:,iEnsmbl),&
                              Ref_AtomNumForVCV(:,iEnsmbl),&
                              N_RmsdAtoms,&  !of Tag 
                              N_VCVAtomsRef(iEnsmbl),&
                              NumberOfAtomPair,&
                              trj_x, trj_y, trj_z,&
                              CoordForRmsd, &
                              TotalProb, NdataPDF(iEnsmbl), enepdf_ens(:,iEnsmbl), pdf_ens(:,iEnsmbl), & 
                              Npairs, q, qq, vcv, & 
                              atm_no_ete, targt_rsidu_atom, npair_for_targt_rsidu_atom)
    enddo
    !@@@This is the sign for pinoint where the end of conformation no. for each is.
    !write(74,"('Reading of ',i2,'th ensenmble ends.')"), iEnsmbl
    !@@@
  enddo

  write(*,'("Total Prob = ", f16.3)') TotalProb

  q(:)    = q(:) / TotalProb  
  qq(:,:) = qq(:,:) / TotalProb
  !******************************************************************************** 

  call MakeVCV(Npairs,q,qq,vcv) !vcv is output 

  !***Output Variance-covariance matrix
  write(62,'(i5)') Npairs
  do i = 1, Npairs
    do j = 1, Npairs
      write(62,'(i6,2x,i6,2x,e16.9)') j,i,vcv(j,i)
    enddo
  enddo
  !************************************ 

  !***Output <q> which is used in projection -->  < q_i-<q> >•v
  do i = 1, Npairs
    write(63,'(i6,2x,e16.9)') i, q(i)
  enddo
  !************************************************************

  call cpu_time(tf)
  write(*,*)
  write(*,'("###########################################")')
  write(*,'("Total CPU time: ",f10.5, 1x, "sec")') tf - ti 
  write(*,'("###########################################")')
end program
