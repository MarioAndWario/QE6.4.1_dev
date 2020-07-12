!
! Copyright (C) 2006-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE write_p_avg(filp, spin_component, firstk, lastk, nvb, ncb, isbare, verbo, output_ppsi, restart_with_ppsi, qshift_, ppsi_v)
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at, bg, ibrav
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : ngm, g, ngm_g, ig_l2g, mill
  USE lsda_mod,             ONLY : nspin
  USE ener,                 ONLY : ef
  USE wvfct,                ONLY : et, nbnd, npwx, wg
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE uspp,                 ONLY : nkb, vkb, okvan
  USE becmod,               ONLY : bec_type, becp, calbec, &
       allocate_bec_type, deallocate_bec_type
  USE noncollin_module,     ONLY : noncolin, npol
  USE ldaU,                 ONLY : lda_plus_u
  USE wavefunctions, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast, mp_sum, mp_barrier, mp_max, mp_get
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_world,             ONLY : world_comm, mpime, nproc  
  USE mp_pools,             ONLY : intra_pool_comm
  USE mp_wave, ONLY : mergewf
  ! USE mp_wave_spinor, ONLY : mergewf_spinor
  use hdf5
  use h5lt

  IMPLICIT NONE
  INCLUDE 'mpif.h'

  INTEGER :: spin_component, nks1, nks2, firstk, lastk, npw
  INTEGER :: iunout, ios, ik, ibnd, jbnd, ipol, nbnd_occ
  COMPLEX(DP) :: zdotc
  COMPLEX(DP), ALLOCATABLE :: ppsi(:,:), ppsi_us(:,:), matp(:,:,:), ppsi_(:,:)
  REAL(DP), ALLOCATABLE :: ppsi_out(:,:,:,:,:), ppsi_us_out(:,:,:,:,:), ppsi_in(:,:,:,:,:), ppsi_us_in(:,:,:,:,:)
  CHARACTER (len=256) :: filp, namefile
  integer, intent(in) :: nvb, ncb
  integer :: icb, ivb
  logical, intent(in) :: isbare, verbo, output_ppsi, restart_with_ppsi, ppsi_v
  real(DP), intent(in) :: qshift_(3)
  logical :: ef_fixed = .TRUE.
  integer :: ik_local, ib, nspin0
  COMPLEX(DP), ALLOCATABLE :: matp_k(:,:,:,:) ! matp_k(iv,ic,ik,ipol)
  integer, allocatable :: g_g ( :, : ), igk_l2g ( :, : ), igk_l2g_ket(:,:), ngk_ket(:)
  integer :: ig, id, npwx_g
  integer, allocatable :: ngk_g ( : )
  real(DP), allocatable :: xk_ket(:,:)
  real(DP) :: umk(3,nkstot)
  integer :: nvnc, nvnc_ket, jbnd_
  real(DP), allocatable :: et_ket(:,:)
  integer :: g_ket(3), g_wfng(3), ig_ppsi_ket, ig_g_g_ket, diff
  integer, allocatable :: map_ppsi2wfng_l(:), map_ik2ikket(:), igwk(:), g_useful(:,:), itmp(:), igwf_l2g ( : )
  integer :: npw_g, npwx_ket, nkstot_ket, ik_ket, jj, nbnd_ket, j, ngg, igwx, ig_wfng, ig_ket
  complex (DP), allocatable, TARGET :: wfng (:)
  complex (DP), POINTER :: wfng2 (:)
  integer, allocatable :: npwx_ket_list(:)
  real(DP) :: delta, qq_temp(3), TOL_Small
  integer :: npwx_max, npwx_max_ket, ispinor, nproc_ket=0
  CHARACTER(LEN=30) :: filename
  INTEGER(HID_T) :: file_id       ! File identifier
  integer(HID_T) :: dset_id
  INTEGER(HID_T) :: plist_id      ! Property list identifier
  integer(HID_T) :: group_id
  INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
  integer(HID_T) :: dataspace
  integer(HID_T) :: memspace
  integer(HSIZE_T) :: dataspace_ngk(1), dataspace_igk(3), dataspace_ppsi(4), dataspace_ppsi_us(4)
  integer, parameter :: rank_npwx = 1, rank_ngk = 2, rank_igk = 3, rank_ppsi = 7, rank_ppsi_us = 7
  integer(HSIZE_T) :: count_npwx(1), offset_npwx(1), dim_npwx(1)
  integer(HSIZE_T) :: count_ngk(2),  offset_ngk(2),  dim_ngk(2)
  integer(HSIZE_T) :: count_igk(3),  offset_igk(3),  dim_igk(3)
  integer(HSIZE_T) :: count_ppsi(7),  offset_ppsi(7),  dim_ppsi(7)
  integer(HSIZE_T) :: count_ppsi_us(7),  offset_ppsi_us(7),  dim_ppsi_us(7)
  integer :: error, info
  integer :: comm
  logical :: flag_use_ik_ket
  
  comm = MPI_COMM_WORLD
  info = MPI_INFO_NULL

  !< Open HDF interface
  call h5open_f(error)

  flag_use_ik_ket = .FALSE.

  TOL_Small = 1.0D-8
  nvnc_ket = -1
  nvnc = nvb + ncb

  IF (lda_plus_u) THEN
     if (ionode) then
        write(*,'(5X,A)') "[WARNING] write_p_avg() not working with LDA+U!"
     endif
     ! CALL errore('write_p_avg', &
     !      'write_p_avg not working with LDA+U',1)
  ENDIF

  CALL allocate_bec_type ( nkb, nbnd, becp)

  IF (nspin==1.or.nspin==4) THEN
     nks1=max(1,firstk) ! = 1
     nks2=min(nkstot, lastk) ! = nkstot = # of total kpoints
     IF (spin_component /= 1)  &
          CALL errore('write_p_avg','incorrect spin_component',1)
  ELSEIF (nspin==2) THEN
     IF (spin_component == 1) THEN
        nks1=max(1,firstk)
        nks2=min(nks/2,lastk)
     ELSEIF (spin_component==2) THEN
        nks1=nks/2 + max(1,firstk)
        nks2=nks/2 + min(nks/2,lastk)
     ELSE
        CALL errore('write_p_avg','incorrect spin_component',1)
     ENDIF
  ENDIF

  ios = 0
  IF ( ionode ) THEN
     iunout=58
     namefile=trim(filp)
     OPEN (unit = iunout, file = namefile, status = 'unknown', form = &
          'formatted', iostat = ios)
     REWIND (iunout)
  ENDIF

  CALL mp_bcast (ios, ionode_id, world_comm)
  IF ( ios/=0 ) CALL errore ('write_p_avg', 'Opening filband file', abs (ios) )

  ! ======
  ! Modules/io_files.f90
  ! iunwfc      = 10

  ! unit_igk_l = 100 + mpime*2
  ! unit_igk_l_formatted = unit_igk_l + 1

  ! unit_ppsi_vc_l = 1000 + mpime*2
  ! unit_ppsi_vc_l_formatted = unit_ppsi_vc_l+1

  !> Find npwx_max=MAX(npwx) across all the proc
  npwx_max=npwx
  call mp_max(npwx_max, world_comm)

  if (output_ppsi) then
     CALL mp_barrier ( world_comm )

     if (ppsi_v) then
        filename="pv.h5"
     else
        filename="pc.h5"
     endif

     if (ionode) then
        !> ROOT create hdf5 file in serial
        call h5fcreate_f(trim(adjustl(filename)), H5F_ACC_TRUNC_F, file_id, error)

        call h5gcreate_f(file_id, '/header', group_id, error)
        call h5gclose_f(group_id, error)
        call h5gcreate_f(file_id, '/kgspace', group_id, error)
        call h5gclose_f(group_id, error)
        call h5gcreate_f(file_id, '/ppsi', group_id, error)
        call h5gclose_f(group_id, error)
        !> we will output all the metadata to /header
        !> meta data for kspace: nbnd, nkstot
        !> meta data for gspace: npwx, nproc
        !> meta data for ppsi: nvnc

        !> need to check nproc each time!
        call hdf5_write_int(file_id,'/header/nbnd',nbnd,error)
        call hdf5_write_int(file_id,'/header/nkstot',nkstot,error)
        call hdf5_write_int(file_id,'/header/npwx_max',npwx_max,error)
        !> npwx is the same within a pool
        call hdf5_write_int(file_id,'/header/nvnc',nvnc,error)
        call hdf5_write_int(file_id,'/header/nproc',nproc,error)
        !> from cartesian to crystal
        call cryst_to_cart ( nkstot, xk, at, - 1 )
        call hdf5_write_double_array(file_id, '/header/xk', (/3,nkstot/), xk, error)
        !> from crystal to cartesian
        call cryst_to_cart ( nkstot, xk, bg, 1 )
        call hdf5_write_double_array(file_id, '/header/et', (/nbnd,nkstot/), et, error)
        !> npwx is different on different proc
        !> npwx_max is the same on different proc
        call hdf5_create_dset(file_id, '/kgspace/npwx', H5T_NATIVE_INTEGER,(/nproc/), error)
        call hdf5_create_dset(file_id, '/kgspace/ngk', H5T_NATIVE_INTEGER,(/nkstot, nproc/), error)
        call hdf5_create_dset(file_id, '/kgspace/igk_l2g', H5T_NATIVE_INTEGER,(/npwx_max, nkstot, nproc/), error)
        call hdf5_create_dset(file_id, '/ppsi/ppsi', H5T_NATIVE_DOUBLE, (/2, npwx_max, npol, nvnc, 3, nkstot, nproc/), error)
        if (okvan) then
           call hdf5_create_dset(file_id, '/ppsi/ppsi_us', H5T_NATIVE_DOUBLE, (/2, npwx_max, npol, nvnc, 3,nkstot, nproc/), error)
        endif

        CALL h5fclose_f(file_id, error)
     endif ! ionode
  endif

  if ((.not. output_ppsi)) then
     if (ionode) then
        allocate ( matp_k(nvb, ncb, nkstot,3) )
        if (ppsi_v) then
           if (isbare) then
              OPEN (unit = 11, file = "vmt.bare.cpv.x.dat", status = 'replace', form = 'unformatted', iostat = ios)
              OPEN (unit = 22, file = "vmt.bare.cpv.y.dat", status = 'replace', form = 'unformatted', iostat = ios)
              OPEN (unit = 33, file = "vmt.bare.cpv.z.dat", status = 'replace', form = 'unformatted', iostat = ios)
           else
              OPEN (unit = 11, file = "vmt.notbare.cpv.x.dat", status = 'replace', form= 'unformatted', iostat = ios)
              OPEN (unit = 22, file = "vmt.notbare.cpv.y.dat", status = 'replace', form= 'unformatted', iostat = ios)
              OPEN (unit = 33, file = "vmt.notbare.cpv.z.dat", status = 'replace', form= 'unformatted', iostat = ios)
           endif
        else
           if (isbare) then
              OPEN (unit = 11, file = "vmt.bare.vpc.x.dat", status = 'replace', form = 'unformatted', iostat = ios)
              OPEN (unit = 22, file = "vmt.bare.vpc.y.dat", status = 'replace', form = 'unformatted', iostat = ios)
              OPEN (unit = 33, file = "vmt.bare.vpc.z.dat", status = 'replace', form = 'unformatted', iostat = ios)
           else
              OPEN (unit = 11, file = "vmt.notbare.vpc.x.dat", status = 'replace', form= 'unformatted', iostat = ios)
              OPEN (unit = 22, file = "vmt.notbare.vpc.y.dat", status = 'replace', form= 'unformatted', iostat = ios)
              OPEN (unit = 33, file = "vmt.notbare.vpc.z.dat", status = 'replace', form= 'unformatted', iostat = ios)
           endif
        endif
        nspin0 = nspin
        if (nspin .eq. 4) then
           write(*,'(A)') "We have spinor, set nspin = 1"
           nspin0 = 1
        elseif (nspin .eq. 2) then
           write(*,'(A)') "We only support nspin = 1 or 4"
           call exit(2)
        endif
        write(11) nspin0, nvb, ncb, nkstot
        write(22) nspin0, nvb, ncb, nkstot
        write(33) nspin0, nvb, ncb, nkstot
     endif

     if (restart_with_ppsi) then
        flag_use_ik_ket = .TRUE.
        CALL mp_barrier ( world_comm )
        if (ppsi_v) then
           filename="pv.h5"
        else
           filename="pc.h5"
        endif

        if (ionode) then
           call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error)
           call hdf5_read_int(file_id, '/header/nbnd', nbnd_ket, error)
           call hdf5_read_int(file_id, '/header/nkstot', nkstot_ket, error)
           call hdf5_read_int(file_id, '/header/npwx_max', npwx_max_ket, error)
           call hdf5_read_int(file_id, '/header/nvnc', nvnc_ket, error)
           call hdf5_read_int(file_id, '/header/nproc', nproc_ket, error)
           call h5fclose_f(file_id, error)

           if (nproc_ket .ne. nproc) then
              write(*,'(A)') "Mismatch of total number of processes between two calculations!"
              call exit(234)
           endif

           if (ppsi_v) then
              if (nvb .ne. nvnc_ket) then
                 write(*,'(A,I5,A,I5,A,I5)') "Mismatch of number of valence bands between two calculations! nvb = ", nvb, " nvnc_ket = ", nvnc_ket, "mpime = ", mpime
                 call exit(234)
              endif
           else
              if (ncb .ne. nvnc_ket) then
                 write(*,'(A,I5,A,I5)') "Mismatch of number of conduction bands between two calculations! nvb = ", nvb, " nvnc_ket = ", nvnc_ket
                 call exit(234)
              endif
           endif

           if (nkstot_ket .ne. nkstot) then
              write(*,'(A)') "Mismatch of number of kpoints between two calculations!"
              call exit(234)
           endif
        endif

        !> Broadcast nbnd_ket, nkstot_ket, npwx_max_ket, nvnc_ket, nproc_ket
        call mp_bcast(nbnd_ket, ionode_id, world_comm)
        call mp_bcast(nkstot_ket, ionode_id, world_comm)
        call mp_bcast(npwx_max_ket, ionode_id, world_comm)
        call mp_bcast(nvnc_ket, ionode_id, world_comm)
        call mp_bcast(nproc_ket, ionode_id, world_comm)
        CALL mp_barrier ( world_comm )

        allocate(npwx_ket_list(nproc_ket))

        if (ionode) then
           ! << read npwx_ket and broadcast to every proc
           call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error)
           call hdf5_read_int_array(file_id, '/kgspace/npwx',(/nproc_ket/), npwx_ket_list, error)
           call h5fclose_f(file_id, error)
        endif

        call mp_bcast(npwx_ket_list, ionode_id, world_comm)

        npwx_ket=npwx_ket_list(mpime+1)

        allocate(xk_ket(3,nkstot_ket))
        allocate(et_ket(nbnd_ket,nkstot_ket))

        !> xk_ket in crystal coordiantes
        !> while original xk in cartesian coordiantes
        if (ionode) then
           call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error)
           call hdf5_read_double_array(file_id, '/header/xk', (/3,nkstot_ket/), xk_ket, error)
           call hdf5_read_double_array(file_id, '/header/et', (/nbnd_ket,nkstot_ket/), et_ket, error)
           call h5fclose_f(file_id, error)
        endif

        call mp_bcast(xk_ket, ionode_id, world_comm)
        call mp_bcast(et_ket, ionode_id, world_comm)
        CALL mp_barrier ( world_comm )

        !> Check that xk and xk_ket match with each other, and find out the Umklapp vector for each pair
        allocate(map_ik2ikket(nkstot))

        !> from cartesian to crystal
        CALL cryst_to_cart ( nkstot, xk, at, - 1 )

        !> [IMPORTANT]
        !> Find mapping between kgq and kg
        do ik = 1, nkstot
           ik_ket = 0
           delta = 0.1d0
           ! For each kc, find the corresponding kv, such that, kv = kc - Q
           do while ((delta .gt. TOL_Small) .and. (ik_ket .lt. nkstot_ket))
              ik_ket = ik_ket+1
              ! ======
              ! kv = kc - Q - umk, kv in FBZ, kc in FBZ, Q in FBZ, kc - Q could be outside of FBZ
              ! kv <-- ik_ket
              ! kv <-- ik

              ! if (ppsi_v) then ! < c k_c | p | v k_v >, umk = kc - Q - kv
              !    qq_temp(:) = xk(:,ik) - qshift_(:) - xk_ket(:,ik_ket)
              ! else ! < v k_v | p | c k_c >, -umk = kv + Q - kc
              !    qq_temp(:) = xk(:,ik) - qshift_(:) - xk_ket(:,ik_ket)
              ! endif

              qq_temp(:) = xk(:,ik) - qshift_(:) - xk_ket(:,ik_ket)

              ! ======
              ! The criterion that a point is in the wigner-seitz FBZ is the minimal distance from origin
              ! ------
              ! ANINT(A [, KIND]) rounds its argument to the nearest whole number.
              ! so -0.5 <= qq_temp(i) <= 0.5, for any i
              do jj=1,3
                 qq_temp(jj) = qq_temp(jj) - anint( qq_temp(jj) )
              enddo

              delta = sqrt((qq_temp(1))**2+(qq_temp(2))**2+(qq_temp(3))**2)
           enddo

           if (delta .gt. TOL_Small) then
              if(ionode) then
                 write(0,*) '  Could not find kpoint #', ik, " qq_temp = ", qq_temp, " xk = ", xk(:,ik), " qshift_ = ", qshift_, " xk_ket = ", xk_ket(:,ik_ket)
                 write(0,*) ' Did you forget to set qshift?'
              endif
              call exit(421)
           endif
           ! ======
           ! kg%f(:,ik) - xct%finiteq(:) = kgq%f(:,xct%indexq(ik)) + xct%umk(:,ik)
           ! ------
           ! [debug]
           ! we will relax this condition!
           ! ------

           ! if ((ik_ket .ne. ik) .or. (ik_ket .eq. 0)) then
           !    write(*,'(A)') "Mismatch of indices of kpoints between two calculations!"
           !    call exit(234)
           ! endif

           map_ik2ikket(ik) = ik_ket
           !> [important]
           !> this umk means Kump when ppsi_v = T and -Kump when ppsi_v = F
           umk(:,ik) = nint( xk(:,ik) - xk_ket(:,ik_ket) - qshift_(:) )
           if (ionode) then
              write(*,'(A,I0,A,3F8.4,A,I0,A,3F8.4,A,I0,A,3F8.4,A)') ">>> k(:,", ik, ")=(", xk(:,ik), ") k'(:,", ik_ket, ")=(", xk_ket(:,ik_ket), ") umk(:,", ik ,")=(", umk(:,ik), ")"
           endif
        enddo

        ! from crystal to cartesian
        CALL cryst_to_cart ( nkstot, xk, bg, 1 )

        !< read npwx_ket for each proc
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
        call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        allocate(ngk_ket(nkstot_ket))
        ! allocate(igk_l2g_ket(npwx_ket,nkstot_ket))
        ! ======
        ! [important] here we use npwx_max_ket instead of npwx_ket
        allocate(igk_l2g_ket(npwx_max_ket,nkstot_ket))

        ngk_ket = 0
        igk_l2g_ket = 0

        ! read (unit_igk_l) ngk_ket(1:nkstot_ket)
        ! read (unit_igk_l) igk_l2g_ket(1:npwx_ket,1:nkstot_ket)

        call hdf5_read_int_hyperslab(file_id, '/kgspace/ngk', (/nkstot_ket,1/), (/0,mpime/), ngk_ket, error)

        call hdf5_read_int_hyperslab(file_id, '/kgspace/igk_l2g', (/npwx_max_ket,nkstot_ket,1/), (/0,0,mpime/), igk_l2g_ket, error)

        ! write(*,*) " ngk_ket(:) = ", ngk_ket(:)
        call h5fclose_f(file_id, error)

     endif
  endif

  if (output_ppsi .or. restart_with_ppsi) then
     ! We will initiate the g_g array for the global Gvectors and igk_l2g index
     !! -----------------------------------
     !! g_g(1:3, 1:ngm_g) gives the fractional coordinates for
     !! all (1:ngm_g) the G vectors, global variable, have a duplicate
     !! at each cpu
     ALLOCATE ( g_g ( 3, ngm_g ) )

     !! ngm_g : global number of Gvectors (summed over procs)
     !! in serial execution, ngm_g = ngm
     DO ig = 1, ngm_g
        DO id = 1, 3
           g_g ( id, ig ) = 0
        ENDDO
     ENDDO

     DO ig = 1, ngm
        g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
        g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
        g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
     ENDDO

     CALL mp_sum ( g_g, intra_pool_comm )

     ! ======
     ! ALLOCATE ( igk_l2g ( npwx, nkstot ) )
     ! ------
     ! [debug]
     ALLOCATE ( igk_l2g ( npwx_max, nkstot ) )

     !! g_g(1:3, igk_l2g(ig, ik)) gives the Gvector for the ig-th Gvector
     !!    of the local ik-th kpoint,
     !! KS(ik) = [\sum_{n=1}^{ik-1} ngk(n)]+1
     !! KE(ik) = [\sum_{n=1}^{ik} ngk(n)]
     !! evc(KS(ik):KE(ik),ibnd) gives all the PW components for (ik,ibnd) indices
     !! g_g(1:3, igk_l2g(KS(ik):KE(ik), ik)) gives all the Gvectors for the local
     !! ik-th kpoint

     DO ik = 1, nkstot
        npw = ngk ( ik )
        DO ig = 1, npw ! this is ig on this proc, which is incomplete
           igk_l2g ( ig, ik ) = ig_l2g ( igk_k (ig, ik) )
        ENDDO
        ! DO ig = npw + 1, npwx
        DO ig = npw + 1, npwx_max
           igk_l2g ( ig, ik ) = 0
        ENDDO
     ENDDO

     !! ---------- Set ngk_g(1:nkstot) --------
     !! compute the global number of GVectors for each k point
     !! ngk(1:nks) contains the number of Gvectors for kpoints in a pool
     !! ngk_g(1:nkstot) contains the number of Gvectors for all the kpoints
     !! ngk_g(iks,ike) contains the number of Gvectors for the kpoints in
     !! current pool
     ALLOCATE ( ngk_g ( nkstot ) )

     !! ngk_g = 0
     !! ngk_g(iks:ike) = ngk(1:nks)
     !! CALL mp_sum( ngk_g, inter_pool_comm )
     !! CALL mp_sum( ngk_g, intra_pool_comm )
     !! ngk_g = ngk_g / nbgrp
     DO ik = 1, nkstot
        ngk_g ( ik ) = 0
     ENDDO

     DO ik = 1, nkstot
        ngk_g ( ik ) = ngk ( ik )
     ENDDO

     CALL mp_sum ( ngk_g, world_comm )

     ! Compute the maximum global Gvector INDEX among all the kpoints
     npw_g = MAXVAL ( igk_l2g ( :, : ) )

     CALL mp_max ( npw_g, world_comm )

     ! Compute the maximum number of Gvector owned by one kpoints among all kpoints
     npwx_g = MAXVAL ( ngk_g ( : ) )

     if (restart_with_ppsi) then
        ALLOCATE ( igwk ( npwx_g ) )
        ALLOCATE (g_useful(3,npwx_g))
     endif

     if (output_ppsi) then
        ! ! output igk_l2g into "igk_l2g.dat"
        ! write (unit_igk_l) ngk (1:nkstot )
        ! write (unit_igk_l) igk_l2g(1:npwx, 1:nkstot)

        ! ! ======
        ! do ik = 1, nkstot
        !    write (unit_igk_l_formatted,'(A, I5, A, I5)') " ik = ", ik, " ngk = ", ngk(ik)
        !    write (unit_igk_l_formatted,'(10I8)') igk_l2g(1:npwx,ik)
        ! enddo
        ! ! ------

        ! < output npwx, ngk and igk_l2g to hdf5 file
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
        call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        ! << output npwx
        !count_npwx(1) = 1
        !offset_npwx(1) = mpime ! starting from 0
        !dim_npwx(1)=nproc
        call hdf5_write_int_hyperslab(file_id, '/kgspace/npwx', (/1/), (/mpime/), (/npwx/),error)

        ! ! this portion of the data in the memory will be sent
        ! ! call h5screate_simple_f(rank_npwx, count_npwx, memspace, error)
        ! ! since local array will be sent in whole, no hyperslab is needed
        ! call h5screate_simple_f(rank_npwx, count_npwx, memspace, error)
        ! ! call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, offset_ngk, count_ngk, error)

        ! call h5dopen_f(file_id, '/kgspace/npwx', dset_id, error)
        ! ! select hyperslab in the file
        ! call h5dget_space_f(dset_id,filespace,error)
        ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_npwx, count_npwx, error)
        ! ! create property list for collective dataset write
        ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! ! write the dataset collectively, the fourth argument is the size of the data in the file
        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, (/npwx/), dim_npwx, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp = plist_id)
        ! call h5pclose_f(plist_id, error)
        ! call h5dclose_f(dset_id, error)
        ! call h5sclose_f(memspace, error)
        ! call h5sclose_f(filespace, error)
        ! >>

        ! << output ngk --> ngk(nkstot,nproc)
        ! count_ngk(1) = nkstot
        ! count_ngk(2) = 1

        ! offset_ngk(1) = 0
        ! offset_ngk(2) = mpime ! starting from 0

        ! dim_ngk(1) = nkstot
        ! dim_ngk(2) = nproc

        call hdf5_write_int_hyperslab(file_id, '/kgspace/ngk', (/nkstot,1/), (/0,mpime/), ngk,error)

        ! call h5screate_simple_f(rank_ngk, count_ngk, memspace, error)

        ! call h5dopen_f(file_id, '/kgspace/ngk', dset_id, error)
        ! ! select hyperslab in the file
        ! call h5dget_space_f(dset_id,filespace,error)
        ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_ngk, count_ngk, error)
        ! ! create property list for collective dataset write
        ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! ! write the dataset collectively, the fourth argument is the size of the data in the file
        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ngk, dim_ngk, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp = plist_id)
        ! call h5pclose_f(plist_id, error)
        ! call h5dclose_f(dset_id, error)
        ! call h5sclose_f(memspace, error)
        ! call h5sclose_f(filespace, error)
        ! >>

        ! << output igk_l2g --> igk_l2g(npwx,nkstot,nproc)
        ! count_igk(1) = npwx_max
        ! count_igk(2) = nkstot
        ! count_igk(3) = 1

        ! offset_igk(1) = 0
        ! offset_igk(2) = 0 ! starting from 0
        ! offset_igk(3) = mpime ! starting from 0

        ! dim_igk(1) = npwx_max
        ! dim_igk(2) = nkstot
        ! dim_igk(3) = nproc

        call hdf5_write_int_hyperslab(file_id, '/kgspace/igk_l2g', (/npwx_max,nkstot,1/), (/0,0,mpime/), igk_l2g, error)

        ! call h5screate_simple_f(rank_igk, count_igk, memspace, error)

        ! call h5dopen_f(file_id, '/kgspace/igk_l2g', dset_id, error)
        ! ! select hyperslab in the file
        ! call h5dget_space_f(dset_id,filespace,error)
        ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_igk, count_igk, error)
        ! ! create property list for collective dataset write
        ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! ! write the dataset collectively, the fourth argument is the size of the data in the file
        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, igk_l2g, dim_igk, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp = plist_id)
        ! call h5pclose_f(plist_id, error)
        ! call h5dclose_f(dset_id, error)
        ! call h5sclose_f(memspace, error)
        ! call h5sclose_f(filespace, error)
        ! >>

        call h5fclose_f(file_id, error)
     endif
  endif

  if ((.not. output_ppsi)) then
     ALLOCATE(matp(nvb,ncb,3))
  endif

  IF (ABS(ef) > 1.0E-6) THEN
     ef_fixed = .FALSE.
     IF ( ionode ) THEN
        write(*,*)
        write(*,'(5X,A)') "Fermi level determined by smearing!"
     ENDIF
  ELSE
     ef_fixed = .TRUE.
     ! ======
     ! Decide nbnd_occ by wg

     if (nspin .eq. 1) then
        nbnd_occ = nint(sum(wg(:,:))/2.0)
     else
        nbnd_occ=nint(sum(wg(:,:)))
     endif

     IF ( ionode ) THEN
        write(*,*)
        write(*,'(5X,A)') "Fermi level determiend by fixed occupation!"
        WRITE(*,'(5X,A,I5)') "nbnd_occ =", nbnd_occ
     ENDIF
  ENDIF

  ! ======
  ! [important]
  ! < Allocate ppsi's
  if ( (.not. restart_with_ppsi) ) then
     ALLOCATE(ppsi(npwx*npol,nvnc))
     ppsi = 0

     if (output_ppsi) then
        allocate(ppsi_out(2,npwx_max,npol,nvnc,3)) ! for each ik and iproc
        ppsi_out = 0
     endif

     IF (okvan) THEN
        ALLOCATE(ppsi_us(npwx*npol,nvnc))
        ppsi_us = 0
        if (output_ppsi) then
           ALLOCATE(ppsi_us_out(2,npwx_max,npol,nvnc,3))
           ppsi_us_out = 0
        endif
     ENDIF
  else ! restart_with_ppsi = T
     ALLOCATE(ppsi(npwx_ket * npol, nvnc_ket))
     ppsi = 0
     allocate(ppsi_in(2,npwx_max_ket,npol,nvnc_ket,3))
     ppsi_in = 0

     ALLOCATE(ppsi_(npwx*npol,nvnc))
     ppsi_ = 0

     IF (okvan) THEN
        ALLOCATE(ppsi_us(npwx_ket * npol, nvnc_ket))
        ppsi_us = 0
        allocate(ppsi_us_in(2,npwx_max_ket,npol,nvnc_ket,3))
        ppsi_us_in = 0
     ENDIF
  endif

  if (verbo .and. (.not. output_ppsi)) then
     if (ppsi_v) then
        write(6,'(A)') "     Write < ik, cb | p | ik_ket, vb >"
        write(6,'(A)') "      ik      cb     ik_ket   vb Pol               Re<p>               Im<p>                |<p>|^2"
        write(6,'(A)')
     else ! <- ppsi_v = F ->
        write(6,'(A)') "     Write < ik, vb | p | ik_ket, cb >"
        write(6,'(A)') "      ik      vb     ik_ket   cb Pol               Re<p>               Im<p>                |<p>|^2"
        write(6,'(A)')        
     endif
  endif

  DO ik = nks1, nks2

     matp(:,:,:) = 0.0D0

     !
     !   Compute the number of occupated bands at this k point
     !
     IF (.NOT. ef_fixed) THEN
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik)<=ef) THEN
              nbnd_occ = ibnd
           ENDIF
        ENDDO
     ENDIF

     IF (nbnd_occ==nbnd) WRITE( stdout, '(5x,/,&
          &"No empty band at point ", i4,3f10.5)') &
          ik,  (xk (ipol, ik) , ipol = 1, 3)


     npw = ngk(ik)

     CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)

     if (restart_with_ppsi) then

        ik_ket = map_ik2ikket(ik)

        ! write(*,*) "ik = ", ik, " ik_ket = ", ik_ket

        ! We should collect evc from each proc to have a evc_tot which includes all the Gvectors, just in case the Gvector needed by ppsi_ket on the i-th proc is not available in evc on the i-th proc
        DO j = 1, npwx_g
           igwk ( j ) = 0
        ENDDO

        ALLOCATE ( itmp ( npw_g ) )
        itmp = 0

        !! ... If the kpoint is in local pool
        DO ig = 1, ngk (ik)
           itmp ( igk_l2g (ig, ik) ) = igk_l2g (ig, ik)
        ENDDO

        CALL mp_sum( itmp, world_comm )

        ! ngg counts the number of useful Gvectors
        ngg = 0

        ! igwk(:) will map the index in g_useful(1:3,:) into the
        ! corresponding index of this useful Gvector in global
        ! list g_g(1:3, 1:ngm_g)
        DO ig = 1, npw_g
           IF ( itmp ( ig ) .EQ. ig ) THEN
              ngg = ngg + 1
              igwk ( ngg ) = ig
           ENDIF
        ENDDO

        DEALLOCATE ( itmp )

        DO ig = 1, ngk_g( ik )
           DO id = 1, 3
              g_useful(id, ig) = g_g (id, igwk(ig))
           ENDDO
        ENDDO

        ALLOCATE ( igwf_l2g ( npw ) )

        DO ig = 1, npw
           igwf_l2g ( ig ) = 0
        ENDDO

        DO ig = 1, npw
           ! ngg stores the global index for each local Gvector
           ngg = igk_l2g ( ig, ik )
           DO j = 1, ngk_g ( ik )
              IF ( ngg .EQ. igwk ( j ) ) THEN
                 ! igwf_l2g(:) is again specific to one kpoint ik_global
                 ! igwf_l2g(:) will map the ig-th local g index into
                 ! its corresponding index in g_useful(:) list
                 igwf_l2g ( ig ) = j
                 EXIT
              ENDIF
           ENDDO
        ENDDO
        ! igwx is a pool-wide variable, that is, it has the same value
        ! in all cpus of one pool        
        igwx = MAXVAL ( igwf_l2g ( 1 : npw ) ) ! === ngk_g(ik_global)
        CALL mp_max ( igwx, intra_pool_comm )

        !> [WORKING]
        ! !> if ik_global is not in current pool, igwx=0*2=0
        ! IF (noncolin) THEN
        !    igwx = igwx * 2
        ! ENDIF

        !> ngk_ket(:) is for ppsi, while ngk(:) is for wfng
        allocate(map_ppsi2wfng_l(ngk_ket(ik_ket)))
        map_ppsi2wfng_l(:) = 0
        ! ======
        ! Find the mapping from ig_ppsi (ig_ket) to ig_wfng with given ik
        ! we also need to consider the umklapp vectors
        ! ------
        ! ig_wfng = [1, igwx], where igwx = nspinor * ngk_g(ik_global)
        ! g_useful(:,1:ngk_g(ik)) = g_g(:,igwk(1:ngk_g(ik)))
        ! igk_l2g(ig_evc=1:npw, ik) = ig_g_g, where npw = ngk(ik)
        ! igwf_l2g(ig_evc=1:npw) = ig_g_useful

        ! so we have igk_l2g_ket for ppsi_ket
        do ig_ppsi_ket = 1, ngk_ket(ik_ket)
           ig_g_g_ket = igk_l2g_ket(ig_ppsi_ket, ik_ket)
           g_ket(:) = g_g(:,ig_g_g_ket)
           do ig_wfng = 1, ngk_g(ik)
              g_wfng(:) = g_useful(:,ig_wfng)
              ! ======
              ! When ppsi_v = T, conj(wfng(g_ket-Gumk)) * ppsi(g_ket)
              ! that is, G = g_ket - Gumk = g_ket - umk
              ! ======
              ! When ppsi_v = F, conj(wfng(g_ket+Gumk)) * ppsi(g_ket)
              ! that is, G = g_ket + G_umk = g_ket - umk
              diff=maxval(abs(g_ket(:) - g_wfng(:) - umk(:,ik)))
              if (diff .eq. 0) then
                 !! map_ppsi2wfng_l(ig_for_ppsi) = ig_for_wfng
                 map_ppsi2wfng_l(ig_ppsi_ket) = ig_wfng
              endif
           enddo
        enddo

        !> read in ppsi_out and ppsi_us_out for the ik_ket-th kpoint and all the polarization directions
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
        call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        call hdf5_read_double_hyperslab(file_id, '/ppsi/ppsi', (/2, npwx_max_ket, npol, nvnc_ket, 3, 1, 1/), (/0,0,0,0,0,ik_ket-1,mpime/), ppsi_in, error)

        if (okvan) then
           allocate(ppsi_us_in(2,npwx_max_ket,npol,nvnc_ket,3))
           call hdf5_read_double_hyperslab(file_id, '/ppsi/ppsi_us', (/2, npwx_max_ket, npol, nvnc_ket, 3, 1, 1/), (/0,0,0,0,0,ik_ket-1,mpime/), ppsi_us_in, error)
        endif

        call h5fclose_f(file_id, error)

     elseif (.not. output_ppsi) then ! <- restart_with_ppsi = F .and. output_ppsi = F -> in this case ppsi_v = T and ik = ik_ket
        ik_ket = ik
     endif

     !> from cartesian to crystal
     CALL cryst_to_cart ( nkstot, xk, at, - 1 )
     if (ionode .and. (.not. output_ppsi)) then
        if (restart_with_ppsi) then
           write(6,'(A,I5,A,"(",3F12.7,")",A,I5,A,"(",3F12.7,")")') ">>>   ik = ", ik, " <k| = ", xk(:,ik), " ik_ket = ", ik_ket, " |k'> = ", xk_ket(:,ik_ket)
        else
           write(6,'(A,I5,A,3F12.7,A,I5,A,3F12.7)') "ik = ", ik, " <k| = ", xk(:,ik), " ik_ket = ", ik_ket, " |k'> = ", xk(:,ik_ket)
        endif
     endif

     !> from crystal to cartesian
     CALL cryst_to_cart ( nkstot, xk, bg, 1 )    
     CALL calbec ( npw, vkb, evc, becp)

     ! The first index of ppsi is arranged in the order of |k+G|^2 and scattered over procs
     DO ipol = 1, 3 ! note that here ipol means 3D cartesian directions, while the npol=2 means spinor wavefunctions.

        CALL mp_barrier ( world_comm )

        ! <- restart_with_ppsi = F ->
        if ( (.not. restart_with_ppsi) ) then
           ! Compute \hat{p} | v > and \hat{p} | c >
           CALL compute_ppsi_vc(ppsi, ppsi_us, ik, ipol, nbnd_occ, nvb, ncb, spin_component)
           ! <- output_ppsi = T ->
           if (output_ppsi) then
              do ispinor = 1, npol
                 do ig = 1, npwx
                    do ib = 1, nvnc
                       ppsi_out(1,ig,ispinor,ib,ipol) = dble(ppsi(ig*ispinor,ib))
                       ppsi_out(2,ig,ispinor,ib,ipol) = DIMAG(ppsi(ig*ispinor,ib))
                    enddo
                 enddo
              enddo
              if (okvan) then
                 do ispinor = 1, npol
                    do ig = 1, npwx
                       do ib = 1, nvnc
                          ppsi_us_out(1,ig,ispinor,ib,ipol) = dble(ppsi_us(ig*ispinor,ib))
                          ppsi_us_out(2,ig,ispinor,ib,ipol) = DIMAG(ppsi_us(ig*ispinor,ib))
                       enddo
                    enddo
                 enddo
              endif
              !> <- output_ppsi = F .and. restart_with_ppsi = F ->
              !> here we consider q=0 (ik_ket = ik) and <ck| p |vk>
           else
              !> Global index for evc
              DO ibnd=nbnd_occ+1,nbnd_occ+ncb
                 DO jbnd=nbnd_occ-nvb+1,nbnd_occ
                    ! [Debugging]
                    ! Relative index for ppsi, ppsi_us
                    jbnd_ = jbnd - (nbnd_occ - nvb)
                    ! Relative band indices for matp
                    ivb = nbnd_occ - jbnd + 1
                    icb = ibnd - nbnd_occ
                    IF (noncolin) THEN
                       ! ======
                       ! LAPACK: ZDOTC forms the dot product of two complex vectors
                       ! ZDOTC = X^H * Y
                       ! zdotc(N, ZX, INCX, ZY, INCY)
                       matp(ivb,icb,ipol) = zdotc(npwx*npol,evc(1,ibnd),1,ppsi(1,jbnd_),1)
                       IF (okvan) THEN
                          matp(ivb,icb,ipol)=                  &
                               matp(ivb,icb,ipol)+             &
                               (0.d0,0.5d0)*(et(ibnd,ik)-et(jbnd,ik))*  &
                               (zdotc(npwx*npol,evc(1,ibnd),1,ppsi_us(1,jbnd_),1) )
                       ENDIF
                    ELSE
                       matp(ivb,icb,ipol) = zdotc(npw,evc(1,ibnd),1,ppsi(1,jbnd_),1)
                       IF (okvan) THEN
                          matp(ivb,icb,ipol)= &
                               matp(ivb,icb,ipol) +  &
                               (0.d0,0.5d0)*zdotc(npw,evc(1,ibnd),1,ppsi_us(1,jbnd_),1)* &
                               (et(ibnd,ik)-et(jbnd,ik))
                       ENDIF
                    ENDIF

                    ! ======
                    ! BGW trunk convention:
                    ! ------
                    ! matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2 / ( et(ibnd,ik)-et(jbnd,ik) )
                    if (isbare) then
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2
                    else
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2 / ( et(ibnd,ik)-et(jbnd,ik) )
                    endif
                    ! ------
                 ENDDO ! jbnd <= ivb
              ENDDO ! ibnd <= icb
           endif
           !> Here we calculate < c k_c | p | v k_v > or < v k_v | p | c k_c > with the ppsi_ket
           !> <- restart_with_ppsi = T ->
        else
           do ispinor = 1, npol
              do ig = 1, npwx_ket
                 do ib = 1, nvnc_ket
                    ppsi(ig*ispinor,ib) = DCMPLX(ppsi_in(1,ig,ispinor,ib,ipol), ppsi_in(2,ig,ispinor,ib,ipol))
                 enddo
              enddo
           enddo

           if (okvan) then
              do ispinor = 1, npol
                 do ig = 1, npwx_ket
                    do ib = 1, nvnc_ket
                       ppsi_us(ig*ispinor,ib) = DCMPLX(ppsi_us_in(1,ig,ispinor,ib,ipol), ppsi_in(2,ig,ispinor,ib,ipol))
                    enddo
                 enddo
              enddo
           endif

           ! <- ppsi_v = T ->
           if (ppsi_v) then
              ! nvb = nvnc_ket
              ! ======
              ! Global index for evc
              ! < c |
              DO ibnd=nbnd_occ+1,nbnd_occ+ncb
                 CALL mp_barrier ( world_comm )
                 ! Collect wavefunctions into wfng
                 !> [WORKING]                 
                 ALLOCATE ( wfng(MAX ( npol*igwx, 1 )) )
                 wfng = ( 0.0D0, 0.0D0 )                 
                 if (npol == 2) wfng2 => wfng(igwx+1:2*igwx)
                 IF (npol == 2) THEN
                    CALL mergewf (evc(1:npwx, ibnd),         wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                    CALL mergewf (evc(npwx+1:2*npwx, ibnd), wfng2, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )                    
                 ELSE
                    CALL mergewf (evc(:, ibnd),              wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                 ENDIF                 

                 ! IF (noncolin) THEN
                 !    CALL mergewf_spinor( evc(:,ibnd), wfng, ngk_g(ik), npw, igwf_l2g, mpime, nproc, ionode_id, world_comm, npwx)
                 !    ! After this mergewf_spinor, each pool-ROOT proc has the PWs for ib-th band and ik_global-th kvector.
                 ! ELSE ! nspinor == 1
                 !    CALL mergewf ( evc (:,ibnd), wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                 ! ENDIF

                 !> Now every proc shares the PWs for ib-th band and ik_global-th kvector
                 call mp_bcast(wfng, ionode_id, world_comm)

                 !> | v >
                 DO jbnd = nbnd_occ-nvb+1, nbnd_occ
                    ! [Debugging]
                    ! Relative index for ppsi, ppsi_us
                    jbnd_ = jbnd - (nbnd_occ - nvb)
                    ! Relative band indices for matp
                    ivb = nbnd_occ - jbnd + 1
                    icb = ibnd - nbnd_occ

                    IF (noncolin) THEN
                       ! ======
                       ! LAPACK: ZDOTC forms the dot product of two complex vectors
                       ! ZDOTC = X^H * Y
                       ! zdotc(N, ZX, INCX, ZY, INCY)
                       ! matp(ivb,icb,ipol) = zdotc(npwx*npol,evc(1,ibnd),1,ppsi(1,jbnd_),1)
                       do ig_ket = 1, ngk_ket(ik_ket)
                          ig_wfng = map_ppsi2wfng_l(ig_ket)
                          if (ig_wfng .ne. 0) then
                             matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + conjg(wfng(ig_wfng)) * ppsi(ig_ket,jbnd_) + conjg(wfng(ig_wfng+igwx)) * ppsi(ig_ket+npwx_ket,jbnd_)
                             IF (okvan) THEN
                                matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + (0.d0,0.5d0)*(et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng)) * ppsi_us(ig_ket,jbnd_) + (0.d0,0.5d0) * (et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng+igwx)) * ppsi_us(ig_ket+npwx_ket,jbnd_)
                             ENDIF
                          endif
                       enddo
                    ELSE
                       do ig_ket = 1, ngk_ket(ik_ket)
                          ig_wfng = map_ppsi2wfng_l(ig_ket)
                          if (ig_wfng .ne. 0) then
                             matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + conjg(wfng(ig_wfng)) * ppsi(ig_ket,jbnd_)
                             IF (okvan) THEN
                                matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + (0.d0,0.5d0) * (et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng)) * ppsi_us(ig_ket,jbnd_)
                             ENDIF
                          endif
                       enddo
                    ENDIF
                    ! ======
                    ! BGW trunk convention:
                    ! ------
                    ! matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2 / ( et(ibnd,ik)-et(jbnd,ik) )
                    if (isbare) then
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2
                    else
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2 / ( et(ibnd,ik)-et_ket(jbnd,ik_ket) )
                    endif
                 enddo ! jbnd
                 DEALLOCATE ( wfng )
              enddo ! ibnd
              ! <- ppsi_v = F ->
           else
              ! ncb = nvnc_ket
              ! ======
              ! Global index for evc
              ! < v |
              DO ibnd=nbnd_occ - nvb + 1, nbnd_occ
                 CALL mp_barrier ( world_comm )
                 !> Collect wavefunctions into wfng
                 !> [WORKING]                 
                 ALLOCATE ( wfng(MAX ( npol*igwx, 1 )) )
                 wfng = ( 0.0D0, 0.0D0 )
                 if (npol == 2) wfng2 => wfng(igwx+1:2*igwx)
                 IF (npol == 2) THEN
                    CALL mergewf (evc(1:npwx, ibnd),         wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                    CALL mergewf (evc(npwx+1:2*npwx, ibnd), wfng2, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )                    
                 ELSE
                    CALL mergewf (evc(:, ibnd),              wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                 ENDIF                                  
                 ! IF (noncolin) THEN
                 !    CALL mergewf_spinor( evc(:,ibnd), wfng, ngk_g(ik), npw, igwf_l2g, mpime, nproc, ionode_id, world_comm, npwx)
                 !    ! After this mergewf_spinor, each pool-ROOT proc has the PWs for ib-th band and ik_global-th kvector.
                 ! ELSE ! nspinor == 1
                 !    CALL mergewf ( evc (:,ibnd), wfng, npw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                 ! ENDIF

                 call mp_bcast(wfng, ionode_id, world_comm)

                 !> | c >
                 DO jbnd = nbnd_occ + 1, nbnd_occ + ncb
                    ! Relative index for ppsi, ppsi_us
                    jbnd_ = jbnd - nbnd_occ
                    ! Relative band indices for matp
                    ivb = nbnd_occ - ibnd + 1
                    icb = jbnd - nbnd_occ

                    IF (noncolin) THEN
                       ! LAPACK: ZDOTC forms the dot product of two complex vectors
                       ! ZDOTC = X^H * Y
                       ! zdotc(N, ZX, INCX, ZY, INCY)
                       !! matp(ivb,icb,ipol) = zdotc(npwx*npol,evc(1,ibnd),1,ppsi(1,jbnd_),1)
                       do ig_ket = 1, ngk_ket(ik_ket)
                          ig_wfng = map_ppsi2wfng_l(ig_ket)
                          if (ig_wfng .ne. 0) then
                             matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + conjg(wfng(ig_wfng)) * ppsi(ig_ket,jbnd_) + conjg(wfng(ig_wfng+igwx)) * ppsi(ig_ket+npwx_ket,jbnd_)
                             IF (okvan) THEN
                                matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + (0.d0,0.5d0) * (et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng))*ppsi_us(ig_ket,jbnd_) + (0.d0,0.5d0) * (et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng+igwx)) * ppsi_us(ig_ket+npwx_ket,jbnd_)
                             ENDIF
                          endif
                       enddo
                    ELSE
                       do ig_ket = 1, ngk_ket(ik_ket)
                          ig_wfng = map_ppsi2wfng_l(ig_ket)
                          if (ig_wfng .ne. 0) then
                             matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + conjg(wfng(ig_wfng)) * ppsi(ig_ket,jbnd_)
                             IF (okvan) THEN
                                matp(ivb,icb,ipol) = matp(ivb,icb,ipol) + (0.d0,0.5d0) * (et(ibnd,ik)-et_ket(jbnd,ik_ket)) * conjg(wfng(ig_wfng)) * ppsi_us(ig_ket,jbnd_)
                             ENDIF
                          endif
                       enddo
                    ENDIF
                    ! ======
                    ! BGW trunk convention:
                    ! ------
                    if (isbare) then
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2
                    else
                       matp(ivb,icb,ipol) = matp(ivb,icb,ipol) * 2 / ( et(ibnd,ik)-et_ket(jbnd,ik_ket) )
                    endif
                 enddo ! jbnd
                 DEALLOCATE ( wfng )
              enddo ! ibnd
           endif
        endif
     ENDDO ! ipol

     ! ======
     ! [debug]
     ! < output ppsi_out and ppsi_us_out
     ! <- when output_ppsi = T, we already set restart_with_ppsi = F ->
     if ( output_ppsi ) then
        ! < output npwx, ngk and igk_l2g to hdf5 file
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, info, error)
        call h5fopen_f(trim(adjustl(filename)), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)
        call hdf5_write_double_hyperslab(file_id, '/ppsi/ppsi', (/2, npwx_max, npol, nvnc, 3, 1, 1/), (/0, 0, 0, 0, 0, ik-1, mpime/), ppsi_out, error)

        ! ! << output ppsi --> ppsi(1:npwx*npol,1:nvnc,1:nproc) --> ppsi_out(2,1:npwx_max,1:npol_spin,1:nvnc, 1:npolarization, 1:nproc)
        ! count_ppsi(1) = 2 ! for complex number
        ! count_ppsi(2) = npwx_max ! instead of npwx
        ! count_ppsi(3) = npol ! for nspinor
        ! count_ppsi(4) = nvnc ! ib
        ! count_ppsi(5) = 3 ! ipol = 1, 3
        ! count_ppsi(6) = 1
        ! count_ppsi(7) = 1

        ! offset_ppsi(1) = 0
        ! offset_ppsi(2) = 0
        ! offset_ppsi(3) = 0
        ! offset_ppsi(4) = 0
        ! offset_ppsi(5) = 0
        ! offset_ppsi(6) = ik-1
        ! offset_ppsi(7) = mpime ! starting from 0

        ! dim_ppsi(1) = 2
        ! dim_ppsi(2) = npwx_max
        ! dim_ppsi(3) = npol
        ! dim_ppsi(4) = nvnc
        ! dim_ppsi(5) = 3
        ! dim_ppsi(6) = nkstot
        ! dim_ppsi(7) = nproc

        ! call h5screate_simple_f(rank_ppsi, count_ppsi, memspace, error)

        ! call h5dopen_f(file_id, '/ppsi/ppsi', dset_id, error)
        ! ! select hyperslab in the file
        ! call h5dget_space_f(dset_id,filespace,error)
        ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_ppsi, count_ppsi, error)
        ! ! create property list for collective dataset write
        ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        ! ! write the dataset collectively, the fourth argument is the size of the data in the file
        ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ppsi_out, dim_ppsi, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp = plist_id)
        ! call h5pclose_f(plist_id, error)
        ! call h5dclose_f(dset_id, error)
        ! call h5sclose_f(memspace, error)
        ! call h5sclose_f(filespace, error)
        ! >>

        IF (okvan) THEN
           ! << output ppsi_us --> ppsi_us(1:npwx*npol,1:nvnc,1:nproc) --> ppsi_us_out(2,1:npwx_max,1:npol_spin,1:nvnc, 1:npolarization, 1:nproc)

           call hdf5_write_double_hyperslab(file_id, '/ppsi/ppsi_us', (/2, npwx_max, npol, nvnc, 3, 1, 1/), (/0, 0, 0, 0, 0, ik-1, mpime/), ppsi_us_out, error)

           ! count_ppsi_us(1) = 2
           ! count_ppsi_us(2) = npwx
           ! count_ppsi_us(3) = npol
           ! count_ppsi_us(4) = nvnc
           ! count_ppsi_us(5) = 3
           ! count_ppsi_us(6) = 1
           ! count_ppsi_us(7) = 1

           ! offset_ppsi_us(1) = 0
           ! offset_ppsi_us(2) = 0
           ! offset_ppsi_us(3) = 0
           ! offset_ppsi_us(4) = 0
           ! offset_ppsi_us(5) = 0
           ! offset_ppsi_us(6) = ik-1
           ! offset_ppsi_us(7) = mpime ! starting from 0

           ! dim_ppsi_us(1) = 2
           ! dim_ppsi_us(2) = npwx_max
           ! dim_ppsi_us(3) = npol
           ! dim_ppsi_us(4) = nvnc
           ! dim_ppsi_us(5) = 3
           ! dim_ppsi_us(6) = nkstot
           ! dim_ppsi_us(7) = nproc

           ! call h5screate_simple_f(rank_ppsi_us, count_ppsi_us, memspace, error)

           ! call h5dopen_f(file_id, '/ppsi/ppsi_us', dset_id, error)
           ! ! select hyperslab in the file
           ! call h5dget_space_f(dset_id,filespace,error)
           ! call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offset_ppsi_us, count_ppsi_us, error)
           ! ! create property list for collective dataset write
           ! call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
           ! call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
           ! ! write the dataset collectively, the fourth argument is the size of the data in the file
           ! call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ppsi_us_out, dim_ppsi_us, error, file_space_id=filespace, mem_space_id=memspace, xfer_prp = plist_id)
           ! call h5pclose_f(plist_id, error)
           ! call h5dclose_f(dset_id, error)
           ! call h5sclose_f(memspace, error)
           ! call h5sclose_f(filespace, error)
        ENDIF

        call h5fclose_f(file_id, error)

     else ! <- output_ppsi = F -> we will output matp

        CALL mp_sum(matp, intra_bgrp_comm)

        IF (ionode) THEN
           ! IF (ik == nks1) &
           !      WRITE (iunout, '(" &p_mat nbnd=",i4,", nks=",i4," /")') &
           !      nbnd, nks2-nks1+1
           ! WRITE (iunout, '(10x,3f10.6,i7)') xk(1,ik),xk(2,ik),xk(3,ik), &
           !      nbnd_occ

           ! DO ipol=1,3
           !    WRITE (iunout, '(i3)') ipol
           !    DO ibnd=nbnd_occ+1,nbnd
           !       WRITE (iunout, '(5f15.8)') &
           !            (abs(matp(ibnd-nbnd_occ,jbnd,ipol))**2, jbnd=1,nbnd_occ)

           !    ENDDO
           ! ENDDO
           ! ======
           ! write to files
           do ipol = 1, 3
              ! < The kpoint index in matp_k always refers to the conduction state!
              ! <- ppsi_v = T .and. output_ppsi = F -> here ik refers to <c k|
              if (ppsi_v) then
                 matp_k(1:nvb,1:ncb,ik,ipol) = matp(1:nvb,1:ncb,ipol)
              else ! <- ppsi_v = F .and. restart_with_ppsi = T -> here ik_ket refers to | c k'>
                 matp_k(1:nvb,1:ncb,ik_ket,ipol) = matp(1:nvb,1:ncb,ipol)
              endif
           enddo ! ipol
           ! write out vmt matrix elements
           DO ibnd = nbnd_occ+1, nbnd_occ+ncb
              DO jbnd = nbnd_occ-nvb+1, nbnd_occ
                 ! Relative band indices
                 ivb = nbnd_occ - jbnd + 1
                 icb = ibnd-nbnd_occ
                 !< formatted output
                 if (verbo) then                    
                    if (ppsi_v) then
                       write(*,'(4I8, " X (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, ibnd, ik_ket, jbnd, matp(ivb,icb,1), abs(matp(ivb,icb,1))**2
                       write(*,'(4I8, " Y (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, ibnd, ik_ket, jbnd, matp(ivb,icb,2), abs(matp(ivb,icb,2))**2
                       write(*,'(4I8, " Z (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, ibnd, ik_ket, jbnd, matp(ivb,icb,3), abs(matp(ivb,icb,3))**2
                    else
                       write(*,'(4I8, " X (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, jbnd, ik_ket, ibnd, matp(ivb,icb,1), abs(matp(ivb,icb,1))**2
                       write(*,'(4I8, " Y (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, jbnd, ik_ket, ibnd, matp(ivb,icb,2), abs(matp(ivb,icb,2))**2
                       write(*,'(4I8, " Z (", ES20.12, ",", ES20.12, ")", ES20.12)') ik, jbnd, ik_ket, ibnd, matp(ivb,icb,3), abs(matp(ivb,icb,3))**2
                    endif
                 endif
              enddo
           enddo
        ENDIF
     endif
     CALL mp_barrier ( world_comm )
     if (restart_with_ppsi) then
        deallocate(map_ppsi2wfng_l)
        deallocate(igwf_l2g)
     endif
  ENDDO ! ik

  if (.not. restart_with_ppsi) then
     DEALLOCATE(ppsi)

     if (output_ppsi) then
        DEALLOCATE(ppsi_out)
     endif

     IF (okvan) THEN
        DEALLOCATE(ppsi_us)
        if (output_ppsi) then
           DEALLOCATE(ppsi_us_out)
        endif
     ENDIF

  else

     DEALLOCATE(ppsi)
     DEALLOCATE(ppsi_)
     deallocate(ppsi_in)

     if (okvan) then
        DEALLOCATE(ppsi_us)
        deallocate(ppsi_us_in)
     endif

  endif

  ! if (output_ppsi .or. restart_with_ppsi) then
  !    close(unit_ppsi_vc_l)
  !    close(unit_ppsi_vc_l_formatted)
  !    close(unit_igk_l)
  !    close(unit_igk_l_formatted)
  ! endif

  if (restart_with_ppsi) then
     deallocate(xk_ket)
     deallocate(igk_l2g_ket)
     deallocate(ngk_ket)
     deallocate(et_ket)
     DEALLOCATE(igwk)
     DEALLOCATE(g_useful)
     deallocate(map_ik2ikket)
  endif

  if ((.not. output_ppsi)) then
     IF (ionode) THEN
        ! ======
        ! unformatted
        ! ------
        ! The kpoint index in matp_k always refers to the conduction state
        write(11) matp_k(1:nvb, 1:ncb, 1:nkstot, 1)
        write(22) matp_k(1:nvb, 1:ncb, 1:nkstot, 2)
        write(33) matp_k(1:nvb, 1:ncb, 1:nkstot, 3)

        deallocate ( matp_k )

        close(11)
        close(22)
        close(33)

        CLOSE(iunout)
     ENDIF

     DEALLOCATE(matp)
  endif

  call h5close_f(error)

CONTAINS
  subroutine hdf5_write_int(loc_id, dset_name, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(in) :: buf !< data buffer
    integer, intent(out) :: errcode !< 1 error code
    integer(HSIZE_T), dimension(1) :: dims_h5type(1)
    logical :: exists_
    integer(HID_T) :: dset_id

    dims_h5type(1) = 0
    !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (exists_) then
       call h5dopen_f(loc_id, dset_name, dset_id, errcode)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, (/buf/), dims_h5type, errcode)
       call h5dclose_f(dset_id, errcode)
    else
       call h5ltmake_dataset_int_f(loc_id, dset_name, 0, dims_h5type, (/buf/), errcode)
    endif

  end subroutine hdf5_write_int

  !/* The following routine reads or writes an array, of type THE_TYPE, from the HDF5 file */
  !> Writes a(n) integer array.
  subroutine hdf5_write_int_array(loc_id, dset_name, dims, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(in), dimension(:) :: dims !< size of the buffer buf
    integer, intent(in), dimension(*) :: buf !< data buffer
    integer, intent(out) :: errcode !< error code
    integer :: rank
    integer(HSIZE_T), dimension(size(dims)) :: dims_h5type
    logical :: exists_
    integer(HID_T) :: dset_id

    rank = size(dims)
    dims_h5type = dims
    !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (exists_) then
       call h5dopen_f(loc_id, dset_name, dset_id, errcode)
       call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, dims_h5type, errcode)
       call h5dclose_f(dset_id, errcode)
    else
       call h5ltmake_dataset_int_f(loc_id, dset_name, rank, dims_h5type, buf, errcode)
    endif

  end subroutine hdf5_write_int_array

  !/* The following routine reads or writes a scalar value, of type THE_TYPE, from the HDF5 file */
  !> Writes a rank-0 real(DP) scalar.
  subroutine hdf5_write_double(loc_id, dset_name, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    real(DP), intent(in) :: buf !< data buffer
    integer, intent(out) :: errcode !< 1 error code
    integer(HSIZE_T), dimension(1) :: dims_h5type(1)
    logical :: exists_
    integer(HID_T) :: dset_id

    dims_h5type(1) = 0
    !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (exists_) then
       call h5dopen_f(loc_id, dset_name, dset_id, errcode)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, (/buf/), dims_h5type, errcode)
       call h5dclose_f(dset_id, errcode)
    else
       call h5ltmake_dataset_double_f(loc_id, dset_name, 0, dims_h5type, (/buf/), errcode)
    endif

  end subroutine hdf5_write_double

  !/* The following routine reads or writes an array, of type THE_TYPE, from the HDF5 file */
  !> Writes a(n) real(DP) array.
  subroutine hdf5_write_double_array(loc_id, dset_name, dims, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(in), dimension(:) :: dims !< size of the buffer buf
    real(DP), intent(in), dimension(*) :: buf !< data buffer
    integer, intent(out) :: errcode !< error code
    integer :: rank
    integer(HSIZE_T), dimension(size(dims)) :: dims_h5type
    logical :: exists_
    integer(HID_T) :: dset_id

    rank = size(dims)
    dims_h5type = dims
    !FHJ: We can`t use H5LTmake_dataset_* if the dataset already exists!
    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (exists_) then
       call h5dopen_f(loc_id, dset_name, dset_id, errcode)
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, dims_h5type, errcode)
       call h5dclose_f(dset_id, errcode)
    else
       call h5ltmake_dataset_double_f(loc_id, dset_name, rank, dims_h5type, buf, errcode)
    endif

  end subroutine hdf5_write_double_array

  !> Creates an empty dataset
  subroutine hdf5_create_dset(loc_id, dset_name, dtype, dims, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer(HID_T), intent(in) :: dtype
    integer, intent(in) :: dims(:)
    integer, intent(out) :: errcode
    integer(HSIZE_T) :: hdims(size(dims))
    integer(HID_T) :: dset_id
    integer(HID_T) :: dspace

    hdims(:) = dims(:)
    call h5screate_simple_f(size(dims), hdims, dspace, errcode)
    call h5dcreate_f(loc_id, dset_name, dtype, dspace, dset_id, errcode)
    call h5dclose_f(dset_id, errcode)
    call h5sclose_f(dspace, errcode)

  end subroutine hdf5_create_dset

  !/* The following routine reads or writes a scalar value, of type THE_TYPE, from the HDF5 file */
  !> Reads a rank-0 integer scalar.
  subroutine hdf5_read_int(loc_id, dset_name, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(inout) :: buf !< data buffer
    integer, intent(out) :: errcode !< 1 error code
    integer(HSIZE_T), dimension(1) :: dims_h5type(1)
    integer, dimension(1) :: buf_1(1)

    dims_h5type(1) = 0
    call h5ltread_dataset_int_f(loc_id, dset_name, buf_1, dims_h5type, errcode)
    buf = buf_1(1)

  end subroutine hdf5_read_int

  !/* The following routine reads or writes an array, of type THE_TYPE, from the HDF5 file */
  !> Reads a(n) integer array.
  subroutine hdf5_read_int_array(loc_id, dset_name, dims, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(in), dimension(:) :: dims !< size of the buffer buf
    integer, intent(inout), dimension(*) :: buf !< data buffer
    integer, intent(out) :: errcode !< error code
    integer :: rank
    integer(HSIZE_T), dimension(size(dims)) :: dims_h5type

    rank = size(dims)
    dims_h5type = dims
    call h5ltread_dataset_int_f(loc_id, dset_name, buf, dims_h5type, errcode)

  end subroutine hdf5_read_int_array

  subroutine hdf5_read_int_hyperslab(loc_id, dset_name, countf, offsetf, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    !> Number of elements to read from the dataset for each dimention
    integer, intent(in) :: countf(:)
    !> Offset when reading dataset from file.
    integer, intent(in) :: offsetf(:)
    !> Data buffer. We treat it as a flat contiguous 1D array.
    integer, intent(inout), dimension(*) :: buf
    integer, intent(out) :: errcode !< error code
    integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
    integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
    integer(HSIZE_T) :: hoffsetf(size(offsetf)) !< Offset for file dataspace
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    logical :: exists_

    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (.not. exists_) call exit(321)
    call h5dopen_f(loc_id, dset_name, dset_id, errcode)
    ! FHJ: Get 2D file dataspace and set selection mask
    call h5dget_space_f(dset_id, dataspace, errcode)
    hcountf(:) = countf(:)
    hoffsetf(:) = offsetf(:)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, hoffsetf, hcountf, errcode)
    ! FHJ: Create flat memory dataspace
    hcountm(1) = product(countf)
    call h5screate_simple_f(1, hcountm, memspace, errcode)
    ! FHJ: Read dataspace
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, buf, hcountm, errcode, memspace, dataspace)
    call h5sclose_f(memspace, errcode)
    call h5sclose_f(dataspace, errcode)
    call h5dclose_f(dset_id, errcode)

  end subroutine hdf5_read_int_hyperslab

  subroutine hdf5_write_int_hyperslab(loc_id, dset_name, countf, offsetf, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    !> Number of elements to read from the dataset for each dimention
    integer, intent(in) :: countf(:)
    !> Offset when reading dataset from file.
    integer, intent(in) :: offsetf(:)
    !> Data buffer. We treat it as a flat contiguous 1D array.
    integer, intent(in), dimension(*) :: buf
    integer, intent(out) :: errcode !< error code
    integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
    integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
    integer(HSIZE_T) :: hoffsetf(size(offsetf)) !< Offset for file dataspace
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    logical :: exists_

    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (.not. exists_) call exit(321)
    call h5dopen_f(loc_id, dset_name, dset_id, errcode)
    ! FHJ: Get 2D file dataspace and set selection mask
    call h5dget_space_f(dset_id, dataspace, errcode)
    hcountf(:) = countf(:)
    hoffsetf(:) = offsetf(:)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, hoffsetf, hcountf, errcode)
    ! FHJ: Create flat memory dataspace
    hcountm(1) = product(countf)
    call h5screate_simple_f(1, hcountm, memspace, errcode)
    ! FHJ: Read dataspace
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buf, hcountm, errcode, memspace, dataspace)
    call h5sclose_f(memspace, errcode)
    call h5sclose_f(dataspace, errcode)
    call h5dclose_f(dset_id, errcode)

  end subroutine hdf5_write_int_hyperslab

  ! ===============

  !/* The following routine reads or writes a scalar value, of type THE_TYPE, from the HDF5 file */
  !> Reads a rank-0 real(DP) scalar.
  subroutine hdf5_read_double(loc_id, dset_name, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    real(DP), intent(inout) :: buf !< data buffer
    integer, intent(out) :: errcode !< 1 error code
    integer(HSIZE_T), dimension(1) :: dims_h5type(1)
    real(DP), dimension(1) :: buf_1(1)

    dims_h5type(1) = 0
    call h5ltread_dataset_double_f(loc_id, dset_name, buf_1, dims_h5type, errcode)
    buf = buf_1(1)

  end subroutine hdf5_read_double

  !/* The following routine reads or writes an array, of type THE_TYPE, from the HDF5 file */
  !> Reads a(n) real(DP) array.
  subroutine hdf5_read_double_array(loc_id, dset_name, dims, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    integer, intent(in), dimension(:) :: dims !< size of the buffer buf
    real(DP), intent(inout), dimension(*) :: buf !< data buffer
    integer, intent(out) :: errcode !< error code
    integer :: rank
    integer(HSIZE_T), dimension(size(dims)) :: dims_h5type

    rank = size(dims)
    dims_h5type = dims
    call h5ltread_dataset_double_f(loc_id, dset_name, buf, dims_h5type, errcode)

  end subroutine hdf5_read_double_array

  subroutine hdf5_read_double_hyperslab(loc_id, dset_name, countf, offsetf, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    !> Number of elements to read from the dataset for each dimention
    integer, intent(in) :: countf(:)
    !> Offset when reading dataset from file.
    integer, intent(in) :: offsetf(:)
    !> Data buffer. We treat it as a flat contiguous 1D array.
    real(DP), intent(inout), dimension(*) :: buf
    integer, intent(out) :: errcode !< error code
    integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
    integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
    integer(HSIZE_T) :: hoffsetf(size(offsetf)) !< Offset for file dataspace
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    logical :: exists_

    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (.not. exists_) call exit(321)
    call h5dopen_f(loc_id, dset_name, dset_id, errcode)
    ! FHJ: Get 2D file dataspace and set selection mask
    call h5dget_space_f(dset_id, dataspace, errcode)
    hcountf(:) = countf(:)
    hoffsetf(:) = offsetf(:)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, hoffsetf, hcountf, errcode)
    ! FHJ: Create flat memory dataspace
    hcountm(1) = product(countf)
    call h5screate_simple_f(1, hcountm, memspace, errcode)
    ! FHJ: Read dataspace
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, buf, hcountm, errcode, memspace, dataspace)
    call h5sclose_f(memspace, errcode)
    call h5sclose_f(dataspace, errcode)
    call h5dclose_f(dset_id, errcode)

  end subroutine hdf5_read_double_hyperslab

  subroutine hdf5_write_double_hyperslab(loc_id, dset_name, countf, offsetf, buf, errcode)
    integer(HID_T), intent(in) :: loc_id !< 1 file id
    character(LEN=*), intent(in) :: dset_name !< 1 dataset name
    !> Number of elements to read from the dataset for each dimention
    integer, intent(in) :: countf(:)
    !> Offset when reading dataset from file.
    integer, intent(in) :: offsetf(:)
    !> Data buffer. We treat it as a flat contiguous 1D array.
    real(DP), intent(in), dimension(*) :: buf
    integer, intent(out) :: errcode !< error code
    integer(HSIZE_T) :: hcountf(size(countf)) !< Count for file dataspace
    integer(HSIZE_T) :: hcountm(1) !< Count for memory dataspace
    integer(HSIZE_T) :: hoffsetf(size(offsetf)) !< Offset for file dataspace
    integer(HID_T) :: dset_id
    integer(HID_T) :: dataspace
    integer(HID_T) :: memspace
    logical :: exists_

    call h5lexists_f(loc_id, dset_name, exists_, errcode)
    if (.not. exists_) call exit(321)
    call h5dopen_f(loc_id, dset_name, dset_id, errcode)
    ! FHJ: Get 2D file dataspace and set selection mask
    call h5dget_space_f(dset_id, dataspace, errcode)
    hcountf(:) = countf(:)
    hoffsetf(:) = offsetf(:)
    call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, hoffsetf, hcountf, errcode)
    ! FHJ: Create flat memory dataspace
    hcountm(1) = product(countf)
    call h5screate_simple_f(1, hcountm, memspace, errcode)
    ! FHJ: Read dataspace
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, buf, hcountm, errcode, memspace, dataspace)
    call h5sclose_f(memspace, errcode)
    call h5sclose_f(dataspace, errcode)
    call h5dclose_f(dset_id, errcode)

  end subroutine hdf5_write_double_hyperslab

  !RETURN
END SUBROUTINE write_p_avg
