!
! Copyright (C) 2008-2012 Georgy Samsonidze
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Converts the output files produced by pw.x to the input files for BerkeleyGW.
!
!-------------------------------------------------------------------------------
!
! BerkeleyGW, Copyright (c) 2011, The Regents of the University of
! California, through Lawrence Berkeley National Laboratory (subject to
! receipt of any required approvals from the U.S. Dept. of Energy).
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! (1) Redistributions of source code must retain the above copyright
! notice, this list of conditions and the following disclaimer.
!
! (2) Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! (3) Neither the name of the University of California, Lawrence
! Berkeley National Laboratory, U.S. Dept. of Energy nor the names of
! its contributors may be used to endorse or promote products derived
! from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! You are under no obligation whatsoever to provide any bug fixes,
! patches, or upgrades to the features, functionality or performance of
! the source code ("Enhancements") to anyone; however, if you choose to
! make your Enhancements available either publicly, or directly to
! Lawrence Berkeley National Laboratory, without imposing a separate
! written license agreement for such Enhancements, then you hereby grant
! the following license: a  non-exclusive, royalty-free perpetual
! license to install, use, modify, prepare derivative works, incorporate
! into other computer software, distribute, and sublicense such
! enhancements or derivative works thereof, in binary and source code
! form.
!
!-------------------------------------------------------------------------------
!
! pw2bgw subroutines:
!
! write_wfng  - generates complex wavefunctions in G-space (normalized to 1)
! real_wfng   - constructs real wavefunctions by applying the Gram-Schmidt
!               process (called from write_wfng)
! write_rhog  - generates real/complex charge density in G-space
!               (units of the number of electronic states per unit cell)
! calc_rhog   - computes charge density by summing over a subset of occupied
!               bands (called from write_rhog), destroys charge density
! write_vxcg  - generates real/complex exchange-correlation potential in
!               G-space (units of Rydberg) [only local part of Vxc]
! write_vxc0  - prints real/complex exchange-correlation potential at G=0
!               (units of eV) [only local part of Vxc]
! write_vxc_r - calculates matrix elements of exchange-correlation potential
!               in R-space (units of eV) [only local part of Vxc]
! write_vxc_g - calculates matrix elements of exchange-correlation potential
!               in G-space (units of eV) [supports non-local Vxc]
! write_vscg  - generates real/complex self-consistent potential in G-space
!               (units of Rydberg) [only local part of Vsc]
! write_vkbg  - generates complex Kleinman-Bylander projectors in G-space
!               (units of Rydberg)
! check_inversion - checks whether real/complex version is appropriate
!               (called from everywhere)
!
! Quantum ESPRESSO stores the wavefunctions in is-ik-ib-ig order
! BerkeleyGW stores the wavefunctions in ik-ib-is-ig order
! the outer loop is over is(QE)/ik(BGW) and the inner loop is over ig
! ik = k-point index, is = spin index, ib = band index, ig = G-vector index
!
! write_wfng reverts the order of is and ik using smap and kmap arrays,
! distributes wavefunctions over processors by ig (either real case or
! spin-polarized case), calls real_wfng that applies the Gram-Schmidt
! process (real case), reverts the order of is and ib (spin-polarized
! case), and writes wavefunctions to disk
!
!-------------------------------------------------------------------------------

!> [important] we only use one kpool, which means nks = nkstot
PROGRAM pw2bgw

  USE constants, ONLY : eps12
  USE control_flags, ONLY : gamma_only
  USE environment, ONLY : environment_start, environment_end
  USE io_files, ONLY : prefix, tmp_dir
  USE io_global, ONLY : ionode, ionode_id
  USE kinds, ONLY : DP
  USE lsda_mod, ONLY : nspin
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE paw_variables, ONLY : okpaw
  USE scf, ONLY : rho_core, rhog_core, rho, v
  USE uspp, ONLY : okvan
  USE ldaU, ONLY : lda_plus_u
  USE exx, ONLY : use_ace, local_thr, nbndproj, ecutfock, aceinit, exxinit
  USE exx_base, ONLY : x_gamma_extrapolation, nq1, nq2, nq3, exxdiv_treatment, yukawa, ecutvcut, exx_grid_init, exx_mp_init, exx_div_check
  USE funct, ONLY: set_exx_fraction, set_screening_parameter, dft_is_hybrid, exx_is_active, stop_exx, dft_is_meta
  USE gvecw, ONLY : ecutwfc
  USE ions_base, ONLY : nat, atm, ityp, tau
  USE symm_base, ONLY : s, ftau, nsym, ft, time_reversal, find_sym
  USE start_k, ONLY : nk1, nk2, nk3
  USE loc_scdm, ONLY : use_scdm, localize_orbitals
  USE loc_scdm_k, ONLY : localize_orbitals_k
  USE loc_scdm,      ONLY : use_scdm, scdm_den, scdm_grd, n_scdm
  USE input_parameters, ONLY : scdmden, scdmgrd, nscdm, n_proj, localization_thr, scdm, ace, nqx1, nqx2, nqx3
  USE fft_rho,  ONLY : rho_g2r  
  USE fft_base,  ONLY : dfftp

  IMPLICIT NONE

  character(len=6) :: codename = 'PW2BGW'

  integer :: real_or_complex
  character ( len = 9 ) :: symm_type
  logical :: wfng_flag
  character ( len = 256 ) :: wfng_file
  logical :: wfng_kgrid
  integer :: wfng_nk1
  integer :: wfng_nk2
  integer :: wfng_nk3
  real (DP) :: wfng_dk1
  real (DP) :: wfng_dk2
  real (DP) :: wfng_dk3
  logical :: wfng_occupation
  integer :: wfng_nvmin
  integer :: wfng_nvmax
  logical :: rhog_flag
  character ( len = 256 ) :: rhog_file
  integer :: rhog_nvmin
  integer :: rhog_nvmax
  logical :: vxcg_flag
  character ( len = 256 ) :: vxcg_file
  logical :: vxc0_flag
  character ( len = 256 ) :: vxc0_file
  character :: vxc_integral
  logical :: vxc_flag
  character ( len = 256 ) :: vxc_file
  integer :: vxc_diag_nmin, vxc_diag_nmax, vxc_offdiag_nmin, vxc_offdiag_nmax
  logical :: vxc_tot_flag
  character ( len = 256 ) :: vxc_tot_file
  integer :: vxc_tot_diag_nmin, vxc_tot_diag_nmax, vxc_tot_offdiag_nmin, vxc_tot_offdiag_nmax
  logical :: vxc_zero_rho_core
  logical :: vscg_flag
  character ( len = 256 ) :: vscg_file
  logical :: vkbg_flag
  character ( len = 256 ) :: vkbg_file
  character ( len = 256 ) :: outdir

  logical :: use_hdf5, use_binary, flag_output_asc_wfn, flag_output_spin, flag_output_mirror
  !> output V_hub matrix elements
  logical :: vhub_flag
  character ( len = 256 ) :: vhub_file
  integer :: vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax
  ! real(DP) :: localization_thr=0.0D0, screening_parameter=-1.0D0, exx_fraction=-1.0D0 !, ecutfock_=-1.0D0
  real(DP) :: screening_parameter=-1.0D0, exx_fraction=-1.0D0 !, ecutfock_=-1.0D0
  ! integer :: n_proj = 0
  logical :: DoLoc = .false.
  logical :: use_ace_ = .true.

  NAMELIST / input_pw2bgw / prefix, outdir, &
       real_or_complex, symm_type, wfng_flag, wfng_file, wfng_kgrid, &
       wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3, &
       wfng_occupation, wfng_nvmin, wfng_nvmax, rhog_flag, rhog_file, &
       rhog_nvmin, rhog_nvmax, vxcg_flag, vxcg_file, vxc0_flag, vxc0_file, &
       vxc_flag, vxc_file, vxc_integral, vxc_diag_nmin, vxc_diag_nmax, &
       vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core, &
       vscg_flag, vscg_file, vkbg_flag, vkbg_file, &
       use_hdf5, use_binary, flag_output_asc_wfn, flag_output_spin, flag_output_mirror, &
       vhub_flag, vhub_file, vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax, &
       use_ace_, vxc_tot_flag, vxc_tot_file, vxc_tot_diag_nmin, vxc_tot_diag_nmax,  vxc_tot_offdiag_nmin, vxc_tot_offdiag_nmax

  integer :: ii, ios
  character ( len = 256 ) :: output_file_name
  character (len=256), external :: trimcheck
  character (len=1), external :: lowercase
  REAL(DP), allocatable :: m_loc(:, :)
  REAL(DP) :: etxc, vtxc  !FZ: change062320
  ! REAL(DP) :: fock3 = 0.0D0

#if defined(__MPI)
  CALL mp_startup ( )
#endif

  CALL environment_start ( codename )

  prefix = 'prefix'
  CALL get_environment_variable ( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM ( outdir ) == ' ' ) outdir = './'
  real_or_complex = 2
  symm_type = 'cubic'
  wfng_flag = .FALSE.
  wfng_file = 'WFN'
  wfng_kgrid = .FALSE.
  wfng_nk1 = 0
  wfng_nk2 = 0
  wfng_nk3 = 0
  wfng_dk1 = 0.0D0
  wfng_dk2 = 0.0D0
  wfng_dk3 = 0.0D0
  wfng_occupation = .FALSE.
  wfng_nvmin = 0
  wfng_nvmax = 0
  rhog_flag = .FALSE.
  rhog_file = 'RHO'
  rhog_nvmin = 0
  rhog_nvmax = 0
  vxcg_flag = .FALSE.
  vxcg_file = 'VXC'
  vxc0_flag = .FALSE.
  vxc0_file = 'vxc0.dat'
  vxc_flag = .FALSE.
  vxc_file = 'vxc.dat'
  vxc_integral = 'g'
  vxc_diag_nmin = 0
  vxc_diag_nmax = 0
  vxc_offdiag_nmin = 0
  vxc_offdiag_nmax = 0
  vxc_tot_file = 'vxc_tot.dat'
  vxc_tot_diag_nmin = 0
  vxc_tot_diag_nmax = 0
  vxc_tot_offdiag_nmin = 0
  vxc_tot_offdiag_nmax = 0

  vxc_zero_rho_core = .TRUE.
  vscg_flag = .FALSE.
  vscg_file = 'VSC'
  vkbg_flag = .FALSE.
  vkbg_file = 'VKB'
  ! (default) output binary format
  use_hdf5 = .false.
  use_binary = .true.
  flag_output_asc_wfn = .false.
  flag_output_spin = .false.
  flag_output_mirror = .false.
  vhub_flag = .false.
  vhub_file = 'vhub.dat'
  vhub_diag_nmin = 0
  vhub_diag_nmax = 0
  vhub_offdiag_nmin = 0
  vhub_offdiag_nmax = 0

  ! nq1 = nqx1
  ! nq2 = nqx2
  ! nq3 = nqx3
  ! exxdiv_treatment_ = trim(exxdiv_treatment)
  ! yukawa_   = yukawa
  ! ecutvcut_ = ecutvcut
  ! use_ace   = ace
  ! nbndproj  = n_proj
  ! local_thr = localization_thr

  use_ace_ = .true.
  use_scdm  = scdm
  scdm_den = scdmden
  scdm_grd = scdmgrd
  n_scdm   = nscdm
  ! write(*,*) "scdmden = ", scdmden, " scdmgrd = ", scdmgrd

  ! IF ( ionode ) THEN
  !    IF ( local_thr > 0.0_dp .AND. .NOT. use_ace ) then
  !       CALL errore('input','localization without ACE not implemented',1)
  !    ENDIF
  ! ENDIF
  ! IF (exx_fraction >= 0.0_DP) CALL set_exx_fraction (exx_fraction)
  ! IF (screening_parameter >= 0.0_DP) &
  !      & CALL set_screening_parameter (screening_parameter)

  ! write(*,*) "x_gamma_extrapolation = ", x_gamma_extrapolation, " nq1 = ", nq1, " nq2 = ", nq2, " nq3 = ", nq3
  ! write(*,*) "exxdiv_treatment = ", exxdiv_treatment, " yukawa = ", yukawa, " ecutvcut = ", ecutvcut, " use_ace = ", use_ace, " nbndproj = ", nbndproj, " local_thr = ", local_thr
  ! write(*,*) "exx_fraction = ", exx_fraction, " ecutfock = ", ecutfock
  ! write(*,*) "screening_parameter = ", screening_parameter

  IF ( ionode ) THEN
     ! flib/inpfile.f90:
     ! open access to input file : pw2bgw.inp
     ! ------
     ! in pw2bgw.inp, we only specify vxc_diag_nmin and vxc_diag_nmax
     ! so vxc_offdiag_nmin = 0, vxc_offdiag_nmax = 0
     ! for G0W0 calculation, the above setting is fine, but for SCGW calculation
     ! we need to include offdiag terms (in single-particle state basis) for Vxc
     CALL input_from_file ( )
     READ ( 5, input_pw2bgw, iostat = ios )

     if (use_hdf5) then
        write(*,'(5X,A)') "We will output in HDF5 format!"
     endif

     if (use_binary) then
        write(*,'(5X,A)') "We will output in binary format!"
     endif

     IF ( ios /= 0 ) CALL errore ( codename, 'input_pw2bgw', abs ( ios ) )

     DO ii = 1, LEN_TRIM (symm_type)
        symm_type(ii:ii) = lowercase (symm_type(ii:ii))
     END DO
     DO ii = 1, LEN_TRIM (vxc_integral)
        vxc_integral(ii:ii) = lowercase (vxc_integral(ii:ii))
     END DO

     IF ( real_or_complex /= 1 .AND. real_or_complex /= 2 ) &
          CALL errore ( codename, 'real_or_complex', 1 )
     IF ( symm_type /= 'cubic' .AND. symm_type /= 'hexagonal' ) &
          CALL errore ( codename, 'symm_type', 1 )
     IF ( vxc_integral /= 'r' .AND. vxc_integral /= 'g' ) &
          CALL errore ( codename, 'vxc_integral', 1 )
  ENDIF

  tmp_dir = trimcheck ( outdir )
  CALL mp_bcast ( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast ( prefix, ionode_id, world_comm )
  CALL mp_bcast ( real_or_complex, ionode_id, world_comm )
  CALL mp_bcast ( symm_type, ionode_id, world_comm )
  CALL mp_bcast ( wfng_flag, ionode_id, world_comm )
  CALL mp_bcast ( wfng_file, ionode_id, world_comm )
  CALL mp_bcast ( wfng_kgrid, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk1, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk2, ionode_id, world_comm )
  CALL mp_bcast ( wfng_dk3, ionode_id, world_comm )
  CALL mp_bcast ( wfng_occupation, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( wfng_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( rhog_flag, ionode_id, world_comm )
  CALL mp_bcast ( rhog_file, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmin, ionode_id, world_comm )
  CALL mp_bcast ( rhog_nvmax, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxcg_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc0_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc_integral, ionode_id, world_comm )
  CALL mp_bcast ( vxc_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_offdiag_nmax, ionode_id, world_comm )

  CALL mp_bcast ( vxc_tot_flag, ionode_id, world_comm )
  CALL mp_bcast ( vxc_tot_file, ionode_id, world_comm )
  CALL mp_bcast ( vxc_tot_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_tot_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vxc_tot_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vxc_tot_offdiag_nmax, ionode_id, world_comm )

  CALL mp_bcast ( vxc_zero_rho_core, ionode_id, world_comm )
  CALL mp_bcast ( vscg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vscg_file, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_flag, ionode_id, world_comm )
  CALL mp_bcast ( vkbg_file, ionode_id, world_comm )
  CALL mp_bcast ( use_binary, ionode_id, world_comm )
  CALL mp_bcast ( use_hdf5, ionode_id, world_comm )
  CALL mp_bcast ( flag_output_asc_wfn, ionode_id, world_comm )
  CALL mp_bcast ( flag_output_spin, ionode_id, world_comm )
  CALL mp_bcast ( flag_output_mirror, ionode_id, world_comm )
  CALL mp_bcast ( vhub_flag, ionode_id, world_comm )
  CALL mp_bcast ( vhub_file, ionode_id, world_comm )
  CALL mp_bcast ( vhub_diag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vhub_diag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( vhub_offdiag_nmin, ionode_id, world_comm )
  CALL mp_bcast ( vhub_offdiag_nmax, ionode_id, world_comm )
  CALL mp_bcast ( use_ace_, ionode_id, world_comm )

  ! ------
  ! PW/src/read_file.f90:
  ! [IMPORTANT]
  ! read PP, allocate potential and WFN, read RHO, initialized FFT grid dfftp and dffts
  ! pw_readfile('rho',ierr)
  ! ------
  ! [Read RHO]:
  ! PW/src/pw_restart.f90: pw_readfile(...)
  ! PW/src/io_rho_xml.f90: read_rho(rho, nspin) ==> read_rho_general(rho, nspin, extension) ==> read_rho_only(rho%of_r, nspin, extension), total, mx, my, mz
  ! ------
  ! [Read WFN]:
  ! PW/src/read_file.f90: read_xml_file_internal(withbs=.true.)
  ! ======
  ! Note that in read_rho_general(rho, nspin, extension), rho is of type TYPE(scf_type)
  ! while in read_rho_only(rho, nspin, extension), rho is of type REAL(DP)
  ! ------
  ! Modules/xml_io_base.f90: read_rho_xml(...)
  ! read_rho_xml():
  ! ======
  ! [IMPORTANT]
  ! We initialize nl(:) in read_file()/ggen(...)
  ! Modules/recvec_subs.f90:
  ! ------
  ! SUBROUTINE ggen ( gamma_only, at, bg, comm, no_global_sort )
  ! ------
  ! ngm : local  number of G vectors (on this processor)
  ! nl(ig_local_g(:)) = ig_global_FFT
  ! ======
  ! [IMPORTANT]
  ! we already have rho%of_g after read_file():
  ! L319@read_file.f90:
  ! DO is = 1, nspin
  !    psic(:) = rho%of_r(:,is)
  !    CALL fwfft ('Dense', psic, dfftp)
  !    rho%of_g(:,is) = psic(nl(:))
  ! END DO
  ! ======
  ! [IMPORTANT]
  ! allocate psic(1:dfftp%nnr)
  ! ------
  ! PW/src/allocate_fft.f90:
  ! ------
  ! ALLOCATE (psic( dfftp%nnr))

  ! ======
  ! pw_readfile('exx')
  ! ------
  ! PW/src/pw_restart.f90: read_exx(ierr) ! exact exchange ==> LDA don't have this term
  !
  CALL read_file ( )

  use_ace = use_ace_
  IF ( dft_is_hybrid() ) THEN
     IF ( local_thr > 0.0_dp .AND. .NOT. use_ace ) then
        CALL errore('input','localization without ACE not implemented',1)
     ENDIF
     IF ( use_scdm ) CALL errore('input','use_scdm not yet implemented',1)
     IF (exx_fraction >= 0.0_DP) CALL set_exx_fraction (exx_fraction)
     IF (screening_parameter >= 0.0_DP) &
          & CALL set_screening_parameter (screening_parameter)
     DoLoc = local_thr.gt.0.0d0

     ! write(*,*) "x_gamma_extrapolation = ", x_gamma_extrapolation, " nq1 = ", nq1, " nq2 = ", nq2, " nq3 = ", nq3
     ! write(*,*) "exxdiv_treatment = ", exxdiv_treatment, " yukawa = ", yukawa, " ecutvcut = ", ecutvcut, " use_ace = ", use_ace, " nbndproj = ", nbndproj, " local_thr = ", local_thr
     ! write(*,*) "exx_fraction = ", exx_fraction, " ecutfock = ", ecutfock
     ! write(*,*) "screening_parameter = ", screening_parameter

     !> [important]
     !> exx_is_active() = T after read_file()
     call stop_exx()
  ENDIF
  !> Force stop exx ==> exx_started = F, such that exxalfa can be initialized in exxinit()

  !! ======
  !! PW/src/read_file.f90:
  !! IF ( lsda ) THEN
  !!    nspin = 2
  !!    npol  = 1
  !! ELSE IF ( noncolin ) THEN
  !!    nspin        = 4
  !!    npol         = 2
  !!    current_spin = 1
  !! ELSE
  !!    nspin        = 1
  !!    npol         = 1
  !!    current_spin = 1
  !! ENDIF
  if (ionode) then
     if (MAX (MAXVAL (ABS (rho_core (:) ) ), MAXVAL (ABS (rhog_core (:) ) ) ) &
          .LT. eps12) then
        WRITE ( 6, '(/,5x,"NLCC is absent")' )
     else
        WRITE ( 6, '(/,5x,"NLCC is present")' )
     endif
  endif
  if (okvan) call errore ( 'pw2bgw', 'BGW cannot use USPP.', 3 )
  if (okpaw) call errore ( 'pw2bgw', 'BGW cannot use PAW.', 4 )
  if (gamma_only) call errore ( 'pw2bgw', 'BGW cannot use gamma-only run.', 5 )
  if (real_or_complex == 1 .AND. vxc_flag .AND. vxc_offdiag_nmax > 0) &
       call errore ( 'pw2bgw', 'Off-diagonal matrix elements of Vxc ' // &
       'with real wavefunctions are not implemented, compute them in ' // &
       'Sigma using VXC.', 7)

  ! PW/src/openfil.f90:
  ! open access to prefix.hub* (ready to write)
  CALL openfil()

  !! PW/src/hinit0.f90:
  !! hinit0() contains gk_sort()
  !! gk_sort() sorts k + G in order of increasing magnitude, up to ecut. And write igk(:) into prefix.igk* files.
  !! PP/src/openfil_pp.f90:
  !! openfil_pp() will open access to ./prefix.wfc* files for reading
  !! set nwordwfc = nbnd * npwx * npol here
  !! twfcollect=.true. ! ==> delete prefix.wfc* files after usage
  CALL hinit0()

  !> [important]
  IF (dft_is_meta()) then                           !FZ:  for metaGGA change062320
     CALL rho_g2r ( dfftp, rho%kin_g, rho%kin_r )   !FZ:  for metaGGA
     CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
  ENDIF

  CALL openfil_pp ( )

  !! Atomic configuration dependent hamiltonian initialization
  !! CALL hinit1()
  IF ( lda_plus_u ) THEN
     IF (ionode) THEN
        write(*,*) "lda_plus_u = T, call orthoUwfc()"
     ENDIF
     CALL orthoUwfc()
  ENDIF

  allocate(m_loc(3,nat))
  m_loc = 0.0D0
  ! write(*,*) "time_reversal = ", time_reversal
  CALL find_sym ( nat, tau, ityp, .not.time_reversal, m_loc)
  deallocate(m_loc)
  !> Initialize hybrid functional calculations

  ! IF (ecutfock <= 0.0_DP) THEN
  !    ! default case
  !    ecutfock_ = 4.0_DP*ecutwfc
  ! ELSE
  !    IF(ecutfock < ecutwfc .OR. ecutfock > ecutrho) CALL errore('iosys', &
  !         'ecutfock can not be < ecutwfc or > ecutrho!', 1)
  !    ecutfock_ = ecutfock
  ! ENDIF

  IF ( dft_is_hybrid() ) THEN
     CALL exx_grid_init()
     CALL exx_mp_init()
     CALL exx_div_check()

     write(*,*) "DoLoc = ", DoLoc
     CALL exxinit(DoLoc)
     IF( DoLoc.and.gamma_only) THEN
        CALL localize_orbitals( )
     ELSE IF (DoLoc) THEN
        CALL localize_orbitals_k( )
     END IF
     IF (use_ace) THEN
        ! CALL aceinit ( DoLoc, fock3)
        CALL aceinit ( DoLoc)
        ! write(*,*) "DoLoc = ", DoLoc, " fock3 = ", fock3
     ENDIF
  ENDIF

  ! write(*,*) "x_gamma_extrapolation = ", x_gamma_extrapolation, " nq1 = ", nq1, " nq2 = ", nq2, " nq3 = ", nq3
  ! write(*,*) "exxdiv_treatment = ", exxdiv_treatment, " yukawa = ", yukawa, " ecutvcut = ", ecutvcut, " use_ace = ", use_ace, " nbndproj = ", nbndproj, " local_thr = ", local_thr
  ! write(*,*) "exx_fraction = ", exx_fraction, " ecutfock = ", ecutfock
  ! write(*,*) "screening_parameter = ", screening_parameter

  if ( ionode ) WRITE ( 6, '("")' )
  ! For Norm-conserving PP calculations, we always use ecutrho = 4 * ecutwfc, which leads to identify dfftp and dffts (doublegrid = F)
  ! So in later FFT functions, dfftp and dffts mean the same thing.
  IF ( wfng_flag ) THEN
     output_file_name = TRIM ( tmp_dir ) // TRIM ( wfng_file )
     IF ( ionode ) WRITE ( 6, '(5x,"call write_wfng")' )
     CALL start_clock ( 'write_wfng' )
     CALL write_wfng ( output_file_name, real_or_complex, symm_type, &
          wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
          wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax, flag_output_asc_wfn, flag_output_spin)
     CALL stop_clock ( 'write_wfng' )
     IF ( ionode ) WRITE ( 6, '(5x,"done write_wfng",/)' )
  ENDIF

  ! IF ( vxcg_flag ) THEN
  !    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxcg_file )
  !    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxcg")' )
  !    CALL start_clock ( 'write_vxcg' )
  !    CALL write_vxcg ( output_file_name, real_or_complex, symm_type, &
  !         vxc_zero_rho_core )
  !    CALL stop_clock ( 'write_vxcg' )
  !    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxcg",/)' )
  ! ENDIF

  ! IF ( vxc0_flag ) THEN
  !    output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc0_file )
  !    IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc0")' )
  !    CALL start_clock ( 'write_vxc0' )
  !    CALL write_vxc0 ( output_file_name, vxc_zero_rho_core )
  !    CALL stop_clock ( 'write_vxc0' )
  !    IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc0",/)' )
  ! ENDIF

  IF ( vxc_flag ) THEN
     output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc_file )
     IF ( vxc_integral .EQ. 'r' ) THEN
        IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_r")' )
        CALL start_clock ( 'write_vxc_r' )
        !CALL write_vxc_r ( output_file_name, &
        !     vxc_diag_nmin, vxc_diag_nmax, &
        !     vxc_offdiag_nmin, vxc_offdiag_nmax, &
        !     vxc_zero_rho_core )
        CALL stop_clock ( 'write_vxc_r' )
        IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_r",/)' )
     ENDIF
     IF ( vxc_integral .EQ. 'g' ) THEN
        IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_g")' )
        CALL start_clock ( 'write_vxc_g' )
        CALL write_vxc_g ( output_file_name, vxc_diag_nmin, vxc_diag_nmax, vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core )
        CALL stop_clock ( 'write_vxc_g' )
        IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_g",/)' )

        ! IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_tot_g")' )
        ! CALL start_clock ( 'write_vxc_tot_g' )
        ! CALL write_vxc_tot_g ( output_file_name, vxc_diag_nmin, vxc_diag_nmax, vxc_offdiag_nmin, vxc_offdiag_nmax, vxc_zero_rho_core )
        ! CALL stop_clock ( 'write_vxc_tot_g' )
        ! IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_tot_g",/)' )
     ENDIF
  ENDIF

  IF ( vxc_tot_flag ) THEN
     output_file_name = TRIM ( tmp_dir ) // TRIM ( vxc_tot_file )
     IF ( ionode ) WRITE ( 6, '(5x,"call write_vxc_tot_g")' )
     CALL start_clock ( 'write_vxc_tot_g' )
     CALL write_vxc_tot_g ( output_file_name, vxc_tot_diag_nmin, vxc_tot_diag_nmax, vxc_tot_offdiag_nmin, vxc_tot_offdiag_nmax, vxc_zero_rho_core )
     CALL stop_clock ( 'write_vxc_tot_g' )
     IF ( ionode ) WRITE ( 6, '(5x,"done write_vxc_tot_g",/)' )
  ENDIF

  ! IF ( vscg_flag ) THEN
  !    output_file_name = TRIM ( tmp_dir ) // TRIM ( vscg_file )
  !    IF ( ionode ) WRITE ( 6, '(5x,"call write_vscg")' )
  !    CALL start_clock ( 'write_vscg' )
  !    CALL write_vscg ( output_file_name, real_or_complex, symm_type )
  !    CALL stop_clock ( 'write_vscg' )
  !    IF ( ionode ) WRITE ( 6, '(5x,"done write_vscg",/)' )
  ! ENDIF

  ! IF ( vkbg_flag ) THEN
  !    output_file_name = TRIM ( tmp_dir ) // TRIM ( vkbg_file )
  !    IF ( ionode ) WRITE ( 6, '(5x,"call write_vkbg")' )
  !    CALL start_clock ( 'write_vkbg' )
  !    CALL write_vkbg ( output_file_name, symm_type, wfng_kgrid, wfng_nk1, &
  !         wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3 )
  !    CALL stop_clock ( 'write_vkbg' )
  !    IF ( ionode ) WRITE ( 6, '(5x,"done write_vkbg",/)' )
  ! ENDIF

  ! since calc_rhog (called from write_rhog) destroys charge density,
  ! it must be called after v_xc (called from write_vxcg, write_vxc0,
  ! write_vxc_r, write_vxc_g)
  IF ( rhog_flag ) THEN
     output_file_name = TRIM ( tmp_dir ) // TRIM ( rhog_file )
     IF ( ionode ) WRITE ( 6, '(5x,"call write_rhog")' )
     CALL start_clock ( 'write_rhog' )
     CALL write_rhog ( output_file_name, real_or_complex, symm_type, &
          rhog_nvmin, rhog_nvmax )
     CALL stop_clock ( 'write_rhog' )
     IF ( ionode ) WRITE ( 6, '(5x,"done write_rhog",/)' )
  ENDIF

  IF (lda_plus_u .and. vhub_flag) THEN
     output_file_name = TRIM ( outdir ) // '/' // TRIM ( vhub_file )
     IF ( ionode ) WRITE ( 6, '(5x,"<BIN> call write_vhub")' )
     CALL start_clock ( 'write_vhub' )
     CALL write_vhub_g ( output_file_name, vhub_diag_nmin, vhub_diag_nmax, vhub_offdiag_nmin, vhub_offdiag_nmax)
     CALL stop_clock ( 'write_vhub' )
     IF ( ionode ) WRITE ( 6, '(5x,"<BIN> done write_vhub",/)' )
  ENDIF

  IF ( ionode ) WRITE ( 6, * )
  IF ( wfng_flag ) CALL print_clock ( 'write_wfng' )
  IF ( rhog_flag ) CALL print_clock ( 'write_rhog' )
  IF ( vxcg_flag ) CALL print_clock ( 'write_vxcg' )
  IF ( vxc0_flag ) CALL print_clock ( 'write_vxc0' )
  IF ( vxc_flag ) THEN
     IF ( vxc_integral .EQ. 'r' ) CALL print_clock ( 'write_vxc_r' )
     IF ( vxc_integral .EQ. 'g' ) CALL print_clock ( 'write_vxc_g' )
  ENDIF
  IF ( vscg_flag ) CALL print_clock ( 'write_vscg' )
  IF ( vkbg_flag ) CALL print_clock ( 'write_vkbg' )
  IF ( wfng_flag .AND. real_or_complex .EQ. 1 ) THEN
     IF ( ionode ) WRITE ( 6, '(/,5x,"Called by write_wfng:")' )
     CALL print_clock ( 'real_wfng' )
  ENDIF

  CALL environment_end ( codename )

  CALL stop_pp ( )

  ! this is needed because openfil is called above
  CALL close_files ( .false. )

  STOP

CONTAINS

  !-------------------------------------------------------------------------------
  !> Refer to Modules/io_base.f90
  SUBROUTINE write_wfng ( output_file_name, real_or_complex, symm_type, &
       wfng_kgrid, wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, &
       wfng_dk3, wfng_occupation, wfng_nvmin, wfng_nvmax, flag_output_asc_wfn, flag_output_spin)

    USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
    USE constants, ONLY : pi, tpi, eps6
    USE fft_base, ONLY : dfftp
    USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
    USE io_files, ONLY : iunwfc, nwordwfc
    USE io_global, ONLY : ionode, ionode_id
    USE ions_base, ONLY : nat, atm, ityp, tau
    USE kinds, ONLY : DP
    USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
    USE lsda_mod, ONLY : nspin, isk
    USE mp, ONLY : mp_sum, mp_max, mp_get, mp_bcast, mp_barrier
    USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
         intra_pool_comm, inter_pool_comm
    USE mp_wave, ONLY : mergewf
    USE mp_world, ONLY : mpime, nproc, world_comm
    USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
    USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
    ! --------------------
    ! s (1:3,1:3,1:nsym): rotation matrix, applying on 1*3 vectors from right
    ! r'= r \cdot S, r=(x,y,z)
    ! --------------------
    ! we will take transpose of s into s_transpose
    ! --------------------
    ! (integer) ftau(1:3,1:nsym): fractional translation vector, in FFT coordinates
    ! for example, for diamond, FFT: nr1=nr2=nr3=64
    ! ftau(1:3,...) = (-16,-16,-16)
    ! \tau_1 \cdot S - ftau = \tau_2, \tau_i's are vectors for basis
    ! --------------------
    ! (real) ft(1:3,1:nsym): fractional translation vector, in crystal coordinates
    ! for diamond, ft(1:3,...) = (-0.25,-0.25,-0.25) = ftau./nr
    ! --------------------
    ! we will take (-1)*ft(:,:) into ft_minus
    USE symm_base, ONLY : s, ftau, nsym, ft
    USE wavefunctions, ONLY : evc
    USE wvfct, ONLY : npwx, nbnd, et, wg
    USE gvecw, ONLY : ecutwfc
    USE matrix_inversion
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_DOUBLE_COMPLEX
#endif
    USE noncollin_module , ONLY : noncolin , npol

    IMPLICIT NONE

    character ( len = 256 ), intent (in) :: output_file_name
    integer, intent (in) :: real_or_complex
    character ( len = 9 ), intent (in) :: symm_type
    logical, intent (in) :: wfng_kgrid
    integer, intent (in) :: wfng_nk1
    integer, intent (in) :: wfng_nk2
    integer, intent (in) :: wfng_nk3
    real (DP), intent (in) :: wfng_dk1
    real (DP), intent (in) :: wfng_dk2
    real (DP), intent (in) :: wfng_dk3
    logical, intent (in) :: wfng_occupation
    integer, intent (in) :: wfng_nvmin
    integer, intent (in) :: wfng_nvmax
    logical, intent (in) :: flag_output_asc_wfn, flag_output_spin

    character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
    logical :: proc_wf, bad_kgrid
    integer :: unit, i, j, k, cell_symmetry, nrecord
    integer :: id, ib, ik, iks, ike, is, ig, ierr
    integer :: nd, ntran, nb, nk_l, nk_g, ns, ng_l, ng_g
    integer :: npw, ngg, npw_g, npwx_g
    integer :: local_pw, ipsour, igwx, ngkdist_g, ngkdist_l
    real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
    real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
    real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
    integer, allocatable :: kmap ( : )
    integer, allocatable :: smap ( : )
    integer, allocatable :: ifmin ( : )
    integer, allocatable :: ifmax ( : )
    integer, allocatable :: itmp ( : )
    integer, allocatable :: ngk_g ( : )
    integer, allocatable :: ipmask ( : )
    integer, allocatable :: igwk ( : )
    integer, allocatable :: igwf_l2g ( : )
    integer, allocatable :: g_g ( :, : )
    integer, allocatable :: igk_l2g ( :, : )
    real (DP), allocatable :: et_g ( :, : )
    real (DP), allocatable :: wg_g ( :, : )
    real (DP), allocatable :: energy ( :, : )
    complex (DP), allocatable, TARGET :: wfng ( : )
    complex (DP), POINTER :: wfng2 ( : )
    complex (DP), allocatable :: wfng_buf ( :, : )
    complex (DP), allocatable :: wfng_dist ( :, :, : )

    INTEGER, EXTERNAL :: atomic_number, global_kpoint_index
    integer, dimension(1:3,1:3,1:nsym) :: s_transpose
    real(DP), dimension(1:3,1:nsym) :: ft_minus
    integer, allocatable :: g_useful( :, : )
    integer, parameter :: unit_spin_x=125, unit_spin_y=126, unit_spin_z=127
    complex (DP) :: sigma_x, sigma_y, sigma_z
    complex (DP), parameter :: imagunit = DCMPLX(0,1)

    integer :: unit_wfn, ascunit=23, ik_global, itran
    character(LEN=30) :: file_wfn
    character(LEN=20) :: ik_global_string
    character(LEN=20) :: ib_string, is_string
    complex (DP) :: norm
    real(DP) :: norm_tol=1.0D-8

    ! In noncollinear case, nspin = 4, npol = 2
    ! [ IMPORTANT ] nspin = 2 ==> proc_wf = .TRUE.

    ! write(*,*) "nspin = ", nspin, " npol = ", npol

    IF ( real_or_complex .EQ. 1 .OR. nspin .EQ. 2 ) THEN
       proc_wf = .TRUE.
    ELSE
       proc_wf = .FALSE.
    ENDIF

    IF (ionode .and. (flag_output_spin)) THEN
       OPEN ( unit = unit_spin_x, file = 'spin.x.dat', form = 'formatted', status = 'replace' )
       WRITE (unit_spin_x, "(A)") '# kpt     band    sigma_x'
       OPEN ( unit = unit_spin_y, file = 'spin.y.dat', form = 'formatted', status = 'replace' )
       WRITE (unit_spin_y, "(A)") '# kpt     band    sigma_y'
       OPEN ( unit = unit_spin_z, file = 'spin.z.dat', form = 'formatted', status = 'replace' )
       WRITE (unit_spin_z, "(A)") '# kpt     band    sigma_z'
    ENDIF

    bad_kgrid = .FALSE.
    IF ( wfng_kgrid ) THEN
       IF ( wfng_nk1 .LE. 0 .OR. wfng_nk2 .LE. 0 .OR. wfng_nk3 .LE. 0 ) &
            bad_kgrid = .TRUE.
    ELSE
       if (ionode) then
          write(*,'(5X,A)') "[ERROR] must set wfng_kgrid"
       endif

       call exit(123)

       IF ( nk1 .LE. 0 .OR. nk2 .LE. 0 .OR. nk3 .LE. 0 ) &
            bad_kgrid = .TRUE.
    ENDIF

    IF ( bad_kgrid .AND. ionode ) THEN
       WRITE ( 6, 101 )
    ENDIF

    CALL date_and_tim ( cdate, ctime )
    WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
    WRITE ( stime, '(A8,24X)' ) ctime(1:8)
    IF ( real_or_complex .EQ. 1 ) THEN
       WRITE ( stitle, '("WFN-Real",24X)' )
    ELSE
       WRITE ( stitle, '("WFN-Complex",21X)' )
    ENDIF

    unit = 4

    nrecord = 1
    nd = 3

    nb = nbnd
    nk_l = nks
    nk_g = nkstot
    ns = nspin
    ng_l = ngm
    ng_g = ngm_g

    iks = global_kpoint_index (nkstot, 1)
    ike = iks + nks - 1

    ! write(*,'(A,I5,A,I5,A,I5,A,I5)') "mpime = ", mpime, " nks = ", nks, " iks = ", iks, " ike = ", ike

    ALLOCATE ( kmap ( nk_g ) )
    ALLOCATE ( smap ( nk_g ) )

    DO i = 1, nk_g
       IF (noncolin) THEN
          kmap ( i ) = i
          smap ( i ) = 1
       ELSE
          j = ( i - 1 ) / ns
          k = i - 1 - j * ns
          kmap ( i ) = j + k * ( nk_g / ns ) + 1
          smap ( i ) = k + 1
       ENDIF
    ENDDO
    ierr = 0
    DO i = 1, nk_g
       !! ik = ik_global
       ik = kmap ( i )
       is = smap ( i )
       IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
            ierr = ierr + 1
    ENDDO
    CALL mp_max ( ierr, world_comm )
    IF ( ierr .GT. 0 ) &
         CALL errore ( 'write_wfng', 'smap', ierr )
    ! -------------------------------
    !              [a1]
    ! ------adot = [a2] \cdot [a1 a2 a3] , in units of Bohr
    !              [a3]
    !> adot(i,j) = a_i \cdot a_j
    ! ai = (a1x, a1y, a1z), i = 1, 2, 3
    ! at = (a1,a2,a3)
    ! omega : volume of unit cell
    alat2 = alat ** 2
    recvol = 8.0D0 * pi**3 / omega

    DO i = 1, nd
       DO j = 1, nd
          adot ( j, i ) = 0.0D0
       ENDDO
    ENDDO
    DO i = 1, nd
       DO j = 1, nd
          DO k = 1, nd
             adot ( j, i ) = adot ( j, i ) + &
                  at ( k, j ) * at ( k, i ) * alat2
          ENDDO
       ENDDO
    ENDDO
    ! -------------------------------
    !               [b1]
    ! ------bdot =  [b2] \cdot [b1 b2 b3], in units of Bohr^-1
    !               [b3]
    !> bdot(i,j) = b_i \cdot b_j
    ! bi = (b1x, b1y, b1z), i = 1, 2, 3
    ! bg = (b1, b2, b3)
    DO i = 1, nd
       DO j = 1, nd
          bdot ( j, i ) = 0.0D0
       ENDDO
    ENDDO
    DO i = 1, nd
       DO j = 1, nd
          DO k = 1, nd
             ! tpiba = 2 \pi / alat
             ! tpiba2 = tpiba^2
             bdot ( j, i ) = bdot ( j, i ) + &
                  bg ( k, j ) * bg ( k, i ) * tpiba2
          ENDDO
       ENDDO
    ENDDO

    ierr = 0
    IF ( ibrav .EQ. 0 ) THEN
       IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
          cell_symmetry = 0
       ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
          cell_symmetry = 1
       ELSE
          ierr = 1
       ENDIF
    ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
       cell_symmetry = 0
    ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
       cell_symmetry = 1
    ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
       cell_symmetry = 0
    ELSE
       ierr = 1
    ENDIF
    IF ( ierr .GT. 0 ) &
         CALL errore ( 'write_wfng', 'cell_symmetry', ierr )

    ntran = nsym
    DO i = 1, ntran
       ! ft_minus(1:3, i) = -ft(1:3, i)
       do id = 1, 3
          !> If an entry of ft is 0.5, ft_minus = ft = 0.5
          if (ABS(ft(id, i) - 0.5D0) < 1.0D-12) then
             ft_minus(id, i) = 0.5D0
          else
             ft_minus(id, i) = - ft(id, i)
          endif
       enddo
       s_transpose(1:3, 1:3, i) = transpose(s(1:3, 1:3, i))
    ENDDO

    ! DO i = 1, ntran
    !    DO j = 1, nd
    !       DO k = 1, nd
    !          r1 ( k, j ) = dble ( s ( k, j, i ) )
    !       ENDDO
    !    ENDDO
    !    CALL invmat ( 3, r1, r2 )
    !    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    !    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    !    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    !    DO j = 1, nd
    !       t2 ( j ) = 0.0D0
    !       DO k = 1, nd
    !          t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
    !       ENDDO
    !       IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
    !            t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
    !       IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
    !            t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    !    ENDDO
    !    DO j = 1, nd
    !       translation ( j, i ) = t2 ( j ) * tpi
    !    ENDDO
    ! ENDDO

    ! CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

    ALLOCATE ( et_g ( nb, nk_g ) )

    DO ik = 1, nk_l
       DO ib = 1, nb
          et_g ( ib, ik ) = et ( ib, ik )
       ENDDO
    ENDDO
#if defined(__MPI)
    CALL poolrecover ( et_g, nb, nk_g, nk_l )
    CALL mp_bcast ( et_g, ionode_id, world_comm )
#endif

    ALLOCATE ( wg_g ( nb, nk_g ) )
    ALLOCATE ( ifmin ( nk_g ) )
    ALLOCATE ( ifmax ( nk_g ) )

    IF ( wfng_occupation ) THEN

       DO ik = 1, nk_g
          DO ib = 1, nb
             IF ( ib .GE. wfng_nvmin .AND. ib .LE. wfng_nvmax ) THEN
                wg_g ( ib, ik ) = 1.0D0
             ELSE
                wg_g ( ib, ik ) = 0.0D0
             ENDIF
          ENDDO
       ENDDO
       DO ik = 1, nk_g
          ifmin ( ik ) = wfng_nvmin
       ENDDO
       DO ik = 1, nk_g
          ifmax ( ik ) = wfng_nvmax
       ENDDO

    ELSE

       DO ik = 1, nk_l
          DO ib = 1, nb
             IF ( wk(ik) == 0.D0 ) THEN
                wg_g(ib,ik) = wg(ib,ik)
             ELSE
                wg_g(ib,ik) = wg(ib,ik) / wk(ik)
             ENDIF
          ENDDO
       ENDDO
#if defined(__MPI)
       CALL poolrecover ( wg_g, nb, nk_g, nk_l )
#endif
       DO ik = 1, nk_g
          ifmin ( ik ) = 0
       ENDDO
       DO ik = 1, nk_g
          ifmax ( ik ) = 0
       ENDDO
       DO ik = 1, nk_g
          DO ib = 1, nb
             IF ( wg_g( ib, ik ) .GT. 0.5D0 ) THEN
                IF ( ifmin ( ik ) .EQ. 0 ) ifmin ( ik ) = ib
                ifmax ( ik ) = ib
             ENDIF
          ENDDO
       ENDDO

    ENDIF

    !! -----------------------------------
    !! g_g(1:3, 1:ngm_g) gives the fractional coordinates for
    !! all (1:ngm_g) the G vectors, global variable, have a duplicate
    !! at each cpu
    ALLOCATE ( g_g ( nd, ng_g ) )

    !! ngm_g : global number of Gvectors (summed over procs)
    !! in serial execution, ngm_g = ngm
    DO ig = 1, ng_g
       DO id = 1, nd
          g_g ( id, ig ) = 0
       ENDDO
    ENDDO

    !! ngm : local number of Gvectors
    !! mill(1:3, ig[1:ngm]) = miller index of G vectors (local to each pool)
    !! G(:) = mill(1)*bg(:,1)+mill(2)*bg(:,2)+mill(3)*bg(:,3)
    !! where bg are the reciprocal lattice basis vectors

    !! Note that this is ig_l2g(:) not igk_l2g(:)
    !! ig_l2g(ig) maps the index of local GVector list (G_i, i=1,ngm) into
    !! global Gvector list g_g(1:3,1:ngm_g)

    !! use local miller indices to recover the global array of Gvectors

    !! g_g (1:3, 1:ngm) is global Gvector list
    !! if global_sort = .true., g_g(1:3,1:ngm) is sorted around
    !! origin point

    !! loop over local Gvectors
    DO ig = 1, ng_l
       g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
       g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
       g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
    ENDDO

    ! ======
    !CALL mp_sum ( g_g, intra_bgrp_comm )
    CALL mp_sum ( g_g, intra_pool_comm )

    !! npwx : maximum number of Gvectors needed by a (ib,ik) wavefunction
    !!        default value npwx = 0
    !!        do nk=1,nks
    !!          npwx = max (npwx, ngk (nk) )
    !!        enddo
    !!
    !! There are nks kpoints shared by this pool
    !! ngk(nk) gives the number of Gvectors used by nk-th kpoint
    !! in this pool
    !!
    !! evc(1:npwx*npol, 1:nbnd) is local wavefunction of ONLY ONE kpoint
    !! and ALL the bands
    !! npol = 1 for scalar-relativistic case
    !! npol = 2 for full-relativistic case
    !!
    !! ALLOCATE( evc( npwx*npol, nbnd ) )
    !!
    !! nks : number of kpoints in this pool
    !! nwordwfc = nbnd*npwx*npol
    !! Note that for a given k point and band index,
    !! evc(:,ib) contains ngk(ik) (k+G) components, and
    !! some void numbers at tail : evc(ngk(ik)+1:npwx,ib) = 0.0

    ALLOCATE ( igk_l2g ( npwx, nk_l ) )

    !! g_g(1:3, igk_l2g(ig, ik)) gives the Gvector for the ig-th Gvector
    !!    of the local ik-th kpoint,
    !! KS(ik) = [\sum_{n=1}^{ik-1} ngk(n)]+1
    !! KE(ik) = [\sum_{n=1}^{ik} ngk(n)]
    !! evc(KS(ik):KE(ik),ibnd) gives all the PW components for (ik,ibnd) indices
    !! g_g(1:3, igk_l2g(KS(ik):KE(ik), ik)) gives all the Gvectors for the local
    !! ik-th kpoint

    DO ik = 1, nk_l
       npw = ngk ( ik )
       DO ig = 1, npw
          igk_l2g ( ig, ik ) = ig_l2g ( igk_k (ig, ik) )
       ENDDO
       DO ig = npw + 1, npwx
          igk_l2g ( ig, ik ) = 0
       ENDDO
    ENDDO

    !! ---------- Set ngk_g(1:nkstot) --------
    !! compute the global number of GVectors for each k point
    !! ngk(1:nks) contains the number of Gvectors for kpoints in a pool
    !! ngk_g(1:nkstot) contains the number of Gvectors for all the kpoints
    !! ngk_g(iks,ike) contains the number of Gvectors for the kpoints in
    !! current pool
    ALLOCATE ( ngk_g ( nk_g ) )

    !! ngk_g = 0
    !! ngk_g(iks:ike) = ngk(1:nks)
    !! CALL mp_sum( ngk_g, inter_pool_comm )
    !! CALL mp_sum( ngk_g, intra_pool_comm )
    !! ngk_g = ngk_g / nbgrp
    DO ik = 1, nk_g
       ngk_g ( ik ) = 0
    ENDDO
    DO ik = 1, nk_l
       ngk_g ( ik + iks - 1 ) = ngk ( ik )
    ENDDO
    !! ... ngk_g(ik_global) = sum of the number of gvectors distributed
    !!                    in every cpu of the same pool
    !! ... Note that in a pool, every cpu has different gvectors, there is
    !!     no common gvector shared by two cpus!!!
    !CALL mp_sum( ngk_g, inter_pool_comm )
    !CALL mp_sum( ngk_g, intra_pool_comm )
    !ngk_g = ngk_g / nbgrp
    CALL mp_sum ( ngk_g, world_comm )

    !! Compute the maximum global Gvector INDEX among all the kpoints
    !! INDEX not number
    npw_g = MAXVAL ( igk_l2g ( :, : ) )

    ! DO NOT USE intra_pool_comm or intra_bgrp_comm here
    ! CALL mp_max( npw_g, intra_pool_comm )
    CALL mp_max ( npw_g, world_comm )

    ! Compute the maximum number of Gvector owned by one kpoints among all kpoints
    ! in general, npwx_g < npw_g
    npwx_g = MAXVAL ( ngk_g ( : ) )

    !! ... For example, suppose there are 2 kpoints in total, kpoint#1
    !!     owns 5 Gvectors(1,3,5,7,9 as in g_g list) , kpoint#2 owns
    !!     4 Gvectors(1,5,10,15 as in g_g list), in this way, npwx_g = 5
    !!     due to kpoint#1,
    !! ... Moreover, kpoint#1's 5-th Gvectors has the index in global
    !!     list g_g as 9, but kpoint#2's 4-th Gvectors has the index in global
    !!     list g_g as 15, in this way, npw_g = 15 due to kpoint#2
    !!
    !! ... also notice that we will need to consider (1,3,5,7,9,10,15) gvectors
    !!     as all the useful gvectors.
    !!
    !! ... It is obvious that npw_g serves as the UPPER bound of INDEX
    !!     for all the useful Gvectors in global g_g list
    !!
    !! ... While npwx_g serves as the LOWER bound of the NUMBER of all the
    !!     useful Gvectors
    !! -----------------------------------------
    !! k vectors xk(1:3, 1:nks=nkstot) from cartesian to crystal

    !! s^{T}_reci = s^{T}_cart \cdot A^{T}
    !! CALL cryst_to_cart ( nk_g / ns, xk, at, - 1 )
    CALL cryst_to_cart ( nkstot, xk, at, - 1 )

    IF ( ionode ) THEN
       OPEN ( unit = unit, file = TRIM ( output_file_name ), &
            form = 'unformatted', status = 'replace' )
       WRITE ( unit ) stitle, sdate, stime
       !! WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho, &
       !!     nk_g / ns, nb, npwx_g, ecutwfc
       if (nspin .eq. 2) then
          WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
               nkstot / nspin, nbnd, npwx_g, ecutwfc
       else
          WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
               nkstot, nbnd, npwx_g, ecutwfc
       endif

       IF ( wfng_kgrid ) THEN
          WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
               wfng_dk1, wfng_dk2, wfng_dk3
       ELSE
          ! WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
          !      0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
          WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
               dble (xk(1,1)), dble (xk(2,1)), dble (xk(3,1))
       ENDIF

       WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
            ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
       WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
            ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
       WRITE ( unit ) ( ( (s_transpose(k, j, i), k = 1, nd ), j = 1, nd ), i = 1, ntran )
       !! WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
       WRITE ( unit ) ( ( ft_minus(j, i), j = 1, nd), i = 1, ntran)
       !! WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
       WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )

       IF (nspin .eq. 2) THEN
          WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot/nspin)
          WRITE ( unit ) ( wk ( ik), ik = 1, nkstot/nspin)
          WRITE ( unit ) ( ( xk ( id, ik), id = 1, nd ), ik = 1, nkstot/nspin)
       ELSE
          WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot)
          IF (noncolin) THEN
             WRITE ( unit ) ( wk ( ik ), ik = 1, nkstot )
          ELSE
             WRITE ( unit ) ( wk ( ik )/2.0D0, ik = 1, nkstot )
          ENDIF
          WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot)
       ENDIF
       !! WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nk_g / ns )
       !! WRITE ( unit ) ( wk ( ik ) * dble ( ns ) / 2.0D0, ik = 1, nk_g / ns )
       !! WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nk_g / ns )
       WRITE ( unit ) ( ifmin ( ik ), ik = 1, nk_g )
       WRITE ( unit ) ( ifmax ( ik ), ik = 1, nk_g )
       WRITE ( unit ) ( ( et_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
       WRITE ( unit ) ( ( wg_g ( ib, ik ), ib = 1, nb ), ik = 1, nk_g )
       WRITE ( unit ) nrecord
       WRITE ( unit ) ng_g
       WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
    ENDIF

    ! output in formatted form, ASCII form
    IF ( ionode ) THEN
       OPEN ( unit = ascunit, file = 'WFN.bin.asc', form = 'formatted', status = 'replace' )
       WRITE ( ascunit,*) ' stitle = ', stitle, ' sdate = ', sdate, ' stime = ', stime
       WRITE ( ascunit,*) ' nspin = ', nspin, ' ngm_g = ', ngm_g, ' ntran = ', ntran, ' cell_symmetry = ', cell_symmetry, ' nat = ', nat, ' ecutrho = ', ecutrho, &
            ' nkstot = ', nkstot, ' nbnd = ', nbnd, ' npwx_g = ', npwx_g, 'npw_g = ', npw_g, ' ecutwfc = ', ecutwfc
       write (ascunit,*) 'FFT grid:'
       IF ( wfng_kgrid ) THEN
          WRITE ( ascunit,*) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
               wfng_dk1, wfng_dk2, wfng_dk3
       ELSE
          WRITE ( ascunit,*) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
               dble (xk(1,1)), dble (xk(2,1)), dble (xk(3,1))
       ENDIF
       WRITE ( ascunit,*) ' omega = ', omega, ' alat = ', alat
       ! -------------------------------------
       write(ascunit,*) 'Lattice vectors = [a1, a2, a3] :'
       do i=1, 3
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') at(i,:)
       enddo
       write(ascunit,*) '-----------------------'
       ! -------------------------------------
       write(ascunit,*) 'adota : '
       do i=1, 3
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') adot(i,:)
       enddo
       write(ascunit,*) '-----------------------'
       ! -------------------------------------
       WRITE ( ascunit,*) ' recvol = ', recvol, 'tpiba = ', tpiba
       ! -------------------------------------
       write(ascunit,*) 'Reciprocal lattice vectors = [b1, b2, b3] :'
       do i=1, 3
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') bg(i,:)
       enddo
       write(ascunit,*) '-----------------------'
       ! -------------------------------------
       write(ascunit,*) 'bdotb : '
       do i = 1, 3
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8,1X)') bdot(i,:)
       enddo
       write(ascunit,*) '-----------------------'
       ! -------------------------------------
       WRITE ( ascunit,*) ' Orignal symmetry matrices (s) : total ', ntran, 'symmetries'
       do itran = 1, ntran
          write(ascunit,'("#" I8)') itran
          do i = 1, 3
             WRITE (ascunit,'(I8,1X,I8,1X,I8)') s(i, 1:3, itran)
          enddo
          write (ascunit,*) '----------------------'
       enddo
       ! -------------------------------------
       WRITE ( ascunit,*) ' Transpose symmetry matrices (s_transpose) : total ', ntran, 'symmetries'
       do itran = 1, ntran
          write(ascunit,'("#" I8)') itran
          do i = 1, 3
             WRITE (ascunit,'(I8,1X,I8,1X,I8)') s_transpose(i, 1:3, itran)
          enddo
          write (ascunit,*) '----------------------'
       enddo
       ! WRITE ( ascunit,*) ' Fractional translation : '
       ! do itran=1,ntran
       !   write(ascunit,'("#" I8)') itran
       !   write(ascunit,'(F15.8,1X,F15.8,1X,F15.8)') translation(1:3,itran)
       ! enddo
       ! write (ascunit,*) '----------------------'
       write(ascunit,*) ' Fractional translation in FFT coordinates (ftau)'
       do itran = 1, ntran
          write(ascunit,'("#" I8)') itran
          write(ascunit,'(I8,1X,I8,1X,I8)') ftau(1:3, itran)
       enddo
       write(ascunit,*) '-------------------------'
       write(ascunit,*) ' Fractional translation in crystal coordinates (ft)'
       do itran = 1, ntran
          write(ascunit,'("#" I8)') itran
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8)') ft(1:3, itran)
       enddo
       write(ascunit,*) '-------------------------'
       write(ascunit,*) ' Fractional translation in crystal coordinates (ft_minus)'
       do itran = 1, ntran
          write(ascunit,'("#" I8)') itran
          write(ascunit,'(F15.8,1X,F15.8,1X,F15.8)') ft_minus(1:3, itran)
       enddo
       write(ascunit,*) '-------------------------'
       !> tau is cartesian coordinates of atoms, in units of alat*Bohr
       WRITE ( ascunit,*) ' Atomic basis : ', ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm( ityp ( i ) ) ), i = 1, nat )
       WRITE ( ascunit,*) ' ngk_g : ', ( ngk_g ( ik_global ), ik_global = 1, nkstot)
       if (noncolin) then
          WRITE ( ascunit,*) ' wk : ', ( wk ( ik_global ), ik_global = 1, nkstot)
       else
          WRITE ( ascunit,*) ' wk : ', ( wk ( ik_global )/2.0D0, ik_global = 1, nkstot)
       endif
       ! kpoints in crystal coordinates
       WRITE ( ascunit,*) ' xk : ', ( ( xk ( id, ik_global ), id = 1, nd ), ik_global = 1, nkstot)
       WRITE ( ascunit,*) ' ifmin : ', ( ifmin ( ik_global ), ik_global = 1, nkstot )
       WRITE ( ascunit,*) ' ifmax : ', ( ifmax ( ik_global ), ik_global = 1, nkstot )
       WRITE ( ascunit,*) ' et_g : ', ( ( et_g ( ib, ik_global ), ib = 1, nbnd ), ik_global = 1,nkstot )
       WRITE ( ascunit,*) ' wg_g : ', ( ( wg_g ( ib, ik_global ), ib = 1, nbnd ), ik_global = 1,nkstot )
       WRITE ( ascunit,*) ' nrecord = ', nrecord
       WRITE ( ascunit,*) ' ngm_g = ', ngm_g
       write(ascunit,*) 'g_g : '
       do i = 1, ngm_g
          WRITE( ascunit,'(I8,1X,I8,1X,I8)') g_g(1:3,i)
       enddo
    ENDIF

    DEALLOCATE ( wg_g )
    DEALLOCATE ( ifmax )
    DEALLOCATE ( ifmin )
    ALLOCATE ( igwk ( npwx_g ) )
    if (ionode) then
       ALLOCATE (g_useful(3,npwx_g))
    endif

    !! <<<<<<<<<<<<< pw_restart.f90 >>>>>>>>>>>>>
    !! ... Define a further local2global map to write gkvectors and
    !!     wfc coherently
    !! ALLOCATE ( igk_l2g_kdip( npwx_g, nks ) )
    !! igk_l2g_kdip = 0
    !! DO ik = iks, ike
    !! CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
    !! igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
    !! END DO
    !! <<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>
    !! ... For real case or nspin = 2, we need to process wf
    !! ... For complex case, we don't need to process wf

    IF ( proc_wf ) THEN
       IF ( MOD ( npwx_g, nproc ) .EQ. 0 ) THEN
          ngkdist_l = npwx_g / nproc
       ELSE
          ngkdist_l = npwx_g / nproc + 1
       ENDIF
       ngkdist_g = ngkdist_l * nproc
       IF ( real_or_complex .EQ. 1 ) &
            ALLOCATE ( energy ( nb, ns ) )
       ALLOCATE ( wfng_buf ( ngkdist_g, ns ) )
       ALLOCATE ( wfng_dist ( ngkdist_l, nb, ns ) )
    ENDIF

    DO i = 1, nk_g
       ik = kmap ( i )
       is = smap ( i )
       ik_global = ik

       IF ( real_or_complex .EQ. 1 ) THEN
          DO ib = 1, nb
             energy ( ib, is ) = et_g ( ib, i )
          ENDDO
       ENDIF

       DO j = 1, npwx_g
          igwk ( j ) = 0
       ENDDO
       !! ... itmp : map global Gvector index into global Gvector index
       !! ... itmp is for each kpoint
       !! ... itmp() will map a useful global Gvector index (in Fortran
       !!     the index is from 1) into itself, which is non-zero
       !!     and map the useless global Gvector index into 0
       !!     because for those useless global Gvectors igk_l2g is 0
       !! ... itmp() serves as a mask to mark the useful global Gvectors
       !!     and itmp() is a global mask, each process has a copy
       !!     of this global mask
       !! ... Compute the maximum global Gvector INDEX among all the kpoints
       ALLOCATE ( itmp ( npw_g ) )
       itmp = 0

       !! ... If the kpoint is in local pool
       IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
          DO ig = 1, ngk ( ik - iks + 1 )
             itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
          ENDDO
       ENDIF

       !! itmp on different cpus in the same pool has different structure
       !! for example itmp(1) = 1, then all the other itmp's on other cpus
       !! in the same pool must have itmp(1) = 0
       !! after sum, itmp is still for one kpoint

       ! ======
       ! DO NOT USE intra_bgrp_comm here
       !CALL mp_sum( itmp, intra_bgrp_comm )

       CALL mp_sum( itmp, world_comm )

       !! ... npw_g gives the upper bound of index for relevant
       !!     global Gvectors in all process/pools.
       !! ... Gvectors with index > npw_g will never be used in all pools
       !! ... itmp(ig) == ig means that g_g(1:3,ig) is used in some pools
       !! ... for example, in global list of Gvectors g_g(1:3, 1:ngm_g),
       !!     there are three consecutive Gvectors g1(1:3), g2(1:3), g3(1:3)
       !!     but only g1 and g3 will be used in the code,
       !! ... ig=1 ==> itmp(ig) == ig ==> ngg = 1 ==> igwk_(1) = ig = 1
       !! ... ig=2 ==> itmp(ig) != ig ==> ngg = 1
       !! ... ig=3 ==> itmp(ig) == ig ==> ngg = 2 ==> igwk_(2) = ig = 3
       !! ... All the useful Gvectors form a list g_useful(1:3,1:npwx_useful)
       !! ... igwk(:) is for each kpoint
       !! ... igwk(:) will map the index in g_useful(1:3,:) into the
       !!     corresponding index of this useful Gvector in global
       !!     list g_g(1:3, 1:ngm_g)
       !! ... ngg counts the number of useful Gvectors.
       ! ---------------------------------------------
       !! ... igwk(1:npwx_g) is a map special to one kpoint
       !!     created purely locally
       !! ... for each ik_global kpoint, igwk(1:ngk_g(ik_global)) are non-zero
       !!     igwk(ngk_g(ik_global)+1:npwx_g) = 0
       ngg = 0

       DO ig = 1, npw_g
          IF ( itmp ( ig ) .EQ. ig ) THEN
             ngg = ngg + 1
             igwk ( ngg ) = ig
          ENDIF
       ENDDO
       DEALLOCATE ( itmp )
       !! ... only done by the root of a pool
       !! ... there is no common kpoint shared by different pools
       !! ... for every pool, ig = 1, ngk_g(ik_global) run over all the
       !!     Gvectors for a kpoint in this pool, which has exactly the
       !!     same number of Gvectors in g_useful(:) list for this
       !!     kpoint, and g_g(1:3,igwk(ig)) will give the corresponding
       !!     Gvectors for all the wavefunctions associated with this kpoint.
       IF ( ionode ) THEN
          !! ... is == 1 for nspin = 1 (scalar relativistic)
          !!             and nspin = 4 (non-collinear SOC)
          !! ... We have written the global g_g(1:3, 1:ngm_g) list into 'WFN'
          IF ( is .EQ. 1 ) THEN
             WRITE ( unit ) nrecord
             WRITE ( unit ) ngk_g ( ik )
             ! for each ik_global, ig = 1:ngk_g(ik_global) will run over
             ! all the Gvectors used in (xk(1:3,ik_global)+G(1:3)) in this
             ! cpu, and has no common Gvectors with other cpus in the same
             ! pool, and the Gvectors outputed into WFN has exactly the
             ! same order in evc(:,ibnd) below
             ! g_useful(:) list
             WRITE ( unit ) ( ( g_g ( id, igwk ( ig ) ), id = 1, nd ), ig = 1, ngk_g ( ik ) )
             DO ig = 1, ngk_g( ik )
                DO id = 1, 3
                   g_useful(id, ig) = g_g (id, igwk(ig))
                ENDDO
             ENDDO
          ENDIF
       ENDIF

       local_pw = 0
       IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
          !! CALL davcio ( evc, 2*nwordwfc, iunwfc, ik - iks + 1, - 1 )
          !! READ the local ik_local-th wavefunction into evc(1:npwx, 1:nbnd)
          !! only ONE kpoints but with all the bands in evc

          !! nwordwfc = nbnd*npwx*npol
          !! for spinor wavefunction : npol = 2
          !! for scalar wavefunction : npol = 1
          !! 2*nwordwfc because we deal with complex numbers
          !! for spinor wavefunction, with ik_global-th kpoint and ib-th band
          !! evc(1:npwx,ib) is upper component
          !! evc(npwx+1:2*npwx,ib) is lower component
          !! - ALLOCATE( evc( npwx*npol, nbnd ) )

          !! after davcio(...), evc(1:npwx*npol,ib) will contain part of the
          !! wavefunction for ib-th band and ik_global-th kpoints stored
          !! in local process, that is, for (ik_global, ib) state, there
          !! should be ngk_g(ik_global) Gvectors, but now evc(1:npwx*npol)
          !! only contains ngk(ik_local) Gvectors.

          !! for spinor, 2*nwordwfc = 2* (nbnd*npwx*npol) = 2* (nbnd*npwx*2)
          !! and therefore evc(1:ngk(ik_local),ib) will contains the upper
          !! component (complex), evc(ngk(ik_local)+1:npwx,ib) = 0.0 + I 0.0
          !! evc(npwx+1:npwx+ngk(ik_local),ib) will contains the lower
          !! component (complex), evc(npwx+ngk(ik_local)+1:2*npwx) = 0 + I 0.0
          CALL davcio ( evc, 2*nwordwfc, iunwfc, ik - iks + 1, - 1 )
          !! local_pw = ng_local : local number of G vectors for a given kpoint = ngk(ik_local)
          local_pw = ngk ( ik - iks + 1 )
       ENDIF

       !! igwf_l2g (1:ng_local) is a local map, the domain is the local
       !! gindex for each kpoint pool
       ALLOCATE ( igwf_l2g ( local_pw ) )
       !! for ( ik_global < iks ) .OR. ( ik_global > ike )
       !! that is, for those kpoints not within current pool, ng_local = 0
       !! and following loop will not be carried out
       DO ig = 1, local_pw
          igwf_l2g ( ig ) = 0
       ENDDO
       !! for ( ik_global < iks ) .OR. ( ik_global > ike )
       !! that is, for those kpoints not within current pool, ng_local = 0
       !! and following loop will not be carried out
       DO ig = 1, local_pw
          !! ... ngg stores the global index for each local Gvector
          ngg = igk_l2g ( ig, ik - iks + 1 )
          !! ... for each ik_global within this pool, j runs over
          !!     all the Gvectors for this kpoint (distributed in all cpus)
          !!     because igwk(:)'s domain is exactly 1:ngk_g(ik_global)
          !!     that is, for each ik_global, igwk(1:ngk_g(ik_global)) is
          !!     non-zero, and igwk(ngk_g(ik_global)+1:npwx_g) = 0
          DO j = 1, ngk_g ( ik )
             !! ... igwk()'s range is global g index
             !!     if ngg == igwk(j)
             !!     then it means that ig-th local g index's image in
             !!     global gindex is actually the j-th useful g vector
             !!     in useful g list g_useful(:).
             !!
             !! ... igwf_l2g(:) is again specific to one kpoint ik_global
             !!     igwf_l2g(:) will map the ig-th local g index into
             !!     its corresponding index in g_useful(:) list
             !! ... igwf_l2g(:)'s domain is just 1:ngk(ik_local)
             IF ( ngg .EQ. igwk ( j ) ) THEN
                igwf_l2g ( ig ) = j
                EXIT
             ENDIF
          ENDDO
       ENDDO

       ALLOCATE ( ipmask ( nproc ) )
       DO j = 1, nproc
          ipmask ( j ) = 0
       ENDDO
       !! ionode_id : index of the i/o node for this image
       !! ionode : ! true if this processor is a i/o node for this image
       ipsour = ionode_id
       !! <<<<<<<<<< mp_pool.f90 >>>>>>>>>>
       !! npool : number of "k-points"-pools
       !! nproc_pool : number of processors within a pool
       !! me_pool : index of the processor within a pool
       !! root_pool : index of the root processor within a pool
       !! my_pool_id : index of my pool

       !! <<<<<<<<<< mp_world.f90 >>>>>>>>>>
       !! nproc : number of processors
       !! mpime : processor index (starts from 0 to nproc-1)
       !! root : index of the root processor
       !! ... If we have more than 1 kpoint pool
       IF ( npool .GT. 1 ) THEN
          !! ... for ik_global not in current pool, we will
          !! enter this IF statement, and therefore for the root process
          !! of other pool, we still have ipmask == 0
          IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
             !! ... All the root_pool cpus are masked as 1
             IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
          ENDIF
          !! ... After the mp_sum command only the root process of
          !!     the pool containing ik_global is masked as 1
          CALL mp_sum ( ipmask, world_comm )
          !! ipsour is a global variable created purely locally.
          DO j = 1, nproc
             !! ... if j == 1, that is, j is the ionode, the root cpu
             !!     of World_Comm, then ipsour = 1 - 1 = 0
             !! ... for the root cpus of other pools, we will not enter
             !!     the IF statement, and ipsour = 0
             !! ... for the root cpu of this pool, and this root cpu is not
             !!     the World_Comm's root, then ipsour /= 0
             !! ... In a word, for each ik_global kpoint, ipsour is only
             !!     nonzero for a root cpu for current pool and this cpu
             !!     is not the World_Comm's root.
             IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
          ENDDO
       ENDIF
       DEALLOCATE ( ipmask )

       igwx = 0
       ierr = 0
       !! ... igwx is a pool-wide variable, that is, it has the same value
       !!     in all cpus of one pool.
       !! ... igwx is specific to one kpoint, ik_global.
       !! ... igwx gives the maximal index of Gvectors in g_useful(:) list
       !!     in current pool, which is exactly ngk_g ( ik_global )
       !! ... intra_pool_comm : cpus in the same kpoint pool
       !! ... for ik_global < iks .OR. ik_global > ike, we have igwx = 0
       IF ( ik .GE. iks .AND. ik .LE. ike ) &
            igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) ) ! === ngk_g(ik_global)

       !! ... if ik_global is not in current pool, igwx = 0
       CALL mp_max ( igwx, intra_pool_comm )

       !> [WORKING]
       ! !! if ik_global is not in current pool, igwx=0*2=0
       ! IF (noncolin) THEN
       !    igwx = igwx * 2
       ! ENDIF

       !! ... if ipsour \= 0, that is, ipsour is a root_pool cpu
       !! ... For each ik_global, it must be in some kpoint pool
       !!     and for this pool, there is a root process where
       !!     ipsour /= 0, that is ipsour /= ionode_id
       !!     for the roots of other pools, we still have ipsour == 0
       IF ( ipsour .NE. ionode_id ) THEN
          !! --------------------------------
          !! mp_get( msg_dest, msg_sour, mpime, dest, sour, ip, gid)
          !! will send msg_sour from sour-th process in gid-th group to
          !! msg_dest from dest-th process in the same gid-th group
          !! ip is the tag of MPI_SEND or MPI_RECV, and is useless herer
          !! --------------------------------
          !! ... send the data igwx from ipsour-th process in
          !!     world_comm to igwx from ionode_id-th (ROOT) process in world_comm
          !! ======
          !! ROOT proc know the size of wavefunction for the ik_global-th kvector, because ROOT proc is responsible to output this wavefunction
          !! ------
          !! ... tag = ip = 1
          !! ... only able to run when current process
          !!     mpime = ipsour or ionode_id
          CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
          !! ... make sure that ipsour and ionode_id has the same igwx
       ENDIF
       ierr = 0
       !! ... for a fixed ik_global, igwx should always be equal to
       !!     ngk_g(ik_global)
       !! ... ik_global < iks .OR. ik_global > ike, that is, ik_global
       !!     is not in current pool, even though igwx = 0 /= ngk_g(ik_global)
       !!     we will not enter this IF statement, and therefore ierr = 0

       ! IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
       !      ierr = 1
       ! CALL mp_max ( ierr, world_comm )
       ! IF ( ierr .GT. 0 ) &
       !      CALL errore ( 'write_wfng', 'igwx ngk_g', ierr )

       !! ------------------------------------------------
       !! ... for ik_global in current pool, we will have
       !!     wfng ( 1 : ngk_g(ik_global) ) stores wavefunction for one band
       !!     and for one kpoint (ik_global) with all the Gvectors.
       !! ... for ik_global not in current pool, igwx = 0
       !!     we have wfng(1:1), and that's why we say wfng(:) stays in the ROOT CPU
       !!     of each pool
       !> [WORKING]
       ALLOCATE ( wfng ( MAX ( npol*igwx, 1 ) ) )
       if (npol == 2) wfng2 => wfng(igwx+1:2*igwx)

       DO ib = 1, nb
          !> for kpoint ik_global not in current pool, igwx = 0, and we will not enter this loop.
          wfng = ( 0.0D0, 0.0D0 )
          !> if we have more than 1 kpoint pools
          IF ( npool .GT. 1 ) THEN
             IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
                !! ... merge the evc(:,ib) for the ib-th band into the total
                !! ... wavefunction wfng(1:ngk_g(ik_global)) [complex]
                !!     for one kpoint, ik_global
                !! -----------------------------------------------------------
                !! mergewf ( pw, pwt, ngwl, ig_l2g, mpime, nproc, root, comm )
                !! -----------------------------------------------------------
                !! pw = evc(:,ib)
                !! pwt = wfng, stored in root process == 0
                !! ngwl = ng_local = local number of G vectors for a given kpoint
                !!        ng_local = ngk(ik_local)
                !! ig_l2g = igwf_l2g = map the ig-th local g index into
                !!          its corresponding index in g_useful(:) list
                !! mpime = me_pool = index of the calling processor within a pool
                !! nproc = nproc_pool = number of processors within a pool
                !! root = root_pool = the process that should receive the data
                !!       for root_pool cpus, ipsour = its mpime - 1 /= 0
                !! comm = intra_pool_comm = comm in a kpoint pool
                !! ... wfng(1:ngk_g(ik_global)) stores the wavefunction
                !!     of a single band for ik_global-th kpoints
                !! ... different cpus in a pool has the same set of kpoints
                !!     but for the same kpoint, different cpus has different
                !!     Gvectors.
                !! ... wfng(1:ngk_g(ik_global)) stays in the root cpu of the
                !!     pool containing ik_global kpoint.
                !!
                !! ... collect wavefunction within a pool into the root_pool
                !!     because for one ik_global kpoint, it could only stay
                !!     in one kpoint pool, with all its Gvectors distributed
                !!     in all the cpus in this pool with no common element.

                !> [WORKING]
                if (npol == 2) then
                   CALL mergewf ( evc(1:npwx, ib),         wfng, local_pw, igwf_l2g, me_pool, nproc_pool, root_pool, intra_pool_comm )
                   CALL mergewf ( evc(npwx+1:2*npwx, ib), wfng2, local_pw, igwf_l2g, me_pool, nproc_pool, root_pool, intra_pool_comm )
                else
                   CALL mergewf ( evc(:, ib),              wfng, local_pw, igwf_l2g, me_pool, nproc_pool, root_pool, intra_pool_comm )
                endif

                ! IF (noncolin) THEN
                !    ! ======
                !    CALL mergewf_spinor( evc(:,ib), wfng, ngk_g(ik), local_pw, igwf_l2g, me_pool, nproc_pool, root_pool, intra_pool_comm, npwx)
                !    ! ------
                !    !! After this mergewf_spinor, each pool-ROOT proc has the PWs for ib-th band and ik_global-th kvector.
                ! ELSE ! nspinor == 1
                !    CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, me_pool, nproc_pool, root_pool, intra_pool_comm )
                ! ENDIF

             ENDIF
             !! ... if ipsour \= 0, that is, ipsour is the root cpu of current
             !!     pool containing ik_global, and this process is not ionode
             !! ... ionode is also the root cpu for World_Comm
             !! ... collect wavefunction from the target pool's
             IF ( ipsour .NE. ionode_id ) THEN
                CALL mp_get ( wfng, wfng, mpime, ionode_id, ipsour, ib, world_comm )
             ENDIF
          ELSE
             !> [WORKING]
             if (npol == 2) then
                CALL mergewf ( evc(1:npwx, ib),         wfng, local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
                CALL mergewf ( evc(npwx+1:2*npwx, ib), wfng2, local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
             else
                CALL mergewf ( evc(:, ib),              wfng, local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
             endif

             ! IF (noncolin) THEN
             !    CALL mergewf_spinor( evc(:,ib), wfng, ngk_g(ik), local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm, npwx)
             !    ! After this mergewf_spinor, each pool-ROOT proc has the PWs for ib-th band and ik_global-th kvector.
             ! ELSE ! nspinor == 1
             !    CALL mergewf ( evc ( :, ib ), wfng, local_pw, igwf_l2g, mpime, nproc, ionode_id, world_comm )
             ! ENDIF

          ENDIF
          ! <- nspin = 2 or (real_wfn, nspin = 1 )->
          IF ( proc_wf ) THEN
             DO ig = 1, igwx
                wfng_buf ( ig, is ) = wfng ( ig )
             ENDDO
             DO ig = igwx + 1, ngkdist_g
                wfng_buf ( ig, is ) = ( 0.0D0, 0.0D0 )
             ENDDO
#if defined(__MPI)
             CALL mp_barrier ( world_comm )
             CALL MPI_Scatter ( wfng_buf ( :, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
                  wfng_dist ( :, ib, is ), ngkdist_l, MPI_DOUBLE_COMPLEX, &
                  ionode_id, world_comm, ierr )
             IF ( ierr .GT. 0 ) &
                  CALL errore ( 'write_wfng', 'mpi_scatter', ierr )
#else
             DO ig = 1, ngkdist_g
                wfng_dist ( ig, ib, is ) = wfng_buf ( ig, is )
             ENDDO
#endif
             !! output wfng into 'WFN'
          ELSE
             !! ... for any ik_global kpoint, it must be in one kpoint pool
             !!     and we have sent wavefunction wfng(:) to wfng(:) in
             !!     ionode, and we have also sent the dimension of wavefunction
             !!     igwx==ngk_g(ik_global) to ionode. In this way, each
             !!     ik_global kpoint's wavefunction is stored in ionode
             !!     and ready to be outputed.
             !! ===============================
             !! (old) For noncolin case, igwx==2*ngk_g(ik_global)
             IF ( ionode ) THEN
                !> Check the norm of each wavefunctions
                norm = DCMPLX(0D0, 0D0)
                !> [WORKING]
                DO ig = 1, npol*igwx
                   norm = norm + CONJG(wfng(ig))*wfng(ig)
                ENDDO

                ! write(*,'(1X,A,I5,A,I5,A,2ES20.12)') "ik #", ik, " ib #", ib, " NORM = ", norm
                IF ( ABS(DBLE(norm)-1.0D0) > norm_tol ) THEN
                   write(*,'(1X,A,I5,A,I5,A,ES20.12)') "ik #", ik, " ib #", ib, " RE[NORM] = ", DBLE(norm)
                   call exit(4321)
                ENDIF

                IF ( ABS(DIMAG(norm)) > norm_tol ) THEN
                   write(*,'(1X,A,I5,A,I5,A,ES20.12)') "ik #", ik, " ib #", ib, " Im[NORM] = ", DIMAG(norm)
                   call exit(1234)
                ENDIF

                WRITE ( unit ) nrecord
                WRITE ( unit ) ngk_g ( ik )
                ! WRITE ( unit ) ( wfng ( ig ), ig = 1, igwx )
                !> [WORKING]
                WRITE ( unit ) ( wfng ( ig ), ig = 1, npol*igwx )

                IF (noncolin .and. (flag_output_spin) ) THEN
                   ! Calculate expectation value of spin
                   sigma_x = 0.0D0
                   sigma_y = 0.0D0
                   sigma_z = 0.0D0
                   DO ig = 1, igwx
                      sigma_x = sigma_x + wfng(ig)*CONJG(wfng(ig+igwx)) + CONJG(wfng(ig))*wfng(ig+igwx)
                      sigma_y = sigma_y + imagunit*(wfng(ig)*CONJG(wfng(ig+igwx)) - CONJG(wfng(ig))*wfng(ig+igwx))
                      sigma_z = sigma_z + CONJG(wfng(ig))*wfng(ig) - CONJG(wfng(ig+igwx))*wfng(ig+igwx)
                   ENDDO
                   ! DO ig = 1, igwx/2
                   !    sigma_x = sigma_x + wfng(ig)*CONJG(wfng(ig+igwx/2)) + CONJG(wfng(ig))*wfng(ig+igwx/2)
                   !    sigma_y = sigma_y + imagunit*(wfng(ig)*CONJG(wfng(ig+igwx/2)) - CONJG(wfng(ig))*wfng(ig+igwx/2))
                   !    sigma_z = sigma_z + CONJG(wfng(ig))*wfng(ig) - CONJG(wfng(ig+igwx/2))*wfng(ig+igwx/2)
                   ! ENDDO

                   WRITE (unit_spin_x, "(I4, 5X, I4, 5X,ES20.12, 5X , ES20.12)") ik, ib, sigma_x
                   WRITE (unit_spin_y, "(I4, 5X, I4, 5X,ES20.12, 5X , ES20.12)") ik, ib, sigma_y
                   WRITE (unit_spin_z, "(I4, 5X, I4, 5X,ES20.12, 5X , ES20.12)") ik, ib, sigma_z
                ENDIF

                if (flag_output_asc_wfn) then
                   unit_wfn = ik_global+1000*ib
                   write(ik_global_string,'(I5)') ik_global
                   write(ib_string,'(I5)') ib
                   file_wfn = 'wfn.bin.k'//trim(adjustl(ik_global_string))//".b"//trim(adjustl(ib_string))//'.dat'

                   open (unit = unit_wfn, file = trim(adjustl(file_wfn)), form = 'formatted', status = 'replace' )

                   write (unit_wfn,"(A,I4,A,I4)") '# ik #', ik_global, ' ib #', ib
                   write (unit_wfn,"(A,'(',G20.12,G20.12,G20.12,5X,')')") '# kvector = ', xk(1:3,ik_global)
                   write (unit_wfn,"(A,I5)") '# ngk = ', ngk_g(ik_global)
                   write (unit_wfn,"(A,G20.12)") '# E = ', et_g(ib,ik_global)

                   if (noncolin) then
                      ! write (unit_wfn, "(3I6,2X,ES20.12,5X,ES20.12,5X,ES20.12,5X,ES20.12)") ( g_useful(1:3,ig), wfng(ig), wfng(ig+igwx/2), ig = 1, igwx/2)
                      write (unit_wfn, "(3I6,2X,ES20.12,5X,ES20.12,5X,ES20.12,5X,ES20.12)") ( g_useful(1:3,ig), wfng(ig), wfng(ig+igwx), ig = 1, igwx)
                   else
                      write (unit_wfn, "(3I6,2X,ES20.12,5X,ES20.12,5X)") ( g_useful(1:3,ig), wfng(ig), ig = 1, igwx)
                   endif

                   close(unit_wfn)
                endif
             ENDIF
          ENDIF
       ENDDO ! band loop

       DEALLOCATE ( wfng )
       DEALLOCATE ( igwf_l2g )

       IF ( proc_wf .AND. is .EQ. ns ) THEN
          IF ( real_or_complex .EQ. 1 ) THEN
             CALL start_clock ( 'real_wfng' )
             CALL real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )
             CALL stop_clock ( 'real_wfng' )
          ENDIF
          DO ib = 1, nb
             DO is = 1, ns
#if defined(__MPI)
                CALL mp_barrier ( world_comm )
                CALL MPI_Gather ( wfng_dist ( :, ib, is ), ngkdist_l, &
                     MPI_DOUBLE_COMPLEX, wfng_buf ( :, is ), ngkdist_l, &
                     MPI_DOUBLE_COMPLEX, ionode_id, world_comm, ierr )
                IF ( ierr .GT. 0 ) &
                     CALL errore ( 'write_wfng', 'mpi_gather', ierr )
#else
                DO ig = 1, ngkdist_g
                   wfng_buf ( ig, is ) = wfng_dist ( ig, ib, is )
                ENDDO
#endif
             ENDDO
             IF ( ionode ) THEN
                WRITE ( unit ) nrecord
                WRITE ( unit ) ngk_g ( ik )
                IF ( real_or_complex .EQ. 1 ) THEN
                   WRITE ( unit ) ( ( dble ( wfng_buf ( ig, is ) ), &
                        ig = 1, igwx ), is = 1, ns )
                ELSE
                   WRITE ( unit ) ( ( wfng_buf ( ig, is ), &
                        ig = 1, igwx ), is = 1, ns )
                ENDIF
             ENDIF
          ENDDO
       ENDIF
    ENDDO ! loop kpoints

    IF (ionode) THEN
       DEALLOCATE ( g_useful )
    ENDIF

    DEALLOCATE ( igwk )
    DEALLOCATE ( ngk_g )
    DEALLOCATE ( igk_l2g )
    DEALLOCATE ( et_g )

    IF ( proc_wf ) THEN
       IF ( real_or_complex .EQ. 1 ) &
            DEALLOCATE ( energy )
       DEALLOCATE ( wfng_buf )
       DEALLOCATE ( wfng_dist )
    ENDIF

    IF ( ionode ) THEN
       CLOSE ( unit = unit, status = 'keep' )
       IF (flag_output_spin) THEN
          CLOSE ( unit = unit_spin_x)
          CLOSE ( unit = unit_spin_y)
          CLOSE ( unit = unit_spin_z)
       ENDIF
    ENDIF

    DEALLOCATE ( g_g )
    DEALLOCATE ( smap )
    DEALLOCATE ( kmap )

    ! nk_g / ns --> nkstot
    CALL cryst_to_cart ( nkstot, xk, bg, 1 )
    CALL mp_barrier ( world_comm )

    RETURN
101 FORMAT ( /, 5X, "WARNING: kgrid is set to zero in the wavefunction file.", &
         /, 14X, "The resulting file will only be usable as the fine grid in inteqp.", / )
  END SUBROUTINE write_wfng

  !-------------------------------------------------------------------------------
  SUBROUTINE real_wfng ( ik, ngkdist_l, nb, ns, energy, wfng_dist )

    USE kinds, ONLY : DP
    USE io_global, ONLY : ionode
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : world_comm

    IMPLICIT NONE

    integer, intent (in) :: ik, ngkdist_l, nb, ns
    real (DP), intent (in) :: energy ( :, : ) ! ( nb, ns )
    complex (DP), intent (inout) :: wfng_dist ( :, :, : ) ! ( ngkdist_l, nb, ns )

    real (DP), PARAMETER :: eps2 = 1.0D-2
    real (DP), PARAMETER :: eps5 = 1.0D-5
    real (DP), PARAMETER :: eps6 = 1.0D-6

    character :: tmpstr*80
    integer :: i, j, k, is, ib, jb, ig, inum, deg, mdeg, inc
    integer :: dimension_span, reduced_span, ierr
    real (DP) :: x
    integer, allocatable :: imap ( :, : )
    integer, allocatable :: inums ( : )
    integer, allocatable :: inull ( : )
    integer, allocatable :: null_map ( :, : )
    real (DP), allocatable :: psi ( :, : )
    real (DP), allocatable :: phi ( :, : )
    real (DP), allocatable :: vec ( : )
    complex (DP), allocatable :: wfc ( : )

    mdeg = 1
    DO is = 1, ns
       DO ib = 1, nb - 1
          deg = 1
          DO jb = ib + 1, nb
             IF ( abs ( energy ( ib, is ) - energy ( jb, is ) ) &
                  .LT. eps5 * dble ( jb - ib + 1 ) ) deg = deg + 1
          ENDDO
          IF ( deg .GT. mdeg ) mdeg = deg
       ENDDO
    ENDDO
    mdeg = mdeg * 2

    ALLOCATE ( imap ( nb, ns ) )
    ALLOCATE ( inums ( ns ) )
    ALLOCATE ( inull ( nb ) )
    ALLOCATE ( null_map ( mdeg, nb ) )

    DO is = 1, ns
       inum = 1
       DO ib = 1, nb
          IF ( ib .EQ. nb ) THEN
             imap ( inum, is ) = ib
             inum = inum + 1
          ELSEIF ( abs ( energy ( ib, is ) - &
               energy ( ib + 1, is ) ) .GT. eps5 ) THEN
             imap ( inum, is ) = ib
             inum = inum + 1
          ENDIF
       ENDDO
       inum = inum - 1
       inums ( is ) = inum
    ENDDO

    ALLOCATE ( wfc ( ngkdist_l ) )
    ALLOCATE ( psi ( ngkdist_l, mdeg ) )
    ALLOCATE ( phi ( ngkdist_l, mdeg ) )
    ALLOCATE ( vec ( ngkdist_l ) )

    DO is = 1, ns
       inc = 1
       inum = inums ( is )
       DO i = 1, inum
          inull ( i ) = 1
          DO ib = inc, imap ( i, is )
             DO ig = 1, ngkdist_l
                wfc ( ig ) = wfng_dist ( ig, ib, is )
             ENDDO
             x = 0.0D0
             DO ig = 1, ngkdist_l
                x = x + dble ( wfc ( ig ) ) **2
             ENDDO
             CALL mp_sum ( x, world_comm )
             IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
             IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
             inull ( i ) = inull ( i ) + 1
             x = 0.0D0
             DO ig = 1, ngkdist_l
                x = x + aimag ( wfc ( ig ) ) **2
             ENDDO
             CALL mp_sum ( x, world_comm )
             IF ( x .LT. eps2 ) null_map ( inull ( i ), i ) = 0
             IF ( x .GT. eps2 ) null_map ( inull ( i ), i ) = 1
             inull ( i ) = inull ( i ) + 1
          ENDDO
          inull ( i ) = inull ( i ) - 1
          inc = imap ( i, is ) + 1
       ENDDO
       inc = 1
       ib = 1
       DO i = 1, inum
          k = 1
          DO j = 1, 2 * ( imap ( i, is ) - inc ) + 1, 2
             IF ( null_map ( j, i ) .EQ. 1 .OR. &
                  null_map ( j + 1, i ) .EQ. 1 ) THEN
                DO ig = 1, ngkdist_l
                   wfc ( ig ) = wfng_dist ( ig, ib, is )
                ENDDO
                IF ( null_map ( j, i ) .EQ. 1 ) THEN
                   DO ig = 1, ngkdist_l
                      phi ( ig, k ) = dble ( wfc ( ig ) )
                   ENDDO
                   k = k + 1
                ENDIF
                IF ( null_map ( j + 1, i ) .EQ. 1 ) THEN
                   DO ig = 1, ngkdist_l
                      phi ( ig, k ) = aimag ( wfc ( ig ) )
                   ENDDO
                   k = k + 1
                ENDIF
                ib = ib + 1
             ENDIF
          ENDDO
          dimension_span = k - 1
          IF ( dimension_span .EQ. 0 ) THEN
             ierr = 201
             WRITE ( tmpstr, 201 ) ik, is, inc
             CALL errore ( 'real_wfng', tmpstr, ierr )
          ENDIF
          DO j = 1, dimension_span
             x = 0.0D0
             DO ig = 1, ngkdist_l
                x = x + phi ( ig, j ) **2
             ENDDO
             CALL mp_sum ( x, world_comm )
             x = sqrt ( x )
             DO ig = 1, ngkdist_l
                phi ( ig, j ) = phi ( ig, j ) / x
             ENDDO
          ENDDO
          !
          ! the Gram-Schmidt process begins
          !
          reduced_span = 1
          DO ig = 1, ngkdist_l
             psi ( ig, 1 ) = phi ( ig, 1 )
          ENDDO
          DO j = 1, dimension_span - 1
             DO ig = 1, ngkdist_l
                vec ( ig ) = phi ( ig, j + 1 )
             ENDDO
             DO k = 1, reduced_span
                x = 0.0D0
                DO ig = 1, ngkdist_l
                   x = x + phi ( ig, j + 1 ) * psi ( ig, k )
                ENDDO
                CALL mp_sum ( x, world_comm )
                DO ig = 1, ngkdist_l
                   vec ( ig ) = vec ( ig ) - psi ( ig, k ) * x
                ENDDO
             ENDDO
             x = 0.0D0
             DO ig = 1, ngkdist_l
                x = x + vec ( ig ) **2
             ENDDO
             CALL mp_sum ( x, world_comm )
             x = sqrt ( x )
             IF ( x .GT. eps6 ) THEN
                reduced_span = reduced_span + 1
                DO ig = 1, ngkdist_l
                   psi ( ig, reduced_span ) = vec ( ig ) / x
                ENDDO
             ENDIF
          ENDDO
          !
          ! the Gram-Schmidt process ends
          !
          IF ( reduced_span .LT. imap ( i, is ) - inc + 1 ) THEN
             ierr = 202
             WRITE ( tmpstr, 202 ) ik, is, inc
             CALL errore ( 'real_wfng', tmpstr, ierr )
          ENDIF
          DO ib = inc, imap ( i, is )
             DO ig = 1, ngkdist_l
                wfng_dist ( ig, ib, is ) = &
                     CMPLX ( psi ( ig, ib - inc + 1 ), 0.0D0, KIND=dp )
             ENDDO
          ENDDO
          inc = imap ( i, is ) + 1
       ENDDO
    ENDDO

    DEALLOCATE ( vec )
    DEALLOCATE ( phi )
    DEALLOCATE ( psi )
    DEALLOCATE ( wfc )
    DEALLOCATE ( null_map )
    DEALLOCATE ( inull )
    DEALLOCATE ( inums )
    DEALLOCATE ( imap )

    RETURN

201 FORMAT("failed Gram-Schmidt dimension span for kpoint =",i6," spin =",i2," band =",i6)
202 FORMAT("failed Gram-Schmidt reduced span for kpoint =",i6," spin =",i2," band =",i6)

  END SUBROUTINE real_wfng

  !-------------------------------------------------------------------------------

  SUBROUTINE write_rhog ( output_file_name, real_or_complex, symm_type, &
       rhog_nvmin, rhog_nvmax )

    USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
    USE constants, ONLY : pi, tpi, eps6
    USE fft_base, ONLY : dfftp
    USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
    USE io_global, ONLY : ionode
    USE ions_base, ONLY : nat, atm, ityp, tau
    USE kinds, ONLY : DP
    USE lsda_mod, ONLY : nspin
    USE mp, ONLY : mp_sum
    ! USE mp_bands, ONLY : intra_bgrp_comm
    USE mp_pools, ONLY : intra_pool_comm
    USE scf, ONLY : rho
    USE symm_base, ONLY : s, ftau, nsym, ft
    USE matrix_inversion

    IMPLICIT NONE

    character ( len = 256 ), intent (in) :: output_file_name
    integer, intent (in) :: real_or_complex
    character ( len = 9 ), intent (in) :: symm_type
    integer, intent (in) :: rhog_nvmin
    integer, intent (in) :: rhog_nvmax

    character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
    integer :: unit, id, is, ig, i, j, k, ierr
    integer :: nd, ns, ng_l, ng_g
    integer :: ntran, cell_symmetry, nrecord
    real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
    real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
    real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
    integer, allocatable :: g_g ( :, : )
    complex (DP), allocatable :: rhog_g ( :, : )

    INTEGER, EXTERNAL :: atomic_number
    INTEGER :: itran
    integer, dimension(1:3,1:3,1:nsym) :: s_transpose
    real(DP), dimension(1:3,1:nsym) :: ft_minus

    CALL date_and_tim ( cdate, ctime )
    WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
    WRITE ( stime, '(A8,24X)' ) ctime(1:8)
    IF ( real_or_complex .EQ. 1 ) THEN
       WRITE ( stitle, '("RHO-Real",24X)' )
    ELSE
       WRITE ( stitle, '("RHO-Complex",21X)' )
    ENDIF

    unit = 4
    nrecord = 1
    nd = 3

    ns = nspin
    ng_l = ngm
    ng_g = ngm_g

    ierr = 0
    IF ( ibrav .EQ. 0 ) THEN
       IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
          cell_symmetry = 0
       ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
          cell_symmetry = 1
       ELSE
          ierr = 1
       ENDIF
    ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
       cell_symmetry = 0
    ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
       cell_symmetry = 1
    ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
       cell_symmetry = 0
    ELSE
       ierr = 1
    ENDIF
    IF ( ierr .GT. 0 ) &
         CALL errore ( 'write_rhog', 'cell_symmetry', ierr )

    ntran = nsym

    ! DO i = 1, ntran
    !    DO j = 1, nd
    !       DO k = 1, nd
    !          r1 ( k, j ) = dble ( s ( k, j, i ) )
    !       ENDDO
    !    ENDDO
    !    CALL invmat ( 3, r1, r2 )
    !    t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
    !    t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
    !    t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
    !    DO j = 1, nd
    !       t2 ( j ) = 0.0D0
    !       DO k = 1, nd
    !          t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
    !       ENDDO
    !       IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
    !            t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
    !       IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
    !            t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
    !    ENDDO
    !    DO j = 1, nd
    !       translation ( j, i ) = t2 ( j ) * tpi
    !    ENDDO
    ! ENDDO

    ! CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

    do itran = 1, ntran
       ! ft_minus(1:3,itran) = -ft(1:3,itran)
       do id = 1, 3
          !> If an entry of ft is 0.5, ft_minus = ft = 0.5
          if (ABS(ft(id, itran) - 0.5D0) < 1.0D-12) then
             ft_minus(id, itran) = 0.5D0
          else
             ft_minus(id, itran) = - ft(id, itran)
          endif
       enddo
       s_transpose(1:3,1:3,itran) = transpose(s(1:3,1:3,itran))
    enddo

    alat2 = alat ** 2
    recvol = 8.0D0 * pi**3 / omega

    DO i = 1, nd
       DO j = 1, nd
          adot ( j, i ) = 0.0D0
       ENDDO
    ENDDO
    DO i = 1, nd
       DO j = 1, nd
          DO k = 1, nd
             adot ( j, i ) = adot ( j, i ) + &
                  at ( k, j ) * at ( k, i ) * alat2
          ENDDO
       ENDDO
    ENDDO

    DO i = 1, nd
       DO j = 1, nd
          bdot ( j, i ) = 0.0D0
       ENDDO
    ENDDO
    DO i = 1, nd
       DO j = 1, nd
          DO k = 1, nd
             bdot ( j, i ) = bdot ( j, i ) + &
                  bg ( k, j ) * bg ( k, i ) * tpiba2
          ENDDO
       ENDDO
    ENDDO
    ! ======
    ! (default) rhog_nvmin = 0, rhog_nvmax = 0
    ! ------
    ! If we skip following IF:
    ! We will use rho%of_g(ig_local,is), which is the charge density of all the valence electrons
    !     write(*,*) "rhog_nvmin = ", rhog_nvmin
    !     write(*,*) "rhog_nvmax = ", rhog_nvmax
    ! ------
    ! After following IF, we have prepared rho%of_g
    IF ( rhog_nvmin .NE. 0 .AND. rhog_nvmax .NE. 0 ) &
         CALL calc_rhog ( rhog_nvmin, rhog_nvmax )

    ALLOCATE ( g_g ( nd, ng_g ) )
    ALLOCATE ( rhog_g ( ng_g, ns ) )

    DO ig = 1, ng_g
       DO id = 1, nd
          g_g ( id, ig ) = 0
       ENDDO
    ENDDO
    DO is = 1, ns
       DO ig = 1, ng_g
          rhog_g ( ig, is ) = ( 0.0D0, 0.0D0 )
       ENDDO
    ENDDO
    ! ======
    ! ig_l2g(ig_local) = ig_global maps the index of local GVector list (G_i, i=1,ngm) into
    ! global Gvector list g_g(1:3,1:ngm_g)
    DO ig = 1, ng_l
       g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
       g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
       g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
    ENDDO
    ! We just rearrange and collect rho%of_g(ig_local,is) in each proc
    ! Each pool has a complete rhog_g(ig_global,is) on ROOT-of-Pool
    DO is = 1, ns
       DO ig = 1, ng_l
          rhog_g ( ig_l2g ( ig ), is ) = rho%of_g ( ig, is )
       ENDDO
    ENDDO

    !! ======
    !! intra_pool_comm = intra_bgrp_comm
    !CALL mp_sum ( g_g, intra_bgrp_comm )
    !CALL mp_sum ( rhog_g, intra_bgrp_comm )
    CALL mp_sum ( g_g, intra_pool_comm )
    CALL mp_sum ( rhog_g, intra_pool_comm )

    ! ======
    ! Multiply volumn of a cell
    DO is = 1, ns
       DO ig = 1, ng_g
          rhog_g ( ig, is ) = rhog_g ( ig, is ) * CMPLX ( omega, 0.0D0, KIND=dp )
       ENDDO
    ENDDO

    IF ( ionode ) THEN
       OPEN ( unit = unit, file = TRIM ( output_file_name ), &
            form = 'unformatted', status = 'replace' )
       WRITE ( unit ) stitle, sdate, stime
       WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
       WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
       WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
            ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
       WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
            ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
       ! WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
       ! WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
       WRITE ( unit ) ( ( ( s_transpose ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
       WRITE ( unit ) ( ( ft_minus ( j, i ), j = 1, nd ), i = 1, ntran )
       WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
       WRITE ( unit ) nrecord
       WRITE ( unit ) ng_g
       WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
       WRITE ( unit ) nrecord
       WRITE ( unit ) ng_g
       IF ( real_or_complex .EQ. 1 ) THEN
          WRITE ( unit ) ( ( dble ( rhog_g ( ig, is ) ), &
               ig = 1, ng_g ), is = 1, ns )
       ELSE
          WRITE ( unit ) ( ( rhog_g ( ig, is ), &
               ig = 1, ng_g ), is = 1, ns )
       ENDIF
       CLOSE ( unit = unit, status = 'keep' )
    ENDIF

    DEALLOCATE ( rhog_g )
    DEALLOCATE ( g_g )

    RETURN

  END SUBROUTINE write_rhog

  !-------------------------------------------------------------------------------
  ! [ IMPORTANT ]
  ! When we calculate DFT using semicore states, but don't want to include them into GPP
  ! calculation of self-energy, we have to enter this subroutine to calculate charge density
  ! in reciprocal space contributed only be outer valence electrons
  ! ------
  ! Refer to PW/src/sum_band.f90: sum_band(), sum_band_k()
  ! ------
  SUBROUTINE calc_rhog (rhog_nvmin, rhog_nvmax)

    ! calc_rhog    Originally By Brad D. Malone    Last Modified (night before his thesis defense)
    ! Computes charge density by summing over a subset of occupied bands

    USE cell_base, ONLY : omega, tpiba2
    ! Modules/fft_base.f90:
    ! dfftp: descriptor for dense grid. Dimensions of the 3D real and reciprocal space FFT grid relative to the charge density and potential ("dense" grid)
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : fwfft, invfft
    USE gvecs, ONLY : doublegrid
    USE io_files, ONLY : nwordwfc, iunwfc
    USE klist, ONLY : xk, nkstot, ngk, nks, igk_k
    ! ======
    ! PW/src/pwcom.f90: isk(npk)
    ! Modules/parameters.f90: npk = 40000
    ! ------
    USE lsda_mod, ONLY : nspin, isk
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_pools, ONLY : inter_pool_comm, kunit, my_pool_id, npool
    USE mp_bands, ONLY : nbgrp, inter_bgrp_comm
    USE noncollin_module, ONLY : nspin_mag, noncolin, npol
    USE scf, ONLY : rho
    USE symme, ONLY : sym_rho, sym_rho_init, sym_rho_deallocate
    USE wavefunctions, ONLY : evc, psic, psic_nc
    USE wvfct, ONLY : wg, npwx, nbnd
    USE spin_orb, ONLY : domag
    USE lsda_mod, ONLY : lsda, nspin, current_spin, isk
    USE uspp,                 ONLY : nkb, vkb
    USE buffers,              ONLY : get_buffer

    IMPLICIT NONE

    integer, intent (in) :: rhog_nvmin
    integer, intent (in) :: rhog_nvmax
    integer, external :: global_kpoint_index
    integer :: npw, ik, is, ib, ig, ir, iks, ike, incr, j, ipol, ibnd
    REAL(DP) :: w1

    iks = global_kpoint_index (nkstot, 1)
    ike = iks + nks - 1

    ! ======
    ! Add warning of number of bands
    if (rhog_nvmin > rhog_nvmax) then
       if (ionode) then
          write(*,'(5X,A)') "[ERROR] rhog_nvmin > rhog_nvmax"
       endif
       call exit(233)
    elseif ( (rhog_nvmin <= 0) .or. (rhog_nvmax <= 0) .or. (rhog_nvmin > nbnd) .or. (rhog_nvmax > nbnd) ) then
       if (ionode) then
          write(*,'(5X,A,I5,A,I5,A,I5)') "[ERROR] invalid band range in (rhog_nvmin, rhog_nvmax) = (", rhog_nvmin, ",", rhog_nvmax, ") with total number of bands = ", nbnd
       endif
       ! call exit(234)
    endif
    ! ======
    ! PW/src/weights.f90:
    ! Calculates weights of Kohn-Sham orbitals used in calculation of rho
    CALL weights ()
    rho%of_r (:, :) = 0.0D0

    ! call MPI_ALLREDUCE(nks,nksmax,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,error)

    ! PW/src/sum_band.f90
    incr = 1
    DO ik = 1, nks
       IF ( lsda ) current_spin = isk(ik)
       npw = ngk (ik)
       ! nwordwfc = nbnd*npwx*npol
       ! 2 for complex
       ! call davcio (evc, 2*nwordwfc, iunwfc, ik_local, -1)
       IF ( nks > 1 ) THEN
          CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
       ENDIF
       IF ( nkb > 0 ) THEN
          CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
       ENDIF
       ! LOOP over selected bands
       DO ibnd = rhog_nvmin, rhog_nvmax, incr
          w1 = wg(ibnd,ik) / omega
          ! <- nspin = 4 ->
          ! PW/src/sum_band.f90: L612
          IF (noncolin) THEN
             psic_nc = (0.D0,0.D0)
             DO ig = 1, npw
                psic_nc(dffts%nl(igk_k(ig,ik)),1)=evc(ig     ,ibnd)
                psic_nc(dffts%nl(igk_k(ig,ik)),2)=evc(ig+npwx,ibnd)
             ENDDO
             CALL invfft ('Wave', psic_nc(:,1), dffts)
             CALL invfft ('Wave', psic_nc(:,2), dffts)
             ! increment the charge density ...
             DO ipol=1,npol
                CALL get_rho(rho%of_r(:,1), dffts%nnr, w1, psic_nc(:,ipol))
             ENDDO
             ! In this case, calculate also the three
             ! components of the magnetization (stored in rho%of_r(ir,2-4))
             IF (domag) THEN
                CALL get_rho_domag(rho%of_r(:,:), dffts%nnr, w1, psic_nc(:,:))
             ELSE
                rho%of_r(:,2:4)=0.0_DP
             ENDIF
             ! <- nspin = 1 or nspin = 2 ->
          ELSE
             !$omp parallel
             CALL threaded_barrier_memset(psic, 0.D0, dffts%nnr*2)
             !$omp do
             DO j = 1, npw
                psic(dffts%nl(igk_k(j,ik))) = evc(j,ibnd)
             ENDDO
             !$omp end do nowait
             !$omp end parallel
             CALL invfft ('Wave', psic, dffts)
             ! ... increment the charge density ...
             ! current_spin = isk(ik_global)
             CALL get_rho(rho%of_r(:,current_spin), dffts%nnr, w1, psic)
          endif ! noncollinear
       enddo ! ibnd = rhog_nvmin, rhog_nvmax
    enddo ! ik = 1, nks

    !! If doublegrid = T, dffts and dfftp differ
    !! Previous block calculates rho%of_r in dffts grid, which is usually coarser than dfftp grid, so we need to interpolate onto a denser dfftp grid
    IF (doublegrid) THEN
       IF (ionode) write(*,*) "[ERROR] doublegrid = T"
       CALL exit(1)
    ENDIF
    !! PW/src/sum_band.f90: L136
    !! rho%of_r(ir,is=1,nspin)
    CALL mp_sum (rho%of_r, inter_pool_comm)
    CALL mp_sum (rho%of_r, inter_bgrp_comm)
    IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
    rho%of_r = rho%of_r / nbgrp
    !! ... bring the unsymmetrized rho(r) to G-space (use psic as work array)
    DO is = 1, nspin
       psic(1:dffts%nnr) = rho%of_r(1:dffts%nnr,is)
       psic(dffts%nnr+1:) = 0.0_dp
       CALL fwfft ('Rho', psic, dffts)
       !! FFTXlib/fft_types.f90:
       !! ngm : my no. of non zero charge density/potential components
       rho%of_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
       rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
    ENDDO
    !! PW/src/symme.f90:
    !! SUBROUTINE sym_rho_init ( gamma_only )
    CALL sym_rho_init (.False.)
    CALL sym_rho (nspin_mag, rho%of_g)
    call sym_rho_deallocate ()

    RETURN

  END SUBROUTINE calc_rhog

  !-------------------------------------------------------------------------------

  SUBROUTINE get_rho(rho_loc, nrxxs_loc, w1_loc, psic_loc)

    IMPLICIT NONE

    INTEGER :: nrxxs_loc
    REAL(DP) :: rho_loc(nrxxs_loc)
    REAL(DP) :: w1_loc
    COMPLEX(DP) :: psic_loc(nrxxs_loc)

    INTEGER :: ir

    !$omp parallel do
    DO ir = 1, nrxxs_loc
       !
       rho_loc(ir) = rho_loc(ir) + &
            w1_loc * ( DBLE( psic_loc(ir) )**2 + &
            AIMAG( psic_loc(ir) )**2 )
       !
    END DO
    !$omp end parallel do

  END SUBROUTINE get_rho

  SUBROUTINE get_rho_gamma(rho_loc, nrxxs_loc, w1_loc, w2_loc, psic_loc)

    IMPLICIT NONE

    INTEGER :: nrxxs_loc
    REAL(DP) :: rho_loc(nrxxs_loc)
    REAL(DP) :: w1_loc, w2_loc
    COMPLEX(DP) :: psic_loc(nrxxs_loc)

    INTEGER :: ir

    !$omp parallel do
    DO ir = 1, nrxxs_loc
       !
       rho_loc(ir) = rho_loc(ir) + &
            w1_loc * DBLE( psic_loc(ir) )**2 + &
            w2_loc * AIMAG( psic_loc(ir) )**2
       !
    END DO
    !$omp end parallel do

  END SUBROUTINE get_rho_gamma

  SUBROUTINE get_rho_domag(rho_loc, nrxxs_loc, w1_loc, psic_loc)

    IMPLICIT NONE

    INTEGER :: nrxxs_loc
    REAL(DP) :: rho_loc(:, :)
    REAL(DP) :: w1_loc
    COMPLEX(DP) :: psic_loc(:, :)

    INTEGER :: ir

    !$omp parallel do
    DO ir = 1, nrxxs_loc
       !
       rho_loc(ir,2) = rho_loc(ir,2) + w1_loc*2.D0* &
            (DBLE(psic_loc(ir,1))* DBLE(psic_loc(ir,2)) + &
            AIMAG(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)))

       rho_loc(ir,3) = rho_loc(ir,3) + w1_loc*2.D0* &
            (DBLE(psic_loc(ir,1))*AIMAG(psic_loc(ir,2)) - &
            DBLE(psic_loc(ir,2))*AIMAG(psic_loc(ir,1)))

       rho_loc(ir,4) = rho_loc(ir,4) + w1_loc* &
            (DBLE(psic_loc(ir,1))**2+AIMAG(psic_loc(ir,1))**2 &
            -DBLE(psic_loc(ir,2))**2-AIMAG(psic_loc(ir,2))**2)
       !
    END DO
    !$omp end parallel do

  END SUBROUTINE get_rho_domag

  ! -------------------------------

  !   SUBROUTINE write_vxcg ( output_file_name, real_or_complex, symm_type, &
  !        vxc_zero_rho_core )

  !     USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  !     USE constants, ONLY : pi, tpi, eps6
  !     USE ener, ONLY : etxc, vtxc
  !     USE fft_base, ONLY : dfftp
  !     USE fft_interfaces, ONLY : fwfft
  !     USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  !     USE io_global, ONLY : ionode
  !     USE ions_base, ONLY : nat, atm, ityp, tau
  !     USE kinds, ONLY : DP
  !     USE lsda_mod, ONLY : nspin
  !     USE mp, ONLY : mp_sum
  !     ! USE mp_bands, ONLY : intra_bgrp_comm
  !     USE mp_pools, ONLY : intra_pool_comm
  !     USE scf, ONLY : rho, rho_core, rhog_core
  !     USE symm_base, ONLY : s, ftau, nsym
  !     USE wavefunctions, ONLY : psic
  !     USE matrix_inversion

  !     IMPLICIT NONE

  !     character ( len = 256 ), intent (in) :: output_file_name
  !     integer, intent (in) :: real_or_complex
  !     character ( len = 9 ), intent (in) :: symm_type
  !     logical, intent (in) :: vxc_zero_rho_core

  !     character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  !     integer :: unit, id, is, ir, ig, i, j, k, ierr
  !     integer :: nd, ns, nr, ng_l, ng_g
  !     integer :: ntran, cell_symmetry, nrecord
  !     real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  !     real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  !     real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  !     integer, allocatable :: g_g ( :, : )
  !     real (DP), allocatable :: vxcr_g ( :, : )
  !     complex (DP), allocatable :: vxcg_g ( :, : )

  !     INTEGER, EXTERNAL :: atomic_number

  !     CALL date_and_tim ( cdate, ctime )
  !     WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  !     WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  !     IF ( real_or_complex .EQ. 1 ) THEN
  !        WRITE ( stitle, '("VXC-Real",24X)' )
  !     ELSE
  !        WRITE ( stitle, '("VXC-Complex",21X)' )
  !     ENDIF

  !     unit = 4
  !     nrecord = 1
  !     nd = 3

  !     ns = nspin
  !     nr = dfftp%nnr
  !     ng_l = ngm
  !     ng_g = ngm_g

  !     ierr = 0
  !     IF ( ibrav .EQ. 0 ) THEN
  !        IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
  !           cell_symmetry = 0
  !        ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
  !           cell_symmetry = 1
  !        ELSE
  !           ierr = 1
  !        ENDIF
  !     ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
  !        cell_symmetry = 0
  !     ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
  !        cell_symmetry = 1
  !     ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
  !        cell_symmetry = 0
  !     ELSE
  !        ierr = 1
  !     ENDIF
  !     IF ( ierr .GT. 0 ) &
  !          CALL errore ( 'write_vxcg', 'cell_symmetry', ierr )

  !     ntran = nsym
  !     DO i = 1, ntran
  !        DO j = 1, nd
  !           DO k = 1, nd
  !              r1 ( k, j ) = dble ( s ( k, j, i ) )
  !           ENDDO
  !        ENDDO
  !        CALL invmat ( 3, r1, r2 )
  !        t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
  !        t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
  !        t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
  !        DO j = 1, nd
  !           t2 ( j ) = 0.0D0
  !           DO k = 1, nd
  !              t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
  !           ENDDO
  !           IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
  !                t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
  !           IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
  !                t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
  !        ENDDO
  !        DO j = 1, nd
  !           translation ( j, i ) = t2 ( j ) * tpi
  !        ENDDO
  !     ENDDO

  !     CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  !     alat2 = alat ** 2
  !     recvol = 8.0D0 * pi**3 / omega

  !     DO i = 1, nd
  !        DO j = 1, nd
  !           adot ( j, i ) = 0.0D0
  !        ENDDO
  !     ENDDO
  !     DO i = 1, nd
  !        DO j = 1, nd
  !           DO k = 1, nd
  !              adot ( j, i ) = adot ( j, i ) + &
  !                   at ( k, j ) * at ( k, i ) * alat2
  !           ENDDO
  !        ENDDO
  !     ENDDO

  !     DO i = 1, nd
  !        DO j = 1, nd
  !           bdot ( j, i ) = 0.0D0
  !        ENDDO
  !     ENDDO
  !     DO i = 1, nd
  !        DO j = 1, nd
  !           DO k = 1, nd
  !              bdot ( j, i ) = bdot ( j, i ) + &
  !                   bg ( k, j ) * bg ( k, i ) * tpiba2
  !           ENDDO
  !        ENDDO
  !     ENDDO

  !     ALLOCATE ( g_g ( nd, ng_g ) )
  !     ALLOCATE ( vxcr_g ( nr, ns ) )
  !     ALLOCATE ( vxcg_g ( ng_g, ns ) )

  !     DO ig = 1, ng_g
  !        DO id = 1, nd
  !           g_g ( id, ig ) = 0
  !        ENDDO
  !     ENDDO
  !     DO is = 1, ns
  !        DO ig = 1, ng_g
  !           vxcg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
  !        ENDDO
  !     ENDDO

  !     DO ig = 1, ng_l
  !        g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
  !        g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
  !        g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  !     ENDDO
  !     vxcr_g ( :, : ) = 0.0D0
  !     IF ( vxc_zero_rho_core ) THEN
  !        rho_core ( : ) = 0.0D0
  !        rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  !     ENDIF
  !     CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !     DO is = 1, ns
  !        DO ir = 1, nr
  !           psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
  !        ENDDO
  !        CALL fwfft ( 'Rho', psic, dfftp )
  !        DO ig = 1, ng_l
  !           vxcg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
  !        ENDDO
  !     ENDDO

  !     ! ======
  !     !! intra_pool_comm = intra_bgrp_comm
  !     ! CALL mp_sum ( g_g, intra_bgrp_comm )
  !     ! CALL mp_sum ( vxcg_g, intra_bgrp_comm )
  !     CALL mp_sum ( g_g, intra_pool_comm )
  !     CALL mp_sum ( vxcg_g, intra_pool_comm )

  !     IF ( ionode ) THEN
  !        OPEN ( unit = unit, file = TRIM ( output_file_name ), &
  !             form = 'unformatted', status = 'replace' )
  !        WRITE ( unit ) stitle, sdate, stime
  !        WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
  !        WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
  !        WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
  !             ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
  !        WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
  !             ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
  !        WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
  !        WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
  !        WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
  !        WRITE ( unit ) nrecord
  !        WRITE ( unit ) ng_g
  !        WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  !        WRITE ( unit ) nrecord
  !        WRITE ( unit ) ng_g
  !        IF ( real_or_complex .EQ. 1 ) THEN
  !           WRITE ( unit ) ( ( dble ( vxcg_g ( ig, is ) ), &
  !                ig = 1, ng_g ), is = 1, ns )
  !        ELSE
  !           WRITE ( unit ) ( ( vxcg_g ( ig, is ), &
  !                ig = 1, ng_g ), is = 1, ns )
  !        ENDIF
  !        CLOSE ( unit = unit, status = 'keep' )
  !     ENDIF

  !     DEALLOCATE ( vxcg_g )
  !     DEALLOCATE ( vxcr_g )
  !     DEALLOCATE ( g_g )

  !     RETURN

  !   END SUBROUTINE write_vxcg

  !   !-------------------------------------------------------------------------------

  !   SUBROUTINE write_vxc0 ( output_file_name, vxc_zero_rho_core )

  !     USE constants, ONLY : RYTOEV
  !     USE ener, ONLY : etxc, vtxc
  !     USE fft_base, ONLY : dfftp
  !     USE fft_interfaces, ONLY : fwfft
  !     USE gvect, ONLY : ngm, mill
  !     USE io_global, ONLY : ionode
  !     USE kinds, ONLY : DP
  !     USE lsda_mod, ONLY : nspin
  !     USE mp, ONLY : mp_sum
  !     ! USE mp_bands, ONLY : intra_bgrp_comm
  !     USE mp_pools, ONLY : intra_pool_comm
  !     USE scf, ONLY : rho, rho_core, rhog_core
  !     USE wavefunctions, ONLY : psic

  !     IMPLICIT NONE

  !     character ( len = 256 ), intent (in) :: output_file_name
  !     logical, intent (in) :: vxc_zero_rho_core

  !     integer :: unit
  !     integer :: is, ir, ig
  !     integer :: nd, ns, nr, ng_l
  !     real (DP), allocatable :: vxcr_g ( :, : )
  !     complex (DP), allocatable :: vxc0_g ( : )

  !     unit = 4
  !     nd = 3

  !     ns = nspin
  !     nr = dfftp%nnr
  !     ng_l = ngm

  !     ALLOCATE ( vxcr_g ( nr, ns ) )
  !     ALLOCATE ( vxc0_g ( ns ) )

  !     DO is = 1, ns
  !        vxc0_g ( is ) = ( 0.0D0, 0.0D0 )
  !     ENDDO

  !     vxcr_g ( :, : ) = 0.0D0
  !     IF ( vxc_zero_rho_core ) THEN
  !        rho_core ( : ) = 0.0D0
  !        rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  !     ENDIF
  !     CALL v_xc ( rho, rho_core, rhog_core, etxc, vtxc, vxcr_g )
  !     DO is = 1, ns
  !        DO ir = 1, nr
  !           psic ( ir ) = CMPLX ( vxcr_g ( ir, is ), 0.0D0, KIND=dp )
  !        ENDDO
  !        CALL fwfft ( 'Rho', psic, dfftp )
  !        DO ig = 1, ng_l
  !           IF ( mill ( 1, ig ) .EQ. 0 .AND. mill ( 2, ig ) .EQ. 0 .AND. &
  !                mill ( 3, ig ) .EQ. 0 ) vxc0_g ( is ) = psic ( dfftp%nl ( ig ) )
  !        ENDDO
  !     ENDDO

  !     ! ======
  !     !! intra_bgrp_comm = intra_pool_comm
  !     ! CALL mp_sum ( vxc0_g, intra_bgrp_comm )
  !     CALL mp_sum ( vxc0_g, intra_pool_comm )

  !     DO is = 1, ns
  !        vxc0_g ( is ) = vxc0_g ( is ) * CMPLX ( RYTOEV, 0.0D0, KIND=dp )
  !     ENDDO

  !     IF ( ionode ) THEN
  !        OPEN (unit = unit, file = TRIM (output_file_name), &
  !             form = 'formatted', status = 'replace')
  !        WRITE ( unit, 101 )
  !        DO is = 1, ns
  !           WRITE ( unit, 102 ) is, vxc0_g ( is )
  !        ENDDO
  !        WRITE ( unit, 103 )
  !        CLOSE (unit = unit, status = 'keep')
  !     ENDIF

  !     DEALLOCATE ( vxcr_g )
  !     DEALLOCATE ( vxc0_g )

  !     RETURN

  ! 101 FORMAT ( /, 5X, "--------------------------------------------", &
  !          /, 5X, "spin    Re Vxc(G=0) (eV)    Im Vxc(G=0) (eV)", &
  !          /, 5X, "--------------------------------------------" )
  ! 102 FORMAT ( 5X, I1, 3X, 2F20.15 )
  ! 103 FORMAT ( 5X, "--------------------------------------------", / )

  !   END SUBROUTINE write_vxc0

  !   !-------------------------------------------------------------------------------

  !   SUBROUTINE write_vxc_r (output_file_name, diag_nmin, diag_nmax, &
  !        offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

  !     USE kinds, ONLY : DP
  !     USE constants, ONLY : rytoev
  !     USE cell_base, ONLY : tpiba2, at, bg
  !     USE ener, ONLY : etxc, vtxc
  !     USE fft_base, ONLY : dfftp
  !     USE fft_interfaces, ONLY : invfft
  !     USE gvect, ONLY : ngm, g
  !     USE io_files, ONLY : nwordwfc, iunwfc
  !     USE io_global, ONLY : ionode
  !     USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
  !     USE lsda_mod, ONLY : nspin, isk
  !     USE mp, ONLY : mp_sum
  !     USE mp_pools, ONLY : inter_pool_comm
  !     ! USE mp_bands, ONLY : intra_bgrp_comm
  !     USE mp_pools, ONLY : intra_pool_comm
  !     USE scf, ONLY : rho, rho_core, rhog_core
  !     USE wavefunctions, ONLY : evc, psic
  !     USE wvfct, ONLY : nbnd

  !     IMPLICIT NONE

  !     character (len = 256), intent (in) :: output_file_name
  !     integer, intent (inout) :: diag_nmin
  !     integer, intent (inout) :: diag_nmax
  !     integer, intent (inout) :: offdiag_nmin
  !     integer, intent (inout) :: offdiag_nmax
  !     logical, intent (in) :: vxc_zero_rho_core

  !     integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2
  !     integer, external :: global_kpoint_index
  !     real (DP) :: dummyr
  !     complex (DP) :: dummyc
  !     real (DP), allocatable :: mtxeld (:, :)
  !     complex (DP), allocatable :: mtxelo (:, :, :)
  !     real (DP), allocatable :: vxcr (:, :)
  !     complex (DP), allocatable :: psic2 (:)

  !     if(diag_nmin > diag_nmax) then
  !        call errore ( 'write_vxc_r', 'diag_nmin > diag_nmax', diag_nmin )
  !     endif
  !     IF (diag_nmin .LT. 1) diag_nmin = 1
  !     IF (diag_nmax .GT. nbnd) then
  !        if (ionode) then
  !           write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
  !        endif
  !        diag_nmax = nbnd
  !     ENDIF
  !     ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

  !     if(offdiag_nmin > offdiag_nmax) then
  !        if (ionode) then
  !           call errore ( 'write_vxc_r', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
  !        endif
  !     endif
  !     IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
  !     IF (offdiag_nmax .GT. nbnd)  then
  !        if (ionode) then
  !           write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
  !        endif
  !        offdiag_nmax = nbnd
  !     ENDIF
  !     noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

  !     IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

  !     unit = 4

  !     iks = global_kpoint_index (nkstot, 1)
  !     ike = iks + nks - 1

  !     IF (ndiag .GT. 0) THEN
  !        ALLOCATE (mtxeld (ndiag, nkstot))
  !        mtxeld (:, :) = 0.0D0
  !     ENDIF
  !     IF (noffdiag .GT. 0) THEN
  !        ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
  !        mtxelo (:, :, :) = (0.0D0, 0.0D0)
  !     ENDIF

  !     ALLOCATE (vxcr (dfftp%nnr, nspin))
  !     IF (noffdiag .GT. 0) ALLOCATE (psic2 (dfftp%nnr))

  !     vxcr (:, :) = 0.0D0
  !     IF ( vxc_zero_rho_core ) THEN
  !        rho_core ( : ) = 0.0D0
  !        rhog_core ( : ) = ( 0.0D0, 0.0D0 )
  !     ENDIF
  !     CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)

  !     DO ik = iks, ike
  !        npw = ngk ( ik - iks + 1 )
  !        CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
  !        IF (ndiag .GT. 0) THEN
  !           DO ib = diag_nmin, diag_nmax
  !              psic (:) = (0.0D0, 0.0D0)
  !              DO ig = 1, npw
  !                 psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
  !              ENDDO
  !              CALL invfft ('Rho', psic, dfftp)
  !              dummyr = 0.0D0
  !              DO ir = 1, dfftp%nnr
  !                 dummyr = dummyr + vxcr (ir, isk (ik)) &
  !                      * (dble (psic (ir)) **2 + aimag (psic (ir)) **2)
  !              ENDDO
  !              dummyr = dummyr * rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x)
  !              ! ======
  !              !! intra_pool_comm = intra_bgrp_comm
  !              ! CALL mp_sum (dummyr, intra_bgrp_comm)
  !              CALL mp_sum (dummyr, intra_pool_comm)
  !              mtxeld (ib - diag_nmin + 1, ik) = dummyr
  !           ENDDO
  !        ENDIF
  !        IF (noffdiag .GT. 0) THEN
  !           DO ib = offdiag_nmin, offdiag_nmax
  !              psic (:) = (0.0D0, 0.0D0)
  !              DO ig = 1, npw
  !                 psic (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib)
  !              ENDDO
  !              CALL invfft ('Rho', psic, dfftp)
  !              DO ib2 = offdiag_nmin, offdiag_nmax
  !                 psic2 (:) = (0.0D0, 0.0D0)
  !                 DO ig = 1, npw
  !                    psic2 (dfftp%nl (igk_k (ig,ik-iks+1))) = evc (ig, ib2)
  !                 ENDDO
  !                 CALL invfft ('Rho', psic2, dfftp)
  !                 dummyc = (0.0D0, 0.0D0)
  !                 DO ir = 1, dfftp%nnr
  !                    dummyc = dummyc + CMPLX (vxcr (ir, isk (ik)), 0.0D0, KIND=dp) &
  !                         * conjg (psic2 (ir)) * psic (ir)
  !                 ENDDO
  !                 dummyc = dummyc &
  !                      * CMPLX (rytoev / dble (dfftp%nr1x * dfftp%nr2x * dfftp%nr3x), &
  !                      0.0D0, KIND=dp)
  !                 ! ======
  !                 !! intra_pool_comm = intra_bgrp_comm
  !                 ! CALL mp_sum (dummyc, intra_bgrp_comm)
  !                 CALL mp_sum (dummyc, intra_pool_comm)
  !                 mtxelo (ib2 - offdiag_nmin + 1, ib - offdiag_nmin &
  !                      + 1, ik) = dummyc
  !              ENDDO
  !           ENDDO
  !        ENDIF
  !     ENDDO

  !     DEALLOCATE (vxcr)
  !     IF (noffdiag .GT. 0) DEALLOCATE (psic2)

  !     IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
  !     IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)

  !     CALL cryst_to_cart (nkstot, xk, at, -1)

  !     IF (ionode) THEN
  !        OPEN (unit = unit, file = TRIM (output_file_name), &
  !             form = 'formatted', status = 'replace')
  !        DO ik = 1, nkstot / nspin
  !           WRITE (unit, 101) xk(:, ik), nspin * ndiag, &
  !                nspin * noffdiag **2
  !           DO is = 1, nspin
  !              IF (ndiag .GT. 0) THEN
  !                 DO ib = diag_nmin, diag_nmax
  !                    WRITE (unit, 102) is, ib, mtxeld &
  !                         (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin), &
  !                         0.0D0
  !                 ENDDO
  !              ENDIF
  !              IF (noffdiag .GT. 0) THEN
  !                 DO ib = offdiag_nmin, offdiag_nmax
  !                    DO ib2 = offdiag_nmin, offdiag_nmax
  !                       WRITE (unit, 103) is, ib2, ib, mtxelo &
  !                            (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
  !                            ik + (is - 1) * nkstot / nspin)
  !                    ENDDO
  !                 ENDDO
  !              ENDIF
  !           ENDDO
  !        ENDDO
  !        CLOSE (unit = unit, status = 'keep')
  !     ENDIF

  !     CALL cryst_to_cart (nkstot, xk, bg, 1)

  !     IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
  !     IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

  !     RETURN

  ! 101 FORMAT (3F13.9, 2I8)
  ! 102 FORMAT (2I8, 2F15.9)
  ! 103 FORMAT (3I8, 2F15.9)

  !   END SUBROUTINE write_vxc_r

  !-------------------------------------------------------------------------------

  SUBROUTINE write_vxc_g (output_file_name, diag_nmin, diag_nmax, offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)

    USE constants, ONLY : rytoev
    USE cell_base, ONLY : tpiba2, at, bg
    USE ener, ONLY : etxc, vtxc
    USE exx, ONLY : vexx
    USE fft_base, ONLY : dfftp
    USE fft_interfaces, ONLY : fwfft, invfft
    USE funct, ONLY : exx_is_active
    USE gvect, ONLY : ngm, g
    USE io_files, ONLY : nwordwfc, iunwfc
    USE io_global, ONLY : ionode
    USE kinds, ONLY : DP
    USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
    USE lsda_mod, ONLY : nspin, isk, current_spin, lsda
    USE mp, ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp_pools, ONLY : intra_pool_comm
    USE scf, ONLY : rho, rho_core, rhog_core
    USE wavefunctions, ONLY : evc, psic
    USE wvfct, ONLY : npwx, nbnd, current_k
    USE noncollin_module , ONLY : noncolin , npol
    USE buffers, ONLY : get_buffer

    IMPLICIT NONE

    character (len = 256), intent (in) :: output_file_name
    integer, intent (inout) :: diag_nmin
    integer, intent (inout) :: diag_nmax
    integer, intent (inout) :: offdiag_nmin
    integer, intent (inout) :: offdiag_nmax
    logical, intent (in) :: vxc_zero_rho_core

    integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
    integer, external :: global_kpoint_index
    complex (DP) :: dummy
    complex (DP), allocatable :: mtxeld (:, :)
    complex (DP), allocatable :: mtxelo (:, :, :)
    real (DP), allocatable :: vxcr (:, :)
    ! ======
    complex (DP) :: hpsi(npwx*npol,nbnd)
    complex (DP), allocatable :: hc(:,:)
    integer :: nspin_
    integer :: kdim, kdmx
    ! integer :: lda, n, m

    ! ------
    if(diag_nmin > diag_nmax) then
       call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
    endif
    IF (diag_nmin .LT. 1) diag_nmin = 1
    IF (diag_nmax .GT. nbnd) then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
       endif
       diag_nmax = nbnd
    ENDIF
    ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

    if(offdiag_nmin > offdiag_nmax) then
       if (ionode) then
          call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
       endif
    endif
    IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
    IF (offdiag_nmax .GT. nbnd)  then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
       endif
       offdiag_nmax = nbnd
    ENDIF
    noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

    IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

    ! ======
    ALLOCATE( hc( nbnd, nbnd) )
    ! ------
    unit = 4

    iks = global_kpoint_index (nkstot, 1)
    ike = iks + nks - 1

    IF (ndiag .GT. 0) THEN
       ALLOCATE (mtxeld (ndiag, nkstot))
       mtxeld (:, :) = (0.0D0, 0.0D0)
    ENDIF
    IF (noffdiag .GT. 0) THEN
       ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
       mtxelo (:, :, :) = (0.0D0, 0.0D0)
    ENDIF

    ALLOCATE (vxcr (dfftp%nnr, nspin))

    vxcr (:, :) = 0.0D0
    IF ( vxc_zero_rho_core ) THEN
       rho_core ( : ) = 0.0D0
       rhog_core ( : ) = ( 0.0D0, 0.0D0 )
    ENDIF
    CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)

    DO ik = iks, ike

       hpsi(:,:) = (0.0D0, 0.0D0)
       evc(:,:) = (0.0D0, 0.0D0)
       hc(:,:) = (0.0D0,0.0D0)
       !! ikk = ik_local
       ikk = ik - iks + 1
       npw = ngk ( ik - iks + 1 )

       IF ( npol == 1 ) THEN
          kdim = npw
          kdmx = npwx
       ELSE
          kdim = npwx*npol
          kdmx = npwx*npol
       ENDIF

       IF ( lsda ) current_spin = isk(ikk)
       ! ! current_k will be used in vloc_psi_nc and vloc_psi_k
       current_k = ikk

       CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
       ! h_psi.f90
       ! PW/src/vloc_psi.f90: we add two more subroutines, vloc_psi_nc_2 and vloc_psi_k_2
       !> lda = npwx, n = npw, m = nbnd, psi = evc
       IF ( noncolin ) THEN
          CALL vloc_psi_nc( npwx, npw, nbnd, evc, vxcr, hpsi)
       ELSE
          CALL vloc_psi_k( npwx, npw, nbnd, evc, vxcr(1, current_spin), hpsi)
       ENDIF

       ! IF ( exx_is_active() ) then
       !    if (ionode) then
       !       call errore ( 'write_vxc_g', 'exx not supported', 1 )
       !    endif
       ! ENDIF

       ! lda = npwx
       ! n = npw
       ! m = nbnd
       ! !> [WORKING]
       ! if (dft_is_meta()) call h_psi_meta (lda, n, m, evc, hpsi)
       ! !
       ! ! ... Here we add the Hubbard potential times psi
       ! !
       ! IF ( lda_plus_u .AND. U_projection.NE."pseudo" ) THEN
       !    !
       !    IF (noncolin) THEN
       !       CALL vhpsi_nc( lda, n, m, evc, hpsi )
       !    ELSE
       !       call vhpsi( lda, n, m, evc, hpsi )
       !    ENDIF
       !    !
       ! ENDIF

       ! !> [WORKING]
       ! ! IF ( exx_is_active() ) CALL vexx( npwx, npw, nbnd, evc, hpsi)
       ! IF ( exx_is_active() ) THEN
       !    IF ( use_ace) THEN
       !       IF (gamma_only) THEN
       !          CALL vexxace_gamma(lda, m, psi, ee, hpsi)
       !       ELSE
       !          CALL vexxace_k(lda, m, psi, ee, hpsi)
       !       END IF
       !    ELSE
       !       CALL vexx( lda, n, m, psi, hpsi, becp )
       !    END IF
       ! END IF

       ! ======
       ! C := alpha*op( A )*op( B ) + beta*C,
       ! zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
       ! hc = (evc).H . hpsi

       CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ( 1.D0, 0.D0 ), evc, kdmx, hpsi, kdmx, ( 0.D0, 0.D0 ), hc, nbnd )

       CALL mp_sum(  hc , intra_bgrp_comm )
       ! CALL mp_sum(  hc , intra_pool_comm )

       ! diagonal
       if (ndiag .gt. 0) then
          do ib = diag_nmin, diag_nmax
             ! mtxeld in units of eV
             mtxeld (ib - diag_nmin + 1, ik) = hc (ib, ib) * CMPLX (rytoev, 0.0D0)
          enddo
       endif

       ! off-diagonal
       if (noffdiag .gt. 0) then
          do ib = offdiag_nmin, offdiag_nmax
             do ib2 = offdiag_nmin, offdiag_nmax
                ! mtxeld in units of eV
                mtxelo (ib2-offdiag_nmin+1, ib-offdiag_nmin+1, ik) = hc (ib2, ib) * CMPLX (rytoev, 0.0D0)
             enddo
          enddo
       endif

    ENDDO ! ik_global

    deallocate (hc)
    DEALLOCATE (vxcr)

    IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
    IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)
    ! xk(1:3) : cartesian ==> crytal
    CALL cryst_to_cart (nkstot, xk, at, -1)

    IF (ionode) THEN
       OPEN (unit = unit, file = TRIM (output_file_name), &
            form = 'formatted', status = 'replace')

       if (nspin .eq. 4) then
          nspin_ = 1
       else
          nspin_ = nspin
       endif

       DO ik = 1, nkstot / nspin_
          WRITE (unit, 101) xk(:, ik), nspin_ * ndiag, &
               nspin_ * noffdiag **2
          DO is = 1, nspin_
             IF (ndiag .GT. 0) THEN
                DO ib = diag_nmin, diag_nmax
                   WRITE (unit, 102) is, ib, mtxeld &
                        (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin_)
                ENDDO
             ENDIF
             IF (noffdiag .GT. 0) THEN
                DO ib = offdiag_nmin, offdiag_nmax
                   DO ib2 = offdiag_nmin, offdiag_nmax
                      WRITE (unit, 103) is, ib2, ib, mtxelo &
                           (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                           ik + (is - 1) * nkstot / nspin_)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       CLOSE (unit = unit, status = 'keep')
    ENDIF

    ! xk(1:3) : crystal ==> cartesian
    CALL cryst_to_cart (nkstot, xk, bg, 1)

    IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
    IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

    RETURN

101 FORMAT (3F13.9, 2I8)
102 FORMAT (2I8, 2F15.9)
103 FORMAT (3I8, 2F15.9)

  END SUBROUTINE write_vxc_g

  !> This subroutine put vxc + vhub + vmetagga + vexx into vxc.dat
  SUBROUTINE write_vxc_tot_g (output_file_name, diag_nmin, diag_nmax, offdiag_nmin, offdiag_nmax, vxc_zero_rho_core)
    USE constants, ONLY : rytoev
    USE cell_base, ONLY : tpiba2, at, bg
    USE ener, ONLY : etxc, vtxc
    USE exx,      ONLY : use_ace, vexx, vexxace_gamma, vexxace_k
    USE fft_base, ONLY : dfftp
    USE fft_interfaces, ONLY : fwfft, invfft
    USE funct, ONLY : exx_is_active, dft_is_meta, get_meta, dft_is_hybrid
    USE gvect, ONLY : ngm, g
    USE io_files, ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
    USE io_global, ONLY : ionode
    USE kinds, ONLY : DP
    USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
    USE lsda_mod, ONLY : nspin, isk, current_spin, lsda
    USE mp, ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp_pools, ONLY : intra_pool_comm
    USE scf, ONLY : rho, rho_core, rhog_core
    USE wavefunctions, ONLY : evc, psic
    USE wvfct, ONLY : npwx, nbnd, current_k
    USE noncollin_module , ONLY : noncolin , npol
    USE realus,   ONLY : real_space
    USE buffers, ONLY : get_buffer
    USE ldaU, ONLY : lda_plus_u, U_projection, wfcU, offsetU, Hubbard_l, nwfcU, Hubbard_lmax, is_hubbard
    USE uspp,     ONLY : vkb, nkb
    USE becmod,   ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
    USE control_flags,    ONLY : gamma_only

    IMPLICIT NONE

    character (len = 256), intent (in) :: output_file_name
    integer, intent (inout) :: diag_nmin
    integer, intent (inout) :: diag_nmax
    integer, intent (inout) :: offdiag_nmin
    integer, intent (inout) :: offdiag_nmax
    logical, intent (in) :: vxc_zero_rho_core

    integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
    integer, external :: global_kpoint_index
    complex (DP) :: dummy
    complex (DP), allocatable :: mtxeld (:, :)
    complex (DP), allocatable :: mtxelo (:, :, :)
    real (DP), allocatable :: vxcr (:, :)
    real (DP), allocatable :: kedtaur (:, :)   !FZ: for metaGGA
    
    complex (DP) :: hpsi(npwx*npol,nbnd)
    complex (DP), allocatable :: hc(:,:)
    integer :: nspin_
    integer :: kdim, kdmx
    integer :: lda, n, m
    real (DP) :: ee

    if (diag_nmin > diag_nmax) then
       call errore ( 'write_vxc_tot_g', 'diag_nmin > diag_nmax', diag_nmin )
    endif
    IF (diag_nmin .LT. 1) diag_nmin = 1
    IF (diag_nmax .GT. nbnd) then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
       endif
       diag_nmax = nbnd
    ENDIF
    ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

    if (offdiag_nmin > offdiag_nmax) then
       if (ionode) then
          call errore ( 'write_vxc_tot_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
       endif
    endif
    IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
    IF (offdiag_nmax .GT. nbnd)  then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
       endif
       offdiag_nmax = nbnd
    ENDIF
    noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

    IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

    ALLOCATE( hc( nbnd, nbnd) )
    unit = 4

    iks = global_kpoint_index (nkstot, 1)
    ike = iks + nks - 1

    IF (ndiag .GT. 0) THEN
       ALLOCATE (mtxeld (ndiag, nkstot))
       mtxeld (:, :) = (0.0D0, 0.0D0)
    ENDIF
    IF (noffdiag .GT. 0) THEN
       ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
       mtxelo (:, :, :) = (0.0D0, 0.0D0)
    ENDIF

    ALLOCATE (vxcr (dfftp%nnr, nspin))
    ALLOCATE (kedtaur (dfftp%nnr, nspin))   !FZ: for meta GGA

    vxcr (:, :) = 0.0D0
    kedtaur (:, :) = 0.0D0    !FZ: for metaGGA

    IF ( vxc_zero_rho_core ) THEN
       rho_core ( : ) = 0.0D0
       rhog_core ( : ) = ( 0.0D0, 0.0D0 )
    ENDIF

    IF (dft_is_meta() .and. (get_meta() /= 4)) then         !FZ: for metaGGA
       CALL v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, vxcr, kedtaur )    !FZ: for metaGGA
    ELSE
       CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxcr)
    ENDIF

    write(*,*) "use_ace = ", use_ace

    ! IF (nkb > 0) THEN
    !    CALL allocate_bec_type (nkb, nbnd, becp, intra_bgrp_comm)
    ! ENDIF

    DO ik = iks, ike

       hpsi(:,:) = (0.0D0, 0.0D0)
       evc(:,:) = (0.0D0, 0.0D0)
       hc(:,:) = (0.0D0,0.0D0)
       !! ikk = ik_local
       ikk = ik - iks + 1
       npw = ngk ( ik - iks + 1 )

       IF ( npol == 1 ) THEN
          kdim = npw
          kdmx = npwx
       ELSE
          kdim = npwx*npol
          kdmx = npwx*npol
       ENDIF

       IF ( lsda ) current_spin = isk(ikk)
       ! ! current_k will be used in vloc_psi_nc and vloc_psi_k
       current_k = ikk

       CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)

       ! IF ( nkb > 0 ) THEN
       !    CALL init_us_2( npw, igk_k(1,ik-iks+1), xk(1,ik-iks+1), vkb )
       ! ENDIF
       ! write(*,*) "SIZE(vkb,1) = ", SIZE(vkb, 1), " SIZE(evc,1) = ", SIZE(evc,1), " SIZE(becp, 1) = ", SIZE(becp, 1)

       IF ( lda_plus_u .AND. (U_projection .NE. 'pseudo') ) THEN
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ikk )
       ENDIF

       !> h_psi.f90
       !> PW/src/vloc_psi.f90: we add two more subroutines, vloc_psi_nc_2 and vloc_psi_k_2
       !> lda = npwx, n = npw, m = nbnd, psi = evc
       IF ( noncolin ) THEN
          CALL vloc_psi_nc( npwx, npw, nbnd, evc, vxcr, hpsi)
       ELSE
          CALL vloc_psi_k( npwx, npw, nbnd, evc, vxcr(1, current_spin), hpsi)
       ENDIF

       lda = npwx
       n = npw
       m = nbnd
       ! IF ( nkb > 0 .AND. .NOT. real_space) THEN
       !    ! CALL start_clock( 'h_psi:calbec' )
       !    CALL calbec ( n, vkb, evc, becp, m )
       !    ! CALL stop_clock( 'h_psi:calbec' )
       !    ! CALL add_vuspsi( lda, n, m, hpsi )
       ! END IF

       !> PW/src/h_psi_meta.f90
       ! This routine computes the specific contribution from the meta-GGA
       ! potential to H*psi; the result is added to hpsi
       if (dft_is_meta()) then
          write(*,*) "AAA"
          call h_psi_meta (lda, n, m, evc, hpsi)
       endif
       !
       ! ... Here we add the Hubbard potential times psi = evc
       !
       IF ( lda_plus_u .AND. U_projection.NE."pseudo" ) THEN
          write(*,*) "BBB"
          IF (noncolin) THEN
             CALL vhpsi_nc( lda, n, m, evc, hpsi )
          ELSE
             call vhpsi( lda, n, m, evc, hpsi )
          ENDIF
       ENDIF

       !> PW/src/exx.f90
       IF ( exx_is_active() ) THEN
          write(*,*) "CCC"
          IF ( use_ace) THEN
             IF (gamma_only) THEN
                CALL vexxace_gamma(lda, m, evc, ee, hpsi)
             ELSE
                write(*,*) "DDD"
                !> [BUGS HERE]
                CALL vexxace_k(lda, m, evc, ee, hpsi)
             END IF
          ELSE
             if (ionode) then
                call errore ( 'write_vxc_tot_g', 'Only support ACE', 1 )
             endif
             CALL vexx( lda, n, m, evc, hpsi, becp )
          END IF
       END IF

       ! ======
       ! C := alpha*op( A )*op( B ) + beta*C,
       ! zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
       ! hc = (evc).H . hpsi
       CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ( 1.D0, 0.D0 ), evc, kdmx, hpsi, kdmx, ( 0.D0, 0.D0 ), hc, nbnd )

       CALL mp_sum(  hc , intra_bgrp_comm )
       ! CALL mp_sum(  hc , intra_pool_comm )

       ! diagonal
       if (ndiag .gt. 0) then
          do ib = diag_nmin, diag_nmax
             ! mtxeld in units of eV
             mtxeld (ib - diag_nmin + 1, ik) = hc (ib, ib) * CMPLX (rytoev, 0.0D0)
          enddo
       endif

       ! off-diagonal
       if (noffdiag .gt. 0) then
          do ib = offdiag_nmin, offdiag_nmax
             do ib2 = offdiag_nmin, offdiag_nmax
                ! mtxeld in units of eV
                mtxelo (ib2-offdiag_nmin+1, ib-offdiag_nmin+1, ik) = hc (ib2, ib) * CMPLX (rytoev, 0.0D0)
             enddo
          enddo
       endif

    ENDDO ! ik_global

    ! IF (nkb > 0) THEN
    !    CALL deallocate_bec_type ( becp )
    ! ENDIF

    deallocate (hc)
    DEALLOCATE (vxcr)
    DEALLOCATE (kedtaur)

    IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
    IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)
    ! xk(1:3) : cartesian ==> crytal
    CALL cryst_to_cart (nkstot, xk, at, -1)

    IF (ionode) THEN
       OPEN (unit = unit, file = TRIM (output_file_name), &
            form = 'formatted', status = 'replace')

       if (nspin .eq. 4) then
          nspin_ = 1
       else
          nspin_ = nspin
       endif

       DO ik = 1, nkstot / nspin_
          WRITE (unit, 101) xk(:, ik), nspin_ * ndiag, &
               nspin_ * noffdiag **2
          DO is = 1, nspin_
             IF (ndiag .GT. 0) THEN
                DO ib = diag_nmin, diag_nmax
                   WRITE (unit, 102) is, ib, mtxeld &
                        (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin_)
                ENDDO
             ENDIF
             IF (noffdiag .GT. 0) THEN
                DO ib = offdiag_nmin, offdiag_nmax
                   DO ib2 = offdiag_nmin, offdiag_nmax
                      WRITE (unit, 103) is, ib2, ib, mtxelo &
                           (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                           ik + (is - 1) * nkstot / nspin_)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       CLOSE (unit = unit, status = 'keep')
    ENDIF

    ! xk(1:3) : crystal ==> cartesian
    CALL cryst_to_cart (nkstot, xk, bg, 1)


    IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
    IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

    RETURN

101 FORMAT (3F13.9, 2I8)
102 FORMAT (2I8, 2F15.9)
103 FORMAT (3I8, 2F15.9)

  END SUBROUTINE write_vxc_tot_g

  !-------------------------------------------------------------------------------

  SUBROUTINE write_vhub_g (output_file_name, diag_nmin, diag_nmax, offdiag_nmin, offdiag_nmax)

    USE constants, ONLY : rytoev
    USE cell_base, ONLY : tpiba2, at, bg
    USE ener, ONLY : etxc, vtxc
    USE exx, ONLY : vexx
    USE fft_base, ONLY : dfftp
    USE fft_interfaces, ONLY : fwfft, invfft
    USE funct, ONLY : exx_is_active
    USE gvect, ONLY : ngm, g
    USE io_files, ONLY : nwordwfc, iunwfc, nwordwfcU, iunhub
    USE io_global, ONLY : ionode
    USE kinds, ONLY : DP
    USE klist, ONLY : xk, nkstot, nks, ngk, igk_k
    USE lsda_mod, ONLY : nspin, isk, current_spin, lsda
    USE ldaU, ONLY : lda_plus_u, U_projection, wfcU, offsetU, Hubbard_l, nwfcU, Hubbard_lmax, is_hubbard
    USE mp, ONLY : mp_sum
    USE mp_pools, ONLY : inter_pool_comm
    USE mp_bands, ONLY : intra_bgrp_comm
    USE mp_pools, ONLY : intra_pool_comm
    USE scf, ONLY : rho, rho_core, rhog_core
    USE wavefunctions, ONLY : evc, psic
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module , ONLY : noncolin , npol
    USE buffers, ONLY : get_buffer
    USE scf, ONLY : v
    USE ions_base,            ONLY : nat, ityp
    ! USE becmod,               ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type

    IMPLICIT NONE

    character (len = 256), intent (in) :: output_file_name
    integer, intent (inout) :: diag_nmin
    integer, intent (inout) :: diag_nmax
    integer, intent (inout) :: offdiag_nmin
    integer, intent (inout) :: offdiag_nmax

    integer :: npw, ik, is, ib, ig, ir, unit, iks, ike, ndiag, noffdiag, ib2, ikk
    integer, external :: global_kpoint_index
    complex (DP) :: dummy
    complex (DP), allocatable :: mtxeld (:, :)
    complex (DP), allocatable :: mtxelo (:, :, :)
    complex (DP) :: hpsi(npwx*npol,nbnd)
    complex (DP), allocatable :: hc(:,:)
    integer :: nspin_
    integer :: kdim, kdmx
    COMPLEX (DP) :: zdotc

    ! ======
    ! COMPLEX(DP) , ALLOCATABLE :: proj(:,:)
    integer :: ldim, is1, ibnd, i, na, m1, nt
    ! integer :: unit_proj = 254
    ! character(LEN=30) :: file_proj
    character(LEN=20) :: ik_string, ib_string, is_string

    ! TYPE (bec_type) :: proj_col     ! proj_col(nwfcU,nbnd)


    if(diag_nmin > diag_nmax) then
       if (ionode) then
          call errore ( 'write_vxc_g', 'diag_nmin > diag_nmax', diag_nmin )
       endif
    endif
    IF (diag_nmin .LT. 1) diag_nmin = 1
    IF (diag_nmax .GT. nbnd) then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting diag_nmax to max number of bands', nbnd
       endif
       diag_nmax = nbnd
    ENDIF
    ndiag = MAX (diag_nmax - diag_nmin + 1, 0)

    if(offdiag_nmin > offdiag_nmax) then
       if (ionode) then
          call errore ( 'write_vxc_g', 'offdiag_nmin > offdiag_nmax', offdiag_nmin )
       endif
    endif
    IF (offdiag_nmin .LT. 1) offdiag_nmin = 1
    IF (offdiag_nmax .GT. nbnd)  then
       if (ionode) then
          write(0,'(a,i6)') 'WARNING: resetting offdiag_nmax to max number of bands', nbnd
       endif
       offdiag_nmax = nbnd
    ENDIF
    noffdiag = MAX (offdiag_nmax - offdiag_nmin + 1, 0)

    IF (ndiag .EQ. 0 .AND. noffdiag .EQ. 0) RETURN

    ALLOCATE( hc( nbnd, nbnd) )

    unit = 4

    iks = global_kpoint_index (nkstot, 1)
    ike = iks + nks - 1

    IF (ndiag .GT. 0) THEN
       ALLOCATE (mtxeld (ndiag, nkstot))
       mtxeld (:, :) = (0.0D0, 0.0D0)
    ENDIF
    IF (noffdiag .GT. 0) THEN
       ALLOCATE (mtxelo (noffdiag, noffdiag, nkstot))
       mtxelo (:, :, :) = (0.0D0, 0.0D0)
    ENDIF

    ! ======
    ldim = 2 * Hubbard_lmax + 1

    ! if (noncolin) then
    !    ALLOCATE( proj(nwfcU,nbnd) )
    ! else
    !    CALL allocate_bec_type ( nwfcU, nbnd, proj_col )
    ! endif
    ! ! ------

    DO ik = iks, ike

       hpsi(:,:) = (0.0D0, 0.0D0)
       evc(:,:) = (0.0D0, 0.0D0)
       hc(:,:) = (0.0D0,0.0D0)
       !! ikk = ik_local
       ikk = ik - iks + 1
       npw = ngk ( ik - iks + 1 )

       IF ( npol == 1 ) THEN
          kdim = npw
          kdmx = npwx
       ELSE
          kdim = npwx*npol
          kdmx = npwx*npol
       ENDIF
       IF ( lsda ) current_spin = isk(ikk)
       CALL davcio (evc, 2*nwordwfc, iunwfc, ik - iks + 1, -1)
       ! if ( nks > 1 .AND. lda_plus_u .AND. (U_projection .NE. 'pseudo') ) then
       IF ( lda_plus_u .AND. (U_projection .NE. 'pseudo') ) THEN
          CALL get_buffer ( wfcU, nwordwfcU, iunhub, ikk )
       ENDIF

       ! ! ======
       ! if (noncolin) then
       !    DO ibnd = 1, nbnd
       !       DO i = 1, nwfcU
       !          proj(i, ibnd) = zdotc (npwx*npol, wfcU (1, i), 1, evc (1, ibnd), 1)
       !       ENDDO
       !    ENDDO

       !    CALL mp_sum ( proj, intra_bgrp_comm )

       ! do ibnd = 1, nbnd
       !    write(ik_string,'(I5)') ik
       !    write(ib_string,'(I5)') ibnd
       !    unit_proj = ik*100+ibnd*1000
       !    file_proj = 'proj.k'//trim(adjustl(ik_string))//".b"//trim(adjustl(ib_string))//'.dat'
       !    open (unit = unit_proj, file = trim(adjustl(file_proj)), form = 'formatted', status = 'replace')
       !    write(unit_proj,'(A)') "#     na    m1    is1     proj(offsetU(na) + m1 + ldim * (is1 - 1), ibnd)"
       ! enddo

       ! do na = 1, nat
       !    nt = ityp (na)
       !    if ( is_hubbard(nt) ) then
       !       ldim = 2 * Hubbard_l(nt) + 1
       !       do m1 = 1, 2 * Hubbard_l(nt) + 1
       !          do is1 = 1, npol
       !             do ibnd = 1, nbnd
       !                unit_proj = ik*100+ibnd*1000
       !                write(unit_proj,'(I5,I5,I5,2F20.11)') na, m1, is1, proj(offsetU(na)+m1+ldim*(is1-1),ibnd)

       !             enddo
       !          enddo
       !       enddo
       !    endif
       ! enddo

       ! do ibnd = 1, nbnd
       !    unit_proj = ik*100+ibnd*10000
       !    close(unit_proj)
       ! enddo

       ! else

       !    CALL calbec ( npw, wfcU, evc, proj_col )

       !    DO ibnd = 1, nbnd
       !       write(ik_string,'(I5)') ik
       !       write(ib_string,'(I5)') ibnd
       !       write(is_string,'(I5)') current_spin
       !       unit_proj = ik*100+ibnd*1000+current_spin*10
       !       file_proj = 'proj.col.k'//trim(adjustl(ik_string))//".b"//trim(adjustl(ib_string))//'.is'//trim(adjustl(is_string))//'.dat'
       !       open (unit = unit_proj, file = trim(adjustl(file_proj)), form = 'formatted', status = 'replace')
       !       write(unit_proj,'(A)') "#     na    m1    proj_col%k(offsetU(na) + m1, ibnd)"
       !    ENDDO

       !    do na = 1, nat
       !       nt = ityp (na)
       !       if ( is_hubbard(nt) ) then
       !          do m1 = 1, 2 * Hubbard_l(nt) + 1
       !             do ibnd = 1, nbnd
       !                unit_proj = ik*100+ibnd*1000+current_spin*10
       !                write(unit_proj,'(I5,I5,2F20.11)') na, m1, proj_col%k(offsetU(na)+m1,ibnd)
       !             enddo
       !          enddo
       !       endif
       !    enddo

       ! endif

       ! ! ======
       ! if (noncolin) then
       !    deallocate(proj)
       ! else
       !    CALL deallocate_bec_type (proj_col)
       ! endif
       ! ! ------

       ! ======
       ! << h_psi.f90 >>
       ! ... Here we add the Hubbard potential times psi
       IF ( lda_plus_u .AND. U_projection.NE."pseudo" ) THEN
          IF (noncolin) THEN
             CALL vhpsi_nc( npwx, npw, nbnd, evc, hpsi )
          ELSE
             CALL vhpsi( npwx, npw, nbnd, evc, hpsi )
          ENDIF
       ENDIF
       ! ======
       ! C := alpha*op( A )*op( B ) + beta*C,
       ! zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
       ! hc = (evc).H . hpsi
       CALL ZGEMM( 'C', 'N', nbnd, nbnd, kdim, ( 1.D0, 0.D0 ), evc, kdmx,  hpsi, kdmx, ( 0.D0, 0.D0 ), hc, nbnd )

       CALL mp_sum(  hc , intra_bgrp_comm )

       ! diagonal
       IF (ndiag .GT. 0) THEN
          DO ib = diag_nmin, diag_nmax
             ! mtxeld in units of eV
             mtxeld (ib - diag_nmin + 1, ik) = hc (ib, ib) * CMPLX (rytoev, 0.0D0)
          ENDDO
       ENDIF
       ! off-diagonal
       IF (noffdiag .gt. 0) THEN
          DO ib = offdiag_nmin, offdiag_nmax
             DO ib2 = offdiag_nmin, offdiag_nmax
                ! mtxeld in units of eV
                mtxelo (ib2-offdiag_nmin+1, ib-offdiag_nmin+1, ik) = hc (ib2, ib) * CMPLX (rytoev, 0.0D0)
             ENDDO
          ENDDO
       ENDIF
    ENDDO ! ik_global

    DEALLOCATE (hc)

    IF (ndiag .GT. 0) CALL mp_sum (mtxeld, inter_pool_comm)
    IF (noffdiag .GT. 0) CALL mp_sum (mtxelo, inter_pool_comm)
    ! xk(1:3) : cartesian ==> crytal
    CALL cryst_to_cart (nkstot, xk, at, -1)

    IF (ionode) THEN
       OPEN (unit = unit, file = TRIM (output_file_name), &
            form = 'formatted', status = 'replace')

       IF (nspin .eq. 4) THEN
          nspin_ = 1
       ELSE
          nspin_ = nspin
       ENDIF

       DO ik = 1, nkstot / nspin_
          WRITE (unit, 101) xk(:, ik), nspin_ * ndiag, &
               nspin_ * noffdiag **2
          DO is = 1, nspin_
             IF (ndiag .GT. 0) THEN
                DO ib = diag_nmin, diag_nmax
                   WRITE (unit, 102) is, ib, mtxeld &
                        (ib - diag_nmin + 1, ik + (is - 1) * nkstot / nspin_)
                ENDDO
             ENDIF
             IF (noffdiag .GT. 0) THEN
                DO ib = offdiag_nmin, offdiag_nmax
                   DO ib2 = offdiag_nmin, offdiag_nmax
                      WRITE (unit, 103) is, ib2, ib, mtxelo &
                           (ib2 - offdiag_nmin + 1, ib - offdiag_nmin + 1, &
                           ik + (is - 1) * nkstot / nspin_)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       CLOSE (unit = unit, status = 'keep')
    ENDIF

    ! ======

    ! ------

    ! ! xk(1:3) : crystal ==> cartesian
    CALL cryst_to_cart (nkstot, xk, bg, 1)

    IF (ndiag .GT. 0) DEALLOCATE (mtxeld)
    IF (noffdiag .GT. 0) DEALLOCATE (mtxelo)

    RETURN

101 FORMAT (3F13.9, 2I8)
102 FORMAT (2I8, 2F15.9)
103 FORMAT (3I8, 2F15.9)

  END SUBROUTINE write_vhub_g

  ! !! --------------------------------------------

  ! SUBROUTINE write_vscg ( output_file_name, real_or_complex, symm_type )

  !   USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  !   USE constants, ONLY : pi, tpi, eps6
  !   USE fft_base, ONLY : dfftp
  !   USE fft_interfaces, ONLY : fwfft
  !   USE gvect, ONLY : ngm, ngm_g, ig_l2g, mill, ecutrho
  !   USE io_global, ONLY : ionode
  !   USE ions_base, ONLY : nat, atm, ityp, tau
  !   USE kinds, ONLY : DP
  !   USE lsda_mod, ONLY : nspin
  !   USE mp, ONLY : mp_sum
  !   USE mp_bands, ONLY : intra_bgrp_comm
  !   USE mp_pools, ONLY : intra_pool_comm
  !   USE scf, ONLY : vltot, v
  !   USE symm_base, ONLY : s, ftau, nsym
  !   USE wavefunctions, ONLY : psic
  !   USE matrix_inversion

  !   IMPLICIT NONE

  !   character ( len = 256 ), intent (in) :: output_file_name
  !   integer, intent (in) :: real_or_complex
  !   character ( len = 9 ), intent (in) :: symm_type

  !   character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  !   integer :: unit, id, is, ir, ig, i, j, k, ierr
  !   integer :: nd, ns, nr, ng_l, ng_g
  !   integer :: ntran, cell_symmetry, nrecord
  !   real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  !   real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  !   real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  !   integer, allocatable :: g_g ( :, : )
  !   real (DP), allocatable :: vscr_g ( :, : )
  !   complex (DP), allocatable :: vscg_g ( :, : )

  !   INTEGER, EXTERNAL :: atomic_number

  !   CALL date_and_tim ( cdate, ctime )
  !   WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  !   WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  !   ! this is supposed to be VSC-Real/Complex but BGW wfn_rho_vxc IO
  !   ! does not recognize VSC header so we are using VXC instead
  !   IF ( real_or_complex .EQ. 1 ) THEN
  !      WRITE ( stitle, '("VXC-Real",24X)' )
  !   ELSE
  !      WRITE ( stitle, '("VXC-Complex",21X)' )
  !   ENDIF

  !   unit = 4
  !   nrecord = 1
  !   nd = 3

  !   ns = nspin
  !   nr = dfftp%nnr
  !   ng_l = ngm
  !   ng_g = ngm_g

  !   ierr = 0
  !   IF ( ibrav .EQ. 0 ) THEN
  !      IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
  !         cell_symmetry = 0
  !      ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
  !         cell_symmetry = 1
  !      ELSE
  !         ierr = 1
  !      ENDIF
  !   ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
  !      cell_symmetry = 0
  !   ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
  !      cell_symmetry = 1
  !   ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
  !      cell_symmetry = 0
  !   ELSE
  !      ierr = 1
  !   ENDIF
  !   IF ( ierr .GT. 0 ) &
  !        CALL errore ( 'write_vscg', 'cell_symmetry', ierr )

  !   ntran = nsym
  !   DO i = 1, ntran
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            r1 ( k, j ) = dble ( s ( k, j, i ) )
  !         ENDDO
  !      ENDDO
  !      CALL invmat ( 3, r1, r2 )
  !      t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
  !      t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
  !      t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
  !      DO j = 1, nd
  !         t2 ( j ) = 0.0D0
  !         DO k = 1, nd
  !            t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
  !         ENDDO
  !         IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
  !              t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
  !         IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
  !              t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
  !      ENDDO
  !      DO j = 1, nd
  !         translation ( j, i ) = t2 ( j ) * tpi
  !      ENDDO
  !   ENDDO

  !   CALL check_inversion ( real_or_complex, nsym, s, nspin, .true., .true., translation )

  !   alat2 = alat ** 2
  !   recvol = 8.0D0 * pi**3 / omega

  !   DO i = 1, nd
  !      DO j = 1, nd
  !         adot ( j, i ) = 0.0D0
  !      ENDDO
  !   ENDDO
  !   DO i = 1, nd
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            adot ( j, i ) = adot ( j, i ) + &
  !                 at ( k, j ) * at ( k, i ) * alat2
  !         ENDDO
  !      ENDDO
  !   ENDDO

  !   DO i = 1, nd
  !      DO j = 1, nd
  !         bdot ( j, i ) = 0.0D0
  !      ENDDO
  !   ENDDO
  !   DO i = 1, nd
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            bdot ( j, i ) = bdot ( j, i ) + &
  !                 bg ( k, j ) * bg ( k, i ) * tpiba2
  !         ENDDO
  !      ENDDO
  !   ENDDO

  !   ALLOCATE ( g_g ( nd, ng_g ) )
  !   ALLOCATE ( vscr_g ( ng_g, ns ) )
  !   ALLOCATE ( vscg_g ( ng_g, ns ) )

  !   DO ig = 1, ng_g
  !      DO id = 1, nd
  !         g_g ( id, ig ) = 0
  !      ENDDO
  !   ENDDO
  !   DO is = 1, ns
  !      DO ig = 1, ng_g
  !         vscg_g ( ig, is ) = ( 0.0D0, 0.0D0 )
  !      ENDDO
  !   ENDDO

  !   DO ig = 1, ng_l
  !      g_g ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
  !      g_g ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
  !      g_g ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  !   ENDDO
  !   vscr_g ( :, : ) = 0.0D0
  !   DO is = 1, ns
  !      DO ir = 1, nr
  !         psic ( ir ) = CMPLX ( v%of_r ( ir, is ) + vltot ( ir ), 0.0D0, KIND=dp )
  !      ENDDO
  !      CALL fwfft ( 'Rho', psic, dfftp )
  !      DO ig = 1, ng_l
  !         vscg_g ( ig_l2g ( ig ), is ) = psic ( dfftp%nl ( ig ) )
  !      ENDDO
  !   ENDDO

  !   CALL mp_sum ( g_g, intra_bgrp_comm )
  !   CALL mp_sum ( vscg_g, intra_bgrp_comm )

  !   IF ( ionode ) THEN
  !      OPEN ( unit = unit, file = TRIM ( output_file_name ), &
  !           form = 'unformatted', status = 'replace' )
  !      WRITE ( unit ) stitle, sdate, stime
  !      WRITE ( unit ) ns, ng_g, ntran, cell_symmetry, nat, ecutrho
  !      WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3
  !      WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
  !           ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
  !      WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
  !           ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
  !      WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
  !      WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
  !      WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
  !      WRITE ( unit ) nrecord
  !      WRITE ( unit ) ng_g
  !      WRITE ( unit ) ( ( g_g ( id, ig ), id = 1, nd ), ig = 1, ng_g )
  !      WRITE ( unit ) nrecord
  !      WRITE ( unit ) ng_g
  !      IF ( real_or_complex .EQ. 1 ) THEN
  !         WRITE ( unit ) ( ( dble ( vscg_g ( ig, is ) ), &
  !              ig = 1, ng_g ), is = 1, ns )
  !      ELSE
  !         WRITE ( unit ) ( ( vscg_g ( ig, is ), &
  !              ig = 1, ng_g ), is = 1, ns )
  !      ENDIF
  !      CLOSE ( unit = unit, status = 'keep' )
  !   ENDIF

  !   DEALLOCATE ( vscg_g )
  !   DEALLOCATE ( vscr_g )
  !   DEALLOCATE ( g_g )

  !   RETURN

  ! END SUBROUTINE write_vscg

  ! !-------------------------------------------------------------------------------

  ! SUBROUTINE write_vkbg (output_file_name, symm_type, wfng_kgrid, &
  !      wfng_nk1, wfng_nk2, wfng_nk3, wfng_dk1, wfng_dk2, wfng_dk3)

  !   USE cell_base, ONLY : omega, alat, tpiba, tpiba2, at, bg, ibrav
  !   USE constants, ONLY : pi, tpi, eps6
  !   USE fft_base, ONLY : dfftp
  !   USE gvect, ONLY : ngm, ngm_g, ig_l2g, g, mill, ecutrho
  !   USE io_global, ONLY : ionode, ionode_id
  !   USE ions_base, ONLY : nat, atm, ityp, tau, nsp
  !   USE kinds, ONLY : DP
  !   USE klist, ONLY : xk, wk, ngk, nks, nkstot, igk_k
  !   USE lsda_mod, ONLY : nspin, isk
  !   USE mp, ONLY : mp_sum, mp_max, mp_get, mp_barrier
  !   USE mp_world, ONLY : mpime, nproc, world_comm
  !   USE mp_bands, ONLY : intra_bgrp_comm, nbgrp
  !   USE mp_pools, ONLY : me_pool, root_pool, npool, nproc_pool, &
  !        intra_pool_comm, inter_pool_comm
  !   USE mp_wave, ONLY : mergewf
  !   USE start_k, ONLY : nk1, nk2, nk3, k1, k2, k3
  !   USE symm_base, ONLY : s, ftau, nsym
  !   USE uspp, ONLY : nkb, vkb, deeq
  !   USE uspp_param, ONLY : nhm, nh
  !   USE wvfct, ONLY : npwx
  !   USE gvecw, ONLY : ecutwfc
  !   USE matrix_inversion

  !   IMPLICIT NONE

  !   character (len = 256), intent (in) :: output_file_name
  !   character ( len = 9 ), intent (in) :: symm_type
  !   logical, intent (in) :: wfng_kgrid
  !   integer, intent (in) :: wfng_nk1
  !   integer, intent (in) :: wfng_nk2
  !   integer, intent (in) :: wfng_nk3
  !   real (DP), intent (in) :: wfng_dk1
  !   real (DP), intent (in) :: wfng_dk2
  !   real (DP), intent (in) :: wfng_dk3

  !   character :: cdate*9, ctime*9, sdate*32, stime*32, stitle*32
  !   integer :: i, j, k, ierr, ik, is, ig, ikb, iat, isp, ih, jh, &
  !        unit, iks, ike, npw, npw_g, npwx_g, ngg, ipsour, &
  !        igwx, local_pw, id, nd, ntran, cell_symmetry, nrecord
  !   real (DP) :: alat2, recvol, t1 ( 3 ), t2 ( 3 )
  !   real (DP) :: r1 ( 3, 3 ), r2 ( 3, 3 ), adot ( 3, 3 )
  !   real (DP) :: bdot ( 3, 3 ), translation ( 3, 48 )
  !   integer, allocatable :: kmap ( : )
  !   integer, allocatable :: smap ( : )
  !   integer, allocatable :: gvec ( :, : )
  !   integer, allocatable :: ngk_g ( : )
  !   integer, allocatable :: igk_l2g ( :, : )
  !   integer, allocatable :: itmp ( : )
  !   integer, allocatable :: igwk ( : )
  !   integer, allocatable :: igwf_l2g ( : )
  !   integer, allocatable :: ipmask ( : )
  !   complex (DP), allocatable :: vkb_g ( : )

  !   INTEGER, EXTERNAL :: atomic_number, global_kpoint_index

  !   IF ( nkb == 0 ) RETURN

  !   CALL date_and_tim ( cdate, ctime )
  !   WRITE ( sdate, '(A2,"-",A3,"-",A4,21X)' ) cdate(1:2), cdate(3:5), cdate(6:9)
  !   WRITE ( stime, '(A8,24X)' ) ctime(1:8)
  !   ! BGW wfn_rho_vxc IO does not recognize VKB header so this file
  !   ! is read directly by SAPO code in BerkeleyGW
  !   WRITE ( stitle, '("VKB-Complex",21X)' )

  !   unit = 4
  !   nrecord = 1
  !   nd = 3

  !   iks = global_kpoint_index (nkstot, 1)
  !   ike = iks + nks - 1

  !   ierr = 0
  !   IF ( ibrav .EQ. 0 ) THEN
  !      IF ( TRIM ( symm_type ) .EQ. 'cubic' ) THEN
  !         cell_symmetry = 0
  !      ELSEIF ( TRIM ( symm_type ) .EQ. 'hexagonal' ) THEN
  !         cell_symmetry = 1
  !      ELSE
  !         ierr = 1
  !      ENDIF
  !   ELSEIF ( abs ( ibrav ) .GE. 1 .AND. abs ( ibrav ) .LE. 3 ) THEN
  !      cell_symmetry = 0
  !   ELSEIF ( abs ( ibrav ) .GE. 4 .AND. abs ( ibrav ) .LE. 5 ) THEN
  !      cell_symmetry = 1
  !   ELSEIF ( abs ( ibrav ) .GE. 6 .AND. abs ( ibrav ) .LE. 14 ) THEN
  !      cell_symmetry = 0
  !   ELSE
  !      ierr = 1
  !   ENDIF
  !   IF ( ierr .GT. 0 ) &
  !        CALL errore ( 'write_vkbg', 'cell_symmetry', ierr )

  !   ntran = nsym
  !   DO i = 1, ntran
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            r1 ( k, j ) = dble ( s ( k, j, i ) )
  !         ENDDO
  !      ENDDO
  !      CALL invmat ( 3, r1, r2 )
  !      t1 ( 1 ) = dble ( ftau ( 1, i ) ) / dble ( dfftp%nr1 )
  !      t1 ( 2 ) = dble ( ftau ( 2, i ) ) / dble ( dfftp%nr2 )
  !      t1 ( 3 ) = dble ( ftau ( 3, i ) ) / dble ( dfftp%nr3 )
  !      DO j = 1, nd
  !         t2 ( j ) = 0.0D0
  !         DO k = 1, nd
  !            t2 ( j ) = t2 ( j ) + r2 ( k, j ) * t1 ( k )
  !         ENDDO
  !         IF ( t2 ( j ) .GE. eps6 + 0.5D0 ) &
  !              t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) + 0.5D0 ) )
  !         IF ( t2 ( j ) .LT. eps6 - 0.5D0 ) &
  !              t2 ( j ) = t2 ( j ) - dble ( int ( t2 ( j ) - 0.5D0 ) )
  !      ENDDO
  !      DO j = 1, nd
  !         translation ( j, i ) = t2 ( j ) * tpi
  !      ENDDO
  !   ENDDO

  !   alat2 = alat ** 2
  !   recvol = 8.0D0 * pi**3 / omega

  !   DO i = 1, nd
  !      DO j = 1, nd
  !         adot ( j, i ) = 0.0D0
  !      ENDDO
  !   ENDDO
  !   DO i = 1, nd
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            adot ( j, i ) = adot ( j, i ) + &
  !                 at ( k, j ) * at ( k, i ) * alat2
  !         ENDDO
  !      ENDDO
  !   ENDDO

  !   DO i = 1, nd
  !      DO j = 1, nd
  !         bdot ( j, i ) = 0.0D0
  !      ENDDO
  !   ENDDO
  !   DO i = 1, nd
  !      DO j = 1, nd
  !         DO k = 1, nd
  !            bdot ( j, i ) = bdot ( j, i ) + &
  !                 bg ( k, j ) * bg ( k, i ) * tpiba2
  !         ENDDO
  !      ENDDO
  !   ENDDO

  !   ALLOCATE ( kmap ( nkstot ) )
  !   ALLOCATE ( smap ( nkstot ) )

  !   DO i = 1, nkstot
  !      j = ( i - 1 ) / nspin
  !      k = i - 1 - j * nspin
  !      kmap ( i ) = j + k * ( nkstot / nspin ) + 1
  !      smap ( i ) = k + 1
  !   ENDDO
  !   ierr = 0
  !   DO i = 1, nkstot
  !      ik = kmap ( i )
  !      is = smap ( i )
  !      IF ( ik .GE. iks .AND. ik .LE. ike .AND. is .NE. isk ( ik ) ) &
  !           ierr = ierr + 1
  !   ENDDO
  !   CALL mp_max ( ierr, world_comm )
  !   IF ( ierr .GT. 0 ) &
  !        CALL errore ( 'write_vkbg', 'smap', ierr )

  !   ALLOCATE ( gvec ( 3, ngm_g ) )
  !   gvec = 0
  !   DO ig = 1, ngm
  !      gvec ( 1, ig_l2g ( ig ) ) = mill ( 1, ig )
  !      gvec ( 2, ig_l2g ( ig ) ) = mill ( 2, ig )
  !      gvec ( 3, ig_l2g ( ig ) ) = mill ( 3, ig )
  !   ENDDO
  !   CALL mp_sum ( gvec, intra_bgrp_comm )

  !   ALLOCATE ( ngk_g ( nkstot ) )
  !   ALLOCATE ( igk_l2g ( npwx, nks ) )
  !   ngk_g = 0
  !   igk_l2g = 0
  !   DO ik = 1, nks
  !      npw = ngk ( ik )
  !      DO ig = 1, npw
  !         igk_l2g ( ig, ik ) = ig_l2g ( igk_k ( ig, ik ) )
  !      ENDDO
  !   ENDDO
  !   DO ik = 1, nks
  !      ngk_g ( ik + iks - 1 ) = ngk ( ik )
  !   ENDDO
  !   CALL mp_sum ( ngk_g, inter_pool_comm )
  !   CALL mp_sum ( ngk_g, intra_pool_comm )
  !   ngk_g = ngk_g / nbgrp
  !   npw_g = MAXVAL ( igk_l2g ( :, : ) )
  !   CALL mp_max ( npw_g, intra_pool_comm )
  !   npwx_g = MAXVAL ( ngk_g ( : ) )

  !   CALL cryst_to_cart (nkstot, xk, at, -1)

  !   IF (ionode) THEN
  !      OPEN (unit = unit, file = TRIM (output_file_name), &
  !           form = 'unformatted', status = 'replace')
  !      WRITE ( unit ) stitle, sdate, stime
  !      WRITE ( unit ) nspin, ngm_g, ntran, cell_symmetry, nat, ecutrho, &
  !           nkstot / nspin, nsp, nkb, nhm, npwx_g, ecutwfc
  !      IF ( wfng_kgrid ) THEN
  !         WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, wfng_nk1, wfng_nk2, wfng_nk3, &
  !              wfng_dk1, wfng_dk2, wfng_dk3
  !      ELSE
  !         WRITE ( unit ) dfftp%nr1, dfftp%nr2, dfftp%nr3, nk1, nk2, nk3, &
  !              0.5D0 * dble ( k1 ), 0.5D0 * dble ( k2 ), 0.5D0 * dble ( k3 )
  !      ENDIF
  !      WRITE ( unit ) omega, alat, ( ( at ( j, i ), j = 1, nd ), i = 1, nd ), &
  !           ( ( adot ( j, i ), j = 1, nd ), i = 1, nd )
  !      WRITE ( unit ) recvol, tpiba, ( ( bg ( j, i ), j = 1, nd ), i = 1, nd ), &
  !           ( ( bdot ( j, i ), j = 1, nd ), i = 1, nd )
  !      WRITE ( unit ) ( ( ( s ( k, j, i ), k = 1, nd ), j = 1, nd ), i = 1, ntran )
  !      WRITE ( unit ) ( ( translation ( j, i ), j = 1, nd ), i = 1, ntran )
  !      WRITE ( unit ) ( ( tau ( j, i ), j = 1, nd ), atomic_number ( atm ( ityp ( i ) ) ), i = 1, nat )
  !      WRITE ( unit ) ( ngk_g ( ik ), ik = 1, nkstot / nspin )
  !      WRITE ( unit ) ( wk ( ik ) * dble ( nspin ) / 2.0D0, ik = 1, nkstot / nspin )
  !      WRITE ( unit ) ( ( xk ( id, ik ), id = 1, nd ), ik = 1, nkstot / nspin )
  !      WRITE ( unit ) ( ityp ( iat ), iat = 1, nat )
  !      WRITE ( unit ) ( nh ( isp ), isp = 1, nsp )
  !      WRITE ( unit ) ( ( ( ( deeq ( jh, ih, iat, is ), &
  !           jh = 1, nhm ), ih = 1, nhm ), iat = 1, nat ), is = 1, nspin )
  !      WRITE ( unit ) nrecord
  !      WRITE ( unit ) ngm_g
  !      WRITE ( unit ) ( ( gvec ( id, ig ), id = 1, nd ), ig = 1, ngm_g )
  !   ENDIF

  !   CALL cryst_to_cart (nkstot, xk, bg, 1)

  !   ALLOCATE ( igwk ( npwx_g ) )

  !   DO i = 1, nkstot

  !      ik = kmap ( i )
  !      is = smap ( i )

  !      igwk = 0

  !      ALLOCATE ( itmp ( npw_g ) )
  !      itmp = 0
  !      IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
  !         DO ig = 1, ngk ( ik - iks + 1 )
  !            itmp ( igk_l2g ( ig, ik - iks + 1 ) ) = igk_l2g ( ig, ik - iks + 1 )
  !         ENDDO
  !      ENDIF

  !      ! CALL mp_sum ( itmp, intra_bgrp_comm )
  !      CALL mp_sum ( itmp, world_comm )

  !      ngg = 0
  !      DO ig = 1, npw_g
  !         IF ( itmp ( ig ) .EQ. ig ) THEN
  !            ngg = ngg + 1
  !            igwk ( ngg ) = ig
  !         ENDIF
  !      ENDDO
  !      DEALLOCATE ( itmp )

  !      IF ( ionode ) THEN
  !         IF ( is .EQ. 1 ) THEN
  !            WRITE ( unit ) nrecord
  !            WRITE ( unit ) ngk_g ( ik )
  !            WRITE ( unit ) ( ( gvec ( j, igwk ( ig ) ), j = 1, 3 ), &
  !                 ig = 1, ngk_g ( ik ) )
  !         ENDIF
  !      ENDIF

  !      local_pw = 0
  !      IF ( ik .GE. iks .AND. ik .LE. ike ) THEN
  !         npw = ngk ( ik - iks + 1 )
  !!         CALL init_us_2 ( npw, igk_k(1, ik-iks+1), xk ( 1, ik ), vkb )
  !         CALL init_us_2 ( npw, igk_k(1, ik-iks+1), xk ( 1, ik-iks+1), vkb )
  !         local_pw = npw
  !      ENDIF

  !      ALLOCATE ( igwf_l2g ( local_pw ) )
  !      igwf_l2g = 0
  !      DO ig = 1, local_pw
  !         ngg = igk_l2g ( ig, ik - iks + 1 )
  !         DO j = 1, ngk_g ( ik )
  !            IF ( ngg .EQ. igwk ( j ) ) THEN
  !               igwf_l2g ( ig ) = j
  !               EXIT
  !            ENDIF
  !         ENDDO
  !      ENDDO

  !      ALLOCATE ( ipmask ( nproc ) )
  !      ipmask = 0
  !      ipsour = ionode_id
  !      IF ( npool .GT. 1 ) THEN
  !         IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
  !            IF ( me_pool .EQ. root_pool ) ipmask ( mpime + 1 ) = 1
  !         ENDIF
  !         CALL mp_sum ( ipmask, world_comm )
  !         DO j = 1, nproc
  !            IF ( ipmask ( j ) .EQ. 1 ) ipsour = j - 1
  !         ENDDO
  !      ENDIF
  !      DEALLOCATE ( ipmask )

  !      igwx = 0
  !      ierr = 0
  !      IF ( ik .GE. iks .AND. ik .LE. ike ) &
  !           igwx = MAXVAL ( igwf_l2g ( 1 : local_pw ) )
  !      CALL mp_max ( igwx, intra_pool_comm )
  !      IF ( ipsour .NE. ionode_id ) &
  !           CALL mp_get ( igwx, igwx, mpime, ionode_id, ipsour, 1, world_comm )
  !      ierr = 0
  !      IF ( ik .GE. iks .AND. ik .LE. ike .AND. igwx .NE. ngk_g ( ik ) ) &
  !           ierr = 1
  !      CALL mp_max ( ierr, world_comm )
  !      IF ( ierr .GT. 0 ) &
  !           CALL errore ( 'write_vkbg', 'igwx ngk_g', ierr )

  !      ALLOCATE ( vkb_g ( MAX ( 1, igwx ) ) )

  !      DO ikb = 1, nkb

  !         vkb_g = ( 0.0D0, 0.0D0 )
  !         IF ( npool .GT. 1 ) THEN
  !            IF ( ( ik .GE. iks ) .AND. ( ik .LE. ike ) ) THEN
  !               CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
  !                    me_pool, nproc_pool, root_pool, intra_pool_comm )
  !            ENDIF
  !            IF ( ipsour .NE. ionode_id ) THEN
  !               CALL mp_get ( vkb_g, vkb_g, mpime, ionode_id, ipsour, ikb, &
  !                    world_comm )
  !            ENDIF
  !         ELSE
  !            CALL mergewf ( vkb ( :, ikb ), vkb_g, local_pw, igwf_l2g, &
  !                 mpime, nproc, ionode_id, world_comm )
  !         ENDIF

  !         IF ( ionode ) THEN
  !            WRITE ( unit ) nrecord
  !            WRITE ( unit ) ngk_g ( ik )
  !            WRITE ( unit ) ( vkb_g ( ig ), ig = 1, igwx )
  !         ENDIF

  !      ENDDO

  !      DEALLOCATE ( vkb_g )
  !      DEALLOCATE ( igwf_l2g )

  !   ENDDO

  !   IF ( ionode ) THEN
  !      CLOSE ( unit = unit, status = 'keep' )
  !   ENDIF

  !   DEALLOCATE ( igwk )
  !   DEALLOCATE ( igk_l2g )
  !   DEALLOCATE ( ngk_g )
  !   DEALLOCATE ( gvec )
  !   DEALLOCATE ( smap )
  !   DEALLOCATE ( kmap )

  !   RETURN

  ! END SUBROUTINE write_vkbg

  !-------------------------------------------------------------------------------

  subroutine check_inversion(real_or_complex, ntran, mtrx, nspin, warn, real_need_inv, tnp)

    ! check_inversion    Originally By David A. Strubbe    Last Modified 11/18/2013
    ! Check whether our choice of real/complex version is appropriate given the
    ! presence or absence of inversion symmetry.

    USE constants, ONLY : eps6
    USE io_global, ONLY : ionode
    USE kinds, ONLY : DP

    implicit none

    integer, intent(in) :: real_or_complex
    integer, intent(in) :: ntran
    integer, intent(in) :: mtrx(3, 3, 48) !< symmetry operations matrices
    integer, intent(in) :: nspin
    logical, intent(in) :: warn !< set to false to suppress warnings, for converters
    logical, intent(in) :: real_need_inv !< use for generating routines to block real without inversion
    !! this is not always true so that it is possible to run real without using symmetries
    real(DP), intent(in) :: tnp(3, 48) !< fractional translations.
    !! optional only to avoid changing external interface for library.

    integer :: invflag, isym, ii, jj, itest
    logical :: origin_inv

    invflag = 0
    origin_inv = .false.
    do isym = 1, ntran
       itest = 0
       do ii = 1, 3
          do jj = 1, 3
             if(ii .eq. jj) then
                itest = itest + (mtrx(ii, jj, isym) + 1)**2
             else
                itest = itest + mtrx(ii, jj, isym)**2
             endif
          enddo
       enddo
       if(itest .eq. 0) then
          invflag = invflag + 1
          !if(present(tnp)) then
          if(sum(abs(tnp(1:3, isym))) < eps6) origin_inv = .true.
          !else
          !  origin_inv = .true.
          !endif
       endif
    enddo
    if(invflag > 0 .and. .not. origin_inv) then
       write(0, '(a)') "WARNING: Inversion symmetry is present only with a fractional translation."
       write(0, '(a)') "Apply the translation so inversion is about the origin, to be able to use the real version."
    endif
    if(invflag .gt. 1) write(0, '(a)') "WARNING: More than one inversion symmetry operation is present."

    !  if(invflag > 0 .and. .not. present(tnp)) then
    !    write(0, '(a)') "WARNING: check_inversion did not receive fractional translations."
    !    write(0, '(a)') "Cannot confirm that inversion symmetry is about the origin for use of real version."
    !  endif

    if(real_or_complex .eq. 2) then
       if(origin_inv .and. warn .and. nspin == 1) then
          if(ionode) &
               write(0, '(a)') "WARNING: Inversion symmetry about the origin is present. The real version would be faster."
       endif
    else
       if(.not. origin_inv) then
          if(real_need_inv) then
             call errore("check_inversion", "The real version cannot be used without inversion symmetry about the origin.", -1)
          endif
          if(ionode) then
             write(0, '(a)') "WARNING: Inversion symmetry about the origin is absent in symmetries used to reduce k-grid."
             write(0, '(a)') "Be sure inversion about the origin is still a spatial symmetry, or you must use complex version instead."
          endif
       endif
       if(nspin > 1) then
          call errore("check_inversion", "Real version may only be used for spin-unpolarized calculations.", nspin)
       endif
    endif

    return

  end subroutine check_inversion

  !-------------------------------------------------------------------------------

END PROGRAM pw2bgw
