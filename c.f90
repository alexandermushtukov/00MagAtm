!======================================================================
!======================================================================
program MagAtm
use omp_lib
use black_body
implicit none

real*8 :: pi=3.141592653589793d0
real*8 :: E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda,T_eff
real*8 :: m_min,m_max,kappa_T,tau_min,tau_max
integer :: n_m,n_mu,n_fi,i,n_E
real*8 :: g14,E_min,E_max,dE,F_E_target,F_target,T_floor
real*8, allocatable :: mas_tau_TkeV(:,:), mas_m_rho(:,:), flux_tot(:,:), flux_E(:,:), J_E(:,:), &
                       kabs_mean(:,:), spec(:,:), num(:,:), den(:,:), B(:), dBdT(:),            &
                       dflux(:), dT(:), dT_sm(:), flux(:), gradF(:), weight_m(:)
integer :: i_iter,jj
character*100 :: file_sp,file_t
real*8 :: eps
real*8 :: sigma_SB_22_keV,T_bot,T_bot_new

! --- Iteration-control variables ---
real*8 :: lambda_T, lambda_bot
real*8 :: max_rel_flux_err, max_rel_dT, rel_flux_err
real*8 :: flux_surf, dm, gradF_scale
real*8 :: tol_flux, tol_dT
real*8 :: m_turn, p_weight, wloc
integer :: iter_max

real(8) :: weights(0:2)

  weights(0) = 0.50d0   ! Central point
  weights(1) = 0.15d0   ! First neighbour on each side
  weights(2) = 0.10d0   ! Second neighbour on each side

  sigma_SB_22_keV = 1.027d2
  eps      = 1.d-2
  T_floor  = 0.1d0

  lambda_T   = 0.25d0   ! Local temperature-correction damping
  lambda_bot = 0.25d0   ! Bottom-temperature damping
  tol_flux   = 1.d-3
  tol_dT     = 1.d-3
  iter_max   = 300

  m_turn   = 100.d0     ! Depth scale where the correction starts to be suppressed
  p_weight = 1.d0       ! Power-law index of depth suppression

  200 format (100(es11.4,"   ") )

  ! === Physical parameters ===
  T_eff     = 2.d0
  m_ns      = 1.4d0
  R6        = 1.d0
  B12       = 1.d2
  theta_B   = 0.1d0
  Z         = 1.d0
  A         = 1.d0
  dot_m_6   = 0.d0
  ln_Lambda = 10.d0
  ! ===========================
  g14 = 1.328d0*m_ns/R6**2

  file_sp ="./res/res_sp_1_"
  open(unit = 25, file = file_sp); close(25)
  file_t ="./res/res_t_1_"
  open(unit = 26, file = file_t); close(26)

  ! === Numerical parameters ===
  m_min = 1.d-2
  m_max = 1.d+4
  n_m   = 50
  n_mu  = 10
  n_fi  = 1
  n_E   = 20
  E_min = 0.1d0
  E_max = 7.d0
  ! ============================

  write(*,*) "# n_m=", n_m

  kappa_T = 0.34d0
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T

  allocate( mas_tau_TkeV(n_m,2), mas_m_rho(n_m,2), flux_tot(n_m,2), flux_E(n_m,2), J_E(n_m,2), &
            kabs_mean(n_m,2), spec(n_E,3), num(n_m,2), den(n_m,2), B(n_m), dBdT(n_m),          &
            dflux(n_m), dT(n_m), dT_sm(n_m), flux(n_m), gradF(n_m), weight_m(n_m) )

  ! === Initial temperature structure ===
  i = 1
  do while (i .le. n_m)
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = T_eff / 1.2d0
    i = i + 1
  end do

  T_bot = T_eff

  ! === Temperature iterations ===
  i_iter = 0
  do while (i_iter .lt. iter_max)

    write(*,*) "# start iter:", i_iter

    ! --- Compute hydrostatic structure for the current temperature profile ---
    call get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)

    F_target = 0.d0
    flux_tot(1:n_m,1:2) = 0.d0
    num(1:n_m,1:2)      = 0.d0
    den(1:n_m,1:2)      = 0.d0

    ! --- Integrate radiative transfer over energy ---
    i = 1
    do while (i .le. n_E)

      E  = E_min * (E_max/E_min)**( dble(i-1) / dble(n_E-1) )
      F_E_target = pi * BB_Intensity_22(E,T_eff)
      dE = E * ( (E_max/E_min)**(1.d0/dble(n_E-1)) - 1.d0 )

      call pol_RT_fixE(flux_E,J_E,kabs_mean,E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda, &
                       mas_m_rho,mas_tau_TkeV,T_bot,n_m,n_mu,n_fi)

      spec(i,1) = E
      spec(i,2) = flux_E(1,1)
      spec(i,3) = flux_E(1,2)

      flux_tot(1:n_m,1:2) = flux_tot(1:n_m,1:2) + flux_E(1:n_m,1:2)*dE

      ! These arrays are kept for diagnostics and possible future hybrid correction
      jj = 1
      do while (jj .le. n_m)
        B(jj) = BB_Intensity_22(E,mas_tau_TkeV(jj,2))/2.d0

        dBdT(jj) = ( BB_Intensity_22(E,mas_tau_TkeV(jj,2)*(1.d0+eps)) - &
                     BB_Intensity_22(E,mas_tau_TkeV(jj,2)*(1.d0-eps)) ) / &
                   ( 2.d0*eps*mas_tau_TkeV(jj,2) )
        dBdT(jj) = dBdT(jj)/2.d0

        num(jj,1:2) = num(jj,1:2) + kabs_mean(jj,1:2) * ( J_E(jj,1:2) - B(jj) ) * dE
        den(jj,1:2) = den(jj,1:2) + kabs_mean(jj,1:2) * dBdT(jj) * dE

        jj = jj + 1
      end do

      F_target = F_target + F_E_target*dE
      i = i + 1
    end do

    ! --- Build total flux profile ---
    flux(1:n_m)  = 0.d0
    dflux(1:n_m) = 0.d0
    gradF(1:n_m) = 0.d0

    i = 1
    do while (i .le. n_m)
      flux(i) = flux_tot(i,1) + flux_tot(i,2)
      i = i + 1
    end do

    i = 2
    do while (i .le. n_m)
      dflux(i) = flux(i) - flux(i-1)
      dm = mas_tau_TkeV(i,1) - mas_tau_TkeV(i-1,1)
      gradF(i) = dflux(i) / max(dm,1.d-30)
      i = i + 1
    end do
    gradF(1) = gradF(2)

    ! --- Build temperature correction in a separate procedure ---
    call build_temperature_correction(dT, dT_sm, max_rel_flux_err, max_rel_dT, &
         mas_tau_TkeV, flux, gradF, F_target, n_m, lambda_T, T_floor, &
         m_turn, p_weight, weights)

    ! --- Print atmospheric structure for the current iteration ---
    write(*,*) '# atmosphere structure, iter =', i_iter
    write(*,*) '# i, m, T_keV, rho, dT, dFdm, F, F_target'

    i = 2
    do while (i .le. n_m)
      write(*,'(I5,1X,3F10.4,1X,ES12.4,1X,F10.4,1X,F12.4,1X,F12.4)') &
            i, mas_tau_TkeV(i,1), mas_tau_TkeV(i,2), mas_m_rho(i,2), &
            dT_sm(i), gradF(i), flux(i), F_target
      i = i + 1
    end do

    ! --- Update the temperature profile ---
    i = 1
    do while (i .le. n_m)
      mas_tau_TkeV(i,2) = max(T_floor, mas_tau_TkeV(i,2) + dT_sm(i))
      i = i + 1
    end do

    ! --- Global correction through the bottom boundary temperature ---
    flux_surf = flux(1)
    !if (flux_surf .gt. 1.d-30) then
    !  T_bot_new = T_bot * (F_target/flux_surf)**0.25d0
    !  T_bot     = (1.d0-lambda_bot)*T_bot + lambda_bot*T_bot_new
    !  T_bot     = max(T_floor, T_bot)
    !end if

    write(*,'(a,i5,2(a,es12.4),a,es12.4,a,es12.4,a,es12.4)') &
         '# iter=', i_iter, &
         '  max_flux_err=', max_rel_flux_err, &
         '  max_rel_dT=',   max_rel_dT, &
         '  Fsurf=',        flux_surf, &
         '  Ftarget=',      F_target, &
         '  T_bot=',        T_bot

    if (max_rel_flux_err .lt. tol_flux .and. max_rel_dT .lt. tol_dT) exit

    i_iter = i_iter + 1
  end do

  write(*,*) "# iteration ended, n_iter =", i_iter

  open(unit = 25, file = file_sp, status = 'old', form='formatted')
  i = 1
  do while (i .le. n_E)
    write(25,*) spec(i,1:3)
    i = i + 1
  end do
  close(25)
  write(*,*) "# end0"

  open(unit = 26, file = file_t, status = 'old', form='formatted')
  i = 1
  do while (i .le. n_m)
    write(26,*) mas_tau_TkeV(i,1:2), flux(i), F_target
    i = i + 1
  end do
  close(26)
  write(*,*) "# end1"

144 return
end program MagAtm


!======================================================================
! Build local temperature correction from the flux gradient
!======================================================================
subroutine build_temperature_correction(dT, dT_sm, max_rel_flux_err, max_rel_dT, &
                                        mas_tau_TkeV, flux, gradF, F_target, n,   &
                                        lambda_T, T_floor, m_turn, p_weight, weights)
implicit none

integer, intent(in) :: n
real(8), intent(in) :: mas_tau_TkeV(n,2), flux(n), gradF(n), F_target
real(8), intent(in) :: lambda_T, T_floor, m_turn, p_weight
real(8), intent(in) :: weights(0:2)

real(8), intent(out) :: dT(n), dT_sm(n)
real(8), intent(out) :: max_rel_flux_err, max_rel_dT

integer :: i
real(8) :: rel_flux_err, gradF_scale, wloc
real(8) :: weight_m(n)

  ! --- Use the current maximum flux-divergence amplitude as normalization ---
  gradF_scale = max( maxval(abs(gradF(2:n))), 1.d-30 )

  ! --- Build depth weights for local temperature correction ---
  i = 1
  do while (i .le. n)
    weight_m(i) = 1.d0 / ( 1.d0 + (mas_tau_TkeV(i,1)/m_turn)**p_weight )
    i = i + 1
  end do

  ! --- Diagnostics and raw local temperature correction ---
  max_rel_flux_err = 0.d0
  max_rel_dT       = 0.d0

  dT(1) = 0.d0
  i = 2
  do while (i .le. n)

    rel_flux_err = abs(flux(i) - F_target) / max(F_target,1.d-30)
    if (rel_flux_err .gt. max_rel_flux_err) max_rel_flux_err = rel_flux_err

    ! Local energy-balance correction:
    ! positive dF/dm -> heating
    ! negative dF/dm -> cooling
    wloc = 1.d0   ! or weight_m(i) if depth suppression is needed
    dT(i) = wloc * lambda_T * gradF(i) / gradF_scale * mas_tau_TkeV(i,2) * 2.d0

    ! Additional suppression in the deepest layers
    if (i .ge. n-5) dT(i) = 0.3d0*dT(i)

    ! Limit relative step to avoid oscillatory behaviour
    if (abs(dT(i)) .gt. 0.05d0*mas_tau_TkeV(i,2)) then
      dT(i) = sign(0.05d0*mas_tau_TkeV(i,2), dT(i))
    end if

    i = i + 1
  end do

  rel_flux_err = abs(flux(1) - F_target) / max(F_target,1.d-30)
  if (rel_flux_err .gt. max_rel_flux_err) max_rel_flux_err = rel_flux_err

  ! --- Smooth the correction ---
  call smooth_array(dT_sm, dT, n, 2, weights)

  ! --- Keep the boundary points unchanged if desired ---
  dT_sm(1) = dT(1)
  dT_sm(n) = dT(n)

  ! --- Compute the maximum relative temperature correction after smoothing ---
  i = 1
  do while (i .le. n)
    if (abs(dT_sm(i)) / max(mas_tau_TkeV(i,2),1.d-30) .gt. max_rel_dT) then
      max_rel_dT = abs(dT_sm(i)) / max(mas_tau_TkeV(i,2),1.d-30)
    end if
    i = i + 1
  end do

return
end subroutine build_temperature_correction


!======================================================================
! Smooth a 1D array using neighbour points and fixed weights
!======================================================================
subroutine smooth_array(arr_out, arr_in, n, nside, weights)
implicit none

integer, intent(in) :: n
integer, intent(in) :: nside
real(8), intent(in) :: arr_in(n)
real(8), intent(in) :: weights(0:nside)
real(8), intent(out) :: arr_out(n)

integer :: i, k, j
real(8) :: tsum, wsum

  do i = 1, n
    tsum = weights(0) * arr_in(i)
    wsum = weights(0)

    do k = 1, nside
      j = i - k
      if (j >= 1) then
        tsum = tsum + weights(k) * arr_in(j)
        wsum = wsum + weights(k)
      end if

      j = i + k
      if (j <= n) then
        tsum = tsum + weights(k) * arr_in(j)
        wsum = wsum + weights(k)
      end if
    end do

    arr_out(i) = tsum / wsum
  end do

return
end subroutine smooth_array


!======================================================================================================================
!======================================================================================================================
subroutine temperature_correction_2(mas_tau_TkeV_new,T_bot_new,n_m,flux_tot,mas_tau_TkeV,T_bot,F_target,alpha,T_floor)
implicit none
real*8,intent(out)::mas_tau_TkeV_new(n_m,2),T_bot_new
integer,intent(in)::n_m
real*8,intent(in)::flux_tot(n_m,2),mas_tau_TkeV(n_m,2),T_bot,alpha,F_target,T_floor
integer::i
real*8::flux(n_m),T_old,T_corr,dT(n_m),dT_(n_m),dflux

  mas_tau_TkeV_new(1:n_m,1:2) = mas_tau_TkeV(1:n_m,1:2)

  i = 1
  do while(i.le.n_m)
    flux(i) = flux_tot(i,1) + flux_tot(i,2)
    i = i+1
  end do

  i = 2
  do while(i.le.n_m)
    dflux = flux(i) - flux(i-1)
    if( dflux.gt.0.d0 )then
      ! == heating ==
      dT(i) = 0.05d0
    else
      ! == cooling ==
      dT(i) = -0.05d0
    end if
    i = i + 1
  end do

  dT_(2)=dT(2)
  i=3
  do while(i.le.n_m-1)
    dT_(i) = (dT(i) + 0.25d0*dT(i-1) + 0.25d0*dT(i+1))/1.5d0
    i = i+1
  end do
  dT_(n_m)=dT(n_m)

  i=2
  do while(i.le.n_m)
    mas_tau_TkeV_new(i,2) = mas_tau_TkeV_new(i,2) + dT_(i)
    i = i+1
  end do

  ! == Correct bottom temperature ==
  if( F_target.gt.flux(n_m) )then
    T_bot_new = T_bot
  else
    T_bot_new = T_bot
  end if

return
end subroutine temperature_correction_2
