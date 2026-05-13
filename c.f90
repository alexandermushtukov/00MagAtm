!======================================================================
!======================================================================
program MagAtm
use omp_lib
use black_body
implicit none
real*8 :: pi=3.141592653589793d0
real*8 :: E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda,T_eff
real*8 :: m_min,m_max,kappa_T,tau_min,tau_max
integer :: n_m,n_mu,n_fi,i,j,n_E,n_threads
real*8 :: g14,E_min,E_max,dE,F_E_target,F_target,T_floor
real*8, allocatable :: mas_tau_TkeV(:,:), mas_m_rho(:,:), flux_tot(:,:), flux_E(:,:), J_E(:,:), J_tot(:,:), &
                       kabs_mean(:,:),k_p(:),k_J(:),k_f(:),dk_p(:),dk_J(:),dk_f(:),du(:),dFz(:),u(:),u_p(:),Fz(:),&
                       spec(:,:),spec_(:,:),&
                       dflux(:), dT(:), dT_sm(:), flux(:), gradF(:), weight_m(:),TkeV_old(:), &
                       flux_E_prev(:,:), J_E_prev(:,:), kabs_mean_prev(:,:), &
                       dk_p_prev(:), dk_J_prev(:), dk_f_prev(:), du_prev(:), dFz_prev(:)

!== arrays for OpenMP over photon energies ==!
real*8, allocatable :: E_grid(:),F_E_target_grid(:)
real*8, allocatable :: flux_E_all(:,:,:),J_E_all(:,:,:)
real*8, allocatable :: dk_p_all(:,:),dk_J_all(:,:),dk_f_all(:,:)
real*8, allocatable :: du_all(:,:),dFz_all(:,:)

integer :: i_iter,jj
character*100 :: file_sp,file_t
real*8 :: eps
real*8 :: sigma_SB_22_keV,help

! --- Iteration-control variables ---
real*8 :: lambda_T, lambda_bot
real*8 :: max_rel_flux_err, max_rel_dT, rel_flux_err
real*8 :: flux_surf, dm, gradF_scale
real*8 :: m_turn, p_weight, wloc
integer :: iter_max
real*8:: delta_surf,eps_Flog,eps_Fmax,eps_rough,delta_surf_,eps_Flog_,eps_Fmax_,eps_rough_

! --- variables for trapezoidal integration over energy ---
real*8 :: E_,E_prev,F_E_target_,F_E_target_prev,dE_
real*8 :: weights(0:2)
logical :: first_E
real*8 :: delta_Flog

  sigma_SB_22_keV = 1.027d2
  eps      = 1.d-2
  T_floor  = 0.1d0

  lambda_T   = 0.25d0
  lambda_bot = 0.15d0
  iter_max   = 1000

  m_turn   = 100.d0
  p_weight = 1.d0

  weights(0) = 1.d0
  weights(1) = 0.6d0
  weights(2) = 0.2d0

  200 format (100(es11.4,"   ") )

  ! === Physical parameters ===
  T_eff     = 0.6d0
  m_ns      = 1.4d0
  R6        = 1.d0
  B12       = 1.d2
  theta_B   = 1.d0
  Z         = 1.d0
  A         = 1.d0
  dot_m_6   = 0.d0
  ln_Lambda = 10.d0
  ! ===========================

  !== OpenMP parameters ==!
  n_threads = 2
  call omp_set_num_threads(n_threads)

  g14 = 1.328d0*m_ns/R6**2

  file_sp ="./res/res_B1e14_theta0_T2_sp"
  open(unit = 25, file = file_sp); close(25)
  file_t ="./res/res_B1e14_theta0_T2_t"
  open(unit = 26, file = file_t); close(26)

  !=== Numerical parameters ===!
  m_min = 1.d-2
  m_max = 1.d+4
  n_m   = 40
  n_mu  = 20
  n_fi  = 16
  n_E   = 80
  E_min = 0.01d0
  E_max = 20.d0
  !============================!

  write(*,*) "# B12=",B12,"  theta_B=",theta_B
  write(*,*) "# n_m=",n_m,"  n_mu=",n_mu,"  n_fi=",n_fi
  write(*,*) "# n_E=",n_E,"  E_min=",E_min,"  E_max=",E_max
  write(*,*) "# n_threads=",n_threads

  kappa_T = 0.4d0
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T

  allocate( mas_tau_TkeV(n_m,2), mas_m_rho(n_m,2), flux_tot(n_m,2), flux_E(n_m,2), J_E(n_m,2), J_tot(n_m,2), &
            kabs_mean(n_m,2),k_p(n_m),k_J(n_m),k_F(n_m),dk_p(n_m),dk_J(n_m),dk_F(n_m),du(n_m),dFz(n_m),u(n_m),u_p(n_m),Fz(n_m), &
            spec(n_E,3),spec_(n_E,3),dflux(n_m), dT(n_m), dT_sm(n_m), flux(n_m), gradF(n_m), weight_m(n_m), TkeV_old(n_m), &
            flux_E_prev(n_m,2), J_E_prev(n_m,2), kabs_mean_prev(n_m,2), &
            dk_p_prev(n_m), dk_J_prev(n_m), dk_f_prev(n_m), du_prev(n_m), dFz_prev(n_m) )

  allocate( E_grid(n_E),F_E_target_grid(n_E), &
            flux_E_all(n_E,n_m,2),J_E_all(n_E,n_m,2), &
            dk_p_all(n_E,n_m),dk_J_all(n_E,n_m),dk_f_all(n_E,n_m), &
            du_all(n_E,n_m),dFz_all(n_E,n_m) )

  ! === Initial temperature structure ===
  i = 1
  do while (i .le. n_m)
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = 0.5d0 + (mas_tau_TkeV(i,1)/ tau_max)**0.25d0  * T_eff / 2.d0
    i = i + 1
  end do
  TkeV_old(1:n_m) = mas_tau_TkeV(1:n_m,2)

  !=== Temperature iterations ===!
  delta_surf_= 1.d3
  eps_Flog_ =  1.d3
  eps_Fmax_ =  1.d3
  eps_rough_=  1.d3
  spec_(1:n_E,1) = 11.11d0; spec_(1:n_E,2) = 11.11d0; spec_(1:n_E,3) = 11.11d0

  i_iter = 0
  do while (i_iter .lt. iter_max)
    write(*,*); write(*,*) "# start iter:", i_iter

    ! --- Compute hydrostatic structure for the current temperature profile ---
    call get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)

    F_target = 0.d0
    flux_tot(1:n_m,1:2) = 0.d0
    J_tot(1:n_m,1:2) = 0.d0

    k_p(1:n_m) = 0.d0
    k_J(1:n_m) = 0.d0
    k_F(1:n_m) = 0.d0
    u(1:n_m) = 0.d0
    Fz(1:n_m) = 0.d0

    !=== Prepare energy grid ===!
    i = 1
    do while(i.le.n_E)
      E_grid(i) = E_min * (E_max/E_min)**( dble(i-1) / dble(n_E-1) )
      F_E_target_grid(i) = pi * BB_Intensity_22(E_grid(i),T_eff)
      i = i+1
    end do

    !=== Parallel radiative transfer over energy ===!
    !$omp parallel do default(shared) private(i,E_,flux_E,J_E,kabs_mean,dk_p,dk_J,dk_f,du,dFz) schedule(dynamic)
    do i = 1, n_E
      E_ = E_grid(i)
      call pol_RT_fixE(flux_E,J_E,kabs_mean,dk_p,dk_J,dk_f,du,dFz,&
                       E_,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda, &
                       mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
      spec(i,1) = E_
      spec(i,2) = flux_E(1,1)
      spec(i,3) = flux_E(1,2)

      flux_E_all(i,1:n_m,1:2) = flux_E(1:n_m,1:2)
      J_E_all(i,1:n_m,1:2)    = J_E(1:n_m,1:2)

      dk_p_all(i,1:n_m) = dk_p(1:n_m)
      dk_J_all(i,1:n_m) = dk_J(1:n_m)
      dk_f_all(i,1:n_m) = dk_f(1:n_m)

      du_all(i,1:n_m)  = du(1:n_m)
      dFz_all(i,1:n_m) = dFz(1:n_m)
    end do
    !$omp end parallel do

    !=== Trapezoidal integration over energy ===!
    i = 2
    do while(i.le.n_E)
      dE_ = E_grid(i) - E_grid(i-1)

      J_tot(1:n_m,1:2) = J_tot(1:n_m,1:2) + &
           0.5d0 * ( J_E_all(i-1,1:n_m,1:2) + J_E_all(i,1:n_m,1:2) ) * dE_
      k_p(1:n_m) = k_p(1:n_m) + &
           0.5d0 * ( dk_p_all(i-1,1:n_m) + dk_p_all(i,1:n_m) ) * dE_
      k_J(1:n_m) = k_J(1:n_m) + &
           0.5d0 * ( dk_J_all(i-1,1:n_m) + dk_J_all(i,1:n_m) ) * dE_
      k_F(1:n_m) = k_F(1:n_m) + &
           0.5d0 * ( dk_f_all(i-1,1:n_m) + dk_f_all(i,1:n_m) ) * dE_
      u(1:n_m) = u(1:n_m) + &
           0.5d0 * ( du_all(i-1,1:n_m) + du_all(i,1:n_m) ) * dE_
      Fz(1:n_m) = Fz(1:n_m) + &
           0.5d0 * ( dFz_all(i-1,1:n_m) + dFz_all(i,1:n_m) ) * dE_
      flux_tot(1:n_m,1:2) = flux_tot(1:n_m,1:2) + &
           0.5d0 * ( flux_E_all(i-1,1:n_m,1:2) + flux_E_all(i,1:n_m,1:2) ) * dE_
      F_target = F_target + &
           0.5d0 * ( F_E_target_grid(i-1) + F_E_target_grid(i) ) * dE_
      i = i+1
    end do

    i = 1
    do while( i.le.n_m )
      k_p(i) = k_p(i)/ (4.d0*100.d0 * BB_Flux24_interval(mas_tau_TkeV(i,2),E_min,E_max,n_E) )
      k_J(i) = k_J(i)/max(u(i),1.d-30)
      k_F(i) = k_F(i)*2.d0/max(Fz(i),1.d-30)
      u_p(i) = 1.333333d-8 * BB_Flux24_interval( mas_tau_TkeV(i,2), E_min, E_max, n_E )
      i = i+1
    end do

    !== now I want real u and u_p ==!
    u(1:n_m) = u(1:n_m)*2.d0/3.d10
    !u_p(1:n_m) = 1.37d-8*( mas_tau_TkeV(1:n_m,2) )**4
    !===============================!

    ! --- Build total flux profile ---
    flux(1:n_m)  = 0.d0
    dflux(1:n_m) = 0.d0
    gradF(1:n_m) = 0.d0

    i = 1
    do while (i .le. n_m)
      flux(i) = flux_tot(i,1) + flux_tot(i,2)
      i = i+1
    end do

    i = 2
    do while (i .le. n_m)
      dflux(i) = flux(i) - flux(i-1)
      dm = mas_tau_TkeV(i,1) - mas_tau_TkeV(i-1,1)
      gradF(i) = dflux(i) / max(dm,1.d-30)
      i = i + 1
    end do
    gradF(1) = gradF(2)

    call get_accur_metrics(delta_surf, eps_Fmax, delta_Flog, eps_Flog, &
                       eps_rough, F_target, Fz, mas_tau_TkeV, n_m)

    TkeV_old(1:n_m) = mas_tau_TkeV(1:n_m,2)

    call temperature_correction(dT,F_target,Fz,mas_tau_TkeV,k_F,k_J,k_P,u,u_p,n_m,E_min,E_max,n_E,delta_surf)

    call smooth_array(dT_sm, dT, n_m, 2, weights)

    do i = 1, n_m
      wloc = 1.d0
      if (i > n_m-6) wloc = 1.d0 - 0.25d0 * ( i-n_m+6 )/6.d0
      dT_sm(i) = sign(1.d0,dT_sm(i)) * &
                 min(abs(dT_sm(i)), wloc*0.05d0*mas_tau_TkeV(i,2))
    end do

    eps_Flog_ = eps_Flog
    spec_(1:n_E,1:3) = spec(1:n_E,1:3)

    !=== Print atmospheric structure for the current iteration ===!
    write(*,*) '# atmosphere structure, iter =', i_iter
    write(*,*) '# i, m, rho, T_keV, dT, F, F_target, J_tot(1:2)'

    i = 2
    do while (i .le. n_m)
      write(*,'(I5,1X,3F10.4,1X,ES12.4,1X,F10.4,1X,F12.4,1X,F12.4,1X,F12.4,1X,F12.4)') &
            i, mas_tau_TkeV(i,1), mas_m_rho(i,2), mas_tau_TkeV(i,2),  &
            dT_sm(i), flux(i), F_target, flux_tot(i,1:2)
      i = i + 1
    end do

    write(*,*)"# :  delta_surf, eps_Flog, eps_Fmax, eps_rough"
    write(*,*)"# ", delta_surf, eps_Flog, eps_Fmax, eps_rough

    !=== Update the temperature profile ===!
    i = 1
    do while (i .le. n_m)
      mas_tau_TkeV(i,2) = max(T_floor, mas_tau_TkeV(i,2) + dT_sm(i))
      i = i + 1
    end do

    open(unit = 25, file = file_sp, status = 'old', form='formatted')
    i = 1
    do while (i .le. n_E)
      write(25,*) spec_(i,1:3)
      i = i + 1
    end do
    close(25)

    open(unit = 26, file = file_t, status = 'old', form='formatted')
    i = 1
    do while (i .le. n_m)
      write(26,*) mas_tau_TkeV(i,1),TkeV_old(i)
      i = i + 1
    end do
    close(26)

    i_iter = i_iter + 1
  end do

144 return
end program MagAtm


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



!========================================================================================
! Fast symmetry test for opposite rays:
!   k(mu,phi)  <->  -k(-mu,phi+pi)
! Checks whether the transformed coordinates in B-frame satisfy:
!   cos(theta_B_opp) = -cos(theta_B)
!   phi_B_opp        = phi_B + pi  (mod 2pi)
! INPUT:
!   m_coord_b(j,k,1) = theta_B
!   m_coord_b(j,k,2) = phi_B
!========================================================================================
subroutine test_opposite_ray_symmetry(m_coord_b,n_mu,n_fi)
implicit none
integer,intent(in) :: n_mu,n_fi
real*8,intent(in) :: m_coord_b(n_mu,n_fi,2)
integer :: j,k,j_,k_
real*8 :: pi
real*8 :: th1,th2,ph1,ph2
real*8 :: err_mu,err_phi
real*8 :: max_err_mu,max_err_phi
real*8 :: dphi

pi = 3.141592653589793d0

max_err_mu  = 0.d0
max_err_phi = 0.d0

write(*,*)
write(*,*) '# ====================================================='
write(*,*) '# TEST: opposite-ray symmetry'
write(*,*) '# ====================================================='

do j = 1, n_mu

  j_ = n_mu + 1 - j

  do k = 1, n_fi

    !== opposite azimuth ==!
    k_ = k + n_fi/2
    k_ = mod(k_ - 1, n_fi) + 1

    th1 = m_coord_b(j ,k ,1)
    ph1 = m_coord_b(j ,k ,2)

    th2 = m_coord_b(j_,k_,1)
    ph2 = m_coord_b(j_,k_,2)

    !== check cos(theta_B) symmetry ==!
    err_mu = abs( cos(th2) + cos(th1) )

    !== check phi_B symmetry modulo 2pi ==!
    dphi = ph2 - ph1 - pi

    !== wrap to [-pi,+pi] ==!
    dphi = dphi - 2.d0*pi*nint(dphi/(2.d0*pi))

    err_phi = abs(dphi)

    max_err_mu  = max(max_err_mu ,err_mu )
    max_err_phi = max(max_err_phi,err_phi)

    if( (err_mu.gt.1.d-10).or.(err_phi.gt.1.d-10) )then
      write(*,'(A,4I5,2ES12.4)') &
      'BAD_PAIR: ',j,k,j_,k_,err_mu,err_phi
    end if

  end do
end do

write(*,*)
write(*,*) '# max error in cos(theta_B): ',max_err_mu
write(*,*) '# max error in phi_B       : ',max_err_phi
write(*,*) '# ====================================================='
write(*,*)

return
end subroutine test_opposite_ray_symmetry
