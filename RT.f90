!============= ToDo ==========================!
! - update temperature correction
! - proton resonance in absorption
!=============================================

!=========================================================================================================
! ...
!=========================================================================================================
subroutine get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)
implicit none
real*8,intent(out)::mas_m_rho(n_m,2)
real*8,intent(in)::g14,dot_m_6,ln_Lambda,tau_min,tau_max,mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m
real*8::x_scale,kappa_T = 0.4d0
real*8::F_tau(n_m,2),mas_x_rho_tau(n_m,5),help
integer::i
  !== get hydrostatical stracure of the atmosphere ==!
  call acc_atm_structure_5_(mas_x_rho_tau,n_m,tau_min,tau_max,x_scale, &
                           g14,mas_tau_TkeV,n_m,dot_m_6,ln_Lambda,0.d0,F_tau)
  mas_m_rho(1:n_m,1) = mas_x_rho_tau(1:n_m,3)/kappa_T         !== colomn density coordinate ==!
  mas_m_rho(1:n_m,2) = mas_x_rho_tau(1:n_m,2)                 !== local mass density ==!

  !== check atmosphere structure ==!
  !i=1; help = 0.d0
  !do while(i.le.n_m)
  !  if( i.gt.1 )then
  !  help = help + (mas_x_rho_tau(i,1)-mas_x_rho_tau(i-1,1))*(mas_m_rho(i,2)+mas_m_rho(i-1,2))/2
  !  end if
  !  write(*,*)mas_m_rho(i,1:2),help
  !  i=i+1
  !end do
  !read(*,*)
return
end subroutine get_hydro_atm_structure



!=======================================================================================================
! ...
!   m_atm_kappa(n_m,n_mu,n_fi,1,1:2) = m_atm_abs_b(i,jj,1,1:2) * kappa_T   !== true absorption ==!
!   m_atm_kappa(n_m,n_mu,n_fi,2,1:2) = m_atm_abs_b(i,jj,2,1:2) * kappa_T   !== absorption due to Compton  ==!
!=======================================================================================================
subroutine pol_RT_fixE(flux_tot,J_E,kabs_mean,dk_p,dk_J,dk_f,du,dFz,&
                       E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,&
                       n_m,n_mu,n_fi)
use black_body
implicit none
real*8,intent(out)::flux_tot(n_m,2),J_E(n_m,2),kabs_mean(n_m,2),dk_p(n_m),dk_J(n_m),dk_f(n_m),du(n_m),dFz(n_m)
real*8::pi=3.141592653589793d0
real*8,intent(in)::E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho(n_m,2),mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m,n_mu,n_fi
real*8::S_therm(n_m,n_mu,n_fi,2),S_0(n_m,n_mu,n_fi,2),S_0_new(n_m,n_mu,n_fi,2),m_atm_kappa(n_m,n_mu,n_fi,2,2)
real*8::S(n_m,n_mu,n_fi,2),I_out(n_mu,n_fi,2),I_out_tot(n_mu,n_fi,2),m_coord_b(n_mu,n_fi,2)
real*8::dmu,dfi,mu,theta,R_b(n_m,n_mu,n_mu,n_fi,2,2),I_e(n_m,n_mu,n_fi,2)
integer::i,k,j,i_max,j_,k_
real*8::eps_S,tol_S,small,help

  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi
  small = 1.d-30
  tol_S = 1.d-3
  i_max = 30   !== maximal number of interractions ==!

  !====================================================================================!
  ! get S_therm - source function due to thermal emission [kappa*B_22]
  !     R_b - redistrubution function [1/ster]
  !     m_atm_kappa - map of opaity [cm^2/g]
  call set_atm_coefficients(S_therm,R_b,m_atm_kappa,m_coord_b,E,B12,g14,theta_B,Z,A,&
                            dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
  !=====================================================================================!

  !== start iterrations ==!
  S_0(1:n_m,1:n_mu,1:n_fi,1:2) = S_therm(1:n_m,1:n_mu,1:n_fi,1:2)    !== [kappa*B_22] : the 1st assumption about source function ==!

  !== iterrations of radiative transfer: start ==!
  i=1
  do while(i.le.i_max)
    I_out_tot(1:n_mu,1:n_fi,1:2) = 0.d0
    flux_tot(1:n_m,1:2) = 0.d0
    J_E(1:n_m,1:2) = 0.d0
    kabs_mean(1:n_m,1:2) = 0.d0

    call RT_iterrations_v2(I_e,S,E,S_0,T_eff,R_b,m_atm_kappa,mas_m_rho,n_m,n_mu,n_fi,m_coord_b,&
                           mas_tau_TkeV(n_m,2),mas_tau_TkeV(n_m-1,2))

    !== integration over (4*pi) ==!
    k = 1
    do while(k.le.n_fi)
      j = 1
      do while(j.le.n_mu)
        mu = -1.d0 + dmu/2 + (j-1)*dmu; theta = acos(mu)
        flux_tot(1:n_m,1:2) = flux_tot(1:n_m,1:2) + I_e(1:n_m,j,k,1:2)*mu * dmu*dfi
        J_E(1:n_m,1:2) = J_E(1:n_m,1:2) + I_e(1:n_m,j,k,1:2) * dmu*dfi
        kabs_mean(1:n_m,1:2) = kabs_mean(1:n_m,1:2) + m_atm_kappa(1:n_m,j,k,1,1:2) * dmu*dfi
        j = j+1
      end do
      k = k+1
    end do
    J_E(1:n_m,1:2) = J_E(1:n_m,1:2) / (4*pi)
    kabs_mean(1:n_m,1:2) = kabs_mean(1:n_m,1:2) / (4*pi)

    !== get new source function and check convergence ==!
    S_0_new(1:n_m,1:n_mu,1:n_fi,1:2) = S_therm(1:n_m,1:n_mu,1:n_fi,1:2) + S(1:n_m,1:n_mu,1:n_fi,1:2)

    eps_S = 0.d0
    help = maxval( abs( S_0_new(1:n_m,1:n_mu,1:n_fi,1:2) - S_0(1:n_m,1:n_mu,1:n_fi,1:2) ) / &
                   max( abs(S_0_new(1:n_m,1:n_mu,1:n_fi,1:2)), small ) )
    eps_S = help

    S_0(1:n_m,1:n_mu,1:n_fi,1:2) = S_0_new(1:n_m,1:n_mu,1:n_fi,1:2)
    I_out_tot(1:n_mu,1:n_fi,1:2) = I_out(1:n_mu,1:n_fi,1:2)
    if( eps_S.lt.tol_S )exit
    i = i+1
  end do
  !== iterrations of radiative transfer: end ==!

  !== get dk_ ==!
  dk_p(1:n_m) = 0.d0; dk_J(1:n_m) = 0.d0; dk_f(1:n_m) = 0.d0
  du(1:n_m) = 0.d0
  dFz(1:n_m) = 0.d0
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      mu = -1.d0 + dmu/2 + (j-1)*dmu
      k = 1
      do while( k.le.n_fi )
        dk_p(i) = dk_p(i) + dmu*dfi* ( m_atm_kappa(i,j,k,1,1) + m_atm_kappa(i,j,k,1,2) ) * BB_Intensity_22(E,mas_tau_TkeV(i,2))!/2
        dk_J(i) = dk_J(i) + dmu*dfi* ( m_atm_kappa(i,j,k,1,1)*I_e(i,j,k,1) + m_atm_kappa(i,j,k,1,2)*I_e(i,j,k,2) )
        !dk_f(i) = dk_f(i) + dmu*dfi/mu/2*sign(1.d0,mu)* ( (m_atm_kappa(i,j,k,1,1)+m_atm_kappa(i,j,k,2,1))*I_e(i,j,k,1) &
        !                                                  + (m_atm_kappa(i,j,k,1,2)+m_atm_kappa(i,j,k,2,2))*I_e(i,j,k,2) )
        !dk_f(i) = dk_f(i) + dmu*dfi/mu/2* ( (m_atm_kappa(i,j,k,1,1)+m_atm_kappa(i,j,k,2,1))*I_e(i,j,k,1) &
        !                                                   + (m_atm_kappa(i,j,k,1,2)+m_atm_kappa(i,j,k,2,2))*I_e(i,j,k,2) )
        if( mu.gt.0.d0 )then
          j_ = n_mu + 1 - j
          k_ = k + n_fi/2; k_ = mod(k_ - 1, n_fi) + 1
          dk_f(i) = dk_f(i) + dmu*dfi/mu/2* &
                   ( (m_atm_kappa(i,j,k,1,1)+m_atm_kappa(i,j,k,2,1))*I_e(i,j,k,1) &
                      + (m_atm_kappa(i,j,k,1,2)+m_atm_kappa(i,j,k,2,2))*I_e(i,j,k,2) &
                      - (m_atm_kappa(i,j_,k_,1,1)+m_atm_kappa(i,j_,k_,2,1))*I_e(i,j_,k_,1) &
                      - (m_atm_kappa(i,j_,k_,1,2)+m_atm_kappa(i,j_,k_,2,2))*I_e(i,j_,k_,2) )
        end if
        du(i)  = du(i) + dmu*dfi*( I_e(i,j,k,1) + I_e(i,j,k,2) )
        dFz(i) = dFz(i) + dmu*dfi* mu*( I_e(i,j,k,1) + I_e(i,j,k,2) )    !== flux in [erg/cm^2] ==!
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
  !===================================!

  !== printing if needed ==!
  k = 1
  do while(k.le.n_fi)
    j = 1
    do while(j.le.n_mu)
      !write(*,*)k*dfi,j*dmu,I_out(j,k,1:2)
      j = j+1
    end do
    !write(*,*)
    k = k+1
  end do
return
end subroutine pol_RT_fixE


!==============================================================================================================================
!  Subroutine pre-calculates souurce function in the atmosphere and scattering redistributoon function.
!    S_0 - source function [kappa*B_22]
!    R - redistrubution function [1/ster]
!    m_atm_kappa - map of opaity [cm^2/g]
!  mas_tau_TkeV(tau,T) - ...
!==============================================================================================================================
subroutine set_atm_coefficients(S_0,R_b,m_atm_kappa,m_coord_b,E,B12,g14,theta_b,Z,A,&
                                dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
use black_body
implicit none
real*8,intent(out)::S_0(n_m,n_mu,n_fi,2),m_atm_kappa(n_m,n_mu,n_fi,2,2),m_coord_b(n_mu,n_fi,2)
real*8,intent(in)::E,B12,g14,theta_b,Z,A,dot_m_6,ln_Lambda,mas_m_rho(n_m,2),mas_tau_TkeV(n_m,2)
real*8::pi=3.141592653589793d0
real*8::m_atm_abs_b(n_m,n_mu,2,2)
real*8::KK_b(n_m,n_mu,2),KK(n_m,n_mu,n_fi,2),R_b(n_m,n_mu,n_mu,n_fi,2,2)
complex*16::m_atm_S(n_m,n_mu,n_mu,n_fi,2,2),m_atm_S_p(n_m,n_mu,n_mu,n_fi,2,2)
real*8::E_cyc
integer::n_m,n_mu,n_fi
real*8::dmu,dfi,mu_i,mu_f,theta_i,theta_f,theta_ib,theta_fb,fi_i,fi_f,fi_ib,fi_fb,n_ib(3),n_fb(3),n_i(3),n_f(3),ksi_i,ksi_f
real*8::help,delta_fi
integer::i,j,k,j1,j2,k1,k2,jj,jj2,kk2
complex*16::Amp_dSigmadOmega_ell_magnitars,Amp_dSigmadOmega_ell_magnitars_p !== functions ==!
real*8::NormWavesEll,NormWavesEll_cvp,NormWavesEll_cvp_,abs_mag_ff_Meszaros_new  !==function==!
real*8::x_scale,kappa_T
real*8 :: xk,frac
integer :: k0

  kappa_T = 0.4d0
  E_cyc = 11.4d0*B12    !== cyclotron energy in keV ==!
 
  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi

  !== array of coordinate transformation ==!
  i = 1
  do while( i.le.n_mu )
    mu_i = -1.d0 + dmu/2 + (i-1)*dmu; theta_i = acos(mu_i)
    j = 1
    do while( j.le.n_fi )
      fi_i = dfi/2 + (j-1)*dfi
      call spherical2cartisian( 1.d0,theta_i,fi_i,n_i(1),n_i(2),n_i(3) )
      call VecRotation3d( n_i(1),n_i(2),n_i(3), 1, -theta_b , n_ib(1),n_ib(2),n_ib(3))
      call cartesian2spherical( n_ib(1),n_ib(2),n_ib(3) , help,theta_ib,fi_ib )
      m_coord_b(i,j,1) = theta_ib
      m_coord_b(i,j,2) = fi_ib
      j = j+1
    end do
    i = i+1
  end do

  !== get ellipticities ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      mu_i = -1.d0 + dmu/2 + (j-1)*dmu
      theta_i = acos(mu_i)
      KK_b(i,j,1) = NormWavesEll_cvp_(2,E,E_cyc,theta_i,mas_m_rho(i,2),Z,A)  !== (+)-mode ==!
      KK_b(i,j,2) = NormWavesEll_cvp_(1,E,E_cyc,theta_i,mas_m_rho(i,2),Z,A)  !== (-)-mode ==!
      !if( isnan(KK_b(i,j,2)) )then
      !  write(*,*)"@",KK_b(i,j,2),E,E_cyc,theta_i,mas_m_rho(i,2),Z,A
      !  read(*,*)
      !end if
      j = j+1
    end do
    i = i+1
  end do

  !== pre-calculate scatterings amps in B-fiesl RF ==!
  i=1
  do while(i.le.n_m)
    j1 = 1
    do while( j1 .le. n_mu )
      mu_i = -1.d0 + dmu/2 + (j1-1)*dmu; theta_i = acos(mu_i)
      j2 = 1
      do while( j2 .le. n_mu )
        mu_f = -1.d0 + dmu/2 + (j2-1)*dmu; theta_f = acos(mu_f)
        ksi_i = KK_b(i,j1,1)
        ksi_f = KK_b(i,j2,1)
        k1 = 1
        fi_i = dfi/2 + (k1-1)*dfi
        theta_ib = m_coord_b(j1,k1,1)
        fi_ib = m_coord_b(j1,k1,2)
        k2 = 1
        do while( k2 .le. n_fi )
          fi_f = dfi/2 + (k2-1)*dfi
          theta_fb = m_coord_b(j2,k2,1)
          fi_fb = m_coord_b(j2,k2,2)
          m_atm_S(i,j1,j2,k2,1,1) = Amp_dSigmadOmega_ell_magnitars(1,1,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S(i,j1,j2,k2,1,2) = Amp_dSigmadOmega_ell_magnitars(1,2,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S(i,j1,j2,k2,2,1) = Amp_dSigmadOmega_ell_magnitars(2,1,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S(i,j1,j2,k2,2,2) = Amp_dSigmadOmega_ell_magnitars(2,2,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)

          m_atm_S_p(i,j1,j2,k2,1,1) = Amp_dSigmadOmega_ell_magnitars_p(1,1,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S_p(i,j1,j2,k2,1,2) = Amp_dSigmadOmega_ell_magnitars_p(1,2,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S_p(i,j1,j2,k2,2,1) = Amp_dSigmadOmega_ell_magnitars_p(2,1,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          m_atm_S_p(i,j1,j2,k2,2,2) = Amp_dSigmadOmega_ell_magnitars_p(2,2,E,E_cyc,mas_m_rho(i,2),&
                                                                   theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
          k2 = k2+1
        end do
        j2 = j2+1
      end do
      j1 = j1+1
    end do
    !write(*,*)"# ",i,n_m
    i = i+1
  end do
  !write(*,*)"#done: scatterings amps"
  !== now we have amplitudes of Compton scattering ==!

  !== get absorption coeffisients in B-field RF ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      mu_i = -1.d0 + dmu/2 + (j-1)*dmu
      theta_i = acos(mu_i)
      !== free-free absorption: cross section in units of sigma_T ==!
      m_atm_abs_b(i,j,1,1) = abs_mag_ff_Meszaros_new( KK_b(i,j,1),E,E_cyc,mas_tau_TkeV(i,2),theta_i,Z,A,mas_m_rho(i,2) )
      m_atm_abs_b(i,j,1,2) = abs_mag_ff_Meszaros_new( KK_b(i,j,2),E,E_cyc,mas_tau_TkeV(i,2),theta_i,Z,A,mas_m_rho(i,2) )  !== ?: does it account for rho? ==!
      !== compton scattering in B-RF: cross section in units of sigma_T ==!
      m_atm_abs_b(i,j,2,1:2) = 0.d0
      j2 = 1
      do while( j2.le.n_mu )
        k2 = 1
        do while( k2.le.n_fi )
          m_atm_abs_b(i,j,2,1:2) = m_atm_abs_b(i,j,2,1:2) &
                               + ( (abs(m_atm_S(i,j,j2,k2,1:2,1)))**2 + (abs(m_atm_S(i,j,j2,k2,1:2,2)))**2 &
                                   + (abs(m_atm_S_p(i,j,j2,k2,1:2,1)))**2 + (abs(m_atm_S_p(i,j,j2,k2,1:2,2)))**2 )*dmu*dfi
          R_b(i,j,j2,k2,1:2,1:2) = (abs(m_atm_S(i,j,j2,k2,1:2,1:2)))**2 * 3/32/pi &
                                 + (abs(m_atm_S_p(i,j,j2,k2,1:2,1:2)))**2 * 3/32/pi !== scattering redistribution function in B-field RF [1/ster] ==!
          k2 = k2+1
        end do
        j2 = j2+1
      end do
      m_atm_abs_b(i,j,2,1:2) = m_atm_abs_b(i,j,2,1:2) * 3/32/pi
      !================================!
      j = j+1
    end do
    i = i+1
  end do
  !write(*,*)"#done: get absorption coeffisients in B-field RF"

  !== get absorption coefficients & sccattering redistribution function in atmospheric RF ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      k = 1
      do while( k.le.n_fi )
        theta_ib = m_coord_b(j,k,1)

        jj = nint( (cos(theta_ib) + 1.d0 - dmu/2) / dmu ) + 1
        jj = max(1, min(n_mu, jj))

        fi_ib = m_coord_b(j,k,2)
        KK(i,j,k,1:2) = KK_b(i,jj,1:2)
        m_atm_kappa(i,j,k,1,1:2) = m_atm_abs_b(i,jj,1,1:2) * kappa_T   !== true absorption ==!
        m_atm_kappa(i,j,k,2,1:2) = m_atm_abs_b(i,jj,2,1:2) * kappa_T   !== absorption due to Compton  ==!

        S_0(i,j,k,1:2) = BB_Intensity_22(E,mas_tau_TkeV(i,2))/2 * m_atm_kappa(i,j,k,1,1:2)   !== initial thermal source function [kappa*B_22] ==!
        !j2 = 1
        !do while(j2.le.n_mu)
        !  k2 = 1
        !  do while( k2.le.n_fi )
        !    theta_fb = m_coord_b(j2,k2,1)
        !    jj2 = nint( (cos(theta_fb) + 1.d0 - dmu/2) / dmu ) + 1
        !    jj2 = max(1, min(n_mu, jj2))!

        !    fi_fb = m_coord_b(j2,k2,2)

        !    delta_fi = fi_fb - fi_ib
        !    ! wrap to [0,2pi)num_tot
        !    delta_fi = delta_fi - 2*pi*floor(delta_fi/(2*pi))
        !    ! (эквивалент if, но работает и для >2pi)
        !    xk   = delta_fi/dfi           ! in [0, n_fi)
        !    k0   = int(xk) + 1            ! 1..n_fi
        !    frac = xk - dble(k0-1)        ! 0..1
        !    k1 = k0 + 1
        !    if (k1.gt.n_fi) k1 = 1        ! periodicity
        !    R(i,j,j2,k,k2,1:2,1:2) = (1.d0-frac)*R_b(i,jj,jj2,k0,1:2,1:2) + frac*R_b(i,jj,jj2,k1,1:2,1:2)!

        !    k2 = k2 + 1
        !  end do
        !  j2 = j2+1
        !end do
        !write(*,*)i,j,k,m_atm_kappa(i,j,k,1:2,1:2)
        !read(*,*)
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
122 return
end subroutine set_atm_coefficients



!==========================================================================================================================
! Solve Radiative Transfer for given source function S_0.
! Output:
!   [I_e] = [B22] - map of intensities
!   [S] = [kappa*B22] - new source function due to scatterings.
! ...
!==========================================================================================================================
subroutine RT_iterrations(I_e,S,E,S_0,T_eff,R_b,m_atm_kappa,mas_m_rho,n_m,n_mu,n_fi,m_coord_b,Tbot,Tprev)
use black_body
implicit none
real*8,intent(out)::S(n_m,n_mu,n_fi,2),I_e(n_m,n_mu,n_fi,2)
integer,intent(in)::n_m,n_mu,n_fi
real*8,intent(in)::E,S_0(n_m,n_mu,n_fi,2),T_eff,m_atm_kappa(n_m,n_mu,n_fi,2,2),mas_m_rho(n_m,2),Tbot,Tprev
real*8,intent(in)::m_coord_b(n_mu,n_fi,2),R_b(n_m,n_mu,n_mu,n_fi,2,2)
integer::i,j,k,i1,i2,ii,i_pol,jj,kk
real*8::pi=3.141592653589793d0
real*8::dmu,dfi,mu,dSigma,tau,dtau,tau_lim,kappa
real*8::kappa_T = 0.4d0, kappa_help
real*8::RR(2,2),frac
integer::q1,q2,qq1,qq2
real*8 :: Bbot, Bprev, dBdm, dBdTauNu, Ibot, dm_bot

  tau_lim = 5.d4   !== the maximal optical distance b/w points ==!
  dmu = 2.d0/n_mu
  dfi = 2*pi/n_fi

  !== get intensity map ==!
  i = 1
  do while(i.le.n_m)
    j = 1
    do while( j.le.n_mu)
      mu = -1.d0 + dmu/2 + (j-1)*dmu
      k = 1
      do while( k.le.n_fi )
        I_e(i,j,k,1:2) = 0.d0
        i_pol = 1
        !== calculate two polarisations separately ==!
        do while(i_pol.le.2)
          tau = 0.d0
          if( mu.ge.0.d0 )then
            !== upward propagation: accounting for underling layers ==!
            ii = i+1       !== accounting for layers below: ii\in [ii+1,n_m] ==!
            do while( (ii.le.n_m).and.(tau.lt.tau_lim) )
              dSigma = mas_m_rho(ii,1) - mas_m_rho(ii-1,1)                       !== colomn density of ii-layer ==!
              kappa = m_atm_kappa(ii,j,k,1,i_pol) + m_atm_kappa(ii,j,k,2,i_pol)  !== total opacity: abs + scattering ==!
              dtau = dSigma * kappa / abs(mu)
              if( dtau.ne.0.d0 )then
                if( ii.ne.(i+1) )then
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) ) * exp(-tau)
                else
                  !== we are at the layer boundary ==!
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) )
                end if
              end if
              tau = tau + dtau
              ii = ii+1
            end do
            !!== add intensity from the lower boundary ==!
            if( ii.ge.n_m )then
              !=== Diffusion-type lower boundary condition, analogous to eq. (2.10) ===!
              Bbot  = BB_Intensity_22(E, Tbot)
              Bprev = BB_Intensity_22(E, Tprev)
              dm_bot = max( mas_m_rho(n_m,1) - mas_m_rho(n_m-1,1), 1.d-30 )
              dBdm = (Bbot - Bprev) / dm_bot
              kappa = m_atm_kappa(n_m,j,k,1,i_pol) + m_atm_kappa(n_m,j,k,2,i_pol)
              kappa = max(kappa, 1.d-30)
              dBdTauNu = dBdm / kappa
              Ibot = 0.5d0 * ( Bbot + mu * dBdTauNu )
              !== optional positivity safeguard ==!
              Ibot = max(0.d0, Ibot)
              I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + Ibot * exp(-tau)
            end if
          else
            !== downward propagation: accounting for upper layers ==!
            ii = i      !== accounting for upper layers starting from the current one ==!
            do while( (ii.ge.1).and.(tau.lt.tau_lim) )
              if( ii.ne.1 )then
                dSigma = mas_m_rho(ii,1) - mas_m_rho(ii-1,1)
              else
                dSigma = mas_m_rho(ii,1)
              end if
              kappa = m_atm_kappa(ii,j,k,1,i_pol) + m_atm_kappa(ii,j,k,2,i_pol)  !== total opacity due to abs and scattering ==!
              dtau = dSigma * kappa / abs(mu)
              if( dtau.ne.0.d0 )then
                if( ii.ne.i )then
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) ) * exp(-tau)
                else
                  !== we are at the boundary of a layer ==!
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) )
                end if
              end if
              !== check coefficient: ( 1.d0 - exp(-dtau) )/dtau it should be fraction of radiation that is created in a layer and leave it ==!
              tau = tau + dtau
              ii = ii-1
            end do
          end if
          i_pol = i_pol+1
        end do
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
  !== getting new source function ==!

  !== get new souse function, in units [kappa*B_22] ==!
  S(1:n_m,1:n_mu,1:n_fi,1:2) = 0.d0
  i = 1
  do while(i.le.n_m)
    j = 1
    do while(j.le.n_mu)
      k = 1
      do while( k.le.n_fi )
        !== integration over (4\pi) ==!
        jj = 1
        do while(jj.le.n_mu)
          kk = 1
          do while(kk.le.n_fi)

            call get_index_for_R(q1,q2,qq1,qq2,frac,jj,j,kk,k,m_coord_b,n_mu,n_fi,dmu,dfi)
            RR(1:2,1:2) = (1.d0-frac)*R_b(i,q1,q2,qq1,1:2,1:2) + frac*R_b(i,q1,q2,qq2,1:2,1:2)
            S(i,j,k,1) = S(i,j,k,1) + ( RR(1,1)*I_e(i,jj,kk,1)*m_atm_kappa(i,jj,kk,2,1) &
                                      + RR(2,1)*I_e(i,jj,kk,2)*m_atm_kappa(i,jj,kk,2,2) ) * dmu*dfi
            S(i,j,k,2) = S(i,j,k,2) + ( RR(1,2)*I_e(i,jj,kk,1)*m_atm_kappa(i,jj,kk,2,1) &
                                      + RR(2,2)*I_e(i,jj,kk,2)*m_atm_kappa(i,jj,kk,2,2) ) * dmu*dfi

            kk = kk+1
          end do
          jj = jj+1
        end do
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
  !== now we have updated source function ==!
return
end subroutine RT_iterrations



!==========================================================================================================================
! Solve Radiative Transfer for given source function S_0.
! Output:
!   [I_e] = [B22] - map of intensities
!   [S] = [kappa*B22] - new source function due to scatterings.
! ...
!   V2: accelerated version
!       1) intensity is built by sweeps over depth
!       2) angular interpolation indices for R are precomputed once
!==========================================================================================================================
subroutine RT_iterrations_v2(I_e,S,E,S_0,T_eff,R_b,m_atm_kappa,mas_m_rho,n_m,n_mu,n_fi,m_coord_b,Tbot,Tprev)
use black_body
implicit none
real*8,intent(out)::S(n_m,n_mu,n_fi,2),I_e(n_m,n_mu,n_fi,2)
integer,intent(in)::n_m,n_mu,n_fi
real*8,intent(in)::E,S_0(n_m,n_mu,n_fi,2),T_eff,m_atm_kappa(n_m,n_mu,n_fi,2,2),mas_m_rho(n_m,2),Tbot,Tprev
real*8,intent(in)::m_coord_b(n_mu,n_fi,2),R_b(n_m,n_mu,n_mu,n_fi,2,2)

integer::i,j,k,ii,i_pol,jj,kk
real*8::pi=3.141592653589793d0
real*8::dmu,dfi,mu,dSigma,dtau,kappa
real*8::kappa_T = 0.4d0
real*8::RR(2,2),frac
integer::q1,q2,qq1,qq2
real*8 :: Bbot, Bprev, dBdm, dBdTauNu, Ibot, dm_bot

!== arrays for accelerated intensity sweeps ==!
real*8::A_up(n_m),tr_up(n_m),I_up(n_m)
real*8::A_dn(n_m),tr_dn(n_m),I_dn(n_m)

!== arrays for precomputed angular mapping ==!
integer::map_q1(n_mu,n_mu,n_fi,n_fi),map_q2(n_mu,n_mu,n_fi,n_fi)
integer::map_qq1(n_mu,n_mu,n_fi,n_fi),map_qq2(n_mu,n_mu,n_fi,n_fi)
real*8 :: map_frac(n_mu,n_mu,n_fi,n_fi)

  dmu = 2.d0/n_mu
  dfi = 2*pi/n_fi

  !== precompute angular map for R interpolation ==!
  jj = 1
  do while(jj.le.n_mu)
    j = 1
    do while(j.le.n_mu)
      kk = 1
      do while(kk.le.n_fi)
        k = 1
        do while(k.le.n_fi)
          call get_index_for_R(q1,q2,qq1,qq2,frac,jj,j,kk,k,m_coord_b,n_mu,n_fi,dmu,dfi)
          map_q1(jj,j,kk,k) = q1
          map_q2(jj,j,kk,k) = q2
          map_qq1(jj,j,kk,k) = qq1
          map_qq2(jj,j,kk,k) = qq2
          map_frac(jj,j,kk,k) = frac
          k = k+1
        end do
        kk = kk+1
      end do
      j = j+1
    end do
    jj = jj+1
  end do

  !== get intensity map ==!
  j = 1
  do while( j.le.n_mu )
    mu = -1.d0 + dmu/2 + (j-1)*dmu
    k = 1
    do while( k.le.n_fi )
      i_pol = 1
      do while(i_pol.le.2)

        if( mu.ge.0.d0 )then
          !== upward propagation: sweep from bottom to top ==!
          A_up(1:n_m) = 0.d0
          tr_up(1:n_m) = 0.d0
          I_up(1:n_m) = 0.d0

          ii = 2
          do while(ii.le.n_m)
            dSigma = mas_m_rho(ii,1) - mas_m_rho(ii-1,1)
            kappa = m_atm_kappa(ii,j,k,1,i_pol) + m_atm_kappa(ii,j,k,2,i_pol)
            dtau = dSigma * kappa / abs(mu)

            if( dtau.ne.0.d0 )then
              A_up(ii) = S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) )
              tr_up(ii) = exp(-dtau)
            else
              A_up(ii) = 0.d0
              tr_up(ii) = 1.d0
            end if

            ii = ii+1
          end do

          !== Diffusion-type lower boundary condition, analogous to (2.10) ==!
          Bbot  = 0.5d0 * BB_Intensity_22(E, Tbot)
          Bprev = 0.5d0 * BB_Intensity_22(E, Tprev)
          dm_bot = max( mas_m_rho(n_m,1) - mas_m_rho(n_m-1,1), 1.d-30 )
          dBdm = (Bbot - Bprev) / dm_bot

          !kappa = m_atm_kappa(n_m,j,k,1,i_pol) + m_atm_kappa(n_m,j,k,2,i_pol)   !== check ==!
          kappa = 0.5d0 * (m_atm_kappa(n_m,j,k,1,i_pol)   + m_atm_kappa(n_m,j,k,2,i_pol) &
                         + m_atm_kappa(n_m-1,j,k,1,i_pol) + m_atm_kappa(n_m-1,j,k,2,i_pol) )
          kappa = max(kappa, 1.d-30)
          dBdTauNu = dBdm / kappa

          Ibot =  Bbot + mu * dBdTauNu
          !Ibot = max(0.d0, Ibot)   !== check: I can be negative ==!

          I_up(n_m) = Ibot
          i = n_m-1
          do while(i.ge.1)
            I_up(i) = A_up(i+1) + tr_up(i+1) * I_up(i+1)
            i = i-1
          end do

          i = 1
          do while(i.le.n_m)
            I_e(i,j,k,i_pol) = I_up(i)
            i = i+1
          end do

        else
          !== downward propagation: sweep from top to bottom ==!
          A_dn(1:n_m) = 0.d0
          tr_dn(1:n_m) = 0.d0
          I_dn(1:n_m) = 0.d0

          ii = 1
          do while(ii.le.n_m)
            if( ii.ne.1 )then
              dSigma = mas_m_rho(ii,1) - mas_m_rho(ii-1,1)
            else
              dSigma = mas_m_rho(ii,1)
            end if

            kappa = m_atm_kappa(ii,j,k,1,i_pol) + m_atm_kappa(ii,j,k,2,i_pol)
            dtau = dSigma * kappa / abs(mu)

            if( dtau.ne.0.d0 )then
              A_dn(ii) = S_0(ii,j,k,i_pol)/kappa * ( 1.d0 - exp(-dtau) )
              tr_dn(ii) = exp(-dtau)
            else
              A_dn(ii) = 0.d0
              tr_dn(ii) = 1.d0
            end if

            ii = ii+1
          end do

          I_dn(1) = A_dn(1)
          i = 2
          do while(i.le.n_m)
            I_dn(i) = A_dn(i) + tr_dn(i) * I_dn(i-1)
            i = i+1
          end do

          i = 1
          do while(i.le.n_m)
            I_e(i,j,k,i_pol) = I_dn(i)
            i = i+1
          end do

        end if

        i_pol = i_pol+1
      end do
      k = k+1
    end do
    j = j+1
  end do
  !== getting new source function ==!

  !== get new souse function, in units [kappa*B_22] ==!
  S(1:n_m,1:n_mu,1:n_fi,1:2) = 0.d0
  i = 1
  do while(i.le.n_m)
    j = 1
    do while(j.le.n_mu)
      k = 1
      do while( k.le.n_fi )
        !== integration over (4\pi) ==!
        jj = 1
        do while(jj.le.n_mu)
          kk = 1
          do while(kk.le.n_fi)

            q1 = map_q1(jj,j,kk,k)
            q2 = map_q2(jj,j,kk,k)
            qq1 = map_qq1(jj,j,kk,k)
            qq2 = map_qq2(jj,j,kk,k)
            frac = map_frac(jj,j,kk,k)

            RR(1:2,1:2) = (1.d0-frac)*R_b(i,q1,q2,qq1,1:2,1:2) + frac*R_b(i,q1,q2,qq2,1:2,1:2)

            S(i,j,k,1) = S(i,j,k,1) + ( RR(1,1)*I_e(i,jj,kk,1)*m_atm_kappa(i,jj,kk,2,1) &
                                      + RR(2,1)*I_e(i,jj,kk,2)*m_atm_kappa(i,jj,kk,2,2) ) * dmu*dfi
            S(i,j,k,2) = S(i,j,k,2) + ( RR(1,2)*I_e(i,jj,kk,1)*m_atm_kappa(i,jj,kk,2,1) &
                                      + RR(2,2)*I_e(i,jj,kk,2)*m_atm_kappa(i,jj,kk,2,2) ) * dmu*dfi

            kk = kk+1
          end do
          jj = jj+1
        end do
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
  !== now we have updated source function ==!
return
end subroutine RT_iterrations_v2


!================================================================================================
! ...
!   R(i,j,j2,k,k2,1:2,1:2) requested
!   R(i,j,j2,k,k2,1:2,1:2) = (1.d0-frac)*R_b(i,jj,jj2,k0,1:2,1:2) + frac*R_b(i,jj,jj2,k1,1:2,1:2)
!================================================================================================
subroutine get_index_for_R(jj,jj2,k0,k1,frac,j,j2,k,k2,m_coord_b,n_mu,n_fi,dmu,dfi)
implicit none
integer,intent(out)::jj,jj2,k0,k1
real*8,intent(out)::frac
real*8,intent(in)::m_coord_b(n_mu,n_fi,2),dmu,dfi
integer,intent(in)::j,j2,k,k2,n_mu,n_fi
real*8::theta_ib,theta_fb,fi_ib,fi_fb,delta_fi,xk
real*8::pi=3.141592653589793d0

  theta_ib = m_coord_b(j,k,1)
  jj = nint( (cos(theta_ib) + 1.d0 - dmu/2) / dmu ) + 1;   jj = max(1, min(n_mu, jj))
  fi_ib = m_coord_b(j,k,2)

  theta_fb = m_coord_b(j2,k2,1)
  jj2 = nint( (cos(theta_fb) + 1.d0 - dmu/2) / dmu ) + 1;  jj2 = max(1, min(n_mu, jj2))

  fi_fb = m_coord_b(j2,k2,2)
  delta_fi = fi_fb - fi_ib
  ! wrap to [0,2pi)
  delta_fi = delta_fi - 2*pi*floor(delta_fi/(2*pi))

  ! (эквивалент if, но работает и для >2pi)
  xk   = delta_fi/dfi           ! in [0, n_fi)
  k0   = int(xk) + 1            ! 1..n_fi
  frac = xk - dble(k0-1)        ! 0..1
  k1 = k0 + 1
  if (k1.gt.n_fi) k1 = 1        ! periodicity

return
end subroutine get_index_for_R



!============================================================================================
! OUTPUT GRID (LOWER-BOUNDARY / NO-DUPLICATE VERSION):
! We output values at the LOWER boundary of each layer, on a log-spaced tau-grid:
!   tau_lower(i) = logspace(tau_min .. tau_max),  i=1..n
! No tau=0 row is written and there is NO duplicated last boundary.
! Therefore:
!   mas_x_rho_tau_2(i,1) = x      at tau = tau_lower(i)
!   mas_x_rho_tau_2(i,2) = rho2   at tau = tau_lower(i)
!   mas_x_rho_tau_2(i,3) = tau_lower(i)
!   mas_x_rho_tau_2(i,4) = eps24  at tau = tau_lower(i)
!   mas_x_rho_tau_2(i,5) = T_keV_ at tau = tau_lower(i)
! Requires tau_min > 0 (if n>1).
!============================================================================================
subroutine acc_atm_structure_5_(mas_x_rho_tau_2,n,tau_min,tau_max,x_scale, &
                              g14,mas_tau_TkeV,n_mas,dot_m_6,ln_Lambda,rho0,F_tau)
implicit none
integer, intent(in)  :: n, n_mas
real*8,  intent(in)  :: tau_min, tau_max, g14, dot_m_6, ln_Lambda, rho0
real*8,  intent(in)  :: mas_tau_TkeV(n_mas,2)

real*8,  intent(out) :: mas_x_rho_tau_2(n,5), x_scale
real*8,  intent(out) :: F_tau(n,2)

! local
integer :: i, steps, max_steps
real*8  :: dx, x, x_old
real*8  :: tau2, tau_old, target_tau, f
real*8  :: rho2, rho_old, drho2
real*8  :: beta, beta_ff, dbeta_dx, Z
real*8  :: T_keV, T_keV_, T_old
real*8  :: eps24, eps_old
real*8  :: rho_floor, beta_floor, T_floor
real*8  :: tau_lower(n), log_ratio, expo

real*8  :: find_H_ordered_inc
external :: find_H_ordered_inc

  rho_floor  = 1.d-30
  beta_floor = 1.d-30
  T_floor    = 1.d-30

  beta_ff = 0.5d0
  Z       = 1.d0

  dx   = 1.d-3
  x    = 0.d0
  tau2 = 0.d0

  rho2 = max(rho0, rho_floor)
  beta = beta_ff

  eps24  = 0.d0
  T_keV_ = 0.d0

  !---- build output tau grid directly: logspace(tau_min..tau_max), size n, no duplicates
  if (n .eq. 1) then
    tau_lower(1) = tau_max
  else
    log_ratio = log(tau_max / tau_min)
    do i = 1, n
      expo = dble(i-1) / dble(n-1)   ! 0..1
      tau_lower(i) = tau_min * exp( expo * log_ratio )
    end do
    tau_lower(n) = tau_max
  end if

  !---- We WRITE rows at tau = tau_lower(i)
  i = 1
  target_tau = tau_lower(i)

  steps = 0
  max_steps = 200000000  ! safety cap

  do while (i .le. n)
    if (steps .ge. max_steps) exit
    steps = steps + 1

    ! save old state for interpolation
    x_old   = x
    tau_old = tau2
    rho_old = rho2
    eps_old = eps24
    T_old   = T_keV_

    ! temperature from tau-table at current tau (before step)
    T_keV = find_H_ordered_inc(tau_old, mas_tau_TkeV, n_mas)
    if (T_keV .le. T_floor) T_keV = T_floor

    ! evolve beta
    if (beta .gt. 0.d0) then
      dbeta_dx = 3.24d-4 * rho2 * Z**2 / max(beta,beta_floor)**3 * ln_Lambda
      beta = max(0.d0, beta - dx*dbeta_dx)
    else
      dbeta_dx = 0.d0
      beta = 0.d0
    end if

    ! evolve rho2
    if (beta .gt. 0.d0) then
      drho2 = ( 0.5d0*0.1d0*g14/T_keV*rho2 + 31.15d0*dot_m_6*0.5d0/T_keV*dbeta_dx ) * dx
    else
      drho2 = ( 0.5d0*0.1d0*g14/T_keV*rho2 ) * dx
    end if
    rho2 = max(rho_floor, rho2 + drho2)

    ! local source + proxy temperature
    if (beta .gt. 0.d0) then
      eps24  = 9.d2 * dot_m_6 * beta * dbeta_dx
      T_keV_ = ( eps24 / 1.75d0 / max(rho2,rho_floor)**2 )**2
    else
      eps24  = 0.d0
      T_keV_ = 0.d0
    end if

    ! evolve optical depth and x
    tau2 = tau2 + 0.4d0 * rho2 * dx
    x    = x + dx

    ! write any crossed targets (may cross several in one dx step)
    do while ( (tau2 .ge. target_tau) .and. (i .le. n) )
      if (tau2 .gt. tau_old) then
        f = (target_tau - tau_old) / (tau2 - tau_old)
      else
        f = 0.d0
      end if

      mas_x_rho_tau_2(i,1) = x_old + f*dx
      mas_x_rho_tau_2(i,2) = rho_old + f*(rho2 - rho_old)
      mas_x_rho_tau_2(i,3) = target_tau
      mas_x_rho_tau_2(i,4) = eps_old + f*(eps24 - eps_old)
      mas_x_rho_tau_2(i,5) = T_old   + f*(T_keV_ - T_old)

      i = i + 1
      if (i .le. n) then
        target_tau = tau_lower(i)
      end if
    end do
  end do

  ! pad if stopped early
  if (i .le. n) then
    do while (i .le. n)
      mas_x_rho_tau_2(i,1) = x
      mas_x_rho_tau_2(i,2) = rho2
      mas_x_rho_tau_2(i,3) = tau_lower(i)
      mas_x_rho_tau_2(i,4) = eps24
      mas_x_rho_tau_2(i,5) = T_keV_
      i = i + 1
    end do
  end if

  x_scale = mas_x_rho_tau_2(n,1) - mas_x_rho_tau_2(1,1)

  ! cumulative source distribution vs tau (normalized)
  F_tau(1,1) = mas_x_rho_tau_2(1,3)
  F_tau(1,2) = 0.d0
  do i = 2, n
    F_tau(i,1) = mas_x_rho_tau_2(i,3)
    F_tau(i,2) = F_tau(i-1,2) + mas_x_rho_tau_2(i,4) * (mas_x_rho_tau_2(i,1) - mas_x_rho_tau_2(i-1,1))
  end do
  if (F_tau(n,2) .gt. 0.d0) F_tau(1:n,2) = F_tau(1:n,2) / F_tau(n,2)

return
end subroutine acc_atm_structure_5_
!=========================================================================================
