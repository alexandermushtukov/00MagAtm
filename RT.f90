!============= ToDo ==========================!
! - flux at the boundaries due to nearest layer emission 
! - make temperature correction
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
real*8::x_scale,kappa_T = 0.34d0
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
!=======================================================================================================
subroutine pol_RT_fixE(flux_tot,E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
implicit none
real*8,intent(out)::flux_tot(n_m,2)
real*8::pi=3.141592653589793d0
real*8,intent(in)::E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho(n_m,2),mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m,n_mu,n_fi
real*8::S_therm(n_m,n_mu,n_fi,2),S_0(n_m,n_mu,n_fi,2),R(n_m,n_mu,n_mu,n_fi,n_fi,2,2),m_atm_kappa(n_m,n_mu,n_fi,2,2)
real*8::S(n_m,n_mu,n_fi,2),I_out(n_mu,n_fi,2),I_out_tot(n_mu,n_fi,2),I_e(n_m,n_mu,n_fi,2)
real*8::dmu,dfi,mu,theta
integer::i,k,j

  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi

  call set_atm_coefficients(S_therm,R,m_atm_kappa,E,B12,g14,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)

  I_out_tot(1:n_mu,1:n_fi,1:2) = 0.d0
  !== start iterrations ==!
  S_0(1:n_m,1:n_mu,1:n_fi,1:2) = S_therm(1:n_m,1:n_mu,1:n_fi,1:2)
  i=1
  do while(i.le.1)
    flux_tot(1:n_m,1:2) = 0.d0
    call RT_iterrations(I_e,S,E,S_0,T_eff,mas_tau_TkeV(n_m,2),R,m_atm_kappa,mas_m_rho,n_m,n_mu,n_fi)
    k = 1
    do while(k.le.n_fi)
     j = 1
      do while(j.le.n_mu)
        mu = -1.d0 + dmu/2 + (j-1)*dmu; theta = acos(mu)
        flux_tot(1:n_m,1:2) = flux_tot(1:n_m,1:2) + I_e(1:n_m,j,k,1:2)*mu * dmu*dfi  !sin(theta)
        j = j+1
      end do
      !write(*,*)
      k = k+1
    end do
    S_0(1:n_m,1:n_mu,1:n_fi,1:2) = S_therm(1:n_m,1:n_mu,1:n_fi,1:2) + S(1:n_m,1:n_mu,1:n_fi,1:2)
    I_out_tot(1:n_mu,1:n_fi,1:2) = I_out(1:n_mu,1:n_fi,1:2)
    i = i+1
  end do

  !== printing ==!
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
subroutine set_atm_coefficients(S_0,R,m_atm_kappa,E,B12,g14,theta_b,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
use black_body
implicit none
real*8,intent(out)::S_0(n_m,n_mu,n_fi,2),R(n_m,n_mu,n_mu,n_fi,n_fi,2,2),m_atm_kappa(n_m,n_mu,n_fi,2,2)
real*8,intent(in)::E,B12,g14,theta_b,Z,A,dot_m_6,ln_Lambda,mas_m_rho(n_m,2),mas_tau_TkeV(n_m,2)
real*8::pi=3.141592653589793d0
real*8::m_atm_abs_b(n_m,n_mu,2,2)
real*8::m_coord_b(n_mu,n_fi,2),KK_b(n_m,n_mu,2),KK(n_m,n_mu,n_fi,2),R_b(n_m,n_mu,n_mu,n_fi,2,2)
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

  kappa_T = 0.34d0
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
        !jj = max(1, min(n_mu, int((cos(theta_ib)+1.d0)/dmu) + 1))

        jj = nint( (cos(theta_ib) + 1.d0 - dmu/2.d0) / dmu ) + 1
        jj = max(1, min(n_mu, jj))

        fi_ib = m_coord_b(j,k,2)
        KK(i,j,k,1:2) = KK_b(i,jj,1:2)
        m_atm_kappa(i,j,k,1,1:2) = m_atm_abs_b(i,jj,1,1:2) * kappa_T   !== true absorption ==!
        m_atm_kappa(i,j,k,2,1:2) = m_atm_abs_b(i,jj,2,1:2) * kappa_T   !== absorption due to Compton  ==!

        S_0(i,j,k,1:2) = BB_Intensity_22(E,mas_tau_TkeV(i,2))/2 * m_atm_kappa(i,j,k,1,1:2)   !== initial thermal source function [kappa*B_22] ==!
        j2 = 1
        do while(j2.le.n_mu)
          k2 = 1
          do while( k2.le.n_fi )
            theta_fb = m_coord_b(j2,k2,1)
            !jj2 = max(1, min(n_mu, int((cos(theta_fb)+1.d0)/dmu) + 1))

            jj2 = nint( (cos(theta_fb) + 1.d0 - dmu/2.d0) / dmu ) + 1
            jj2 = max(1, min(n_mu, jj2))

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
            R(i,j,j2,k,k2,1:2,1:2) = (1.d0-frac)*R_b(i,jj,jj2,k0,1:2,1:2) + frac*R_b(i,jj,jj2,k1,1:2,1:2)

            k2 = k2 + 1
          end do
          j2 = j2+1
        end do
        !write(*,*)i,j,k,m_atm_kappa(i,j,k,1:2,1:2)
        !read(*,*)
        k = k+1
      end do
      j = j+1
    end do
    !m_atm_kappa(i,1:n_mu,1:n_fi,1:2,1:2) = m_atm_kappa(i,1:n_mu,1:n_fi,1:2,1:2)  !== it is opcity [cm^2/g] ==!
    !write(*,*)"## ",i,mas_m_rho(i,1:2),KK_b(i,4,1:2)
    i = i+1
  end do
122 return
end subroutine set_atm_coefficients



!==========================================================================================================================
! ...
!==========================================================================================================================
subroutine RT_iterrations(I_e,S,E,S_0,T_eff,T_bottom,R,m_atm_kappa,mas_m_rho,n_m,n_mu,n_fi)
use black_body
implicit none
real*8,intent(out)::S(n_m,n_mu,n_fi,2),I_e(n_m,n_mu,n_fi,2)
integer,intent(in)::n_m,n_mu,n_fi
real*8,intent(in)::E,S_0(n_m,n_mu,n_fi,2),T_eff,T_bottom,R(n_m,n_mu,n_mu,n_fi,n_fi,2,2),m_atm_kappa(n_m,n_mu,n_fi,2,2),mas_m_rho(n_m,2)
integer::i,j,k,i1,i2,ii,i_pol,jj,kk
real*8::pi=3.141592653589793d0
real*8::dmu,dfi,mu,dSigma,tau,dtau,tau_lim,kappa
real*8::kappa_T = 0.34d0
real*8::F_E_target

  F_E_target = pi * BB_Intensity_22(E,T_eff)   !== used later to calculate emission from the bottom ==!

  tau_lim = 50.d0   !== the maximal optical distance b/w points ==!
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
            ii = i+1       !== accounting for layers below ==!
            do while( (ii.le.n_m).and.(tau.lt.tau_lim) )
              dSigma = mas_m_rho(ii,1) - mas_m_rho(ii-1,1)                       !== colomn density of ii-layer ==!
              kappa = m_atm_kappa(ii,j,k,1,i_pol) + m_atm_kappa(ii,j,k,2,i_pol)  !== total opacity due to abs and scattering ==!
              dtau = dSigma * kappa / abs(mu)
              if( dtau.ne.0.d0 )then
                if( ii.ne.(i+1) )then
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol) * dSigma /abs(mu) * ( 1.d0 - exp(-dtau) )/dtau * exp(-tau)
                else
                  !== we are at the layer boundary ==!
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol) * dSigma /abs(mu) * ( 1.d0 - exp(-dtau) )/dtau
                end if
              end if
              !== check coefficient: ( 1.d0 - exp(-dtau) )/dtau it should be fraction of radiation that is created in a layer and leave it ==!
              tau = tau + dtau
              ii = ii+1
            end do
            !== add intensity from the lower boundary ==!
            if( ii.ge.n_m )then
              !== we still see the bottom of the atmosphere ==!
              I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + BB_Intensity_22(E,T_bottom)/2 * exp(-tau) + 3*mu*(F_E_target/2)/4/pi  * exp(-tau)
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
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol) * dSigma /abs(mu) * ( 1.d0 - exp(-dtau) )/dtau * exp(-tau)
                else
                  !== we are at the boundary of a layer ==!
                  I_e(i,j,k,i_pol) = I_e(i,j,k,i_pol) + S_0(ii,j,k,i_pol) * dSigma /abs(mu) * ( 1.d0 - exp(-dtau) )/dtau
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

  !== get new souse function, in units kappa*B_22 ==!
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
            S(i,j,k,1) = S(i,j,k,1) + ( R(i,jj,j,kk,k,1,1)*I_e(i,jj,kk,1) + R(i,jj,j,kk,k,2,1)*I_e(i,jj,kk,2) )*kappa_T * dmu*dfi
            S(i,j,k,2) = S(i,j,k,2) + ( R(i,jj,j,kk,k,1,2)*I_e(i,jj,kk,1) + R(i,jj,j,kk,k,2,2)*I_e(i,jj,kk,2) )*kappa_T * dmu*dfi
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
    tau2 = tau2 + 0.34d0 * rho2 * dx
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
!============================================================================================

