!========================================================================================
!  m_atm_S(:,:,:,:,:,:) - scattering matrix, (n_m,n_mu,n_mu,2*n_fi,2,2)
!  mas_tau_TkeV(tau,T) - ...
!========================================================================================
subroutine set_atm_structure()
use abs_mag_ff_Meszaros
implicit none
real*8::pi=3.141592653589793d0
real*8,allocatable::m_atm_I(:,:,:,:),m_atm_T_rho(:,:),m_atm_sigma(:,:,:,:,:),mas_x_rho_tau(:,:),mas_tau_TkeV(:,:),F_tau(:,:),mas_m_rho(:,:)
complex*16,allocatable::m_atm_S(:,:,:,:,:,:)
real*8::B12,m_ns,R6,g14,theta_b,m_min,m_max,E,E_cyc,dot_m_6,ln_Lambda
integer::n_m,n_mu,n_fi
real*8::dmu,dfi,mu_i,mu_f,theta_i,theta_f,theta_ib,theta_fb,fi_i,fi_f,fi_ib,fi_fb,n_ib(3),n_fb(3),n_i(3),n_f(3)
real*8::help
integer::i,j,j1,j2,k1,k2
complex*16::Amp_dSigmadOmega !== functions ==!

real*8::x_scale,kappa_T,tau_max,tau_min

  kappa_T = 0.34d0
  E = 2.d0

  !== physical paramters ==!
  m_ns = 1.4d0; R6 = 1.d0
  g14 = 1.328d0*m_ns/R6**2
  B12 = 1.d2            !== surface B-field strength ==!
  E_cyc = 11.4d0*B12    !== cyclotron energy in keV ==!
  theta_b = 1.d0        !== B-field inclination ==!
  dot_m_6 = 0.d0
  ln_Lambda = 10.d0
  !==========================!
 
  !== numerical paramters ==!
  m_min = 1.d-2        !== the maximal column dencity [g/cm^2] ==!
  m_max = 1.d+2        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 1000
  n_mu = 20;  n_fi = 20
  !=========================!

  allocate( m_atm_I(n_m,n_mu,n_fi,2),m_atm_T_rho(n_m,2),m_atm_S(n_mu,n_mu,n_fi,n_fi,2,2) )  !== intensity; T & rho; S-matrix ==!
  allocate( m_atm_sigma(n_m,n_mu,n_fi,2,2),mas_x_rho_tau(n_m,5),mas_tau_TkeV(n_m,2),F_tau(n_m,2),mas_m_rho(n_m,2) )   !== free-free absorption and Compton ==!
  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi

  !== get hydro stracure of atmosphere ==!
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T
  call acc_atm_structure_4(mas_x_rho_tau,n_m,tau_min,tau_max,x_scale, &
                           g14,mas_tau_TkeV,n_m,dot_m_6,ln_Lambda,0.d0,F_tau)
  mas_m_rho(1:n_m,1) = mas_x_rho_tau(1:n_m,3)/kappa_T
  mas_m_rho(1:n_m,2) = mas_x_rho_tau(1:n_m,2)
  !i=1
  !do while(i.le.n_m)
  !  write(*,*)mas_m_rho(i,1:2); i=i+1
  !end do
  !== now we have hydro structure of atmosphere ==!

  j1 = 1
  do while( j1 .le. n_mu )
    mu_i = dmu/2 + (j1-1)*dmu; theta_i = acos(mu_i)
    j2 = 1
    do while( j2 .le. n_mu )
      mu_f = dmu/2 + (j2-1)*dmu; theta_f = acos(mu_f)
      k1 = 1
      do while( k1 .le. n_fi )
        fi_i = dfi/2 + (k1-1)*dfi
        call spherical2cartisian( 1.d0,theta_i,fi_i,n_i(1),n_i(2),n_i(3) )
        call VecRotation3d( n_i(1),n_i(2),n_i(3), 1, -theta_b , n_ib(1),n_ib(2),n_ib(3))
        call cartesian2spherical( n_ib(1),n_ib(2),n_ib(3) , help,theta_ib,fi_ib )
        k2 = 1
        do while( k2 .le. n_fi )
          fi_f = dfi/2 + (k2-1)*dfi
          call spherical2cartisian(1.d0,theta_f,fi_f,n_f(1),n_f(2),n_f(3))
          call VecRotation3d( n_f(1),n_f(2),n_f(3), 1, -theta_b , n_fb(1),n_fb(2),n_fb(3))
          call cartesian2spherical( n_fb(1),n_fb(2),n_fb(3) , help,theta_fb,fi_fb )
          m_atm_S(j1,j2,k1,k2,1,1) = Amp_dSigmadOmega(1,1,E,E_cyc,theta_ib,theta_fb,fi_ib,fi_fb)
          m_atm_S(j1,j2,k1,k2,1,2) = Amp_dSigmadOmega(1,2,E,E_cyc,theta_ib,theta_fb,fi_ib,fi_fb)
          m_atm_S(j1,j2,k1,k2,2,1) = Amp_dSigmadOmega(2,1,E,E_cyc,theta_ib,theta_fb,fi_ib,fi_fb)
          m_atm_S(j1,j2,k1,k2,2,2) = Amp_dSigmadOmega(2,2,E,E_cyc,theta_ib,theta_fb,fi_ib,fi_fb)
          !write(*,*)m_atm_S(4,j1,j2,k1,k2,1:2,1:2); read(*,*)
          k2 = k2+1
        end do
        k1 = k1+1
      end do
      j2 = j2+1
    end do
    write(*,*)"# ",j1,n_mu
    j1 = j1+1
  end do
  !== now we have amplitudes of Compton scattering ==!
return
end subroutine set_atm_structure
