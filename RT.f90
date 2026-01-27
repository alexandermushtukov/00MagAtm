!========================================================================================
!  m_atm_S(:,:,:,:,:,:) - scattering matrix, (n_m,n_mu,n_mu,2*n_fi,2,2)
!  mas_tau_TkeV(tau,T) - ...
!========================================================================================
subroutine set_atm_structure()
use abs_mag_ff_Meszaros
implicit none
real*8::pi=3.141592653589793d0
real*8,allocatable::m_atm_I(:,:,:,:),m_atm_T_rho(:,:),m_atm_sigma(:,:,:,:,:),m_atm_abs_b(:,:,:,:),mas_x_rho_tau(:,:),mas_tau_TkeV(:,:),F_tau(:,:),mas_m_rho(:,:)
real*8,allocatable::m_coord_b(:,:,:),KK_b(:,:,:),KK(:,:,:,:)
complex*16,allocatable::m_atm_S(:,:,:,:,:,:,:)
real*8::B12,m_ns,R6,g14,theta_b,m_min,m_max,E,E_cyc,dot_m_6,ln_Lambda
integer::n_m,n_mu,n_fi
real*8::dmu,dfi,mu_i,mu_f,theta_i,theta_f,theta_ib,theta_fb,fi_i,fi_f,fi_ib,fi_fb,n_ib(3),n_fb(3),n_i(3),n_f(3),Z,A,ksi_i,ksi_f
real*8::help
integer::i,j,k,j1,j2,k1,k2,jj
complex*16::Amp_dSigmadOmega,Amp_dSigmadOmega_ell_magnitars !== functions ==!
real*8::NormWavesEll,NormWavesEll_cvp  !==function==!

real*8::x_scale,kappa_T,tau_max,tau_min

  kappa_T = 0.34d0
  E = 2.d0

  !== physical paramters ==!
  m_ns = 1.4d0; R6 = 1.d0
  g14 = 1.328d0*m_ns/R6**2
  B12 = 1.d2            !== surface B-field strength ==!
  E_cyc = 11.4d0*B12    !== cyclotron energy in keV ==!
  theta_b = 1.d0        !== B-field inclination with respect to the normal [rad] ==!
  dot_m_6 = 0.d0
  ln_Lambda = 10.d0
  Z = 1.d0
  A = 1.d0
  !==========================!
 
  !== numerical paramters ==!
  m_min = 1.d-2        !== the maximal column dencity [g/cm^2] ==!
  m_max = 1.d+3        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 30
  n_mu = 15;  n_fi = 15
  !=========================!

  allocate( m_atm_I(n_m,n_mu,n_fi,2),m_atm_T_rho(n_m,2),m_atm_S(n_m,n_mu,n_mu,n_fi,n_fi,2,2),m_coord_b(n_mu,n_fi,2) )
         !== intensity; T & rho; S-matrix ==!
  allocate( m_atm_sigma(n_m,n_mu,n_fi,2,2),m_atm_abs_b(n_m,n_mu,2,2),mas_x_rho_tau(n_m,5),mas_tau_TkeV(n_m,2),F_tau(n_m,2),mas_m_rho(n_m,2) )
         !== free-free absorption and Compton ==!
  allocate( KK_b(n_m,n_mu,2),KK(n_m,n_mu,n_fi,2) )
  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi


  !== set up some temperature sctructure ==!
  i = 1
  do while( i.le. n_m )
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = 2.d0
    i = i+1
  end do

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
      !write(*,*)i,j,mu_i,theta_i,m_coord_b(i,j,1)
      j = j+1
    end do
    i = i+1
  end do

  !== get hydrostatical stracure of the atmosphere ==!
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T
  call acc_atm_structure_4(mas_x_rho_tau,n_m,tau_min,tau_max,x_scale, &
                           g14,mas_tau_TkeV,n_m,dot_m_6,ln_Lambda,0.d0,F_tau)
  mas_m_rho(1:n_m,1) = mas_x_rho_tau(1:n_m,3)/kappa_T
  mas_m_rho(1:n_m,2) = mas_x_rho_tau(1:n_m,2)
  !== now we have hydro structure of atmosphere ==!

  !== get ellipticities ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      mu_i = -1.d0 + dmu/2 + (i-1)*dmu
      theta_i = acos(mu_i)
      KK_b(i,j,1) = NormWavesEll_cvp(1,E,E_cyc,theta_i,mas_m_rho(i,2),Z,A)  !== X-mode ==!
      KK_b(i,j,2) = NormWavesEll_cvp(2,E,E_cyc,theta_i,mas_m_rho(i,2),Z,A)  !== O-mode ==!
      j = j+1
    end do
    i = i+1
  end do

  !== precalculate scatterings amps ==!
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
        do while( k1 .le. n_fi )
          fi_i = dfi/2 + (k1-1)*dfi
          theta_ib = m_coord_b(j1,k1,1)
          fi_ib = m_coord_b(j1,k1,2)
          k2 = 1
          do while( k2 .le. n_fi )
            fi_f = dfi/2 + (k2-1)*dfi
            theta_fb = m_coord_b(j2,k2,1)
            fi_fb = m_coord_b(j2,k2,2)
            m_atm_S(i,j1,j2,k1,k2,1,1) = Amp_dSigmadOmega_ell_magnitars(1,1,E,E_cyc,mas_m_rho(i,2),theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
            m_atm_S(i,j1,j2,k1,k2,1,2) = Amp_dSigmadOmega_ell_magnitars(1,2,E,E_cyc,mas_m_rho(i,2),theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
            m_atm_S(i,j1,j2,k1,k2,2,1) = Amp_dSigmadOmega_ell_magnitars(2,1,E,E_cyc,mas_m_rho(i,2),theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
            m_atm_S(i,j1,j2,k1,k2,2,2) = Amp_dSigmadOmega_ell_magnitars(2,2,E,E_cyc,mas_m_rho(i,2),theta_ib,theta_fb,fi_ib,fi_fb,Z,A,ksi_i,ksi_f)
            k2 = k2+1
          end do
          k1 = k1+1
        end do
        j2 = j2+1
      end do
      j1 = j1+1
    end do
    write(*,*)"# ",i,n_m
    i = i+1
  end do
  !== now we have amplitudes of Compton scattering ==!


  !== get absorption coeffisients in B-field RF ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      mu_i = -1.d0 + dmu/2 + (i-1)*dmu
      theta_i = acos(mu_i)
      m_atm_abs_b(i,j,1,1) = absorption_mag_ff_Meszaros_term( 1,E,E_cyc,mas_tau_TkeV(i,2),theta_i,Z,A,mas_m_rho(i,2) )
      m_atm_abs_b(i,j,1,2) = absorption_mag_ff_Meszaros_term( 2,E,E_cyc,mas_tau_TkeV(i,2),theta_i,Z,A,mas_m_rho(i,2) )  !== ?: does it account for rho? ==!
  
      !== compton scattering in B-RF ==!
      m_atm_abs_b(i,j,2,1:2) = 0.d0
      j2 = 1
      do while( j2.le.n_mu )
        k2 = 1
        do while( k2.le.n_fi )
          m_atm_abs_b(i,j,2,1:2) = m_atm_abs_b(i,j,2,1:2) + (m_atm_S(i,j,j2,1,k2,1:2,1) + m_atm_S(i,j,j2,1,k2,1:2,2))*dmu*dfi
          k2 = k2+1
        end do
        j2 = j2+1
      end do
      !================================!

      j = j+1
    end do
    i = i+1
  end do

  !== get absorption coefficients in atmospheric RF ==!
  i = 1
  do while( i.le.n_m )
    j = 1
    do while( j.le.n_mu )
      k = 1
      do while( k.le.n_fi )
        theta_ib = m_coord_b(j,k,1)
        jj = ( cos(theta_ib) + 1.d0 )/dmu + 1    !== fix if: improve adding interpolation ==!
        KK(i,j,k,1:2) = KK_b(i,jj,1:2)
        m_atm_sigma(i,j,k,1,1:2) = m_atm_abs_b(i,jj,1,1:2)   !== true absorption ==!
        m_atm_sigma(i,j,k,2,1:2) = m_atm_abs_b(i,jj,2,1:2)   !== absorption due to Compton pol 1  ==!
        !write(*,*)i,j,k,m_atm_sigma(i,j,k,1:2,1:2)
        !read(*,*)
        k = k+1
      end do
      j = j+1
    end do
    i = i+1
  end do
  !== now we have absorption coefficients ==!
122 return
end subroutine set_atm_structure




!============================================================================================
! ...
! We assume coherent scattering.
! See appendix B in 2204.12271
!============================================================================================
complex*16 function Amp_dSigmadOmega_ell_magnitars(s1,s2,E,Ecyc,rho,theta1,theta2,fi1,fi2,Z,A,ksi_i,ksi_f)
implicit none
integer,intent(in)::s1,s2
real*8,intent(in)::E,Ecyc,rho,theta1,theta2,fi1,fi2,Z,A,ksi_i,ksi_f
complex*16::ii=(0.d0,1.d0)
complex*16::g,f,amp
real*8::c_mod,f_gamma_n_v2_aver,NormWavesEll,NormWavesEll_cvp  !==function==!
real*8::E2,dkfdzi,y_f,z_f,Gamma,c_i,c_f,s_i,s_f,coeff
integer::det

  !!Gamma = f_gamma_n_v2_aver(0.d0,1,Ecyc/511)   !== Landau level width according to Pavlov [mc^2] ==!
  !!Gamma = Gamma*511      !== [mc^2]->[keV] ==!

  !ksi_i = NormWavesEll_cvp(1,E,Ecyc,theta1,rho,Z,A)
  !ksi_f = NormWavesEll_cvp(1,E,Ecyc,theta2,rho,Z,A)

  !g = E/(E+Ecyc+ii*Gamma)*exp(+ii*(fi1-fi2))  !== we exclude consideration of cyclotron line ==!   
  !f = E/(E-Ecyc+ii*Gamma)*exp(-ii*(fi1-fi2))
  g = E/(E+Ecyc)*exp(+ii*(fi1-fi2))
  f = E/(E-Ecyc)*exp(-ii*(fi1-fi2))

  c_i=cos(theta1);  c_f=cos(theta2)
  s_i=sin(theta1);  s_f=sin(theta2)
  coeff=1.d0/sqrt(1.d0+ksi_i**2)/sqrt(1.d0+ksi_f**2)

  if((s1.eq.1).and.(s2.eq.1))then
    amp=2*ksi_f*ksi_i*s_f*s_i + g*(1.d0-ksi_i*c_i)*(1.d0-ksi_f*c_f) + f*(1.d0+ksi_i*c_i)*(1.d0+ksi_f*c_f)
  end if
  if((s1.eq.1).and.(s2.eq.2))then
    amp=2*ksi_i*s_f*s_i + g*(c_f+ksi_f)*(ksi_i*c_i-1.d0) + f*(c_f-ksi_f)*(ksi_i*c_i+1.d0)
  end if
  if((s1.eq.2).and.(s2.eq.1))then
    amp=2*ksi_f*s_f*s_i + g*(c_i+ksi_i)*(ksi_f*c_f-1.d0) + f*(c_i-ksi_i)*(ksi_f*c_f+1.d0)
  end if
  if((s1.eq.2).and.(s2.eq.2))then
    amp=2*s_i*s_f + g*(ksi_i+c_i)*(ksi_f+c_f) + f*(ksi_i-c_i)*(ksi_f-c_f)
  end if
  Amp_dSigmadOmega_ell_magnitars=amp*coeff

return
end function Amp_dSigmadOmega_ell_magnitars
