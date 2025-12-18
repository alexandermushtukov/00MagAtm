!========================================================================================
!  m_atm_S(:,:,:,:,:,:) - scattering matrix, (n_m,n_mu,n_mu,2*n_fi,2,2)
!========================================================================================
subroutine set_atm_structure()
implicit none
real*8::pi=3.141592653589793d0
real*8,allocatable::m_atm_I(:,:,:,:),m_atm_T_rho(:,:)
complex*16,allocatable::m_atm_S(:,:,:,:,:,:)
real*8::B12,m_ns,R6,g14,theta_b,m_max,E,E_cyc
integer::n_m,n_mu,n_fi
real*8::dmu,dfi,mu_i,mu_f,delta_fi
integer::i,j,j1,j2,k
complex*16::Amp_dSigmadOmega !== functions ==!
  
  E = 2.d0

  m_ns = 1.4d0; R6 = 1.d0
  g14 = 1.328d0*m_ns/R6**2
  B12 = 1.d2          !== surface B-field strength ==!
  E_cyc = 11.4d0*B12  !== cyclotron energy in keV ==!
  theta_b = 1.d0      !== B-field inclination ==!
  m_max = 1.d3        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 1000
  n_mu = 40
  n_fi = 40
  allocate( m_atm_I(n_m,n_mu,n_fi,2),m_atm_T_rho(n_m,2),m_atm_S(n_m,n_mu,n_mu,2*n_fi,2,2) )  !== intensity; T & rho; S-matrix ==!
  dmu = 2.d0/n_mu; dfi = 2*pi/n_fi
  j1 = 1
  do while( j1 .le. n_mu )
    mu_i = dmu/2 + (j1-1)*dmu
    j2 = 1
    do while( j2 .le. n_mu )
      mu_f = dmu/2 + (j2-1)*dmu
      k = 1
      do while( k .le. 2*n_fi )
        delta_fi = dfi/2 + (k-1)*dfi
        m_atm_S(1:n_m,j1,j2,k,1,1) = Amp_dSigmadOmega(1,1,E,E_cyc,acos(mu_i),acos(mu_f),0.d0,delta_fi)
        m_atm_S(1:n_m,j1,j2,k,1,2) = Amp_dSigmadOmega(1,2,E,E_cyc,acos(mu_i),acos(mu_f),0.d0,delta_fi)
        m_atm_S(1:n_m,j1,j2,k,2,1) = Amp_dSigmadOmega(2,1,E,E_cyc,acos(mu_i),acos(mu_f),0.d0,delta_fi)
        m_atm_S(1:n_m,j1,j2,k,2,2) = Amp_dSigmadOmega(2,2,E,E_cyc,acos(mu_i),acos(mu_f),0.d0,delta_fi)
        k = k+1
      end do
      j2 = j2+1
    end do
    write(*,*)j1,n_mu
    j1 = j1+1
  end do
return
end subroutine set_atm_structure
