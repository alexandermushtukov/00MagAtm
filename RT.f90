!========================================================================================
!========================================================================================
subroutine set_atm_structure()
implicit none
real*8::pi=3.141592653589793d0
real*8,allocatable::m_atm_I(:,:,:,:),m_atm_T_rho(:,:)
real*8::B12,g14,theta_b,m_max
integer::n_m,n_mu,n_fi
real*8::dmu,dfi
  g14 = 1.d0
  B12 = 1.d2
  theta_b = 1.d0
  m_max = 1.d3
  n_m = 1000
  n_mu = 40
  n_fi = 100
  allocate( m_atm_I(n_m,n_mu,n_fi,1),m_atm_T_rho(n_m,2) )  !== intensity; T & rho ==!
  dmu = 2.d0/n_mu; dfi = pi/n_fi
return
end subroutine set_atm_structure
