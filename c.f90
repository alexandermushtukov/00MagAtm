!======================================================================
!======================================================================
program MagAtm
implicit none
real*8::pi=3.141592653589793d0
real*8::E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda    !== physical parameters ==!
real*8::m_min,m_max,kappa_T,tau_min,tau_max
integer::n_m,n_mu,n_fi,i
real*8::g14
real*8,allocatable::mas_tau_TkeV(:,:),mas_m_rho(:,:)
  write(*,*)"#1"
  !== physical paramters ==!
  E = 2.d0
  m_ns = 1.4d0
  R6 = 1.d0
  B12 = 1.d2
  theta_B = 0.5d0
  Z = 1.d0
  A = 1.d0
  dot_m_6 = 0.d0
  ln_Lambda = 10.d0
  !==========================!
  g14 = 1.328d0*m_ns/R6**2

  !== numerical paramters ==!
  m_min = 1.d-2        !== the maximal column dencity [g/cm^2] ==!
  m_max = 1.d+3        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 40
  n_mu = 20;  n_fi = 20
  !=========================!

  kappa_T = 0.34d0
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T

  allocate( mas_tau_TkeV(n_m,2),mas_m_rho(n_m,2) )

  !== set up some temperature sctructure ==!
  i = 1
  do while( i.le. n_m )
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = 2.d0
    i = i+1
  end do

  call get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)
  call pol_RT_fixE(E,B12,g14,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
  write(*,*)"#2"
return
end program MagAtm
