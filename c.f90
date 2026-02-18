!======================================================================
!======================================================================
program MagAtm
use omp_lib
implicit none
real*8::pi=3.141592653589793d0
real*8::E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda    !== physical parameters ==!
real*8::m_min,m_max,kappa_T,tau_min,tau_max
integer::n_m,n_mu,n_fi,i
real*8::g14,flux_tot(2)
real*8,allocatable::mas_tau_TkeV(:,:),mas_m_rho(:,:)
integer::n_stream,i_task
  
  !!call test_absorption_mag_ff(); goto 144

  n_stream = 2

  !== physical paramters ==!
  E = 2.d0
  m_ns = 1.4d0
  R6 = 1.d0
  B12 = 2.d1
  theta_B = 0.d0
  Z = 1.d0
  A = 1.d0
  dot_m_6 = 0.d0
  ln_Lambda = 10.d0
  !==========================!
  g14 = 1.328d0*m_ns/R6**2

  !== numerical paramters ==!
  m_min = 1.d-2        !== the maximal column dencity [g/cm^2] ==!
  m_max = 1.d+4        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 200 !400
  n_mu = 8  !18
  n_fi = 1!2 !18
  !=======================!

  kappa_T = 0.34d0
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T

  allocate( mas_tau_TkeV(n_m,2),mas_m_rho(n_m,2) )

  !== set up some initial temperature sctructure ==!
  i = 1
  do while( i.le. n_m )
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = 0.5d0
    i = i+1
  end do

  call get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)

  E = 0.34d0
  do while(E.le.10.d0)
    call pol_RT_fixE(flux_tot,E,B12,g14,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
    write(*,*)E,flux_tot(1:2)
    E = E*1.2
  end do

  goto 144

  call omp_set_num_threads(n_stream)
  !$omp parallel default(none)&
  !$omp shared(B12,g14,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)&
  !$omp private(i_task,E,flux_tot)
    i_task = omp_get_thread_num()+1
    write(*,*)i_task
    E = 0.1d0 + 1.d0*(i_task-1)
    call pol_RT_fixE(flux_tot,E,B12,g14,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
    write(*,*)E,flux_tot(1:2)
  !$omp end parallel
144 return
end program MagAtm
