!======================================================================
!======================================================================
program MagAtm
use omp_lib
use black_body
implicit none
real*8::pi=3.141592653589793d0
real*8::E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda,T_eff    !== physical parameters ==!
real*8::m_min,m_max,kappa_T,tau_min,tau_max
integer::n_m,n_mu,n_fi,i,n_E
real*8::g14,E_min,E_max,dE,F_E_target,F_target,flux,alpha,T_old,T_corr
real*8,allocatable::mas_tau_TkeV(:,:),mas_m_rho(:,:),flux_tot(:,:),flux_E(:,:),spec(:,:)
integer::n_stream,i_task,i_iter
character*100::file_sp,file_t

  200 format (100(es11.4,"   ") )
  !!call test_absorption_mag_ff(); goto 144
  !!call test_cross_sections(); goto 144

  n_stream = 2

  !== physical paramters ==!
  T_eff = 1.d0
  m_ns = 1.4d0
  R6 = 1.d0
  B12 = 1.d2
  theta_B = 0.1d0
  Z = 1.d0
  A = 1.d0
  dot_m_6 = 0.d0
  ln_Lambda = 10.d0
  !==========================!
  g14 = 1.328d0*m_ns/R6**2

  file_sp ="./res/res_sp_6"
  open(unit = 25, file = file_sp); close(25)
  file_t ="./res/res_t_6"
  open(unit = 26, file = file_t); close(26)


  !== numerical paramters ==!
  m_min = 1.d-2        !== the maximal column dencity [g/cm^2] ==!
  m_max = 1.d+4        !== the maximal column dencity [g/cm^2] ==!
  n_m  = 200  !400
  n_mu = 20   !18
  n_fi = 1   !2 !18
  n_E = 20
  E_min = 0.1d0
  E_max = 7.d0
  !=======================!

  write(*,*)"# n_m=",n_m

  kappa_T = 0.34d0
  tau_min = m_min*kappa_T
  tau_max = m_max*kappa_T

  allocate( mas_tau_TkeV(n_m,2),mas_m_rho(n_m,2),flux_tot(n_m,2),flux_E(n_m,2),spec(n_E,3) )

  !== set up some initial temperature sctructure ==!
  i = 1
  do while( i.le. n_m )
    mas_tau_TkeV(i,1) = m_min * (m_max/m_min)**( dble(i-1) / dble(n_m-1) )
    mas_tau_TkeV(i,2) = T_eff
    i = i+1
  end do

  !== teperature iterrrations ==!
  i_iter = 0
  do while(i_iter.le.5)
    !== get hydro structure of the atmosphere ==!
    write(*,*)"# iter:",i_iter
    call get_hydro_atm_structure(mas_m_rho,g14,dot_m_6,ln_Lambda,tau_min,tau_max,n_m,mas_tau_TkeV)

    F_target = 0.d0
    flux_tot(1:n_m,1:2) = 0.d0
    i = 1
    do while(i.le.n_E)
      E = E_min * (E_max/E_min)**( real(i-1) / real(n_E-1) )
      F_E_target = pi * BB_Intensity_22(E,T_eff)
      dE  = E*( (E_max/E_min)**(1.0/(n_E-1)) - 1.0 )
      call pol_RT_fixE(flux_E,E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
      !write(*,200)E,flux_E(1,1:2)
      spec(i,1) = E
      spec(i,2) = flux_E(1,1)
      spec(i,3) = flux_E(1,2)
      flux_tot(1:n_m,1:2) = flux_tot(1:n_m,1:2) + flux_E(1:n_m,1:2)*dE
      F_target = F_target + F_E_target*dE
      !write(*,200)E,flux_tot(1:14,1)
      write(*,*)"# ",i,flux_E(1,1:2)
      i = i+1
    end do

    i = 1
    do while(i.le.n_m)
      flux = flux_tot(i,1)+flux_tot(i,2)
      if( flux.le.0.d0 )then
        flux = F_E_target/2
      end if

      alpha = 0.2d0   ! damping
      T_old  = mas_tau_TkeV(i,2)
      T_corr = mas_tau_TkeV(i,2) * ( F_target/ flux )**0.25
      !T_corr = T_old * (F_goal/flux)**0.25
       mas_tau_TkeV(i,2)  = (1d0-alpha)*T_old + alpha*T_corr

      !mas_tau_TkeV(i,2) = mas_tau_TkeV(i,2) * ( F_target/ flux )**0.25


      !write(*,*)i,flux_tot(i,1:2),F_E_target,mas_tau_TkeV(i,2)
      i = i+1
    end do
    i_iter = i_iter+1
  end do
  write(*,*)"# iteration ended"

  open(unit = 25, file = file_sp, status = 'old',form='formatted')
  i = 1
  do while(i.le.n_E)
    write(25,*)spec(i,1:3)
    i = i+1
  end do
  close(25)
  write(*,*)"# end0"
 
  open(unit = 26, file = file_t, status = 'old',form='formatted')
  i = 1
  do while(i.le.n_m)
    write(26,*)mas_tau_TkeV(i,1:2),flux_tot(i,1)+flux_tot(i,2),F_target
    i = i+1
  end do
  close(26)
  write(*,*)"# end1"

  goto 144

  call omp_set_num_threads(n_stream)
  !$omp parallel default(none)&
  !$omp shared(B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)&
  !$omp private(i_task,E,flux_tot)
    i_task = omp_get_thread_num()+1
    write(*,*)i_task
    E = 0.1d0 + 1.d0*(i_task-1)
    call pol_RT_fixE(flux_tot,E,B12,g14,T_eff,theta_B,Z,A,dot_m_6,ln_Lambda,mas_m_rho,mas_tau_TkeV,n_m,n_mu,n_fi)
    write(*,*)E,flux_tot(1,1:2)
  !$omp end parallel
144 return
end program MagAtm
