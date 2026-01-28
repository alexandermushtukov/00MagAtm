!======================================================================
!======================================================================
program MagAtm
implicit none
real*8::pi=3.141592653589793d0
real*8::E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda
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

  call pol_RT_fixE(E,B12,m_ns,R6,theta_B,Z,A,dot_m_6,ln_Lambda)
  write(*,*)"#2"
return
end program MagAtm
