!========================================================================================
!========================================================================================
subroutine set_atm_structure()
implicit none
real*8::B12,theta_b,m_max
integer::n_m,n_mu,n_fi
  B12 = 1.d2
  theta_b = 1.d0
  m_max = 1.d3
  n_m = 1000
  n_mu = 40
  n_fi = 100
return
end subroutine set_atm_structure
