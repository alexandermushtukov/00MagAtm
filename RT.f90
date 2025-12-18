!========================================================================================
!========================================================================================
subroutine set_atm_structure()
implicit none
real*8::pi=3.141592653589793d0
real*8,allocatable::m_atm_I(:,:,:,:),m_atm_T_rho(:,:)
complex*16,allocatable::m_atm_S(:,:,:,:,:)
real*8::B12,m_ns,R6,g14,theta_b,m_max
integer::n_m,n_mu,n_fi
real*8::dmu,dfi
integer::i,j,j1,j2,k
  m_ns = 1.4d0; R6 = 1.d0
  g14 = 1.328d0*m_ns/R6**2
  B12 = 1.d2       !== surface B-field strength ==!
  theta_b = 1.d0   !== B-field inclination ==!
  m_max = 1.d3     !== the maximal column dencity [g/cm^2] ==!
  n_m  = 1000
  n_mu = 40
  n_fi = 40
  allocate( m_atm_I(n_m,n_mu,n_fi,2),m_atm_T_rho(n_m,2),m_atm_S(n_m,n_mu,n_mu,2*n_fi,4) )  !== intensity; T & rho; S-matrix ==!
  dmu = 2.d0/n_mu; dfi = pi/n_fi
  i = 1
  do while( i.le.n_m )
    j1 = 1
    do while( j1 .le. n_mu )
      j2 = 1
      do while( j2 .le. n_mu )
        k = 1
        do while( k .le. 2*n_fi )
          m_atm_S(i,j1,j2,k,1:4) = (0.d0,1.d0)
          k = k+1
        end do
        j2 = j2+1
      end do
      j1 = j1+1
    end do
    write(*,*)i
    i = i+1
  end do
return
end subroutine set_atm_structure
