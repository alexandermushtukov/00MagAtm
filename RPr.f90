!========================================================================================================
! Function calculates complex scattering amplitudes for electrons.
!   ksi_i, ksi_f - ellipticities of photons before and after scattering event.
!   Note: we assume coherent scattering.
!   See appendix B in 2204.12271
!========================================================================================================
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


!========================================================================================================
! Function calculates complex scattering amplitudes for protons.
!   ksi_i, ksi_f - ellipticities of photons before and after scattering event.
!   Note: we assume coherent scattering.
!   See appendix B in 2204.12271
!========================================================================================================
complex*16 function Amp_dSigmadOmega_ell_magnitars_p(s1,s2,E,Ecyc,rho,theta1,theta2,fi1,fi2,Z,A,ksi_i,ksi_f)
implicit none
integer,intent(in)::s1,s2
real*8,intent(in)::E,Ecyc,rho,theta1,theta2,fi1,fi2,Z,A,ksi_i,ksi_f
complex*16::Amp_dSigmadOmega_ell_magnitars  !== function ==!
real*8::pi=3.141592653589793d0
   Amp_dSigmadOmega_ell_magnitars_p = Amp_dSigmadOmega_ell_magnitars(s1,s2,E,Ecyc/1838,rho,pi-theta1,pi-theta2,fi1,fi2,Z,A,ksi_i,ksi_f)
   Amp_dSigmadOmega_ell_magnitars_p = Amp_dSigmadOmega_ell_magnitars_p/1838
return
end function Amp_dSigmadOmega_ell_magnitars_p


!=====================================================================================================
!=====================================================================================================
subroutine test_cross_sections()
implicit none
real*8::pi=3.141592653589793d0
real*8::B12,theta_i,theta_f,fi_i,fi_f,E_cyc,E,ksi_i,ksi_f
complex*16::Amp_dSigmadOmega_ell_magnitars_p,Amp_dSigmadOmega_ell_magnitars  !== functions ==!
real*8:: NormWavesEll_cvp_ !== function ==!
complex*16::amp
real*8::sc_e,sc_p,singma,dtheta,dfi,sigma_e,sigma_p ,rho,Z,A
integer::n_theta,n_fi,i,j
  n_theta = 50; n_fi = 30
  dtheta = pi/n_theta
  dfi = 2*pi/n_fi

  Z = 1.d0; A = 1.d0
  rho = 1.d0
  B12 = 1.d2
  E_cyc = 11.6d0*B12
  theta_i = 0.d0
  fi_i = 1.d0

  E = 0.1d0
  do while(E.le.1000.d0)
    ksi_i = NormWavesEll_cvp_(1,E,E_cyc,theta_i,rho,Z,A)  !== (+)-mode ==!
    sigma_e = 0.d0
    sigma_p = 0.d0
    i=1
    do while( i.le.n_theta )
      theta_f = dtheta/2 + (i-1)*dtheta
      j = 1
      do while(j.le.n_fi)
        fi_f = dfi/2 + (j-1)*dfi

        ksi_f = NormWavesEll_cvp_(2,E,E_cyc,theta_i,rho,Z,A)
        amp = Amp_dSigmadOmega_ell_magnitars(1,1,E,E_cyc,rho,theta_i,theta_f,fi_i,fi_f,Z,A,ksi_i,ksi_f)
        sigma_e = sigma_e + (abs(amp))**2 * 3/32/pi *dfi*dtheta*sin(theta_f)
        amp = Amp_dSigmadOmega_ell_magnitars_p(1,1,E,E_cyc,rho,theta_i,theta_f,fi_i,fi_f,Z,A,ksi_i,ksi_f)
        sigma_p = sigma_p + (abs(amp))**2 * 3/32/pi *dfi*dtheta*sin(theta_f)

        ksi_f = NormWavesEll_cvp_(1,E,E_cyc,theta_i,rho,Z,A)
        amp = Amp_dSigmadOmega_ell_magnitars(1,1,E,E_cyc,rho,theta_i,theta_f,fi_i,fi_f,Z,A,ksi_i,ksi_f)
        sigma_e = sigma_e + (abs(amp))**2 * 3/32/pi *dfi*dtheta*sin(theta_f)
        amp = Amp_dSigmadOmega_ell_magnitars_p(1,1,E,E_cyc,rho,theta_i,theta_f,fi_i,fi_f,Z,A,ksi_i,ksi_f)
        sigma_p = sigma_p + (abs(amp))**2 * 3/32/pi *dfi*dtheta*sin(theta_f)

        j = j+1
      end do
      i = i+1
    end do
    write(*,*)E,sigma_e,sigma_p
    E = E*1.1d0
  end do
return
end subroutine test_cross_sections


!====================================================================================================================
! Elliptisity of normal waves: (-i)(Ey/Ex), where Oy \in (k,B), Ox \perp (k,B).
! Cold plasma + vacuum effects are taken into account.
!   alpha - defines the type of photon: X-mode (alpha=1) or O-mode (alpha=0 or 2)
! This definition of elliptisity of related to the one from Lai(2003): NormWavesEll = -K .
! [E]=[E_cyc]=[keV].
!====================================================================================================================
real*8 function NormWavesEll_cvp_(alpha,E,Ecyc,theta,rho,Z,A)
implicit none
integer,intent(in)::alpha
real*8,intent(in)::E,Ecyc,theta,rho,Z,A
real*8::B12
real*8::K_1,K_2,Kz_1,Kz_2
complex*16::beta,epsilon_pv,c_K_1,c_K_2,c_Kz_1,c_K_z2,K_plus,K_min
dimension epsilon_pv(3,3)
real*8::M_Stokes_rot
dimension M_Stokes_rot(4,4)
real*8::aa,n_1,n_2
  
  B12=Ecyc/11.6d0
  call Dielectric_Tensor_Plasma_Vac(epsilon_pv,M_Stokes_rot,beta,aa,c_K_1,c_K_2,K_plus,K_min,c_Kz_1,c_K_z2,&
                                    n_1,n_2,E,theta,B12,rho,Z,A)
  if(alpha.eq.1)then
    NormWavesEll_cvp_ = -real(K_plus)
  else
    NormWavesEll_cvp_ = -real(K_min)
  end if
return
end function NormWavesEll_cvp_



!================================================================================================================
! The function calculates the cross section of free-free absorption in units of Thomson scattering cross-section.
! Note that the absorption cross secton is dependent in local density (proportional).
! Temperature participates in sigma_0 only.
! (4.4.37) Meszaros 1992.
! We provide this function with the ellipticity of initial photons.
!================================================================================================================
real*8 function abs_mag_ff_Meszaros_new(K_j,E,E_cyc,T_keV,theta,Z,A,rho)
use abs_mag_ff_Meszaros
implicit none
real*8,intent(in)::K_j,E,E_cyc,T_keV,theta,Z,A,rho
real*8::koef_e_0,koef_e_min,koef_e_plus,g_perp,g_par,g_R,res,K_jz,K_1,K_2,Kz_1,Kz_2
real*8::NormWavesEll  !==function==!
real*8::sigma_0

  sigma_0 = 2.86d0/E**3 * Z**2/A * rho/sqrt(T_keV)   !==non-magnetic absorption cross section in units if Thomson scattering CS==!

  K_jz=0.d0
  !koef_e_0    = (K_j*sin(theta)-K_jz*cos(theta))**2 / (1.d0+K_j**2+K_jz**2)
  !koef_e_plus = (1.d0+(K_j*cos(theta)+K_jz*sin(theta)))**2/2/ (1.d0+K_j**2+K_jz**2)  !==2003MNRAS338_233H_(3.1)==!
  !koef_e_min  = (1.d0-(K_j*cos(theta)+K_jz*sin(theta)))**2/2/ (1.d0+K_j**2+K_jz**2)  !==2003MNRAS338_233H_(3.1)==!

  koef_e_0    = ( K_j*sin(theta) )**2 / ( 1.d0 + K_j**2 )
  koef_e_plus = ( 1.d0 + K_j*cos(theta) )**2 / 2/ ( 1.d0 + K_j**2 )  !==2003MNRAS338_233H_(3.1)==!
  koef_e_min  = ( 1.d0 - K_j*cos(theta) )**2 / 2/ ( 1.d0 + K_j**2 )  !==2003MNRAS338_233H_(3.1)==!
  call g( E/511, T_keV/511, E_cyc/511, theta, g_perp, g_par, g_R)
  res = E**2/(E+E_cyc)**2*koef_e_plus*g_perp + koef_e_min*g_R + koef_e_0*g_par
  abs_mag_ff_Meszaros_new = res*sigma_0
return
end function abs_mag_ff_Meszaros_new
!============================================================


