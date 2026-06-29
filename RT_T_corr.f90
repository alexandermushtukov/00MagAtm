!===================================================================================
! Getting metrics for temperature profile.
!===================================================================================
subroutine get_accur_metrics(delta_surf,eps_Fmax,delta_Flog,eps_Flog,eps_rough,F_target,Fz,mas_tau_TkeV,n_m)
implicit none
real*8,intent(out)::delta_surf,eps_Fmax,delta_Flog,eps_Flog,eps_rough
real*8,intent(in)::F_target,Fz(n_m),mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m
real*8::help,help2,w,rel_err,drel
integer::i,n_atm

  n_atm = n_m - 1   !== last point is lower boundary, not atmosphere ==!

  delta_surf = ( Fz(1) - F_target ) / max(F_target,1.d-30)

  eps_Fmax = 0.d0
  delta_Flog = 0.d0
  eps_Flog = 0.d0
  eps_rough = 0.d0
  help = 0.d0;  help2 = 0.d0

  i=2
  do while(i.le.n_atm)
    w = log( mas_tau_TkeV(i,1) / mas_tau_TkeV(i-1,1) )

    rel_err = ( Fz(i) - F_target ) / max(F_target,1.d-30)

    eps_Fmax = max( eps_Fmax, abs(rel_err) )
    delta_Flog = delta_Flog + w * rel_err
    eps_Flog = eps_Flog + w * rel_err**2
    help = help + w

    drel = ( Fz(i) - Fz(i-1) ) / max(F_target,1.d-30)
    eps_rough = eps_rough + w * drel**2
    help2 = help2 + w

    i = i+1
  end do

  if( help.gt.0.d0 )then
    delta_Flog = delta_Flog / help
    eps_Flog = sqrt( eps_Flog / help )
  else
    delta_Flog = 0.d0
    eps_Flog = 0.d0
  end if

  if( help2.gt.0.d0 )then
    eps_rough = sqrt( eps_rough / help2 )
  else
    eps_rough = 0.d0
  end if

return
end subroutine get_accur_metrics



!==================================================================================
!==================================================================================
subroutine temperature_correction(dT,F_target,Fz,mas_tau_TkeV,k_F,k_J,k_P,u,u_p,n_m,E_min,E_max,n_E,delta_surf)
use black_body
implicit none
real*8,intent(out)::dT(n_m)
real*8,intent(in)::F_target,Fz(n_m),mas_tau_TkeV(n_m,2),k_F(n_m),k_J(n_m),k_p(n_m),u(n_m),u_p(n_m),E_min,E_max,delta_surf
integer,intent(in)::n_m,n_E
real*8::help(n_m),f_int(n_m),damp_fact,damp_fact_2
integer::i,n_atm
real*8::kappa_T=0.4d0

  n_atm = n_m - 1   !== last point is lower boundary, not atmosphere ==!

  damp_fact   = 0.03d0  !0.3d0
  damp_fact_2 = 0.005d0 !0.05d0

  !if (abs(delta_surf) .lt. 0.3d0) then
  !  damp_fact = 0.10d0
  !  damp_fact_2 = 0.02d0
  !endif

  !if (abs(delta_surf) .lt. 0.1d0) then
  !  damp_fact = 0.03d0
  !  damp_fact_2 = 0.005d0
  !endif

  dT(1:n_m) = 0.d0
  help(1:n_m) = 0.d0
  f_int(1:n_m) = 0.d0

  !== integrand for the non-local term; exclude lower boundary ==!
  i = 1
  do while( i.le.n_atm )
    f_int(i) = k_F(i)/kappa_T * ( F_target - Fz(i) )
    i = i+1
  end do

  !== cumulative integral from surface to active atmosphere point ==!
  help(1) = 0.d0

  i = 2
  do while( i.le.n_atm )
    if( i.eq.2 )then
      help(i) = mas_tau_TkeV(1,1) * f_int(1)
    else
      help(i) = help(i-1) + ( mas_tau_TkeV(i-1,1) - mas_tau_TkeV(i-2,1) ) &
                            * ( f_int(i-1) + f_int(i-2) ) / 2.d0
    end if
    i = i+1
  end do

  !== temperature correction only for real atmosphere points ==!
  i = 1
  do while( i.le.n_atm )
    dT(i) = mas_tau_TkeV(i,2)/16.d0/( BB_Flux24_interval(mas_tau_TkeV(i,2),E_min,E_max,n_E) *100.d0 ) &
            * ( 3.d10/k_p(i) * ( k_J(i)*u(i) - k_p(i)*u_p(i) ) &
            + k_J(i)/k_p(i) * ( help(i) + 2.d0*( F_target - Fz(1) ) ) )

    dT(i) = dT(i)*damp_fact
    dT(i) = sign(1.d0,dT(i)) * min( abs(dT(i)), damp_fact_2 * mas_tau_TkeV(i,2) )
    i = i+1
  end do

  !== lower boundary is fixed; do not correct it ==!
  dT(n_m) = 0.d0

return
end subroutine temperature_correction
