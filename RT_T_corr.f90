!===================================================================================
! Getting metrics for temperature profile.
!===================================================================================
subroutine get_accur_metrics(delta_surf,eps_Fmax,delta_Flog,eps_Flog,eps_rough,F_target,Fz,mas_tau_TkeV,n_m)
implicit none
real*8,intent(out)::delta_surf,eps_Fmax,delta_Flog,eps_Flog,eps_rough
real*8,intent(in)::F_target,Fz(n_m),mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m
real*8::help,help2,w,rel_err,drel
integer::i
  delta_surf = ( Fz(1) - F_target ) / max(F_target,1.d-30)

  eps_Fmax = 0.d0
  delta_Flog = 0.d0
  eps_Flog = 0.d0
  eps_rough = 0.d0
  help = 0.d0;  help2 = 0.d0
  i=2
  do while(i.le.n_m)
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
    delta_Flog = 0.d0;  eps_Flog = 0.d0
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
subroutine temperature_correction(dT,F_target,Fz,mas_tau_TkeV,k_F,k_J,k_P,u,u_p,n_m)
use black_body
implicit none
real*8,intent(out)::dT(n_m)
real*8,intent(in)::F_target,Fz(n_m),mas_tau_TkeV(n_m,2),k_F(n_m),k_J(n_m),k_p(n_m),u(n_m),u_p(n_m)
integer,intent(in)::n_m
real*8::help
integer::i,j
real*8::kappa_T=0.4d0
  i = 1
  do while( i.le.n_m )
    help = 0.d0
    !== integration from 0 to given point ==!
    j = 1
    do while(j.le.(i-1))
      if(j.eq.1)then
        help = help + ( mas_tau_TkeV(j,1) - 0.d0 ) *k_F(j)/kappa_T* ( F_target - Fz(j) )
      else
        help = help + ( mas_tau_TkeV(j,1) - mas_tau_TkeV(j-1,1) ) &
                      * (k_F(j) + k_F(j-1))/2  /kappa_T* ( F_target - Fz(j) )
      end if
      j = j+1
    end do
    !== end integration ==!
    help = help + 2.d0 * ( F_target - Fz(1) )
    dT(i) = mas_tau_TkeV(i,2)/ 16 / ( BB_Flux24(mas_tau_TkeV(i,2))*100 ) &
             * ( 3.d10/k_p(i)* ( k_J(i)*u(i) - k_p(i)*u_p(i) ) + k_J(i)/k_p(i)*help )
    dT(i) = dT(i)*0.2d0
    dT(i) = sign(1.d0,dT(i)) * min( abs(dT(i)), mas_tau_TkeV(i,2)/10 )
    i = i+1
  end do
return
end subroutine temperature_correction


