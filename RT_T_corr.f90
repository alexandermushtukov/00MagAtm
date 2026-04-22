!===================================================================================
!===================================================================================
subroutine get_accur_metrics(eps_surf,eps_Fmax,eps_Flog,F_target,Fz,mas_tau_TkeV,n_m)
implicit none
real*8,intent(out)::eps_surf,eps_Fmax,eps_Flog
real*8,intent(in)::F_target,Fz(n_m),mas_tau_TkeV(n_m,2)
integer,intent(in)::n_m
real*8::help
integer::i
  eps_surf = abs(F_target - Fz(1))/F_target
  eps_Fmax = 0.d0
  eps_Flog = 0.d0; help = 0.d0
  i=2
  do while(i.le.n_m)
    eps_Fmax = max( eps_Fmax, abs(F_target - Fz(1))/F_target )
    eps_Flog = eps_Flog + log( mas_tau_TkeV(i,1) - mas_tau_TkeV(i-1,1) ) * ( (F_target - Fz(1))/F_target )**2
    help = help + log( mas_tau_TkeV(i,1) - mas_tau_TkeV(i-1,1) )
    i = i+1
  end do
  eps_Flog = sqrt(eps_Flog/help)
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


