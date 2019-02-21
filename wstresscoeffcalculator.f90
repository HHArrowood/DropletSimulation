!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Uses velocities computed on surface of drop to find the 
!>coefficients of the Gegenbauer expansion of the stresses.
!>Run this subroutine right after computing w flow. 
!> @brief
!> Numerically integrates w.r.t. theta (trapezoidal method for now)
!> Computes \f$I_n,G_n\f$ from series of the form
!> \f$\Sum_{n=2}^{\infty}I_n(r)C_n^{-1/2}(\cos(\theta)) \f$
!> for velocities and stresses
!
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation - H. Arrowood
! 28_Nov_2017 - Updated to working version - H. Arrowood
!
!> @param[in] wu, wv  Cart components of pertvel near surface of drop
!> @param[in] thetavecforstress  theta-coords of points
!> @param[in] nsurfacepoints   number of integration pts near surface    
!> @param[out] I2, G2, G4, G6, G8, G2HO    
!
!----------------------------------------------------------------------- 

subroutine stresscoeffs 
	use globalinfo
	real(kind=8),external                  :: GegenbauerC !>Gegen fcn
	real(kind=8), dimension(nsurfpoints+1) :: e21int, h21int, h22int, &
		&h41int, h61int, h81int, d2h2int, h21intHO
	real(kind=8)		               :: e21, h21, h22, h41, h61, h81,&
		&d2dh2
		!> gegen coeffs of vel. e21 is coeff of radial velocity for n=2;
		!> hn1, hn2, are coeffs of tangential velocity,
		!>hn1 at r=a+epsilon,and  hn2 at r=a+2*epsilon
		
!first compute coefficients of wtheta and wr at n evenly spaced points. 
!I need two values, for r=a+epsilon and r=a+2*epsilon, in order to 
!compute the first and second derivatives of H2, and the rest just at 
!r=a+epsilon

!evaluation of wu and wv happens inside pertvelT

!-------Find integrands for stress coeffs-------------------------------
		do istress=1, nsurfpoints+1
	!r derivative of radial velocity for n=2 at r=a+epsilon
		e21int(istress) = -3.0/2.0*(sin(thetavecforstress(istress))*1.0&
			&/eps*wU1(istress)+cos(thetavecforstress(istress))*1.0/eps*&
			&wV1(istress))*sin(thetavecforstress(istress))*cos(thetavec&
			&forstress(istress))
		
	!r derivative of tangential velocity at r=a+epsilon
		h21int(istress) =  3.0*(cos(thetavecforstress(istress))*1.0/eps&
			&*wU1(istress)-sin(thetavecforstress(istress))*1.0/eps*wV1&
			&(istress))*GegenbauerC(2,thetavecforstress(istress))
			!print*, "h21int", h21int
		h21intHO(istress) =  3.*(3.*(cos(thetavecforstress(istress))*1.&
		&/eps*wU1(istress)-sin(thetavecforstress(istress))*1.0/eps*wV1&
			&(istress))-1.5*(cos(thetavecforstress(istress))*1.0/eps&
			&*wU2(istress)-sin(thetavecforstress(istress))*1.0/eps*wV2&
			&(istress))+1./3.*(cos(thetavecforstress(istress))*1.0/eps&
			&*wU3(istress)-sin(thetavecforstress(istress))*1.0/eps*wV3&
			&(istress)))*GegenbauerC(2,thetavecforstress(istress))
			!print*, "h21intHO", h21intHO
			!print*, "WU1", WV1
			!print*, "WU2", WV2
			!print*, "WU3", WV3
		h41int(istress) = 42.0*(cos(thetavecforstress(istress))*1.0/eps&
			&*wU1(istress)-sin(thetavecforstress(istress))*1.0/eps*wV1(&
			&istress))*GegenbauerC(4,thetavecforstress(istress))
		h61int(istress) = 165.0*(cos(thetavecforstress(istress))*1.0/ep&
			&s*wU1(istress)-sin(thetavecforstress(istress))*1.0/eps*wV1&
			&(istress))*GegenbauerC(6,thetavecforstress(istress))
		h81int(istress) = 420.0*(cos(thetavecforstress(istress))*1.0/ep&
			&s*wU1(istress)-sin(thetavecforstress(istress))*1.0/eps*wV1&
			&(istress))*GegenbauerC(8,thetavecforstress(istress))
		!tangential velocity at r=a+2*epsilon, divided by eps already in
		!preparation for the second derivative w.r.t. r
		h22int(istress) = 3.0*(cos(thetavecforstress(istress))*1.0/eps*&
			&wU2(istress)-sin(thetavecforstress(istress))*1.0/eps*wV2(i&
			&stress))*GegenbauerC(2,thetavecforstress(istress))
			
		!Second r derivative of tan velocity
		d2h2int = 1/eps*h22int-2/eps*h21int
	end do
	
!--------Integrate to obtain r derivs of velocity coeffs----------------
!@fixme maybe derive a Gauss quadrature for this instead of using trap. 
	e21 =  (2.0*sum(e21int)-e21int(1)-e21int(nsurfpoints+1))* &
		&(pi/(2.0*nsurfpoints))

	h21 = (2.0*sum(h21int)-h21int(1)-h21int(nsurfpoints+1))* &
		&(pi/(2.0*nsurfpoints))		
	G2HO = (2.0*sum(h21intHO)-h21intHO(1)-h21intHO(nsurfpoints+1))* &
		&(pi/(2.0*nsurfpoints))
	h41 = (2.0*sum(h41int)-h41int(1)-h41int(nsurfpoints+1))* &
		&(pi/(2.0*nsurfpoints))
	h61 = (pi/(2.0*nsurfpoints))*(2.0*sum(h61int)-h61int(1)- &
		&h61int(nsurfpoints+1))
	h81 = (pi/(2.0*nsurfpoints))*(2.0*sum(h81int)-h81int(1)- &
		&h81int(nsurfpoints+1))

	
	h22 = (pi/(2.0*nsurfpoints))*(2.0*sum(h22int)-h22int(1)- &
		&h22int(nsurfpoints+1))
 

	d2dh2 = (pi/(2.0*nsurfpoints))*(2.0*sum(d2h2int)-d2h2int(1)-&
		&d2h2int(nsurfpoints+1))
		
	

!---now assemble into stress coeffs------------------------------------- 

!tangential stress coefficients are just r-derivatives
	G2 = h21	
	G4 = h41
	G6 = h61
	G8 = h81

!normal stress coeff for n=2
	I2 = G2+1.0/2.0*R*d2dh2-2.0*e21

end subroutine stresscoeffs

