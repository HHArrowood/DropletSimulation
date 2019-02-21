!
!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH; adapted from work of C. Falcon, UNC-CH
!
! DESCRIPTION: 
!> Computes perturbation velocity using trapezoidal integration
!> @brief
!> Computes perturbation velocities, using a combo of 2D and 3D 
!> integrals. Uses far field approx if drop is far away from interface,
!> full expression if close
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation - H. Arrowood
! 04_Feb_2017 - Added r-derivative integral for stress coeff
! 
!
!-----------------------------------------------------------------------
subroutine wStressTN () !(wu, wv)
	use globalinfo
	implicit none

	real (kind = 8)    :: wstbackflow(nsurfpoints+1), wstsidesphere(nsu&
		&rfpoints+1), wstbelowsphere(nsurfpoints+1), wstabovesphere(nsu&
		&rfpoints+1),tempbackflowst, tempsidespherest, tempbelowsphere&
		&st, tempabovespherest, G2int(nsurfpoints+1), dwthetadrintegran&
		&d(nsurfpoints+1)
	real (kind=8), external :: w2StIntegrandBackflow,GegenbauerC, &
		&w2StIntegrandBelowSphere, w2StIntegrandPartialSphere,&
		&w2StIntegrandZetaSphere, w2StIntegrandZetaVert, w2StIntegrand&
		&R, w2StIntegrandZ, dwthetadr
	integer (kind=4) wsti, wsti2, wsti3
 

	integer  ierr,i, ist, istressnew

print*, "XN inside wst", XN

!--------For stress calculator------------------------------------------
!first initialize n evenly-spaced integration points on drop surface

	thetavecforstress=(/(Real(i),i=0,nsurfpoints)/)*pi/Real(nsurfpoints)

!just to see if this matches the interface used by w
!print*, "sx inside Stress Calculator", sx
!----Initialize the velocities at each region of the perturbed interface

	wstbackflow	    = real(0.0, kind=8)
	wstsidesphere 	= real(0.0, kind=8)
	wstbelowsphere 	= real(0.0, kind=8)
	wstabovesphere 	= real(0.0, kind=8)
	AreaReflux      = real(0.0, kind=8)
	AreaEntrain     = real(0.0, kind=8)



	if (minval(sy) < maxval(sy)) then !check that interface is perturbed
		!print*,"newstress is happening!"
!----------------Initialize Interpolated Interface----------------------
		
		if (allocated(cwinterpx)) deallocate(cwinterpx, stat=ierr)
		if (allocated(cwinterpy)) deallocate(cwinterpy, stat=ierr)
		
		allocate(cwinterpx(1000000), stat=ierr)
		if (ierr /= 0) print*, "cwinterpx : Allocation failed" 
		!> @fixme  make this abort 

		allocate(cwinterpy(1000000), stat=ierr)
		if (ierr /= 0) print*, "cwinterpy : Allocation failed"

!----Do loop over dwthetadr stress evaluation pts-----------------------
		do wsti = 1, nsurfpoints+1!+2+XN+
			
			px=thetavecforstress(wsti)
!--Initialize these integration vars to 0------------------------------
			cwinternalcount=0.0;
			cwinterpx = 0;
			cwinterpy = 0;

			tempbackflowst	= real(0.0, kind=8)
			tempsidespherest 	= real(0.0, kind=8)
			tempbelowspherest	= real(0.0, kind=8)
			tempabovespherest	= real(0.0, kind=8)

!-----Now actually evaluate final integrals over each region------------
!print*, "line 125"
!print*, "flagb", flagb
			if(flagb <XN+1 ) then
				!print *, "w0"
				call trapz1(w2StIntegrandBackflow, xflagb, sx(XN+1), &        !trapz1(f,a,b,h1,r) f the function to integrate, (a, b) the interval of 
					&numtrapz, tempbackflowst)                                !integration, h1 the subinterval length, r the result
				!call simp2(w2StIntegrandBackflow, xflagb, sx(XN+1), integthres, tempbackflowst)
				!print*, "after integration"
			endif
!print*, "line 131"
			if (flagu /= 0) then
				!print *, "w1"
				call trapz1(w2StIntegrandPartialSphere,real(0.0, kind=8&
				&), xflagl, numtrapz, tempbelowspherest)
				!call simp2(w2StIntegrandPartialSphere, real(0.0, kind=8), &
				!	&xflagl, integthres, tempbelowspherest)
				print *, "w2"
				call trapz1(w2StIntegrandZetaSphere, max(-R, sy(1)), R,&
					&numtrapz, tempsidespherest)
				!call simp2(w2StIntegrandZetaSphere, max(-R, sy(1)), R, &
				!	&integthres, tempsidespherest)
				!print *, "w3"
				call trapz1(w2StIntegrandZetaVert, R, yend, numtrapz, &
					&tempabovespherest)
				!call simp2(w2StIntegrandZetaVert, R, yend, integthres, &
				!	&tempabovespherest)
				!print *, "w3.5"
			elseif (flagl /= 0) then
				!print *, "w4"
				call trapz1(w2StIntegrandPartialSphere, real(0.0, kind=&
				&8), xflagl, numtrapz, tempbelowspherest)
				!call simp2(w2StIntegrandPartialSphere, real(0.0, kind=8), &
				!	&xflagl, integthres, tempbelowspherest)
				print *, "w5"
				call trapz1(w2StIntegrandZetaSphere, max(-R, sy(1)), ye&
				&nd, numtrapz, tempsidespherest)
				!call simp2(w2StIntegrandZetaSphere, max(-R, sy(1)), yend, &
				!	&integthres, tempsidespherest)
				!print *, "w5.5"
			else
				!print *, "w6" 

				!print *, "below sphere limits", real(0.0, kind=8), xflagb
				call trapz1(w2StIntegrandBelowSphere, max(sx(1), real(0&
				&.0, kind=8)), xflagb, numtrapz, tempbelowspherest)
				!call simp2(w2StIntegrandBelowSphere, max(sx(1), real(0.0, &
				!	&kind=8)), xflagb, integthres, tempbelowspherest)

				!print *, "belowsphere", tempbelowsphere

			endif

			wstbackflow(wsti) 		= tempbackflowst
			wstbelowsphere(wsti) 	= tempbelowspherest
			wstsidesphere(wsti) 	= tempsidespherest
			wstabovesphere(wsti) 	= tempabovespherest
			
		end do	
	endif !the if that checks if the interface is perturbed
	
!----Make vector of dwthetadr values------------------------------------

		do wsti2=1, nsurfpoints+1
			if(wsti2==nsurfpoints+1) then
				print*,"actual dwthetadrintegrand final element", (-wst&
				&backflow(wsti2)+wstbelowsphere(wsti2)+wstsidesphere &!this is some really hacky code. think of something better later.
				&(wsti2)+ wstabovesphere(wsti2))*drhogover8mu
					dwthetadrintegrand(wsti2) = 0
			elseif(isnan(wstsidesphere(wsti2)))then !Detect NaNs due to integrator issues maybe also add a threshold for "singularities"
				wstsidesphere(wsti2)=0
				dwthetadrintegrand(wsti2) = (-wstbackflow(wsti2)+wstbel&
				&owsphere(wsti2)+wstsidesphere(wsti2)+ wstabovesphere(w&
				&sti2))*drhogover8mu
				
			else		
				dwthetadrintegrand(wsti2) = (-wstbackflow(wsti2)+wstbel&
				&owsphere(wsti2)+wstsidesphere(wsti2)+ wstabovesphere(w&
				&sti2))*drhogover8mu
			endif
			
		end do
	


!print*, "wstbackflow", wstbackflow	
!print*,	"wstbelowsphere", wstbelowsphere 	
!print*,"wstsidesphere", wstsidesphere	
!print*,	"wstabovesphere", wstabovesphere 	
	!if(abs(wbelowsphere(1)) > 1000.0) then
	!stop
	!endif
	do istressnew = 1, nsurfpoints+1
		G2int(istressnew) =  3.0*(dwthetadrintegrand(istressnew))*Gegen&
		&bauerC(2,thetavecforstress(istressnew))
	end do
	!print*, "G2Newint", G2int
	G2new = (2.0*sum(G2int)-G2int(1)-G2int(nsurfpoints+1))*(pi/(2.0*nsu&
	&rfpoints))
end subroutine wStressTN


!========== For Velocity Integrands =========!

function w2StIntegrandBackflow(rho)
	use globalinfo
	implicit none

	real (kind=8) rho, zcoord, w2StIntegrandBackflow
	real (kind=8), external :: w2StIntegrandZeta

	!print *, "w1"
	!print*, "flagb", flagb
	!print*, "sx", sx
	call interpbridge( XN+2-max(flagb-2, 1), sx(max(flagb-2,1):XN+1), &
	&sy(max(flagb-2, 1):XN+1), rho, zcoord)
	!print *, "afterw1"
	myrho   = rho !this is where the value of myrho is assigned...gets used further down
	call trapz1(w2StIntegrandZeta, yend, zcoord, numtrapz, w2StIntegran&
	&dBackflow)
	!call simp2(w2StIntegrandZeta, yend, zcoord, integthres, w2StIntegrandBackflow)

	cwinternalcount = cwinternalcount +1;
	cwinterpx(cwinternalcount) = rho;
	cwinterpy(cwinternalcount) = zcoord;

end

function w2StIntegrandZetaSphere(zeta)
	use globalinfo
	implicit none

	real (kind=8) zeta, xupper, xlower, w2StIntegrandZetaSphere
	real (kind=8), external :: w2StIntegrandRho

	myzeta  = zeta
	!print*, "myzeta", myzeta
	!print *, "w2"
	call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(&
	&1:min(flagb+2, XN+1)), zeta, xupper)
	!print *, "w2"
	print*, "zeta", zeta
	print*, "R", R
	xlower  = sqrt(R**2 - zeta**2)
	if(isnan(xlower))then !to fix the case due to numerical error where xlower became a NaN because zeta was off by 1E-16
		xlower=0.0
	endif
	
	!Now to check that xupper isn't accidentally inside the drop due to interpolation error
	if(xupper**2+zeta**2<R**2)then
		print*, "upper bound was inside!!"
		xupper=xlower !So basically, if the interpolation accidentally puts it inside, just let the layer thickness go to 0
	endif
	call trapz1(w2StIntegrandRho, xlower, xupper, numtrapz, w2StIntegra&
	&ndZetaSphere)
	!call simp2(w2StIntegrandRho, xlower, xupper, integthres, w2StIntegrandZetaSphere)

	cwinternalcount = cwinternalcount +1;
	cwinterpx(cwinternalcount) = xupper;
	cwinterpy(cwinternalcount) = zeta;
end

function w2StIntegrandPartialSphere(rho)
use globalinfo
implicit none

	real (kind=8) rho, zcoord, w2StIntegrandPartialSphere
	real (kind=8), external :: w2StIntegrandZeta

	myrho = rho
	!print *, "w3"
	call interpbridge(min(flagl+3, XN+1), sx(1:min(flagl+3, XN+1)), sy(&
	&1:min(flagl+3, XN+1)), rho, zcoord)
	!print *, "w3"
	call trapz1(w2StIntegrandZeta, zcoord, -R, numtrapz, w2StIntegrandP&
	&artialSphere)
	!call simp2(w2StIntegrandZeta, zcoord, -R, integthres, w2StIntegrandPartialSphere)
	cwinternalcount = cwinternalcount +1;
	cwinterpx(cwinternalcount) = rho;
	cwinterpy(cwinternalcount) = zcoord;
end

function w2StIntegrandZetaVert(zeta)
	use globalinfo
	implicit none

	real (kind=8) zeta, xupper, w2StIntegrandZetaVert
	real (kind=8), external :: w2StIntegrandRho
	integer (kind=4) temp(1), tempmini

	myzeta  = zeta

	call interpbridge(min(flagb+2, XN+1), sy(1:min(flagb+2, XN+1)), sx(&
	&1:min(flagb+2, XN+1)), zeta, xupper)

	call trapz1(w2StIntegrandRho, real(0.0, kind=8), xupper, numtrapz,w&
	&2StIntegrandZetaVert)
	!call simp2(w2StIntegrandRho, real(0.0, kind=8), xupper, integthres, w2StIntegrandZetaVert)
	
	cwinternalcount = cwinternalcount +1;
	cwinterpx(cwinternalcount) = xupper;
	cwinterpy(cwinternalcount) = zeta;
end


function w2StIntegrandBelowSphere(rho)
	use globalinfo
	implicit none

	real (kind=8) rho, zcoord, w2StIntegrandBelowSphere
	real (kind=8), external :: w2StIntegrandZeta


	!print *, "w5p2"
	call interpbridge(max(flagl+2, XN+1), sx(1:max(flagl+2, XN+1)), sy(&
	&1:max(flagl+2, XN+1)), rho, zcoord)
	!print *, "w5p2"
	myrho   = rho
	!print *, zcoord, yend, myrho
	call trapz1(w2StIntegrandZeta, zcoord, yend, numtrapz, w2StIntegran&
	&dBelowSphere)
	!call simp2(w2StIntegrandZeta, zcoord, yend, integthres, w2StIntegrandBelowSphere)

end                                           




!=====================================================
!Integrands
!=====================================================

function w2StIntegrandZeta(zeta)
	use globalinfo
	implicit none

	real (kind=8) zeta, w2StIntegrandZeta
	real (kind=8), external :: w2StIntegrandR!, w2StIntegrandZ!, w2IntegrandRFF, w2IntegrandZFF 


	myzeta = zeta

		w2StIntegrandZeta =w2StIntegrandR()
	
end



function w2StIntegrandRho(rho)
	use globalinfo
	implicit none

	real (kind=8) rho, w2StIntegrandRho
	real (kind=8), external :: w2StIntegrandR!, w2IntegrandZ, w2IntegrandRFF, w2IntegrandZFF


	myrho = rho
	!if ( ((py-myzeta)**2+(px-myrho)**2)>singthres) then



	w2StIntegrandRho = w2StIntegrandR()

end




function w2StIntegrandR()!(ellipticK,ellipticE,ellipticK1, ellipticE1)
! this is the integrand after theta integration
	use globalinfo
	implicit none

	real (kind=8) w2StIntegrandR,ILogR!,ellipticE, ellipticK,ellipticE1, ellipticK1
	real (kind=8), external :: dWthetadr! I1R,I2R,I3R,I4R,I5R,I6R,I7R,I8R,

	w2StIntegrandR=0.0;
	ILogR = 0.0;


	if (px > 0.0) then

!for now, I'm going to assume the pseudosingularities are the same as in the velocities

		if (abs (myzeta*px + py * myrho ) < logsing) then
		!call simp2(LogTermThetaR, real(0.0,kind=8), real(2.0*pi,kind=8),integthres, ILogR)

		call trapz1(dWthetadr,  real(0.0,kind=8),real(2.0*pi,kind=8), l&
		&ogtrapzbig, ILogR)
		!call simp2(dWthetadr,  real(0.0,kind=8),real(2.0*pi,kind=8), integthres, ILogR)
		else
		call trapz1(dWthetadr,  real(0.0,kind=8),real(2.0*pi,kind=8), l&
		&ogtrapz, ILogR)
		!call simp2(dWthetadr,  real(0.0,kind=8),real(2.0*pi,kind=8), integthres, ILogR)
		endif




  !print*, "ILogR", ILogR
		

		w2StIntegrandR = ILogR*myrho!/pi
	endif
	!print *, "w2IntegrandR", w2IntegrandR

end





function dWthetadr(theta)
	use globalinfo
	implicit none
	
	real (kind=8), intent(in) ::  theta
	real (kind=8) :: dWthetadr,thetax, y1, y2, y3
	real (kind=8) :: ry, dWU1, dWU2, phicoeffderiv, dphidx1term1, dphid&
	&x1term2,dphidx1term3, dphidx1term4, dphidx1term5, dWV1, dWV2, dph&
	&idx3term1,dphidx3term2, dphidx3term3, dphidx3term4, dphidx3term5, &
	&dWU, dWV
	!First, we need the wu and wv terms, differentiated w.r.t. r. I broke 
	!these each up into several terms to hopefully streamline the debugging 
	!a little.
	!x-values on drop surface will be input as R (drop radius) from globalinfo and a vector
	!of thetax values, so I can use these to greatly reduce the number of 
	!of operations. 

	thetax = px
	y1 = myrho *cos(theta)
	y2= myrho* sin(theta)
	y3 = myzeta
	ry = sqrt(y1**2+y2**2+y3**2)
	!print*, y1**2+y2**2+y3**2
	!print*, "ry", ry
	!print*, "sqrt(75)", sqrt(75.)

	dWU1 =  -((R**4*(-(R*y3) + (y1**2 + y2**2 + y3**2)*Cos(thetax))*Sin&
	&(thetax))/((y1**2 + y2**2 + y3**2)**2.5*((R**2*(R**2 + y1**2 + y2*&
	&*2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2&
	&**2 + y3**2))**1.5)) - (R**4*Cos(thetax)*(-(R*y1) + (y1**2 + y2**2&
	& + y3**2)*Sin(thetax)))/((y1**2 + y2**2 + y3**2)**2.5*((R**2*&
	& (R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(t&
	&hetax)))/(y1**2 + y2**2 + y3**2))**1.5) + (3*(-(R*y3) + (y1**2 + y&
	&2**2 + y3**2)*Cos(thetax))*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**&
	&2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3*&
	&*2))*(y1**2 + y2**2 + y3**2 - R*y3*Cos(thetax) - R*y1*Sin(thetax))&
	&*(-(R*y1) + (y1**2 + y2**2 + y3**2)*Sin(thetax)))/((y1**2 + y2**2 &
	&+ y3**2)**1.5*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) -&
	& 2*R*y1*Sin(thetax))**3) -  (3*(y3 - R*Cos(thetax))*(-y1 + R*Sin(t&
	&hetax))*(-R + y3*Cos(thetax) + y1*Sin(thetax)))/(y2**2 + (y3 - R*C&
	&os(thetax))**2 + (y1 - R*Sin(thetax))**2)**2.5 + ((-y3 + R*Cos(the&
	&tax))*Sin(thetax))/(y2**2 + (y3 - R*Cos(thetax))**2 + (y1 - R*Sin(&
	&thetax))**2)**1.5 + (Cos(thetax)*(-y1 + R*Sin(thetax)))/(y2**2 + (&
	&y3 - R*Cos(thetax))**2 + (y1 - R*Sin(thetax))**2)**1.5

		 
	dWU2 = ((-R**2 + y1**2 + y2**2 + y3**2)*(R*Cos(thetax)*(-(y1*(2*y1*&
	&*4 + 2*y2**4 + 3*y2**2*y3**2 + y3**4 - R**2*(y1**2 + y2**2 - 2*y3*&
	&*2) + y1**2*(4*y2**2 + 3*y3**2))) + R*(y1**4 + y1**2*y2**2 + y2**2&
	&*y3**2 + y3**4)*Sin(thetax)) + y3*(-(R*(y1**4 + R**2*(2*y1**2 - y2&
	&**2 - y3**2) + 3*y1**2*(y2**2 + y3**2) + 2*(y2**2 + y3**2)**2)*Sin&
	&(thetax)) + y1*((y1**2 + y2**2 + y3**2)*(2*R**2 + y1**2 + y2**2 + &
	&y3**2) + R**2*y1*y3*Sin(2*thetax)))))/((y1**2 + y2**2 + y3**2)**2.&
	&5*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(&
	&thetax))**2*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(&
	&thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)))
		
	!Note that the coefficient of the dphi term is zero on the drop surface
	!So there is no need to take a derivative of the phi term, just the coeff
	!when applying the product rule to this term.     
	phicoeffderiv = -((R*(-R**2+ry**2))/ry**3)
		
		
	dphidx1term1 =    (4*R**3*y1*y3 + 7*R*y1**3*y3 + 7*R*y1*y2**2*y3 + &
	&7*R*y1*y3**3 - 6*R*y1*y3*(y1**2 + y2**2 + y3**2)*Cos(2*thetax) - 1&
	&4*R**2*y1**2*y3*Sin(thetax) - 3*y1**4*y3*Sin(thetax) - 6*R**2*y2**&
	&2*y3*Sin(thetax) - 6*y1**2*y2**2*y3*Sin(thetax) - 3*y2**4*y3*Sin(t&
	&hetax) - 6*R**2*y3**3*Sin(thetax) - 6*y1**2*y3**3*Sin(thetax) - 6*&
	&y2**2*y3**3*Sin(thetax) - 3*y3**5*Sin(thetax) + R*Cos(thetax)*(R*y&
	&1*(3*y1**2 + 3*y2**2 - 5*y3**2) - 3*(y1**4 + y2**4 - 3*y3**4)*Sin(&
	&thetax)) - 3*R*y1**2*y2**2*Sin(2*thetax) + 3*R*y1**2*y3**2*Sin(2*t&
	&hetax) + 3*R*y2**2*y3**2*Sin(2*thetax))/(R**2*(R**2 + y1**2 + y2**&
	&2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax))**2*Sqrt((R**2&
	&*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(t&
	&hetax)))/(y1**2 + y2**2 + y3**2)))
	   
	   
	dphidx1term2 = (3*R**2*(-(R*y1) + (y1**2 + y2**2 + y3**2)*Sin(theta&
	&x))*(R*(y1**2 + y2**2 - y3**2)*Cos(thetax) + y3*(y1**2 + y2**2 + y&
	&3**2 - 2*R*y1*Sin(thetax))))/((y1**2 + y2**2 + y3**2)*((R**2*(R**2&
	& + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)&
	&))/(y1**2 + y2**2 + y3**2))**1.5*(-R**4 + R**3*y3*Cos(thetax) + R*&
	&*3*y1*Sin(thetax) + y1**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*Sqrt(&
	&(R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*&
	&Sin(thetax)))/(y1**2 + y2**2 + y3**2)) + y2**2*Sqrt(R**4/(y1**2 + &
	&y2**2 + y3**2))*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*&
	&Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)) + y3**&
	&2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*Sqrt((R**2*(R**2 + y1**2 + y2&
	&**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y&
	&2**2 + y3**2))))

		
	dphidx1term3 =   (-3*(y1**2 + y2**2 + y3**2)*((R**4*Sin(thetax))/Sq&
	&rt(R**4/(y1**2 + y2**2 + y3**2)) + R*y1*(-Sqrt(R**4/(y1**2 + y2**2&
	& + y3**2)) + Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos&
	&(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2))))*(R**3*y&
	&3*(2*R**2 + y1**2 + y2**2 + y3**2 - 2*R*y1*Sin(thetax) - (2*R**2*S&
	&qrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R&
	&*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)))/Sqrt(R**4/(y1**2 + y2*&
	&*2 + y3**2))) + Cos(thetax)*(-(R**4*(y1**2 + y2**2 + 3*y3**2)) + &
	&Sqrt(R**4/(y1**2 + y2**2 + y3**2))*(y1**2 + y2**2 + y3**2)**2*&
	&Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*&
	&R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)))))/(R**3*(R**2 + y1**2&
	& + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax))*&
	&(-R**4 + R**3*y3*Cos(thetax) + R**3*y1*Sin(thetax) + y1**2*Sqrt(R*&
	&*4/(y1**2 + y2**2 + y3**2))*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3*&
	&*2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3&
	&**2)) + y2**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*Sqrt((R**2*(R**2 &
	&+ y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax))&
	&)/(y1**2 + y2**2 + y3**2)) + y3**2*Sqrt(R**4/(y1**2 + y2**2 + y3**&
	&2))*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) &
	&- 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)))**2)
		
			   
	dphidx1term4 =    (-3*y3*(y1**2 + y2**2 + y3**2)*Sin(thetax))/(R**2&
	&*(R**4/Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + R*Sqrt(R**2)*y3*Cos(th&
	&etax) + R*Sqrt(R**2)*y1*Sin(thetax)))

	 
	dphidx1term5 = (3*(y1**2 + y2**2 + y3**2)*(R*Sqrt(R**2)*y3 + (R**4*&
	&Cos(thetax))/Sqrt(R**4/(y1**2 + y2**2 + y3**2)))*(R*Sqrt(R**2)*y1 &
	&+ (R**4*Sin(thetax))/Sqrt(R**4/(y1**2 + y2**2 + y3**2))))/(R**7*Sq&
	&rt(R**2)*((R*Sqrt(R**2))/Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + y3*C&
	&os(thetax) + y1*Sin(thetax))**2)


	WV1 = (-3*(-y3 + R*Cos(thetax))**2*(2*Cos(thetax)*(-y3 + R*Cos(thet&
	&ax)) + 2*Sin(thetax)*(-y1 + R*Sin(thetax))))/(2.*(y2**2 + (-y3 + R&
	&*Cos(thetax))**2 + (-y1 + R*Sin(thetax))**2)**2.5) +(2*Cos(thetax)&
	&*(-y3 + R*Cos(thetax)))/(y2**2 + (-y3 + R*Cos(thetax))**2 + (-y1 +&
	& R*Sin(thetax))**2)**1.5 - (2*Cos(thetax)*(-y3 + R*Cos(thetax)) + &
	&2*Sin(thetax)*(-y1 + R*Sin(thetax)))/(2.*(y2**2 + (-y3 + R*Cos(the&
	&tax))**2 + (-y1 + R*Sin(thetax))**2)**1.5) + (3*R**3*(-((R**2*y3)/&
	&(y1**2 + y2**2 + y3**2)) + R*Cos(thetax))**2*(2*Cos(thetax)*(-((R*&
	&*2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(thetax)) + 2*Sin(thetax)*(&
	&-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(thetax))))/(2.*(y1**2&
	& + y2**2 + y3**2)**1.5*((R**4*y2**2)/(y1**2 + y2**2 + y3**2)**2 + &
	&(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(thetax))**2 + (-((R*&
	&*2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(thetax))**2)**2.5) - (2*R*&
	&*3*Cos(thetax)*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta&
	&x)))/((y1**2 + y2**2 + y3**2)**1.5*((R**4*y2**2)/(y1**2 + y2**2 + &
	&y3**2)**2 + (-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(thetax))&
	&**2 + (-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(thetax))**2)**&
	&1.5) + (R*(2*Cos(thetax)*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R&
	&*Cos(thetax)) + 2*Sin(thetax)*(-((R**2*y1)/(y1**2 + y2**2 + y3**2)&
	&) + R*Sin(thetax))))/(2.*Sqrt(y1**2 + y2**2 + y3**2)*((R**4*y2**2)&
	&/(y1**2 + y2**2 + y3**2)**2 + (-((R**2*y3)/(y1**2 + y2**2 + y3**2)&
	&) + R*Cos(thetax))**2 + (-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*&
	&Sin(thetax))**2)**1.5)
		
		 
	WV2 = (y3*(-R**2+ry**2)*(2.*R**2*y1**2*y3+y1**4*y3+2.*R**2*y2**2*y&
	&3+2.*y1**2*y2**2*y3+y2**4*y3+2.*R**2*y3**3+2.*y1**2*y3**3+2.*y2**&
	&2*y3**3+y3**5+R*(-4.*y1**4-4.*y2**4-7.*y2**2*y3**2-3.*y3**4+R**2*&
	&(2.*y1**2+2.*y2**2-y3**2)-y1**2*(8.*y2**2+7.*y3**2))*cos(thetax)+&
	&R**2*y3*ry**2*cos(2.*thetax)-3.*R**3*y1*y3*sin(thetax)+R*y1**3*y3*&
	&sin(thetax)+R*y1*y2**2*y3*sin(thetax)+R*y1*y3**3*sin(thetax)+R**2*&
	&y1**3*sin(2.*thetax)+R**2*y1*y2**2*sin(2.*thetax)+R**2*y1*y3**2*&
	&sin(2.*thetax)))/(ry**5*(R**2+ry**2-2.*R*y3*cos(thetax)-2.*R*y1*&
	&sin(thetax))**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(thetax)-&
	&2.*R*y1*sin(thetax)))/ry**2))

		  
	dphidx3term1 = (-y3*(3.*ry**4+R**2*(5.*y1**2+5.*y2**2+13.*y3**2))*&
	&cos(thetax)+1./2.*R*(2.*R**2*y1**2-y1**4+2.*R**2*y2**2-2.*y1**2*&
	&y2**2-y2**4+10.*R**2*y3**2+12.*y1**2*y3**2+12.*y2**2*y3**2+13.*&
	&y3**4-3.*(y1**4+y2**4-2.*y2**2*y3**2-3.*y3**4+2.*y1**2*(y2**2-&
	&y3**2))*cos(2.*thetax)-4.*R*y1*(y1**2+y2**2+5.*y3**2)*sin(thetax)+&
	&12.*y1**3*y3*sin(2.*thetax)+12.*y1*y2**2*y3*sin(2.*thetax)+12.*y1*&
	&y3**3*sin(2.*thetax)))/(R**2*(R**2+ry**2-2.*R*y3*cos(thetax)-2.*R*&
	&y1*sin(thetax))**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(thetax)-&
	&2.*R*y1*sin(thetax)))/ry**2))
		
		
	dphidx3term2 =  (3*(y1**2 + y2**2 + y3**2)**2*((R**6*(R*y3 - (y1**2&
	& + y2**2 + y3**2)*Cos(thetax))**2)/  (y1**2 + y2**2 + y3**2)**3 - &
	&(R**5*y3*(R*y3 - (y1**2 + y2**2 + y3**2)*Cos(thetax))*(R**2 + y1**&
	&2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1*&
	&*2 + y2**2 + y3**2)**3 - ((R**4/(y1**2 + y2**2 + y3**2))**1.5*(R**&
	&2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax&
	&))*(Sqrt(R**4/(y1**2 + y2**2 + y3**2)) - Sqrt((R**2*(R**2 + y1**2 &
	&+ y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2&
	& + y2**2 + y3**2))))/R**2))/(R**3*((R**2*(R**2 + y1**2 + y2**2 + y&
	&3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + &
	&y3**2))**1.5*(-R**4 + R**3*y3*Cos(thetax) + R**3*y1*Sin(thetax) + &
	& y1**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*Sqrt((R**2*(R**2 + y1**2&
	& + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**&
	&2 + y2**2 + y3**2)) + y2**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*Sqr&
	&t((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(thetax) - 2*R*y&
	&1*Sin(thetax)))/(y1**2 + y2**2 + y3**2)) + y3**2*Sqrt(R**4/(y1**2 &
	&+ y2**2 + y3**2))*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y&
	&3*Cos(thetax) - 2*R*y1*Sin(thetax)))/(y1**2 + y2**2 + y3**2))))
		
		
	dphidx3term3 = -((3.*ry**2*((R**4*cos(thetax))/(R**2/ry)+R*y3*&
	&(-(R**2/ry)+sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(thetax)-&
	&2.*R*y1*sin(thetax)))/(ry**2))))*(R**3*y3*(2.*R**2+ry**2-&
	&2.*R*y1*sin(thetax)-(2.*R**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*&
	&cos(thetax)-2.*R*y1*sin(thetax)))/ry**2))/(R**2/ry))+cos(thetax)*&
	&(-R**4*(y1**2+y2**2+3.*y3**2)+(R**2/ry)*ry**4*sqrt((R**2*(R**2+&
	&ry**2-2.*R*y3*cos(thetax)-2.*R*y1*sin(thetax)))/(ry**2)))))/(R**3*&
	&(R**2+ry**2-2.*R*y3*cos(thetax)-2.*R*y1*sin(thetax))*(-R**4+R**3*y&
	&3*cos(thetax)+R**3*y1*sin(thetax)+y1**2*(R**2/ry)*sqrt((R**2*&
	&(R**2+ry**2-2.*R*y3*cos(thetax)-2.*R*y1*sin(thetax)))/ry**2)+&
	&y2**2*(R**2/ry)*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(thetax)-2.*R*&
	&y1*sin(thetax)))/ry**2)+y3**2*(R**2/ry)*sqrt((R**2*(R**2+ry**2-2.*&
	&R*y3*cos(thetax)-2.*R*y1*sin(thetax)))/ry**2))**2))
		
		
	dphidx3term4 = -((3.*ry**2*(ry+y3*cos(thetax)))/(R**2*(R**2*ry+&
		&R**2*y3*cos(thetax)+R**2*y1*sin(thetax))))
		
		
	dphidx3term5 = (3.*ry**2*(R**2*y3+(R**4*cos(thetax))/(R**2/ry))**2)&
	&/(R**7*R*(ry+y3*cos(thetax)+y1*sin(thetax))**2)
		
if(isnan(dWU)) then
print*, "dWU isNan"
	if(isnan(dWU1))then
	print*, "dWU1 isNaN"
	endif
	if(isnan(dWU2))then
	print*, "dWU2 isNaN"
	endif
	if(isnan(phicoeffderiv*(dphidx1term1+dphidx1term2+&
		&dphidx1term3+dphidx1term4+dphidx1term5)))then
	print*, "phi isNaN"
	endif
endif
	
		
if(isnan(dWV)) then
print*, "dWV isNan"
endif	
		dWU = dWU1+dWU2+phicoeffderiv*(dphidx1term1+dphidx1term2+&
		&dphidx1term3+dphidx1term4+dphidx1term5)
		dWV = dWV1+dWV2+phicoeffderiv*(dphidx3term1+dphidx3term2+&
		&dphidx3term3+dphidx3term4+dphidx3term5)
		
		!print*, "WU", WU
		!print*, "WV", WV!verified
		
		dWthetadr = cos(thetax)*dWU - sin(thetax)*dWV
		return 
end 



