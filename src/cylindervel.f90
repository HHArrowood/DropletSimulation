!For Non-Uniform Interface !Full2d & FF

!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH; adapted from work of C. Falcon, UNC-CH
!
! DESCRIPTION: 
!> Initializes spacial comp of stokes flow in a cylinder, used for 
!> interpolation later rather than recomputing every timestep 
!> @brief
!> I reworked Claudia's code, replacing her reflections with the drop 
!> reflection, and the third reflection with my crude first approx
!> to the third reflection
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation - H. Arrowood
! 02_Jan_2018 - Implemented reflections - H. Arrowood   
!
!-----------------------------------------------------------------------
subroutine cylindervelinit()
	use globalinfo
	implicit none
	integer index
	
!------Set up discretization pts----------------------------------------
	do index = 1, NZ
		WZ(index) = real(index-1.0, kind=8)*(LZ-epsilon)/NZ
		myk(index) = real(index-1.0, kind =8)*2.0*pi/(LZ-epsilon)
	end do
	
	do index = 1, NR
		WR(index) = real(index-1.0, kind=8)*LR/NR
	end do
!------Now call H and G functions to compute reflections----------------	
	
	call Hfunc(NZ, WZ+epsilon, myHZ)
	call Gfunc(NZ, WZ+epsilon, myGZ)
	call Hfunc(NR, WR, myHR)
	call Gfunc(NR, WR, myGR)

	firstpartZ = (WZ+epsilon)/2.0*(myHZ+myGZ);
	firstpartR = WR/2.0*(myHR+myGR);

	!print *, "firstpartZ", firstpartZ(10)
	!print *, "firstpartR", firstpartR(10)
		
	do index = 1, upperRangeZ
		zcoordinateZ(index) = (index-1.0)*2.0*pi/(LZ-epsilon)
	end do
	
	do index = 1, upperRangeR
		zcoordinateR(index) = (index-1.0)*2.0*pi/LR
	end do
	
	call cylindervelgrid()
	
!	print *, cylindervelZ
!	stop
end subroutine cylindervelinit

!-----------------------------------------------------------------------  
!> @author 
!> C. Falcon, UNC-CH; later edited by H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Used in cylindervelinit, which initializes the spacial component of 
!> Stokes flow in a cylinder 
!> @brief
!> 
! REVISION HISTORY:
! 16_Jan_2018 - Added documentation - H. Arrowood  
!
!-----------------------------------------------------------------------
subroutine cylindervelgrid()
use globalinfo
	implicit none
	
	real (kind=8) :: FZ(NZ), FR(NR), BESSI
	!real (kind=8), external :: sign
	integer myi, WRi, WZi
  	real    ( kind = 4 ) wsavez(4*NZ+15), wsaver(4*NR+15)
  	complex ( kind = 4 ) tempFR(NR), tempFZ(NZ)
  	real (kind = 8) :: besseli0R(NR), besseli1R(NR), besseli0Z(NZ), &
  	&besseli1Z(NZ)
	
    
	do myi = 1, XN+1
		!r-component of velocity
		do WRi = 1, NR
			besseli0R(WRi) = BESSI(0,WR(WRi)*x(myi))
			besseli1R(WRi) = BESSI(1,WR(WRi)*x(myi))
		end do
		
		tempFR = real((x(myi)*firstpartR*besseli0R-myGR*besseli1R)*0.5,&
			& kind=4)
		tempFR(1) = 0.0
!		print *, firstpartR
!		
!		print *, myGR
!    
!		print *, tempFR(16384)
!		stop
		call cffti ( NR, wsaver )
  		call cfftb ( NR, tempFR, wsaver )
  		FR = real(aimag(tempFR)/NR*LR/pi, kind=8)
  		cylindervelR(1:upperRangeR, myi) = FR(1:upperRangeR)
  		
		!z-component of velocity
		do WZi = 1, NZ
			besseli0Z(WZi) = BESSI(0,(WZ(WZi)+epsilon)*x(myi))
			besseli1Z(WZi) = BESSI(1,(WZ(WZi)+epsilon)*x(myi))
		end do
		
		tempFZ = real((x(myi)*firstpartZ*besseli1Z+myHZ*besseli0Z)*0.5,&
			& kind=4)
		call cffti ( NZ, wsavez )
		call cfftb ( NZ, tempFZ, wsavez )
		FZ = real(exp(imagi*epsilon*abs(myk))*real(tempFZ/NZ, kind=8)*&
			&(LZ-epsilon)/pi+3.0*R/pi*epsilon*((-2.0*R**2/(3.0*R0**2)+&
			&1.0)*x(myi)**2/R0**2+log(epsilon*0.5*R0)-1.0), kind=8)
		cylindervelZ(1:upperRangeZ, myi) = FZ(1:upperRangeZ)
  		
    end do

end subroutine cylindervelgrid




!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Uses stress coeffs and calculation based on the work of Haberman
!> to compute Stokes velocity  
!
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation - H. Arrowood
!
!> @param[in] G2,G4, G6, G8 pert stress coeffs    
!> @param[out] su, sv components of stokes velocity on interface pts      
!
!-----------------------------------------------------------------------
subroutine stokes () !(su, sv)
	use globalinfo
	implicit none
	
	real (kind=8) :: k1(XN+1), k2(XN+1),u3r(XN+1),u3z(XN+1),myr(XN+1), &
		mytheta(XN+1), tempx(4), tempy(4), myinterp1, myinterp2, &
		myinterp3, myinterp4 
	real(kind=8)	:: B2, B4, B6, B8, D2, D4, D6, D8!, u3constapprox
	real(kind=8) 	:: stokesvelocityr(XN+1), stokesvelocitytheta(XN+1)&
		&, stokesvelocity3r(XN+1), stokesvelocity3theta(XN+1), stokeste&
		&stx(XN+1), stokestesty(XN+1)
 	real(kind=8), external   :: LegendreP, GegenbauerC
	real (kind=8), external :: sign
	integer (kind=4) :: starti, sizecyl, i, j, k, ii
	
	k1 = 0
	k2 = 0

	sizecyl = size(cylindervelR, 2)
	
    do i=1, XN+1
        myr(i)  = dsqrt(sx(i)**2+sy(i)**2) !just doing this here rather than outside the do loop
	    mytheta(i) = dacos(sy(i)/myr(i))!atan(sx(i)/sy(i)) doesn't map to [0,pi] 
	   ! print*, "arctan", atan(sx(i)/sy(i))
	    !print*, "adcos", adcos(sy(i)/myr(i))

		!for non-uniform interface
		if (sx(i)<=R0cl) then
		starti = floor((sx(i)*(XNclose**2.0)/R0cl)**(1.0/2.0)+1.0)
		else
		starti = floor((sx(i)-R0cl)*XNfar/(R0-R0cl)+XNclose +1.0)
		endif

		starti = max(starti, 1);
		starti = min(starti, sizecyl-3);
	 
!---------horizontal velocity component (second reflection)-------------
!for reference,  interpbridge(N, interpx, interpy, xval, yval)----------
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upp&
		&erRangeR, starti), abs(sy(i)), myinterp1)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upp&
		&erRangeR, starti+1), abs(sy(i)), myinterp2)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upp&
		&erRangeR, starti+2), abs(sy(i)), myinterp3)
		call interpbridge(upperRangeR, zcoordinateR, cylindervelR(1:upp&
		&erRangeR, starti+3), abs(sy(i)), myinterp4)
		

        tempx = (/(startx(starti+j), j=0,3)/) 
		tempy = (/ myinterp1, myinterp2, myinterp3, myinterp4 /)
	 	tempy=tempy*sign(sy(i))
		call interpbridge(4, tempx, tempy, max(min(sx(i), R0), real(0.0&
		&, kind=8)), k1(i))
		
!----------vertical velocity component (second reflection)--------------
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upp&
		&erRangeZ, starti), abs(sy(i)), myinterp1)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upp&
		&erRangeZ, starti+1), abs(sy(i)), myinterp2)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upp&
		&erRangeZ, starti+2), abs(sy(i)), myinterp3)
		call interpbridge(upperRangeZ, zcoordinateZ, cylindervelZ(1:upp&
		&erRangeZ, starti+3), abs(sy(i)), myinterp4)
		 
		tempy = (/ myinterp1, myinterp2, myinterp3, myinterp4 /)
		call interpbridge(4, tempx, tempy, max(min(sx(i), R0), real(0.0&
		&, kind=8)), k2(i))
		

    end do
!print*, "mytheta", mytheta
!print*, "myr", myr
!print*, "k1*U", k1*U
!print*, "k2*U", k2*U
	!r and theta components of stokes velocities (first reflection)

	B2 = -R**3*(3*muin*U+R*mu*G2)/(6*(mu+muin))
	D2 = -(-6*R*mu*U-9*R*muin*U-R**2*mu*G2)/(6*(mu+muin))

	B4 = 0.0!R**6*mu*G4/(14*(mu+muin))
	D4 = 0.0!R**4*mu*G4/(14*(mu+muin))

	B6 = 0.0!R**8*mu*G6/(22*(mu+muin))
	D6 = 0.0!R**6*mu*G6/(22*(mu+muin))

	B8 = 0.0!R**10*mu*G8/(30*(mu+muin))
	D8 = 0.0!R**8*mu*G8/(30*(mu+muin))
!print*, "U in cylindervel before addition", U
v2atorigin = 1./(2.*Pi)*(1./(R0**3*(1.+ (muin/mu)))*R*U*(R0**2*(-8.566&
	&226429811636-12.849339644717457*(muin/mu))+6.515546450454206*R**2&
	&*(muin/mu))) !this has been checked! -H.
!	print*, "U in cylindervel after addition", U
	do k=1, XN+1
		stokesvelocityr(k)=-(LegendreP(1,mytheta(k))*(-(U)+B2*1/(myr(k)&
			&**3)+D2*(1/myr(k)))+LegendreP(3,mytheta(k))*(B4*1/(myr(k)&
			&**5)+D4*(1/myr(k)**3))+LegendreP(5,mytheta(k))*(B6*1/&
			&(myr(k)**7)+D6*(1/myr(k)**5))+LegendreP(7,mytheta(k))*&
			&(B8*1/(myr(k)**9)+D8*(1/myr(k)**7)))

		if (k==1) then
			stokesvelocitytheta(k)=0.0
		else
			stokesvelocitytheta(k)=1/dsin(mytheta(k))*(GegenbauerC(2,&
				&mytheta(k))*(-2*(U)-B2*1/(myr(k)**3)+D2*1/myr(k))+Gege&
				&nbauerC(4,mytheta(k))*(-3*B4*1/(myr(k)**5)-D4*1/(myr(k&
				&)**3))+GegenbauerC(6,mytheta(k))*(-5*B6*1/(myr(k)**7)-&
				&3*D6*1/(myr(k)**5))+GegenbauerC(8,mytheta(k))*(-7*B8*1&
				&/(myr(k)**9)-5*D8*1/(myr(k)**7)))
		endif
	end do

!------Approx Third Reflection/Leading Order Error on Drop Surface------
!------This causes the code to violate boundary conditions on both the 
!surface of drop and cylinder walls, but is the approximation I used to 
!compute the drag-------------------------------------------------------


	
do k=1, XN+1
		stokesvelocity3r(k)=-LegendreP(1,mytheta(k))*(-(v2atorigin)+(-R&
		&**3*(3*muin*v2atorigin+R*mu*G2)/(6*(mu+muin)))*1./(myr(k)&
		&**3)+((6*R*mu*v2atorigin+9*R*muin*v2atorigin+R**2*mu*G2)/(6*(m&
		&u+muin)))*(1/myr(k)))!-(LegendreP(1,mytheta(k))*(v2atorigin+(-R**&
		!&3*(-3*muin*v2atorigin)/(6*(mu+muin)))*1/(myr(k)**3)+(-6*R*mu*&
		!&v2atorigin+9*R*muin*v2atorigin)/(6*(mu+muin))*(1/myr(k))))
!(-R**3*(-3*muin*U)/(6*(mu+muin))) !b2
!(-6*R*mu*U+9*R*muin*U)/(6*(mu+muin))!d2

		if (k==1) then
			stokesvelocity3theta(k)=0
		else
			stokesvelocity3theta(k)=1/dsin(mytheta(k))*(GegenbauerC(2,&
			&mytheta(k))*(-2*(v2atorigin)-(-(R**3*( 3.*muin*v2atorigi&
			&n))/(6.*(mu + muin)))*1/(myr(k)**3)+((6*R*mu*v2atorigi&
			&n+9*R*muin*v2atorigin)/(6*(mu+muin)))*1/myr(k)))
				!1/dsin(mytheta(k))*(GegenbauerC(2,&
				!&mytheta(k))*(2*v2atorigin-(-R**3*(-3*muin*v2atorigin)/&
				!&(6*(mu+muin)))*1/(myr(k)**3)+(-6*R*mu*v2atorigin+9*R*&
				!&muin*v2atorigin)/(6*(mu+muin))*1/myr(k)))
		endif
	end do	
	
!------Stokes Velocity Full Assembly-------------------------- 
	!print*, "stokesvelocitytheta", stokesvelocitytheta
	do ii=1,XN+1
		if(ii==1)then !force first point horizontal vel to be 0. numerical error in sin(pi) was making it off
			if (myr(ii)>1.2*R)then 
				su(ii) = 0.0!-(dsin(mytheta(ii))*stokesvelocityr(ii)+dcos(mytheta(ii))*stokesvelocitytheta(ii))!+ k1(ii)*U!&
					!&-(dsin(mytheta(ii))*stokesvelocity3r(ii)+dcos(mytheta(ii))*stokesvelocity3theta(ii))
				sv(ii) = -(dcos(mytheta(ii))*stokesvelocityr(ii)-dsin(m&
				&ytheta(ii))*stokesvelocitytheta(ii))+ k2(ii)*U!&
					!&-(dcos(mytheta(ii))*stokesvelocity3r(ii)-dsin(mytheta(ii))*stokesvelocity3theta(ii)-v2atorigin)
			else
				su(ii) = 0.0!-((dsin(mytheta(ii))*stokesvelocityr(ii)+dcos(mytheta(ii))*stokesvelocitytheta(ii)))!&
				!&+(dsin(mytheta(ii))*stokesvelocity3r(ii)+dcos(mytheta(ii))*stokesvelocity3theta(ii)))
				sv(ii) = -((dcos(mytheta(ii))*stokesvelocityr(ii)-dsin(&
				&mytheta(ii))*stokesvelocitytheta(ii))&
				&+(dcos(mytheta(ii))*stokesvelocity3r(ii)-dsin(mytheta(&
				&ii))*stokesvelocity3theta(ii)))
			endif
		else
			if (myr(ii)>1.2*R)then
				su(ii) = -(dsin(mytheta(ii))*stokesvelocityr(ii)+dcos(m&
				&ytheta(ii))*stokesvelocitytheta(ii))+ k1(ii)*U!&
					!&-(dsin(mytheta(ii))*stokesvelocity3r(ii)+dcos(mytheta(ii))*stokesvelocity3theta(ii))
				sv(ii) = -(dcos(mytheta(ii))*stokesvelocityr(ii)-dsin(&
				&mytheta(ii))*stokesvelocitytheta(ii))+ k2(ii)*U!&
					!&-(dcos(mytheta(ii))*stokesvelocity3r(ii)-dsin(mytheta(ii))*stokesvelocity3theta(ii)-v2atorigin)
			else
				su(ii) = -((dsin(mytheta(ii))*stokesvelocityr(ii)+dcos(&
				&mytheta(ii))*stokesvelocitytheta(ii))&
				&+(dsin(mytheta(ii))*stokesvelocity3r(ii)+dcos(mytheta(&
				&ii))*stokesvelocity3theta(ii)))
				sv(ii) = -((dcos(mytheta(ii))*stokesvelocityr(ii)-dsin(&
				&mytheta(ii))*stokesvelocitytheta(ii))&
				&+(dcos(mytheta(ii))*stokesvelocity3r(ii)-dsin(mytheta(&
				&ii))*stokesvelocity3theta(ii)))
			endif
		endif
	end do	
print*, "U", U
print*, "v2 at origin", v2atorigin
print*, "sin(pi)", sin(pi)
print*, "svr(1)", stokesvelocityr(1)
print*, "su at bottom of drop", su(1)

end subroutine stokes


!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH; adapted from work of C. Falcon, UNC-CH
!
! DESCRIPTION: 
!> Evaluates the H function of the second reflection for a drop in a 
!> cylinder 
!> @brief
!> I reworked Claudia's code, adding in the reflections for the drop 
!> rather than the sphere case. Note that this is the spacial part only,
!> all V(t)-dependence has been factored out. 
! REVISION HISTORY:
!> 12_Dec_2017 - adapted H expression - H. Arrowood
!> 16_Jan_2018 - added documentation - H. Arrowood
!> 23_Jan_2018 - verified H - H. Arrowood   
!
! @param Num - the number of points at which to evaluate the function
! @param lambda - the variable of integration in Hfunc
! @return ReturnH - the H-vals at each point
!-----------------------------------------------------------------------
subroutine Hfunc(Num, lambda, ReturnH)
use globalinfo
	implicit none
	
	integer index, Num
	real (kind=8) :: lambda(Num), besselk0(Num), besselk1(Num),&
	  	besseli1(Num), besseli2(Num), besseli0(Num), ReturnH(Num),&
		BESSK, BESSI
	
	do index = 1, Num
		besselk0(index) = BESSK(0,R0*lambda(index))
		besselk1(index) = BESSK(1,R0*lambda(index))
		besseli0(index) = BESSI(0,R0*lambda(index))
		besseli1(index) = BESSI(1,R0*lambda(index))
		besseli2(index) = BESSI(2,R0*lambda(index))
	end do
		
	ReturnH = (R*(R0*lambda*besseli0*(-(4. + 6.*(muin/mu) + &
		&R**2*(muin/mu)*lambda**2)*real(besselk0) + R0*(2. + 3.*&
		&(muin/mu))*lambda*real(besselk1)) + besseli1*((8. + 12.*&
		&(muin/mu) + (2.*R**2*(muin/mu) + R0**2*(2. + 3.*(muin/mu)))*&
		&lambda**2)*real(besselk0) - R0*lambda*(4. + 6.*(muin/mu) + &
		&R**2*(muin/mu)*lambda**2)*real(besselk1))))/((1. + (muin/mu))*& 
		&(R0*lambda*besseli0**2 - 2.* besseli0*besseli1 - R0*lambda*&
		&besseli1**2))


end subroutine

!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH; adapted from work of C. Falcon, UNC-CH
!
! DESCRIPTION: 
!> Evaluates the G function of the second reflection for a drop in a 
!> cylinder 
!> @brief
!> I reworked Claudia's code, adding in the reflections for the drop 
!> rather than the sphere case. Note that this is the spacial part only,
!> all V(t)-dependence has been factored out. 
! REVISION HISTORY:
!> 12_Dec_2017 - adapted G expression - H. Arrowood
!> 16_Jan_2018 - added documentation - H. Arrowood
!> 23_Jan_2018 - final verification - H. Arrowood   
!
! @param Num - the number of points at which to evaluate the function
! @param lambda - the variable of integration in Gfunc
! @return ReturnH - the G-vals at each point
!-----------------------------------------------------------------------
subroutine Gfunc(Num, lambda, ReturnG)
	use globalinfo
	implicit none
	
	integer index, Num
	real (kind=8) :: lambda(Num), besselk0(Num), besselk1(Num), &
		& besselk2(Num), besseli1(Num), besseli2(Num), besseli0(Num), &
		& ReturnG(Num), BESSK, BESSI

	do index = 1, Num
		besselk0(index) = BESSK(0,R0*lambda(index))		
		besselk1(index) = BESSK(1,R0*lambda(index))
		besselk2(index) = BESSK(2,R0*lambda(index))
		besseli0(index) = BESSI(0,R0*lambda(index))
		besseli1(index) = BESSI(1,R0*lambda(index))
		besseli2(index) = BESSI(2,R0*lambda(index))
	end do
		
	 ReturnG = (lambda**2*(-R**3*R0*(muin/mu)*lambda*besseli1*&
		&real(besselk0) + 2.*R*R0**2*besseli1*real(besselk0) + 3.*R*&
		&R0**2*(muin/mu)*besseli1*real(besselk0) + 2.*R*R0**2*besseli0*&
		&real(besselk1) - 2.*R**3*(muin/mu)*besseli0*real(besselk1) + &
		&3.*R*R0**2*(muin/mu)*besseli0*real(besselk1) - R**3*R0*&
		&(muin/mu)*lambda*besseli1*real(besselk1)))/((1. + (muin/mu))*&
		&(-R0*lambda*besseli0**2 + 2.*besseli0*besseli1 + R0*lambda*&
		&besseli1**2))
		
end

function sign(val)
	implicit none
	
	real(kind=8) val, sign
	
	if (val< 0) then
		sign = -1.0
	else
		sign = 1.0
	endif
end function


