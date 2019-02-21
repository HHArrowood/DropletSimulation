
!----------------------------------------------------------------------- 
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Uses Newton's method to find a new interface point on a trajectory
!> between the trajectories of the preceding and following points. 
!> @brief
!> Interpolation routine designed with geometry of our problem in mind. 
!
! REVISION HISTORY:
! 11_May_2018 - Implemented - H. Arrowood
!
!> @param[in] x1in,y1in,x2in,y2in                
!> @return xout, yout
!----------------------------------------------------------------------- 
subroutine NewtonInterpolation(x1in,y1in,x2in,y2in,xout,yout)
	use globalinfo
	implicit none

	real(kind=8), intent(in) :: x1in,y1in,x2in,y2in !Cartesian coords of input points
	real(kind=8)   :: r1, r2, theta1, theta2,thet!Polar coords of input points
	real(kind=8) :: B2, D2,K, l, ri
	real(kind=8)  :: rnew, rnew1, rnew2, r1result, r2result, rguess, th&
	&etanew, elar(2), lout,rout, left, right
	integer(kind=4) :: i
	real(kind=8), external :: NewtonIter1,NewtonIter2, GegenbauerC!, Bisection
	real(kind=8) :: xout, yout !! final points

!----To begin, convert input x and y vals into r and theta vals----------

	print*, "x1 inside routine", x1in
	
	r1 = dsqrt(x1in**2+y1in**2)
	theta1 = dacos(y1in/r1)
	r2 = dsqrt(x2in**2+y2in**2)
	theta2 = dacos(y2in/r2)
	print*, "theta1", theta1
	print*, "theta2", theta2
	thet=(theta1+theta2)/2.
	print*, "dsqrt(x2in**2+y2in**2)diff", dsqrt(x2in**2+y2in**2)-sqrt(x&
	&2in**2+y2in**2)
	!U=.1
	!G2 = -.01
	U=Unewforthirdref !to get through third ref, U=U+v2atorigin
	!G2=-.01
!Find K value, the average of the streamfunction evaluated at x1 and at x2
	B2 = -(R**3*(G2*mu*R + 3.*muin*U))/(6.*(mu + muin))
	D2 = -(-(G2*mu*R**2) - 6.*mu*R*U - 9.*muin*R*U)/(6.*(mu + muin))
	k = (GegenbauerC(2,theta1)*(-U*r1**2+B2*1/r1+D2*r1)+GegenbauerC(2, &
	&theta2)*(-U*r2**2+B2*1./r2+D2*r2))/2.
print*, "k", k
print*, "r1", r1
print*, "r2", r2
print*,"theta", thet
print*, "B2", B2
print*, "D2", D2
!Use Newton's iteration, applied four times to get machine precision under
!reasonable circumstances, to find rnew

Print*, "r1 result", NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(Ne&
&wtonIter1(r1,k,B2, D2, thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)&
&,k,B2,D2,thet)
print*, "r2 result", NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(Ne&
&wtonIter1(r2,k,B2, D2, thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)&
&,k,B2,D2,thet)

r1result=NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(r1&
&,k,B2, D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)
r2result= NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(r&
&2,k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)

	if(r1result>=R)then
		rnew=r1result
	elseif(r2result>=R)then
		rnew=r2result
	else

!		elar(1)=R
!		elar(2)=max(r1,r2, 1.1*R)
	
!		do i=1,50
!		    print*, "elar", elar

			
!		    print*, "elar(1)", elar(1)
!		    l = elar(1)
!		    ri = elar(2)
!		    print*, "elar(2)", elar(2)
		   
!			call Bisection(l, ri, k, B2, D2, thet, lout, rout)
!			elar(1)=lout
!			elar(2)=rout
!		print*, "LandR", lout, rout
!		end do

		rnew=max(r1,r2)!R*1.0000000000001!
	endif
	
!rnew1 = NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(r1,k,B2, D2, thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)
!print*, "rnew1", rnew1
!now with same guess, but deflated polynomial
!rnew2 = NewtonIter2(NewtonIter2(NewtonIter2(NewtonIter2(NewtonIter2(r1,k,B2, D2, thet,rnew1),k,B2,D2,thet,rnew1),k,B2,D2,thet,rnew1),k,B2,D2,thet,rnew1),k,B2,D2,thet,rnew1)
!print*, "rnew2", rnew2

!now do what Wilkinson(1963) recommends, and iterate with the appropriate root from deflated polynomial
!rguess= max(rnew1,rnew2)

!rnew = NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(NewtonIter1(max(r1,r2),k,B2, D2, thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet),k,B2,D2,thet)

!Let's give the bisection method a try. 30 iterations should get us pretty close to machine precision

!	elar(1)=R
!	elar(2)=max(r1,r2, 1.1*R)
	
!	do i=1,50
	 !   print*, "elar", elar

	    
	 !   print*, "elar(1)", elar(1)
!	    l = elar(1)
!	    ri = elar(2)
	 !   print*, "elar(2)", elar(2)
	   
!		call Bisection(l, ri, k, B2, D2, thet, lout, rout)
!		elar(1)=lout
!		elar(2)=rout
		!print*, "LandR", lout, rout
!	end do

!	rnew=elar(2)
	print*, "rnew", rnew
!Convert rnew,thetanew into xout and yout
	xout = rnew*dsin(thet)
	yout = rnew*dcos(thet)
end subroutine NewtonInterpolation 


!----------------------------------------------------------------------- 
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Newton's method for the cubic function with roots equal to the desired rnew 
!> for the above subroutine to find a new interface point on a trajectory
!> between the trajectories of the preceding and following points. 
!> @brief
!> Interpolation routine designed with geometry of our problem in mind. 
!
! REVISION HISTORY:
! 11_May_2017 - Implemented - H. Arrowood
!
!> @param[in]  r               
!> @return rout
!----------------------------------------------------------------------- 


function NewtonIter1(rr, k, B2, D2,theta)

	use globalinfo
	implicit none
	
	real(kind=8), intent(in) :: rr, k, B2, D2, theta
	real(kind=8) :: f, fp
	real(kind=8) :: NewtonIter1
	real(kind=8), external :: GegenBauerC
	
	!U=Unewforthirdref

	f=-U*rr**3+D2*rr**2-K/GegenbauerC(2,theta)*rr+B2
	fp=-3.*U*rr**2+2.*D2*rr-K/GegenbauerC(2,theta)
	
	NewtonIter1 = rr-f/fp
	return
end function NewtonIter1


function NewtonIter2(rr, k, B2, D2,theta, rt1)
!for the deflated polyomial
	use globalinfo
	implicit none
	
	real(kind=8), intent(in) :: rr, k, B2, D2, theta, rt1
	real(kind=8) :: f, fp
	real(kind=8) :: NewtonIter2
	real(kind=8), external :: GegenBauerC
	
	!U=Unewforthirdref

	f=-U*rr**3+D2*rr**2-K/GegenbauerC(2,theta)*rr+B2
	fp=-3.*U*rr**2+2.*D2*rr-K/GegenbauerC(2,theta)
	
	NewtonIter2 = rr-f/((fp-f)*(rr-rt1))
	return
end function NewtonIter2
	

subroutine Bisection(l, ri, k, B2, D2, theta, lout, rout)
	use globalinfo
	implicit none
	
	real(kind=8), intent(in) :: l, ri, k, B2, D2, theta
	real(kind=8) ::fl, fr, fm, m
	real(kind=8), external ::GegenbauerC
	real(kind=8) :: lout, rout
	print*, "l", l
	print*, "ri", ri
	
	fl=-U*l**3+D2*l**2-K/GegenbauerC(2,theta)*l+B2
	fr=-U*ri**3+D2*ri**2-K/GegenbauerC(2,theta)*ri+B2
	m = (l+ri)/2
	fm = -U*m**3+D2*m**2-K/GegenbauerC(2,theta)*m+B2
	print*, "fl", fl
	print*, "fr", fr
	print*, "fm", fm
	!Bisection(1)=0.0
	!Bisection(2)=0.0
	if(fl>0.0 .and. fm<0.0)then
			lout=l
			rout=m
		!	print*, "bisect1", lout, rout
			!return
	elseif(fl<0.0 .and. fm>0.0)then
			lout=l
			rout=m
		!	print*, "bisect2", lout, rout
			!return
	elseif(fr<0.0 .and. fm>0.0)then
			lout=m
			rout=ri
			!print*, "bisect3", lout, rout
			!return
	elseif(fr>0.0 .and. fm<0.0)then
			lout=m
			rout=ri
			!print*, "bisect4", lout, rout
	else
		print*, "interval not big enough (already converged to drop sur&
		&face)"
		lout=l
		rout=ri
		print*, "bisect5", lout, rout
		!return
	endif
end subroutine Bisection
	
	
	

