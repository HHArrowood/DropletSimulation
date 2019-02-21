!-----------------------------------------------------------------------
! Main Routine for Oil Droplet Rising Through Sharply Stratified Fluid
!-----------------------------------------------------------------------
!
!MODULE: Global Information
!
!> @author
!> H. Arrowood, UNC-CH. Adapted from code by C. Falcon, UNC-CH
!
!DESCRIPTION:
!>Defines global parameters and initializes global variables. 
!
!REVISION HISTORY:
!06 Nov 2017 - Final Version documented
!-----------------------------------------------------------------------

module globalinfo
	implicit none
	
!========Actual Experimental Parameters================================
	real (kind=8), parameter :: rho0pctKI = 1.24017
	real (kind=8), parameter :: rho2pctKI = 1.25335
	real (kind=8), parameter :: rho4pctKI = 1.26521
    real (kind=8), parameter :: rhotop 	  = rho0pctKI!<top density g/cc
	real (kind=8), parameter :: rhobottom = rho4pctKI!1.38!<bot density g/cc
	real (kind=8), parameter :: rhodrop	  = 0.87372 !<drop density g/cc
	real (kind=8), parameter :: mu        = .5*(2.782011+2.712486)!2.74822!(2.782011+2.712486+2.74822)/3.!(2.712486+2.74822)/2.!(2.782011+2.712486+2.74822)/3.!<dyn visc am mPa.s/100
	real (kind=8), parameter :: muin	  = .021672  !<dyn visc drop "
	real (kind=8), parameter :: R         = .05075!-2*.00093!.05564-2.*.00202!.0508!.05075!.053!.07!<radius of drop cm
	real (kind=8), parameter :: R0        = 0.7874 !<radius cylinder cm
	real (kind=8), parameter :: y0        = -1.5!-7.89601/5.0*.78+.4!R0+.4 !<init interface cm
	
!========Experimental parameters of equivalent falling drop problem=====
	!(since Claudia's w code depends on the falling geometry)

	real (kind=8), parameter :: rhot            = 2*rhodrop-rhobottom
		!density of top fluid
	real (kind=8), parameter :: rhob            = 2*rhodrop-rhotop
		!density of bottom fluid
	real (kind=8), parameter :: rhos            = rhodrop
		!density of the sphere
	real (kind=8) :: U                          = 0.0
		!initial velocity of sphere
	real (kind=8), parameter :: g               = 981.0 !gravity cm/s**2
	real (kind=8), parameter :: pi 	            = 2.0*dasin(1.0d0)

!========Numerical Parameters===========================================
	real (kind=8), parameter :: maxTime     = 50. !time to run to
	real (kind=8), parameter :: dt          = .1 !time step
	real (kind=8), parameter :: numtrapz    = 0.005 
		!h for trapezoidal trapz1trapz1
	real (kind=8), parameter :: integthres  = 0.1E-5!simpson integration
	real (kind=8), parameter :: singthres   = 0.1E-6 
		!threshold singularity x=y
	real (kind=8), parameter :: logsing     = 0.1 
		!threshold log "sing" x=y
	real (kind=8), parameter :: logtrapzbig = 0.5
		!h trapz for log around "sing"
	real (kind=8), parameter :: logtrapz    = 0.1 !htrapz for log
	real (kind= 8)           :: FFR         = 2.*R !rad of full soln cm
	real (kind=8) ,parameter :: dx          = 0.1*R 


!=======================Non-Uniform Interface===========================

	integer (kind=4)           ::  RorZ !, XN=ceiling(R0/dx)
	real (kind=8 )             ::  R0cl  	= 2.*R
	integer (kind=4),parameter ::  XNfar 	= 30, XNclose=40
	integer (kind=4)           ::  XN 		= XNclose+XNfar
	!real (kind=8)             ::	 dx = R0/(XNclose+XNfar) 
		!initial grid size interface
	!integer (kind=4)          ::  RorZ, XN=ceiling(R0/dx)




!====================Dependent Parameters===============================
 	real (kind=8) 	:: ms               = 4.0/3.0*pi*R**3*rhos 				
		!mass of the sphere
	real (kind=8)	:: oneoversixpiamuK = (1-2.10444*(R/R0) +2.08877* &
			((R/R0)**3))* 1.0/(6.0*pi*R*mu) !CHANGE THIS TO DROP STUFF
		!reflection correction
	
	real (kind=8)	:: stresspertcoeff = -0.25*g*(rhot-rhob)*R*2.0*pi    	
		!coefficient of the perturbatio stress
	!real (kind=8)	:: stresspert      = 0.0;   !perturbation stress
	real (kind=8) 	:: buoyancytop     = -4.0/3.0*pi*R**3*g*rhot    
	      !buoyant force when sphere is above the interface
	real (kind=8) 	:: buoyancybottom  = -4.0/3.0*pi*R**3*g*rhob   
	      !buoyant force when sphere is below the interface
	real (kind=8) 	:: buoyancyCoeff1  = -pi*g/3.0*(rhob-rhot)           	
		!first coeff for the buoyant force when sphere is in interface
	real (kind=8) 	:: buoyancyCoeff2  = -2.0*pi*g/3.0*R**3*(rhob+rhot)
		!second coeff for the buoyant force when sphere is in interface
	real (kind=8) 	:: drhogover8mu    = (rhob-rhot)*g/(8.0*mu)
		!coefficient for the perturbation flow
	real (kind=8) 	:: myt			   = 0.					


!===============Interface Interpolation Stuff===========================
	real (kind=8), dimension(:), allocatable :: x, y, sx, sy, su ,sv, & 
						&wu ,wv, cinterpx, cinterpy,cwinterpx, cwinterpy
	real (kind=8)                            :: yend, xflagb, xflagl, &
												&px, py, myrho,myzeta, &
												&xvalG2, yvalG2,xvalI2,&
												&yvalI2
	integer (kind=4)                         :: flagb, flagu, flagl, &
								 &flagt, cinternalcount, cwinternalcount
	!flagb is the index of the interface point where the backflow begins 
!===============Fourier Components======================================
	real (kind=8), parameter    :: epsilon		= 0.001
	real (kind=8), parameter    :: LZ    	    = 10.0         
		!domain window of lambda for z-component of cylinder vel
	real (kind=8), parameter    :: LR     	  	= 10.0    
		!domain window of lambda for R-component of cylinder vel
	integer (kind=4), parameter :: NZ      		= 2**14                        
		!number of discretization pts, z direction
	integer (kind=4), parameter :: NR           = 2**16!14
		!number of discretization pts, r direction
	real (kind=8), parameter    :: u3coeff      = -2.1044428*R/R0 + &
												  &2.1800173*R**3/R0**3
		!@fixme I need to change this to the correct droplet expression 
	real (kind=8), parameter    :: upperZ 		= 40.0;
		!> IDK what this is, read it some more
	integer (kind=4), parameter :: upperRangeZ  = ceiling((LZ-epsilon) &
												  *upperZ/(2.0*pi)+1.0)
	integer (kind=4), parameter ::	upperRangeR = ceiling(LR*upperZ/ &
													(2.0*pi)+1.0)
	real (kind=8)               :: cylindervelR(upperRangeR, XNclose + &
				   XNfar+1), cylindervelZ(upperRangeZ, XNclose+XNfar+1)
	real (kind=8)               :: zcoordinateZ(upperRangeZ), &
									&zcoordinateR(upperRangeR)
	complex (kind=8), parameter          :: MINUSONE      = -1.0d0	
	complex (kind=8)                   :: imagi		= zSQRT(MINUSONE)

	real (kind=8)               :: WZ(NZ), WR(NR), myHZ(NZ), myGZ(NZ), &
		 &myHR(NR), myGR(NR), firstpartZ(NZ), firstpartR(NR), myk(NZ)

	real (kind=8)               :: AreaReflux,AreaSpherePortion, &
		&AreaEntrain, startx(XNclose+XNfar+1), starty(XNclose+XNfar+1),&
		&wforceE, wforceR, wforce, ArchBouyancy,ArchBE, stresspert, &
		&stresspertA,ArchBR, stresspertE,stresspertReflux,stresspertG2,&
		&stresspertAG2, stresspertEG2,stresspertRefluxG2,stresspertI2, &
		&stresspertAI2, stresspertEI2,stresspertRefluxI2
		
!===============Stress Coefficient Calculator===========================
    integer(kind=4), parameter  :: nsurfpoints  = 30!5!0 !number of points 
		!to use in integration over surface of drop for stress coeffs.
	real (kind=8), dimension(nsurfpoints+1) :: thetavecforstress, &
		surfx1, surfy1, surfx2,surfx3, surfy2, surfy3, wU1, wV1, wU2, w&
		&V2, wU3, wV3
		!to build the points on drop surface where w will be evaluated.
	real (kind=8)               :: eps          = .00000001 
		!epsilon for stress calculation
    real (kind=8)               :: I2, G2, G4, G6, G8, G2New, G2HO,Unew&
    &forthirdref, v2atorigin, G2reciprocal, I2reciprocal
		!Stress Coeffs
end module globalinfo


!-----------------------------------------------------------------------
!Main Routine for Droplet
!-----------------------------------------------------------------------
!
!PROGRAM: Full Simulation
!
!> @author 
!> H. Arrowood, UNC-CH; C. Falcon, UNC-CH
!
!Description: 
!> Main routine; finds interface, drop vel, etc at each time step
!
!REVISION HISTORY:
!06 Nov 2017 - Initial documented/nicely formatted version
!
!----------------------------------------------------------------------- 
program fullsimulation
    use globalinfo
    implicit none
	
	integer (kind=4)            :: i, ierr, flag, ix, iy, index, iv
		!indices for various do loops
	real (kind=8)               :: velocity(ceiling(maxTime/dt)+1), &
		stresspertvect(ceiling(maxTime/dt)+1),  abserr=0.001, &
		relerr=0.001, xout, yout, xin1, yin1, xin2, yin2
		!store sphere velocity and stress due to perturbation vel at 
		!each timestep
	real (kind=8), dimension(:), allocatable :: V, VP ! V stores 
		!position of interface [xvals, yvals], VP stores velocities 
		![uvals, vvals]
	Character(len=256)           :: filename
	external rhoode !ODE solver for advecting interface


!----Initialize Interface-----------------------------------------------
 	allocate(x(XNfar+XNclose+1), stat=ierr)
		if (ierr /= 0) print*, "x : Allocation failed"
	
	allocate(y(XNfar +XNclose+1), stat=ierr)
		if (ierr /= 0) print*, "y : Allocation failed"

	do i=1,XNclose + XNfar+1 !This actually defines the values of the x 
		!values. More points closer to center where they'll be advected
		if (i<XNclose+2) then
			x (i) = (i-1)**(2) *(R0cl/XNclose**(2))
		else 
			x( i ) = (i - XNclose-1)* (R0 -R0cl ) / XNfar + R0cl
		endif 

	end do 
	y = y0  !initialize all y-values undisturbed -Holly. 

	startx = x;     !save beginning interface
	starty = y;
	
	
!Here's a code block to test out my interpolation routine. 	
	
	!xin1 = (R+.0001)*sin(pi/2+.5)
	!yin1 = (R+.0001)*cos(pi/2+.5)
	!xin2 = (R+.00001)*sin(pi/2+.4)
	!yin2 = (R+.00001)*cos(pi/2+.4)
	!print*, "xin1, yin1", xin1, yin1
	!print*, "xin2,yin2", xin2, yin2
	!call NewtonInterpolation(xin1,yin1,xin2,yin2,xout,yout)

	!print*, "xout", xout
	!print*, "yout", yout
	
!initialize for the second ref spacial component -Holly, Jan 19 2018
	allocate(cinterpx(1000000), stat=ierr)
 		if (ierr /= 0) print*, "cintertpx : Allocation failed"

	allocate(cinterpy(1000000), stat=ierr)
		if (ierr /= 0) print*, "cintertpy : Allocation failed"

!initialize interface interpolated from w (I think this takes care of 
!the values of w on the interface points that are found by interpolation
!rather than from the full integ of the Oseen stuff -Holly,01/19/18)
	allocate(cwinterpx(1000000), stat=ierr)
		if (ierr /= 0) print*, "cintertpx : Allocation failed"

	allocate(cwinterpy(1000000), stat=ierr)
		if (ierr /= 0) print*, "cintertpy : Allocation failed"
!---------------keep track of how many interface points there are, for 
!use in defining various arrays later -Holly----------------------------


	cinternalcount=XN+1;  
	cwinternalcount=XN+1;

    
!------Initialize the interpolated interface points---------------------
	do ix=1, max(cinternalcount,cwinternalcount) 
		cinterpx(ix)  = x(ix);
		cinterpy(ix)  = y(ix);
		cwinterpx(ix) = x(ix);
		cwinterpy(ix)= y(ix);
	end do

!--------Set up files to save entrainment volume and WForce-------------
	open (unit =9,file = 'VolumeTrack.dat')
		write(9,*) "Volume of Entrainment, Volume of Reflux, Volume of &
			&Portion of Sphere"

	open (unit =8,file = 'WForce.dat')
		write(8,*) " ArchBER, Wforce = 6pimuAstresscoeff(wFE - wFR), &
			&SphARchBoyancy, ArchBR,wforceR,wforceE"

	
!----------To Do: set up reflections for the following------------------	
!----------Compute spacial component of second reflection on a grid, to 
!be used to interpolate second reflection values later on---------------
	  
	call cylindervelinit() 
	
!----------Initialize Stress Coefficients-------------------------------	
    I2 = 0.0
    G2 = 0.0
    G4 = 0.0
    G6 = 0.0
    G8 = 0.0
    G2New = 0.0
    G2HO = 0.0

 
 !call NewtonInterpolation(.5001*sin(Pi/3),.5001*cos(Pi/3), .5002*sin(pi/3-.01), .5002*cos(pi/3-.01),xout, yout)
!print*, "xout", xout
!print*, "yout", yout
!-----------------------------------------------------------------------
!**********************  TIME LOOP *************************************
!-----------------------------------------------------------------------

	
	do index = 0,ceiling(maxTime/dt)

!allocate all the positions and velocities. XN changes at each timestep,
!so this has to be inside the time loop 
	allocate(V(2*(XN+1)), stat=ierr)
		if (ierr /= 0) print*, "V : Allocation failed"
	allocate(VP(2*(XN+1)), stat=ierr)
		if (ierr /= 0) print*, "VP : Allocation failed"
		
!initialize stokes flow
	allocate(su(XN+1), stat=ierr)
		if (ierr /= 0) print*, "x : Allocation failed"
	allocate(sv(XN+1), stat=ierr)
		if (ierr /= 0) print*, "y : Allocation failed"
		su = 0 
		sv = 0 


!initialize w flow
	allocate(wu(XN+1), stat=ierr)
		if (ierr /= 0) print*, "x : Allocation failed"
	allocate(wv(XN+1), stat=ierr)
		if (ierr /= 0) print*, "y : Allocation failed"
		wu = 0
		wv = 0

!Stores the interface x and y coordinates in a single vector -Holly		
		V(1:XN+1)        = x        
		V(XN+2:2*(XN+1)) = y     


	print*, "current number of interface points", XN+1
	print*, index
	
!===================== Write  interface ================================

	write (filename, fmt='(a,f10.2,a)') 'interface',real(index+1),'.dat'

	open (unit =2,file = filename,form='formatted')
		write(2,*) "x,y at time=", myt, "rhos", rhos

		do ix=1, XN+1
			write(2,*), x(ix), ",", y(ix),","
		end do


!===================== Write interpolated interface ====================


	write (filename, fmt='(a,f10.2,a)') 'interpolation',real(index +1),&
	&'.dat'
	open (unit =6,file = filename,form='formatted')
!open (unit =6,file = 'interpolation.dat',form='formatted')
	write(6,*) "interpolated interface x,y at time=", myt

		do ix=1, max(cinternalcount,cwinternalcount)
			write(6,*), cinterpx(ix), ",", cinterpy(ix),",", &
			&cwinterpx(ix), ",", cwinterpy(ix),","
		end do



!====================== ODE SOLVER =====================================

allocate(sx(XN+1), stat=ierr)
		if (ierr /= 0) print*, "sx : Allocation failed"
	
	allocate(sy(XN+1), stat=ierr)
		if (ierr /= 0) print*, "sy : Allocation failed"

!print*, "iteration #", index, "right before ode solver"
	flag = 1
	print*, "about to call solver"
	call r8_rkf45 (rhoode, 2*(XN+1), V, VP, index*dt, (index+1.0)*dt, &
		&relerr, abserr, flag ) !this is solves the ode to advect the 
		!interface, with the Runge Kutta method found in rk4.f90 -Holly
	print*, "solver done"
	x = V(1:XN+1)   !unwrap interface positions from V
	y = V(XN+2:2*(XN+1))



!================Write data files=======================================


	write (filename, fmt='(a,f10.2,a)') 'stokes',real(index+1),'.dat'

!print *, filename

	open (unit =4,file = filename,form='formatted')
		write(4,*) "us,sv at time=", myt

		do ix=1, XN+1
			write(4,*), su(ix), ",", sv(ix), ","
		end do


	write (filename, fmt='(a,f10.2,a)') 'wpert',real(index+1),'.dat'

!print *, filename

	open (unit =5,file = filename,form='formatted')
		write(5,*) "wu,wv at time=", myt

		do ix=1, XN+1
			write(5,*), wu(ix), ",", wv(ix), ","
		end do

!print *, filename

	open (unit =9,file = 'VolumeTrack.dat')
		write(9,*), AreaEntrain, ",", AreaReflux, ",", AreaSpherePortion


	open (unit =8,file = 'WForce.dat')
		write(8,*), ArchBE, ",", wforce, ",", ArchBouyancy, ",", Arch&
		&BR,",",wforceR,",", wforceE


	open (unit =11,file = 'sphereVel.dat')
		write(11,*), U, ",", myt, ",", yend

print*, "right before fillgaps in do loop"
    call fillgaps() !add more interface points as needed.
    
    velocity(index+1) = U
	stresspertvect(index+1)=stresspert

	print *, U, ",", myt, ",", yend
	
!=====+==Deallocate velocities and positions============================
    if (allocated(V)) deallocate(V,stat=ierr)
    if (allocated(VP)) deallocate(VP,stat=ierr)


    if (allocated(wu)) deallocate(wu,stat=ierr)
	if (allocated(wv)) deallocate(wv,stat=ierr)

	if (allocated(su)) deallocate(su,stat=ierr)
	if (allocated(sv)) deallocate(sv,stat=ierr)
		
	if (allocated(sx)) deallocate(sx, stat=ierr)
	if (allocated(sy)) deallocate(sy, stat=ierr)

	end do
!======================End of Time Loop=================================
	
!=======Print velocities and stresses at the end========================
	print *, "velocity"
		do iv = 1, ceiling(maxTime/dt)+1
			print *, velocity(iv), ","
		end do
	
	print *, "stress"
		do iv = 1, ceiling(maxTime/dt)+1
			print *, stresspertvect(iv), ","
		end do
!=======Print interface points at the end===============================
	print *, " "
	print *, "*********************************************************"
	print *, "x"
		do ix=1, XN+1
			print *, x(ix), ","
		end do
	print *, "y"
		do iy=1, XN+1
			print *, y(iy), ","
		end do
	print *, "*********************************************************"
	print *, " "
end program fullsimulation


!-----------------------------------------------------------------------
!> @author
!> H. Arrowood, UNC-CH; based on work of C. Falcon, UNC-CH
!
! DESCRIPTION:
!> Computes all necessary velocities on interface to feed into the Runge
!!Kutta Fehlberg solver that advects the interface

! REVISION HISTORY:
! 06 Nov 2017 - Added documentation
!
!>@return V   Points on new interface
!>@return VP  Velocities at interface points
!-----------------------------------------------------------------------


subroutine rhoode(T, V, VP)
	use globalinfo
	implicit none
	
	real (kind=8)     :: T, sr(XN+1), V(2*(XN+1)), VP(2*(XN+1)), sunew,&
	& svnew 
		!time, r-values of points, positions, velocities
	integer (kind=4)  :: ierr, i, ix, iy
		! I don't know what ierr does exactly, but i, ix, and iy are 
		!just various indices
	myt = T
	print*, "time", T
	
!-----Unwrap interface points from V------------------------------------
	sx = V(1:XN+1)  
	sy = V(XN+2:2*(XN+1))
    
!------This sets up the integrands for w calculation; determines where 
!the perturbed interface intersects the original interface I think------

	call setSpecialPositions(sx, sy, XN, R, flagb, flagu, flagl, flagt,&
	& xflagb, xflagl, yend)
!------Find perturbation velocity---------------------------------------
   print*,"XN", XN
	!call wStressTN()	
	call wTN() !wu,wv
!------computes  coefficients of stresses, using velocities computed on
! the drop surface in wTN()---------------------------------------------
	call stresscoeffs() 
	
	!print*, "done with new stress"
!----Find velocity of sphere based on w, ustokes of previous timestep---
	!print*, "about to call sphvel"
	call sphvel() !This finds U -Holly
	!print*, "finished calling sphvel"
    print*,"U outside sphvel", U
    
!-----------Find Stokes velocity in a cylinder, uses stress coeffs------
	print*, "about to call stokes"
	call stokes() !su, sv
	print*, "finished calling stokes"
!------Set velocity to zero if interface inside the drop----------------
	! @fixme: Shouldn't this be an error message instead?
	sr = dsqrt(sx**2+sy**2)
	do i = 1, XN+1
		if (sr(i) <= R*1.000000001) then 
			wu(i) = 0.0
			wv(i) = 0.0
			!sunew = 0.0
			!svnew = 0.0
			!sunew = cos(atan(sx(i)/sy(i)))*(cos(atan(sx(i)/sy(i)))*su(i)-sin(atan(sx(i)/sy(i)))*sv(i)) !kill radial velocity, keep tangential velocity
			!svnew = sin(atan(sx(i)/sy(i)))*(cos(atan(sx(i)/sy(i)))*su(i)-sin(atan(sx(i)/sy(i)))*sv(i))
			!su(i) = sunew!0.0
			!sv(i) = svnew!0.0
		elseif(i==XN+1) then
			wu(i) = 0.0
			wv(i) = 0.0
		print*, "fake velocity stuff happening"
		endif
	end do

!----Put together full velocities--------------------------------------
	VP(1:XN+1) = su+wu 
	VP(XN+2:2*(XN+1)) = (sv+wv)
!print*, "x vel", VP(1:XN+1)
!print*, "y vel", VP(XN+2:2*(XN+1))
!	print*,"end of rhoode"
end subroutine rhoode



!-----------------------------------------------------------------------
!> @author
!>H. Arrowood, UNC-CH, C. Falcon, UNC-CH
!
! DESCRIPTION:
!> Balances stokes and perturbation forces with buoyancy force to find
!! velocity of drop
!
! REVISION HISTORY
! 7 Nov 2017 - added documentation
! 2 Jan 2018 - added indirect computation of I2


subroutine sphvel()
	use globalinfo
	implicit none
	
	real (kind=8) :: buoyancy, stressbelowsphere, stresssidesphere, &
		stressabovesphere, stressbackflow, stressbelowsphereA, &
		stresssidesphereA, stressabovesphereA, stressbackflowA, &
		stressbelowsphereE, stresssidesphereE, stressabovesphereE, &
		stressbackflowE,cindex, ArchBuoyancyDrop, IndirectI2
	real (kind=8) :: HabermanStressPert, stressdiff, Ureflections !to test out my scheme	
	real (kind=8), external :: stresstail1D, stressIntegrandFlat1D,&
		stressIntegrandsphere1D, stressIntegrand1D,stresstail1DA, &
		stressIntegrandFlat1DA, stressIntegrandsphere1DA, &
		stressIntegrand1DA, stresstail1DE, stressIntegrandFlat1DE,&
		stressIntegrandsphere1DE, stressIntegrand1DE
	real (kind=8) :: stressbelowsphereG2, stresssidesphereG2, &
	&stressabovesphereG2, stressbackflowG2, stressbelowsphereAG2, &
	&stresssidesphereAG2, stressabovesphereAG2, stressbackflowAG2, &
	&stressbelowsphereEG2, stresssidesphereEG2, stressabovesphereEG2, &
	&stressbackflowEG2
	real (kind=8), external :: stresstail1DG2, stressIntegrandFlat1DG2,&
	&stressIntegrandsphere1DG2, stressIntegrand1DG2,stresstail1DAG2, &
	&stressIntegrandFlat1DAG2, stressIntegrandsphere1DAG2, &
	&stressIntegrand1DAG2, stresstail1DEG2, stressIntegrandFlat1DEG2,&
	&stressIntegrandsphere1DEG2, stressIntegrand1DEG2
	Character(len=256):: filename
	integer i,ix, ierr	!I'm pretty sure these are taken care of in 
		!rhoode, this is redundant
	
!-----buoyancy for a two layer fluid; taking care of cases where drop is
!in the interface------------------------------------------------------
	if (sy(XN+1)>=R) then
		buoyancy = buoyancybottom;
	elseif (abs(sy(XN+1))<R) then
		buoyancy = buoyancyCoeff1*(3.0*R**2*sy(XN+1)-sy(XN+1)**3)+ &
		buoyancyCoeff2;
	else
		buoyancy = buoyancytop;
	endif

!-----initialize interpolated interface points--------------------------

	if (allocated(cinterpx)) deallocate(cinterpx, stat=ierr)
	if (allocated(cinterpy)) deallocate(cinterpy, stat=ierr)


	allocate(cinterpx(1000000), stat=ierr)
	if (ierr /= 0) print*, "cintertpx : Allocation failed"

	allocate(cinterpy(1000000), stat=ierr)
	if (ierr /= 0) print*, "cintertpy : Allocation failed"

	cinternalcount=0.0;
	cinterpx =0;
	cinterpy = 0;


	!calculate stress force
	stresspert = 0.0;
	if (maxval(sy) > minval(sy)) then
		call setSpecialPositions(sx, sy, XN, R, flagb, flagu, flagl, fl&
		&agt, xflagb, xflagl, yend)
			
		stressbelowsphere = 0.0;
		stresssidesphere = 0.0;
		stressabovesphere = 0.0;
		stressbackflow = 0.0;
		
		if  (flagu /= 0) then
			!call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), &
			!&integthres,stressbelowsphere)
			!call simp(stressIntegrandsphere1D, max(-R, sy(1)), R, &
			!&integthres, stresssidesphere)
			!call simp(stressIntegrand1D, R, yend, integthres, &
			!&stressabovesphere)
			call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), numtr&
			&apz, stressbelowsphere)
			call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), R, num&
			&trapz, stresssidesphere)
			call trapz1(stressIntegrand1D, R, yend, numtrapz, stressabo&
			&vesphere)

		elseif  (flagl /= 0) then
			!call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), &
			!&integthres, stressbelowsphere)
			!call simp(stressIntegrandsphere1D, max(-R, sy(1)), yend, &
			!&integthres, stresssidesphere)

			call trapz1(stressIntegrand1D, sy(1), max(-R, sy(1)), numtr&
			&apz, stressbelowsphere)
			call trapz1(stressIntegrandsphere1D, max(-R, sy(1)), yend, &
			&numtrapz, stresssidesphere)

		else
			!call simp(stressIntegrandFlat1D, sx(1), xflagb, integthres,&
			!&stressbelowsphere)
			call trapz1(stressIntegrandFlat1D, sx(1), xflagb, numtrapz,&
			& stressbelowsphere)
			

		endif

		if (flagb /= XN+1) then
			call trapz1(stresstail1D, xflagb, sx(XN+1),real(0.01,kind&
			&=8), stressbackflow)
		endif
	stresspert  = stresspertcoeff*(stressbelowsphere + stresssidesphere&
	& + stressabovesphere -stressbackflow)


!------------Calculate G2-----------------------------------------------
!Uses same integration routine as Claudia's above, but only part of the integrand


		stresspertG2 = 0.0;
			
		stressbelowsphereG2 = 0.0;
		stresssidesphereG2 = 0.0;
		stressabovesphereG2 = 0.0;
		stressbackflowG2 = 0.0;
		
		if  (flagu /= 0) then
			!call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), &
			!&integthres,stressbelowsphere)
			!call simp(stressIntegrandsphere1D, max(-R, sy(1)), R, &
			!&integthres, stresssidesphere)
			!call simp(stressIntegrand1D, R, yend, integthres, &
			!&stressabovesphere)
			call trapz1(stressIntegrand1DG2, sy(1), max(-R, sy(1)), num&
			&trapz, stressbelowsphereG2)
			call trapz1(stressIntegrandsphere1DG2, max(-R, sy(1)), R, n&
			&umtrapz, stresssidesphereG2)
			call trapz1(stressIntegrand1DG2, R, yend, numtrapz, stressa&
			&bovesphereG2)

		elseif  (flagl /= 0) then
			!call simp(stressIntegrand1D, sy(1), max(-R, sy(1)), &
			!&integthres, stressbelowsphere)
			!call simp(stressIntegrandsphere1D, max(-R, sy(1)), yend, &
			!&integthres, stresssidesphere)

			call trapz1(stressIntegrand1DG2, sy(1), max(-R, sy(1)), &
			&numtrapz, stressbelowsphereG2)
			call trapz1(stressIntegrandsphere1DG2, max(-R, sy(1)), yend&
			&, numtrapz, stresssidesphereG2)

		else
			!call simp(stressIntegrandFlat1D, sx(1), xflagb, integthres,&
			!&stressbelowsphere)
			call trapz1(stressIntegrandFlat1DG2, sx(1), xflagb, numtrap&
			&z, stressbelowsphereG2)
			

		endif

		if (flagb /= XN+1) then
			call trapz1(stresstail1DG2, xflagb, sx(XN+1),real(0.01,kind&
			&=8), stressbackflowG2)
		endif

		stresspertG2  = 2.*pi*(rhotop-rhobottom)*(stressbelowsphereG2 +&
		& stresssidesphereG2 + stressabovesphereG2 +stressbackflowG2)!*g

		print*, "below", stressbelowsphereG2
		print*, "beside", stresssidesphereG2
		print*, "above", stressabovesphereG2
		print*, "backflow", stressbackflowG2

		
	endif

G2reciprocal = (3./(4.*R**2*pi*mu))*stresspertG2!*3./(4.*R**2*pi*mu)



!-------calculate archimedean force of fluid----------------------------
	stresspertA = 0.0;
	if (maxval(sy) > minval(sy)) then
		call setSpecialPositions(sx, sy, XN, R, flagb, flagu, flagl, &
		&flagt, xflagb, xflagl, yend)
		stressbelowsphereA = 0.0;
		stresssidesphereA = 0.0;
		stressabovesphereA = 0.0;
		stressbackflowA = 0.0;
		if  (flagu /= 0) then
			call simp(stressIntegrand1DA, sy(1), max(-R, sy(1)), &
			&integthres, stressbelowsphereA)
			call simp(stressIntegrandsphere1DA, max(-R, sy(1)), R,&
			& integthres, stresssidesphereA)
			call simp(stressIntegrand1DA, R, yend, integthres, &
			&stressabovesphereA)

		elseif  (flagl /= 0) then
			call simp(stressIntegrand1DA, sy(1), max(-R, sy(1)), &
			&integthres, stressbelowsphereA)
			call simp(stressIntegrandsphere1DA, max(-R, sy(1)), yend, &
			&integthres, stresssidesphereA)
		else
			call simp(stressIntegrandFlat1DA, sx(1), xflagb, integth&
			&res, stressbelowsphereA)
		endif

		if (flagb /= XN+1) then
			call simp(stresstail1DA, xflagb, sx(XN+1), integthres, &
			&stressbackflowA)
		endif
	stresspertA  = -g*(rhob-rhot)*(stressbelowsphereA + stresssidesphe&
	&reA + stressabovesphereA-stressbackflowA)
	endif


!-----------calculate archimedean force of fluid------------------------
	stresspertE = 0.0;
	if (maxval(sy) > minval(sy)) then
		stressbelowsphereE = 0.0;
		stresssidesphereE = 0.0;
		stressabovesphereE = 0.0;
		stressbackflowE = 0.0;

		if  (flagu /= 0) then
			call simp(stressIntegrand1DE, sy(1), max(-R, sy(1)), &
			&integthres, stressbelowsphereE)
			call simp(stressIntegrandsphere1DE, max(-R, sy(1)), R, &
			&integthres, stresssidesphereE)
			call simp(stressIntegrand1DE, R, yend, integthres, stress&
			&abovesphereE)

		elseif  (flagl /= 0) then
			call simp(stressIntegrand1DE, sy(1), max(-R, sy(1)), integ&
			&thres, stressbelowsphereE)
			call simp(stressIntegrandsphere1DE, max(-R, sy(1)), yend,&
			& integthres, stresssidesphereE)
		else
		call simp(stressIntegrandFlat1DE, sx(1), xflagb, integthres, &
		&stressbelowsphereE)
		endif

		if (flagb /= XN+1) then
			call simp(stresstail1DE, xflagb, sx(XN+1), integthres, &
			&stressbackflowE)
		endif

	stresspertE = g*(rhob-rhot)*(stressbelowsphereE + stresssidesphereE&
	& + stressabovesphereE-stressbackflowE)
	endif



!------Now put everything together--------------------------------------

	wforceE= stressbelowsphere + stresssidesphere + stressabovesphere
	wforceR=stressbackflow
	wforce = oneoversixpiamuK*stresspert

	ArchBouyancy =oneoversixpiamuK*(g*ms + buoyancy)
	ArchBE=-stresspertA *oneoversixpiamuK
	ArchBR =-g*(rhob-rhot)*(-stressbackflowA)*oneoversixpiamuK
	stresspertReflux= -g* ( rhob-rhot)*stressbackflowE
	

	
!------stresspert in the case of the sphere, for checking computation of
!G2 and I2. The stress due purely to the perturbation velocity should 
!be the same in both cases, though of course the stokes stress is also 
!dependent of the pert vel stress coeffs. So I can use this setup to 
!indirectly compute I2 just as in the sphere case, then end by 
!assembling all the pieces for the drop case---------------------------- 
   HabermanStressPert = 4.0*pi*R**2*(mu/3.0)*(G2-I2)
	print*, "*************************************************************"
	print*,"Values for Force on Sphere"
	!print*, "New Force:        ", stresspertG2+stresspertI2
    print*, "Claudia's Force:  ", stresspert 
    print*, "Haberman Force:   ", HabermanStressPert
	print*, "*************************************************************"
	print*, "Tangential Components of Force"
	print*, "New:    ",stresspertG2
	print*, "Hab:    ", 4./3.*R**2*pi*mu*G2
	print*, "Normal Components of Force"
	!print*, "New:    ", stresspertI2
	print*, "Hab:    ", -4./3.*R**2*pi*mu*I2
	print*, "*************************************************************"
	print*, "G2 Values"
	print*, "G2 direct:   ", G2
    print*, "G2 HO:       ", G2HO
!	print*, "G2 analytic: ", G2New
	print*, "G2reciprocal:", G2reciprocal


	  
    IndirectI2 = G2reciprocal-3/(4*pi*R**2*mu)*stresspert!G2New-3/(4*pi*R**2*mu)*stresspert
 
    !print*, "stresspertcoeff", stresspertcoeff
    print*, "Indirectly computed I2", IndirectI2
	print*, "I2", I2
    
	ArchBuoyancyDrop = (g*ms+buoyancy)!oneoversixpiamuK!buoyancy takes care 
		!of the density variation. So basically g*ms+buoyancy is what 
		!I'm calling "buoyancy force" in my notes oneoversixpiamuK is 
		!the reflection force correction stuff, I'll need to fix it 

!---terminal velocity for free space drop stuff, see p. 99 of notes-----
	!U = (mu+muin)/(2*mu+3*muin)*1/(2*R*pi*mu)*(ArchBuoyancyDrop-2*R*pi*mu/&
	!&(3*(mu+muin))*(R*(mu*(-2*IndirectI2+G2)+2*muin*(-IndirectI2+G2))))
!---velocity for stratified drop, see p. 62 and the Mathematica Notebook
! "ApproximateThirdReflection.nb"-------------------------------------- 

	G2 = -G2reciprocal!G2!New!-G2 !made this positive just to see how it looks...
	I2 = -IndirectI2!"
	U = -((-ArchBuoyancyDrop+(2.0943951023931953d0*R**2.0*mu*((G2-2.0*I&
	&2)*mu+2.0*(G2-1.0*I2)*muin))/(mu+muin))/(-((2.0943951023931953d0*&
	&R*mu*(6.0*mu+9.0*muin))/(mu+muin))+(2.0*R**2*mu*(mu+1.5*muin)*&
	&((6.515546450454206d0*R**2*muin)/mu+R0**2*(-8.566226429811636d0&
	&-(12.849339644717457d0*muin)/mu)))/(R0**3*(mu+muin)*(1.0+muin&
	&/mu))))
		

	UReflections = -((4.1887902047863905d0*R**3*(-981.0d0)*(-rhobottom+&
	&rhos))/(-((12.566370614359172d0*R*mu*(mu+1.5*muin))/(mu +muin))+(&
	&2.*R**2*mu*(mu+1.5*muin)*((6.515546450454206d0*R**2*muin)/mu+R0**2&
	&*(-8.566226429811636d0-(12.849339644717457d0*muin)/mu)))/(R0**3*(m&
	&u+muin)*(1.+muin/mu))))
		
	v2atorigin = 1./(2.*Pi)*(1./(R0**3*(1.+(muin/mu)))*R*U*(R0**2*(-8.5&
	&66226429811636-12.849339644717457*(muin/mu))+6.515546450454206*R**&
	&2*(muin/mu))) !this has been checked! -H.
	
	print*, "Uneater", (ArchBuoyancyDrop+stresspert+2.*R**2*mu**2*G2/&
	&(3.*(mu+muin)))*(mu+muin)/(mu+1.5*muin)/(4*pi*R*mu)/(1+ 1./&
	&(2.*Pi)*(1./(R0**3*(1.+ (muin/mu)))*R*(R0**2*(-8.566&
	&226429811636-12.849339644717457*(muin/mu))+6.515546450454206*R**2&
	&*(muin/mu))))
	
	Unewforthirdref=U+ v2atorigin
	!UClaudia = oneoversixpiamuK*(g*ms + buoyancy + stresspert)
	print*,"U", U
	print*,"URefBottom", UReflections
	print*, "URefTop", -((4.1887902047863905d0*R**3*(-981.0d0)*(-rhotop&
	&+rhos))/(-((12.566370614359172d0*R*mu*(mu+1.5*muin))/(mu +muin))+(&
	&2.0*R**2*mu*(mu+1.5*muin)*((6.515546450454206d0*R**2*muin)/mu+R0**&
	&2*(-8.566226429811636d0-(12.849339644717457d0*muin)/mu)))/(R0**3*(&
	&mu+muin)*(1.+muin/mu))))
	!print*, "UClaudia", UClaudia
end subroutine sphvel

!-----------------------------------------------------------------------
!> @author
!> C. Falcon, UNC-CH
!
!DESCRIPTION:
! Determines a number of special positions on the interface that define
!! the intersections of the current interface with the unperturbed
!! interface. See Claudia's notebook for details.
!
! @param myx the x positions of the interface
! @param myy the y positions of the interface
!
! Does not return anything, but sets the following variables that
! define the aforementioned intersection cases and points that are
! contained in globalinfo.

! Flags that define which cases are active. All integers with 0=inactive,
! 1=active:
! - flagb
! - flagu
! - flagl
! - flagt
! Flags that define where the specific points are in the case that they
! are active:
! - xflagb : x value of point of intersection of perturbed interface with original interface
! - xflagl: The point on the interface w/ same y val as bottom of sphere
! - yend
subroutine setSpecialPositions(myx, myy, XN, R, flagb, flagu, flagl,&
& flagt, xflagb, xflagl, yend)
	!determine special positions on the interface
	implicit none
	
	! In 
	integer (kind=4), intent(in) 	:: XN 
	real (kind=8)				 	:: myx(XN+1), myy(XN+1)
	real (kind=8), intent(in) 	 	:: R
	
	! Out
	integer (kind=4), intent(out)	:: flagb, flagu, flagl,flagt
	real (kind=8), intent(out) 		::  xflagl, yend, xflagb
	
	! Local
	integer (kind=4)   				:: temp(1), tempmaxi
	
	flagl = 0
	flagu = 0
	 
!----find position of backflow; that is, bit that's higher than starting
	yend = myy(XN+1)
	
	if(maxval(myy) > yend) then
		!flagb should be the index of the last point below the original interface
		flagb = XN+1  
		do tempmaxi = 1, XN
			if (myy(tempmaxi+1) > yend .and. myy(tempmaxi) <= yend) then
				flagb = tempmaxi
				exit
			endif
		end do
		!tempmaxi will now hold the starting index of the window of length
		!5 around the intersection of the new interface with the starting position
		tempmaxi = min(flagb-2, XN-3)
		tempmaxi = max(tempmaxi, 1)
		!print *, 'sp 1'
		!finds the x value where the interface crosses the original interface 
		!within the above window
		call interpbridge(5, myy(tempmaxi:tempmaxi+4), myx(tempmaxi: &
			&tempmaxi+4),  yend, xflagb)
		!print *, 'sp 2'
	else
	   flagb = XN+1;
	   xflagb = myx(XN+1);
	endif 
	 
!-----find x position of bottom of sphere-------------------------------
	
	if (yend >= -R) then	
		temp = minloc(abs(myy+R))
		flagl = temp(1)
		if (myy(1) >= -R) then
			xflagl = 0
		else
			!print *, 'sp 3'
             print*, min(flagb+2,XN+1), "N test1"
			call interpbridge(min(flagb+2, XN+1), myy(1:min(flagb+2, &
				&XN+1)), myx(1:min(flagb+2, XN+1)), -R, xflagl)
			!print *, 'sp 4'
		endif
	endif
	 
!------find x position of top of sphere---------------------------------
	
	if (yend >= R) then	
    
    		!flagu=is 1 if interface is past sphere top
		flagu = 1;
        
		temp = minloc(abs(myy-R))
		flagt = temp(1)
	endif
end subroutine setSpecialPositions


!-----------------------------------------------------------------------
!> @author
!> C. Falcon, UNC-CH
!
! DESCRIPTION:
!> Adds in more interface points afer each time step to fill gaps. Takes 
!! in x and y values of interface points, returns new set of x and y 
!! values with the gaps filled in using cubic spline interpolation.
!
! REVISION HISTORY:
! 07 Nov 2017 - Added documentation and standardized - H. Arrowood
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @author
!> C. Falcon, UNC-CH
!
! DESCRIPTION:
!> Adds in more interface points afer each time step to fill gaps. Takes 
!! in x and y values of interface points, returns new set of x and y 
!! values with the gaps filled in using cubic spline interpolation.
!
! REVISION HISTORY:
! 07 Nov 2017 - Added documentation and standardized - H. Arrowood
!-----------------------------------------------------------------------

subroutine fillgaps()
	use globalinfo
	implicit none
	
	real (kind=8), dimension(:), allocatable :: newx, newy, newx1, newy1 !working x
		!and y arrays
	real (kind=8)                            :: theta, dist, newpt,nx,ny !distance 
		!between consecutive interface points, new x/y val from spline
	integer (kind=4)                         :: internalcount, ierr, &
		&xi, posi, hi, ii, counter1, counter2 !various indices
	

!----Kill points within 1E-8 of the interface-------------------------
!-----Count number of points to remove
!	counter1 = 0
!	do hi=1,XN+1
!		if(sqrt(x(hi)**2+y(hi)**2) < (R+1E-8)) then
!			counter1=counter1+1
!		elseif((sqrt((x(hi+1)-x(hi))**2+(y(hi+1)-y(hi))**2)<.01 .and. hi>3))then
!			counter1=counter1+1
!		endif
!	end do
!---Allocate newx1, newy1 interface points-------------------------------
	
	!initialize new interface
!	allfillgapsocate(newx1(XN+1-counter1), stat=ierr)
!		if (ierr /= 0) print*, "newx1 : Allocation failed"
	
!	allocate(newy1(XN+1-counter1), stat=ierr)
!		if (ierr /= 0) print*, "newy1 : Allocation failed"
		
!-----Put all x, y points satisfying our criteria into the newx1 and 
!newy1 arrays-----------------------------------------------------------
!	counter2=0
!	do ii=1,XN+1
!		if(sqrt(x(ii)**2+y(ii)**2) < (R+1E-8)) then
!		  counter2=counter2+1
!		  print*, "killing points inside drop"
!		elseif(((sqrt((x(ii+1)-x(ii))**2+(y(ii+1)-y(ii))**2)<.01).and. ii>2))then
!			counter2=counter2+1
!			print*,"killing points that are too close"
!		else
!			newx1(ii-counter2)=x(ii)
!			newy1(ii-counter2)=y(ii)
!		endif
!	end do
	
!-----Make x and y without the bad points-------------------------------
 !   XN = XN-counter1
	
!	x   = newx1
!	y   = newy1
	
	!initialize new interface
	allocate(newx(2*(XN+1)), stat=ierr)
		if (ierr /= 0) print*, "newx : Allocation failed"
	
	allocate(newy(2*(XN+1)), stat=ierr)
		if (ierr /= 0) print*, "newy : Allocation failed"

	call setSpecialPositions(x, y, XN, R, flagb, flagu, flagl, flagt,&  
	& xflagb, xflagl, yend)
    
	internalcount   = 0.0
	
!---Check distances between consecutive points and add a point if larger
!than dx----------------------------------------------------------------	
   	do xi=1,XN
		
			internalcount       = internalcount+1.0
			newx(internalcount) = x(xi)
			newy(internalcount) = y(xi)
		    
			dist = sqrt((x(xi+1)-x(xi))**2+(y(xi+1)-y(xi))**2)
			theta = atan(x(xi)/y(xi))

			if  ((sqrt(x(xi)**2+y(xi)**2) < (2*R) .and. dist > dx) .or.&!Two cases in which interpolation happens, either close to drop and dist>dx or far away and dist>R/2. 
				&dist >R/2) then
				if (sqrt(x(xi)**2+y(xi)**2) < R*1.05 .and. sqrt(x(xi+1)&
				&**2+y(xi+1)**2) < R*1.05) then
					print*, "interpolation very close to drop happening"
!Case1: 				
		               internalcount           = internalcount+1
 						!Call my interpolation routine instead of hers
						print*, "Calling Newton Interpolation"
						call NewtonInterpolation(x(xi), y(xi), x(xi+1),&
						& y(xi+1), nx, ny)
		             
							newx(internalcount) = nx!newpt
							newy(internalcount) = ny
						
					
				else
						print*, "Claudia's interpolation happening"
						!cubic interpolation to find point to fill gap.
						internalcount           = internalcount+1
					if(x(xi) > 2.0*R) then
						newx(internalcount) = 0.5*(x(xi+1)+x(xi))
						posi = min(XN+1.0, xi+3.0)
							!print *, "fillgap1"
						call interpbridge( 7, x(posi-6:posi), y(posi-6:&
						&posi), 0.5*(x(xi+1)+x(xi)), newpt)
							!print *, "fillgap2"
						newy(internalcount) = newpt
					else
						print*, "interpolation close to drop happening"
							posi = max(xi-3, 1)
							if(flagl /= 0.0 .and. xi > flagl) then
								newy(internalcount)=0.5*(y(xi+1)+y(xi))

								!print *, "fillgap3"
					if (y(xi+1) >= y(xi)) then 
						call interpbridge( 7, y(posi:posi+6), x(posi:po&
						&si+6), 0.5*(y(xi+1)+y(xi)), newpt)
					else
						call interpbridge( 7, y(posi+6:posi:-1), x(posi&
						&+6: posi:-1), 0.5*(y(xi+1)+y(xi)), newpt)
					endif
								!print *, "fillgap4"
					if (newpt <= max(x(xi+1), x(xi)) .and. newpt >= min&
					&(x(xi+1), x(xi))) then
						newx(internalcount) = newpt
					elseif (((0.5*(y(xi+1)+y(xi)))**2 + (0.5*(x(xi+1)+x&
					&(xi)))**2) > R**2) then
						newx(internalcount) = 0.5*(x(xi+1)+x(xi))
					else
						newx(internalcount) = sqrt((R+1E-5)**2 - (0.5*&
						&(y(xi+1)+y(xi)))**2)
					endif
					else
						newx(internalcount) = 0.5*(x(xi+1)+x(xi))
						call interpbridge( 7, x(posi:posi+6), y(posi:po&
						&si+6), 0.5*(x(xi+1)+x(xi)), newpt)

						if (newpt <= max(y(xi+1), y(xi)) .and. newpt >=&
						& min(y(xi+1), y(xi))) then
							newy(internalcount) = newpt
						elseif (((0.5*(y(xi+1)+y(xi)))**2 + (0.5*(x(xi+&
						&1)+x(xi)))**2) > R**2) then
									newy(internalcount) = 0.5*(y(xi+1)+&
									&y(xi))
						else
									newy(internalcount) = sqrt((R+1E-5)&
									&**2 - (0.5*(x(xi+1)+x(xi)))**2)
						endif
					endif
				endif
			endif
		endif
	
	end do
	newx(internalcount+1) = x(XN+1)
	newy(internalcount+1) = y(XN+1)

	if (allocated(x)) deallocate(x,stat=ierr)
	if (allocated(y)) deallocate(y,stat=ierr)
	
	XN = internalcount
	allocate(x(XN+1), stat=ierr)
		if (ierr /= 0) print*, "fillgap - x : Allocation failed"
	
	allocate(y(XN+1), stat=ierr)
		if (ierr /= 0) print*, "fillgap - y : Allocation failed"
	
	x   = newx(1:internalcount+1)
	y   = newy(1:internalcount+1)
    
	if (allocated(newx)) deallocate(newx,stat=ierr)
	if (allocated(newy)) deallocate(newy,stat=ierr)
end subroutine fillgaps

!-----------------------------------------------------------------------
!> @author
!> C. Falcon, UNC-CH
!
!DESCRIPTION:
! Takes in existing interface points
! Takes in an xval for a point to be interpolated, and returns the 
! interpolated y value
!
! @param N the number of points
! @param interpx x coordinates of existing interface points
! @param interpy y coordinates of existing interface points
! @param xval the x value of the point to be interpolated
! @return yval the y value of the point to be interpolated
!
! REVISION HISTORY:
! 07 Nov 2017 - Added documentation and standardized - H. Arrowood
!-----------------------------------------------------------------------

subroutine interpbridge(N, interpx, interpy, xval, yval)
	use globalinfo
	implicit none
	! In/out variables
	integer (kind=4), intent(in) :: N
	real (kind=8) :: interpx(N), interpy(N), xval
	real (kind=8), intent(out) :: yval
	! Working variables
	integer (kind=4) :: setmin(1), mini, maxi, tempi, tempj, &
		interpchecki=1
	real (kind=8) :: d(N), checkorder(N-1)
	real (kind=8) :: checkmin, checkmax
	
	
!-------make first element of interpolated y vals the same as the first 
! element of original values--------------------------------------------	
	if (N == 1) then
		yval = interpy(1)
!-------otherwise, determine whether data are increasing or decreasing--
	else
		!checkmin will be positive if strictly increasing, negative otherwise
		checkorder = interpx(2:N) - interpx(1:N-1) 
		checkmin = minval(checkorder)
		interpchecki = 1
		if (checkmin <= 0) then !decreasing input x vals
			mini = 1      !I guess just initializing these dudes?
			maxi = 1
			
			do interpchecki=1, N
				do tempi = maxi, N-1 !this could really just be 1, right?
					if (checkorder(tempi) > 0) then !this would be the case if the input x vals increased somewhere
						mini = tempi  !mini replaced with the index where x vals start to be increasing
						maxi = N   !maxi bumped up to the last index value
						do tempj = tempi, N-1 !check through all x vals once they start increasing
							if (checkorder(tempj) < 0) then !if it starts decreasing somewhere
								maxi = tempj !make the spot where it starts decreasing again the new maxi
								exit !quit the do loop
							endif
						end do
						exit !quit once we find the end of where things are increasing
					endif !end the case where there's an increasing region
				end do !end the do loop checking if there's a region where things are increasing
				if (xval<=interpx(maxi) .and. xval>=interpx(mini)) then
					exit !check that the value at which we're interpolating falls within the little increasing region we've just set; if so, move on
				endif
			end do
		else !so x vals strictly increasing
			mini = 1
			maxi = N
		endif
		
		if (xval > interpx(maxi) .or. xval < interpx(mini) .or. &
			&interpchecki == N) then
		    if (xval>interpx(maxi)) then
				print *, "FATAL ERROR: x value greater than largest interpolation point"
			elseif (xval < interpx(mini)) then
				print *, "FATAL ERROR: x value less than smallest interpolation point"
			else
				print *, "FATAL ERROR: interpchecki", interpchecki, N
			endif
			print *, "time", myt, "xval", xval, "flagl", flagl, &
				&"flagb", flagb, "flagu", flagu
			print *, "interpx"
			do interpchecki = 1, N
				print *, interpx(interpchecki)
			end do
			print *, "interpy"
			do interpchecki = 1, N
				print *, interpy(interpchecki)
			end do
			print *, "x"
			do interpchecki = 1, XN+1
				print *, sx(interpchecki)
			end do
			print *, "y"
			do interpchecki = 1, XN+1
				print *, sy(interpchecki)
			end do

			print *, interpx(maxi), xval, interpx(mini)
			print*, "FATAL ERROR: interpolation failure, see diagnostics above"
			stop
		endif		
		call spline_pchip_set (maxi-mini+1, interpx(mini:maxi), &
			&interpy(mini:maxi), d)
		call spline_pchip_val (maxi-mini+1, interpx(mini:maxi), &
			&interpy(mini:maxi), d, 1, xval, yval)
	endif
end subroutine interpbridge
