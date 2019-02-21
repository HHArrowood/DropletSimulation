
!-----------------------------------------------------------------------  
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Computes perturbation velocity using trapezoidal integration
!> @brief
!> A function to evaluate the analytic r-derivative of the theta 
!> component of the Oseen integrand in order to compute the stresses at 
!> each point on the drop surface with greater accuracy than allowed by
!> numerical differentiation. Will be integrated over 3d region, just 
!> like the log term of W. 
!
! REVISION HISTORY:
! 03_Feb_2018 - Wrote and fully verified expressions - H. Arrowood
!    
!@fixme Make this look nice, finish documentation
!-----------------------------------------------------------------------

function dWthetadr(R, theta, y1, y2, y3)
real (kind=16), intent(in) :: r, theta, y1, y2, y3
real (kind=16) :: dWthetadr
real (kind=16) :: ry,WU1, WU2, phicoeffderiv, dphidx1term1, dphidx1term2,&
	&dphidx1term3, dphidx1term4, dphidx1term5, WV1, WV2, dphidx3term1,&
	&dphidx3term2, dphidx3term3, dphidx3term4, dphidx3term5, WU, WV
!First, we need the wu and wv terms, differentiated w.r.t. r. I broke 
!these each up into several terms to hopefully streamline the debugging 
!a little.
!x-values on drop surface will be input as R (drop radius) and a vector
!of theta values, so I can use these to greatly reduce the number of 
!of operations. 


ry = sqrt(y1**2+y2**2+y3**2)
print*, y1**2+y2**2+y3**2
print*, "ry", ry
print*, "sqrt(75)", sqrt(75.)

WU1 =  -((R**4*(-(R*y3) + (y1**2 + y2**2 + y3**2)*Cos(theta))*Sin(theta))/&
     & ((y1**2 + y2**2 + y3**2)**2.5*((R**2*&
     & (R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/(y1**2 + y2**2 + y3**2))**&
     &    1.5)) - (R**4*Cos(theta)*(-(R*y1) + (y1**2 + y2**2 + y3**2)*Sin(theta)))/&
     &   ((y1**2 + y2**2 + y3**2)**2.5*((R**2*&
     & (R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/(y1**2 + y2**2 + y3**2))**1.5&
     & ) + (3*(-(R*y3) + (y1**2 + y2**2 + y3**2)*Cos(theta))*&
     & Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &   (y1**2 + y2**2 + y3**2))*(y1**2 + y2**2 + y3**2 - R*y3*Cos(theta) - R*y1*Sin(theta))*&
     & (-(R*y1) + (y1**2 + y2**2 + y3**2)*Sin(theta)))/&
     &   ((y1**2 + y2**2 + y3**2)**1.5*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta))**3) - &
     &  (3*(y3 - R*Cos(theta))*(-y1 + R*Sin(theta))*(-R + y3*Cos(theta) + y1*Sin(theta)))/&
     &   (y2**2 + (y3 - R*Cos(theta))**2 + (y1 - R*Sin(theta))**2)**2.5 + &
     &  ((-y3 + R*Cos(theta))*Sin(theta))/(y2**2 + (y3 - R*Cos(theta))**2 + (y1 - R*Sin(theta))**2)**1.5 + &
     &  (Cos(theta)*(-y1 + R*Sin(theta)))/(y2**2 + (y3 - R*Cos(theta))**2 + (y1 - R*Sin(theta))**2)**1.5

     
WU2 =         ((-R**2 + y1**2 + y2**2 + y3**2)*(R*Cos(theta)*&
     &   (-(y1*(2*y1**4 + 2*y2**4 + 3*y2**2*y3**2 + y3**4 - R**2*(y1**2 + y2**2 - 2*y3**2) +& 
     &          y1**2*(4*y2**2 + 3*y3**2))) + R*(y1**4 + y1**2*y2**2 + y2**2*y3**2 + y3**4)*Sin(theta)) + &
     &  y3*(-(R*(y1**4 + R**2*(2*y1**2 - y2**2 - y3**2) + 3*y1**2*(y2**2 + y3**2) + &
     &          2*(y2**2 + y3**2)**2)*Sin(theta)) +& 
     &     y1*((y1**2 + y2**2 + y3**2)*(2*R**2 + y1**2 + y2**2 + y3**2) + R**2*y1*y3*Sin(2*theta)))))/&
     &  ((y1**2 + y2**2 + y3**2)**2.5*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta))**2*&
     &Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/(y1**2 + y2**2 + y3**2))&
     &)
	
!Note that the coefficient of the dphi term is zero on the drop surface
!So there is no need to take a derivative of the phi term, just the coeff
!when applying the product rule to this term.     
phicoeffderiv = -((R*(-R**2+ry**2))/ry**3)
	
	
dphidx1term1 =    (4*R**3*y1*y3 + 7*R*y1**3*y3 + 7*R*y1*y2**2*y3 + 7*R*y1*y3**3 -& 
     &6*R*y1*y3*(y1**2 + y2**2 + y3**2)*Cos(2*theta) - 14*R**2*y1**2*y3*Sin(theta) - 3*y1**4*y3*Sin(theta) - &
     &6*R**2*y2**2*y3*Sin(theta) - 6*y1**2*y2**2*y3*Sin(theta) - 3*y2**4*y3*Sin(theta) - 6*R**2*y3**3*Sin(theta) - &
     &6*y1**2*y3**3*Sin(theta) - 6*y2**2*y3**3*Sin(theta) - 3*y3**5*Sin(theta) + &
     &R*Cos(theta)*(R*y1*(3*y1**2 + 3*y2**2 - 5*y3**2) - 3*(y1**4 + y2**4 - 3*y3**4)*Sin(theta)) - &
     &3*R*y1**2*y2**2*Sin(2*theta) + 3*R*y1**2*y3**2*Sin(2*theta) + 3*R*y2**2*y3**2*Sin(2*theta))/&
     &(R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta))**2*&
     &Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/(y1**2 + y2**2 + y3**2))&
     &)
   
   
dphidx1term2 =    (3*R**2*(-(R*y1) + (y1**2 + y2**2 + y3**2)*Sin(theta))*&
     &(R*(y1**2 + y2**2 - y3**2)*Cos(theta) + y3*(y1**2 + y2**2 + y3**2 - 2*R*y1*Sin(theta))))/&
     &  ((y1**2 + y2**2 + y3**2)*((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &   (y1**2 + y2**2 + y3**2))**1.5*(-R**4 + R**3*y3*Cos(theta) + R**3*y1*Sin(theta) +& 
     &  y1**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &   Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &     (y1**2 + y2**2 + y3**2)) + y2**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &   Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &     (y1**2 + y2**2 + y3**2)) + y3**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &   Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &     (y1**2 + y2**2 + y3**2))))

	
dphidx1term3 =   (-3*(y1**2 + y2**2 + y3**2)*((R**4*Sin(theta))/Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + &
     &  R*y1*(-Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + &
     &     Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &       (y1**2 + y2**2 + y3**2))))*&
     &(R**3*y3*(2*R**2 + y1**2 + y2**2 + y3**2 - 2*R*y1*Sin(theta) - &
     &     (2*R**2*Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &          (y1**2 + y2**2 + y3**2)))/Sqrt(R**4/(y1**2 + y2**2 + y3**2))) + &
     &  Cos(theta)*(-(R**4*(y1**2 + y2**2 + 3*y3**2)) + &
     &     Sqrt(R**4/(y1**2 + y2**2 + y3**2))*(y1**2 + y2**2 + y3**2)**2*&
     &      Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &        (y1**2 + y2**2 + y3**2)))))/&
     &  (R**3*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta))*&
     &(-R**4 + R**3*y3*Cos(theta) + R**3*y1*Sin(theta) + &
     &   y1**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &    Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &      (y1**2 + y2**2 + y3**2)) + y2**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &    Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &      (y1**2 + y2**2 + y3**2)) + y3**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &    Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &      (y1**2 + y2**2 + y3**2)))**2)
	
           
dphidx1term4 =    (-3*y3*(y1**2 + y2**2 + y3**2)*Sin(theta))/&
     &  (R**2*(R**4/Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + R*Sqrt(R**2)*y3*Cos(theta) + R*Sqrt(R**2)*y1*Sin(theta)))

 
dphidx1term5 = (3*(y1**2 + y2**2 + y3**2)*(R*Sqrt(R**2)*y3 + (R**4*Cos(theta))/Sqrt(R**4/(y1**2 + y2**2 + y3**2)))*&
     &(R*Sqrt(R**2)*y1 + (R**4*Sin(theta))/Sqrt(R**4/(y1**2 + y2**2 + y3**2))))/&
     & (R**7*Sqrt(R**2)*((R*Sqrt(R**2))/Sqrt(R**4/(y1**2 + y2**2 + y3**2)) + y3*Cos(theta) + y1*Sin(theta))**2)


WV1 =         (-3*(-y3 + R*Cos(theta))**2*(2*Cos(theta)*(-y3 + R*Cos(theta)) + 2*Sin(theta)*(-y1 + R*Sin(theta))))/&
     & (2.*(y2**2 + (-y3 + R*Cos(theta))**2 + (-y1 + R*Sin(theta))**2)**2.5) +& 
     &(2*Cos(theta)*(-y3 + R*Cos(theta)))/(y2**2 + (-y3 + R*Cos(theta))**2 + (-y1 + R*Sin(theta))**2)**1.5 - &
     &(2*Cos(theta)*(-y3 + R*Cos(theta)) + 2*Sin(theta)*(-y1 + R*Sin(theta)))/&
     & (2.*(y2**2 + (-y3 + R*Cos(theta))**2 + (-y1 + R*Sin(theta))**2)**1.5) +& 
     &(3*R**3*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta))**2*&
     &   (2*Cos(theta)*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta)) + &
     &     2*Sin(theta)*(-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(theta))))/&
     & (2.*(y1**2 + y2**2 + y3**2)**1.5*((R**4*y2**2)/(y1**2 + y2**2 + y3**2)**2 + &
     &      (-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta))**2 + &
     &      (-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(theta))**2)**2.5) -& 
     &(2*R**3*Cos(theta)*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta)))/&
     & ((y1**2 + y2**2 + y3**2)**1.5*((R**4*y2**2)/(y1**2 + y2**2 + y3**2)**2 +& 
     &      (-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta))**2 + &
     &      (-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(theta))**2)**1.5) +& 
     &(R*(2*Cos(theta)*(-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta)) +& 
     &     2*Sin(theta)*(-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(theta))))/&
     & (2.*Sqrt(y1**2 + y2**2 + y3**2)*((R**4*y2**2)/(y1**2 + y2**2 + y3**2)**2 +& 
     &      (-((R**2*y3)/(y1**2 + y2**2 + y3**2)) + R*Cos(theta))**2 + &
     &      (-((R**2*y1)/(y1**2 + y2**2 + y3**2)) + R*Sin(theta))**2)**1.5)
	
     
 WV2 = (y3*(-R**2+ry**2)*(2.*R**2*y1**2*y3+y1**4*y3+2.*R**2*y2**2*y3+&
	&2.*y1**2*y2**2*y3+y2**4*y3+2.*R**2*y3**3+2.*y1**2*y3**3+2.*y2**2*&
	&y3**3+y3**5+R*(-4.*y1**4-4.*y2**4-7.*y2**2*y3**2-3.*y3**4+R**2*&
	&(2.*y1**2+2.*y2**2-y3**2)-y1**2*(8.*y2**2+7.*y3**2))*cos(theta)+&
	&R**2*y3*ry**2*cos(2.*theta)-3.*R**3*y1*y3*sin(theta)+R*y1**3*y3*&
	&sin(theta)+R*y1*y2**2*y3*sin(theta)+R*y1*y3**3*sin(theta)+&
	&R**2*y1**3*sin(2.*theta)+R**2*y1*y2**2*sin(2.*theta)+R**2*y1*y3**2*&
	&sin(2.*theta)))/(ry**5*(R**2+ry**2-2.*R*y3*cos(theta)-2.*R*y1*&
	&sin(theta))**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(theta)-&
	&2.*R*y1*sin(theta)))/ry**2))

      
 dphidx3term1 = (-y3*(3.*ry**4+R**2*(5.*y1**2+5.*y2**2+13.*y3**2))*&
	&cos(theta)+1./2.*R*(2.*R**2*y1**2-y1**4+2.*R**2*y2**2-2.*y1**2*&
	&y2**2-y2**4+10.*R**2*y3**2+12.*y1**2*y3**2+12.*y2**2*y3**2+13.*&
	&y3**4-3.*(y1**4+y2**4-2.*y2**2*y3**2-3.*y3**4+2.*y1**2*(y2**2-&
	&y3**2))*cos(2.*theta)-4.*R*y1*(y1**2+y2**2+5.*y3**2)*sin(theta)+&
	&12.*y1**3*y3*sin(2.*theta)+12.*y1*y2**2*y3*sin(2.*theta)+&
	&12.*y1*y3**3*sin(2.*theta)))/(R**2*(R**2+ry**2-2.*R*y3*cos(theta)-&
	&2.*R*y1*sin(theta))**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(theta)-&
	&2.*R*y1*sin(theta)))/ry**2))
	
	
 dphidx3term2 =  (3*(y1**2 + y2**2 + y3**2)**2*((R**6*(R*y3 - (y1**2 + y2**2 + y3**2)*Cos(theta))**2)/&
     &  (y1**2 + y2**2 + y3**2)**3 - (R**5*y3*(R*y3 - (y1**2 + y2**2 + y3**2)*Cos(theta))*&
     &    (R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/(y1**2 + y2**2 + y3**2)**3 - &
     & ((R**4/(y1**2 + y2**2 + y3**2))**1.5*&
     &    (R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta))*&
     &    (Sqrt(R**4/(y1**2 + y2**2 + y3**2)) -& 
     &      Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &        (y1**2 + y2**2 + y3**2))))/R**2))/&
     &  (R**3*((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &  (y1**2 + y2**2 + y3**2))**1.5*(-R**4 + R**3*y3*Cos(theta) + R**3*y1*Sin(theta) + &
     & y1**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &  Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &    (y1**2 + y2**2 + y3**2)) + y2**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &  Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &    (y1**2 + y2**2 + y3**2)) + y3**2*Sqrt(R**4/(y1**2 + y2**2 + y3**2))*&
     &  Sqrt((R**2*(R**2 + y1**2 + y2**2 + y3**2 - 2*R*y3*Cos(theta) - 2*R*y1*Sin(theta)))/&
     &    (y1**2 + y2**2 + y3**2))))
	
	
 dphidx3term3 = -((3.*ry**2*((R**4*cos(theta))/(R**2/ry)+R*y3*&
	&(-(R**2/ry)+sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(theta)-&
	&2.*R*y1*sin(theta)))/(ry**2))))*(R**3*y3*(2.*R**2+ry**2-&
	&2.*R*y1*sin(theta)-(2.*R**2*sqrt((R**2*(R**2+ry**2-2.*R*y3*&
	&cos(theta)-2.*R*y1*sin(theta)))/ry**2))/(R**2/ry))+cos(theta)*&
	&(-R**4*(y1**2+y2**2+3.*y3**2)+(R**2/ry)*ry**4*sqrt((R**2*(R**2+&
	&ry**2-2.*R*y3*cos(theta)-2.*R*y1*sin(theta)))/(ry**2)))))/(R**3*&
	&(R**2+ry**2-2.*R*y3*cos(theta)-2.*R*y1*sin(theta))*(-R**4+&
	&R**3*y3*cos(theta)+R**3*y1*sin(theta)+y1**2*(R**2/ry)*sqrt((R**2*&
	&(R**2+ry**2-2.*R*y3*cos(theta)-2.*R*y1*sin(theta)))/ry**2)+&
	&y2**2*(R**2/ry)*sqrt((R**2*(R**2+ry**2-2.*R*y3*cos(theta)-2.*R*&
	&y1*sin(theta)))/ry**2)+y3**2*(R**2/ry)*sqrt((R**2*(R**2+ry**2-2.*&
	&R*y3*cos(theta)-2.*R*y1*sin(theta)))/ry**2))**2))
	
	
 dphidx3term4 = -((3.*ry**2*(ry+y3*cos(theta)))/(R**2*(R**2*ry+&
	&R**2*y3*cos(theta)+R**2*y1*sin(theta))))
	
    
 dphidx3term5 = (3.*ry**2*(R**2*y3+(R**4*cos(theta))/(R**2/ry))**2)/&
	&(R**7*R*(ry+y3*cos(theta)+y1*sin(theta))**2)
	
	
	WU = WU1+WU2+phicoeffderiv*(dphidx1term1+dphidx1term2+&
	&dphidx1term3+dphidx1term4+dphidx1term5)
	WV = WV1+WV2+phicoeffderiv*(dphidx3term1+dphidx3term2+&
	&dphidx3term3+dphidx3term4+dphidx3term5)
	
	print*, "WU", WU
	print*, "WV", WV!verified
	
	dWthetadr = cos(theta)*WU - sin(theta)*WV
	return 
end function dWthetadr
