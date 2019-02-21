subroutine trapz1(f,a,b,h1,r)
!===========================================================
! int_trap.f: integration by trapezoid rule of f(x) on [a,b]
!-------------------------------------------------
! f     - Function to integrate (supplied above)
! a	- Lower limit of integration
! b	- Upper limit of integration
! R	- Result of integration (out)
! n	- number of intervals
!=================================================
Real*8 a, b, f, r ,dx, x,h,h1
!double precision a,b,f,r,dx,x
Integer*4 n, i
! interger*4 i


h =  min ((b-a)/4, h1)
!correct to even

n= (b-a)/h + 1

n=2*ceiling(((b-a)/h + 1)/2)


r = 0.d0

dx = (b-a)/(n-1)

Do i=2,n-1
x = a+(i-1)*dx
r = r + f(x)
End do

r = (r + (f(a)+f(b))/2.d0)*dx

Return
end subroutine trapz1
