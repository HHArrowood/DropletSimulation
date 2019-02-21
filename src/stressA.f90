
!================= force at leading order ====================

function stresstail1DA ( xval )
	use globalinfo
	implicit none
	
	real (kind=8) eta, xval, stresstail1DA
	integer (kind=4) startingi
	
	startingi = max(flagb-2, 1)
	call interpbridge(XN+2-startingi, sx(startingi:XN+1), sy(startingi:&
	&XN+1), xval, eta)
	
 	stresstail1DA =xval*2*pi*(eta-yend)
end
 
function stressIntegrandFlat1DA ( xval )
	use globalinfo
	implicit none

	real (kind=8) xval, eta, stressIntegrandFlat1DA
	integer (kind=4) endingi

	endingi = min(flagb+2, XN+1)
	call interpbridge( endingi, sx(1:endingi), sy(1:endingi), xval, eta)
 	stressIntegrandFlat1DA = xval*2*pi*(yend-eta)
end

function stressIntegrandsphere1DA(yval)
	use globalinfo
	implicit none
	
	real (kind=8) yval, eta, stressIntegrandsphere1DA
	integer (kind=4) tempflagb

	tempflagb = min(flagb+2, XN+1)

	call interpbridge( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yva&
	&l, eta )
	stressIntegrandsphere1DA = pi*(  eta**2-(R**2-yval**2))
end


function stressIntegrand1DA(yval)
	use globalinfo
	implicit none	
	
	integer (kind=4) tempflagb
	real (kind=8) yval, eta, stressIntegrand1DA
	
	tempflagb = min(flagb+2, XN+1)
	call interpbridge ( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yv&
	&al, eta )
!print *, "tempflagb",tempflagb, "yval", yval,"eta",eta
	stressIntegrand1DA =pi*(eta)**2;
end


!================= force order epsilon ====================

function stresstail1DE ( xval )
	use globalinfo
	implicit none

	real (kind=8) eta, xval, stresstail1DE
	integer (kind=4) startingi

	startingi = max(flagb-2, 1)
	call interpbridge(XN+2-startingi, sx(startingi:XN+1), sy(startingi:&
	&XN+1), xval, eta)

	stresstail1DE =2*pi*(xval*(eta - yend)*(2*R**3 - 3*R**2*(-R - sy(1)&
	&) + (eta**2 + eta*yend+ yend**2)*(-R - sy(1))))/(2.*R**3)
end

function stressIntegrandFlat1DE ( xval )
	use globalinfo
	implicit none

	real (kind=8) xval, eta, stressIntegrandFlat1DE
	integer (kind=4) endingi

	endingi = min(flagb+2, XN+1)
	call interpbridge( endingi, sx(1:endingi), sy(1:endingi), xval, eta)

	stressIntegrandFlat1DE = - 2*pi* (xval*(eta - yend)*(2*R**3 &
	- 3*R**2*(-R - sy(1)) + (eta**2 + eta*yend + yend**2)*(-R - sy(1)))&
	&)/(2.*R**3)
end

function stressIntegrandsphere1DE(yval)
	use globalinfo
	implicit none

	real (kind=8) yval, eta, stressIntegrandsphere1DE
	integer (kind=4) tempflagb

	tempflagb = min(flagb+2, XN+1)

	call interpbridge( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yva&
	&l, eta )

	stressIntegrandsphere1DE = - 2*pi*((-eta**2 + R**2 - yval**2)*(2*R*&
	&*3 - 3*R**2*(-R - sy(1)) + 3*yval**2*(-R - sy(1))))/(4.*R**3)
end


function stressIntegrand1DE(yval)
	use globalinfo
	implicit none

	integer (kind=4) tempflagb
	real (kind=8) yval, eta, stressIntegrand1DE

	tempflagb = min(flagb+2, XN+1)
	call interpbridge ( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yv&
	&al, eta )

	!print *, "tempflagb",tempflagb, "yval", yval,"eta",eta
	stressIntegrand1DE = (eta**2*pi*(2*R**3 - 3*R**2*(-R - sy(1)) + 3*y&
	&val**2*(-R - sy(1))))/(2.*R**3);
end


