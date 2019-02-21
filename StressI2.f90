
!outputs interpolated interfaces

function stresstail1DI2 ( xvalI2temp)
	use globalinfo
	implicit none
	
	real (kind=8) :: eta, stresstail1DI2, xvalI2temp
	real (kind=8), external :: I2IntegrandY
	integer (kind=4) :: startingi
	
	xvalI2=xvalI2temp
	startingi = max(flagb-2, 1)
	call interpbridge(XN+2-startingi, sx(startingi:XN+1), sy(startingi:XN+1), xvalI2, eta)
	
 	!stresstail1DI2 =-xval*(eta*(R**2-3.0*(eta**2+xval**2))/sqrt(eta**2+xval**2)**3)&
 	!	+xval*(yend*(R**2-3.0*(xval**2+yend**2)))/sqrt(xval**2+yend**2)**3&
 	!	-xval*6.0*log(eta+sqrt(xval**2+eta**2))+xval*6.0*log(yend+sqrt(xval**2+yend**2))
	call trapz1(I2IntegrandY, yend, eta, numtrapz, stresstail1DI2)

	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = eta;
	cinterpx(cinternalcount) = xvalI2;

end
 
function stressIntegrandFlat1DI2 ( xvalI2temp ) !already integrated in y, altered to be done numerically.
	use globalinfo
	implicit none

	real (kind=8) ::  eta, xvalI2temp, stressIntegrandFlat1DI2
	real (kind=8), external :: I2IntegrandY
	integer (kind=4) :: endingi

	xvalI2=xvalI2temp
	endingi = min(flagb+2, XN+1)
	call interpbridge( endingi, sx(1:endingi), sy(1:endingi), xvalI2, eta)

 !	stressIntegrandFlat1D = -xval*(yend*(R**2-3.0*(yend**2+xval**2))/sqrt(yend**2+xval**2)**3)&
 !		+xval*(eta*(R**2-3.0*(xval**2+eta**2)))/sqrt(xval**2+eta**2)**3&
 	!	-6.0*log((yend+sqrt(xval**2+yend**2))**xval)+6.0*log((eta+sqrt(xval**2+eta**2))**xval)
	call trapz1(I2IntegrandY, yend, eta, numtrapz, stressIntegrandFlat1DI2)
	
	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = eta;
	cinterpx(cinternalcount) = xvalI2;
end

function stressIntegrandsphere1DI2(yvalI2temp) !this time already integrated once in x.
	use globalinfo
	implicit none
	
	real (kind=8) eta, yvalI2temp, stressIntegrandsphere1DI2, xmin
	real (kind=8), external :: I2IntegrandX	
	integer (kind=4) tempflagb
	
	yvalI2=yvalI2temp
	tempflagb = min(flagb+2, XN+1)

	call interpbridge( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yvalI2, eta )
	xmin = Sqrt(R**2-yvalI2**2)
	if(isnan(xmin))then
		xmin = 0
	endif

	!stressIntegrandsphere1D = 2.0/R*(R**2-yval**2)-eta**2*(-R**2+3.0*(yval**2+eta**2))/(yval**2+eta**2)**(1.5)
	call trapz1(I2IntegrandX, xmin, eta, numtrapz, stressIntegrandsphere1DI2)

	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = yvalI2;
	cinterpx(cinternalcount) = eta;
end


function stressIntegrand1DI2(yvalI2temp) !Again, first integration is in X
	use globalinfo
	implicit none	
	
	integer (kind=4) tempflagb
	real (kind=8) eta, yvalI2temp, stressIntegrand1DI2
	real (kind=8), external :: I2IntegrandX
	
	yvalI2 = yvalI2temp
	tempflagb = min(flagb+2, XN+1)
	call interpbridge ( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yvalI2, eta )
!print *, "tempflagb",tempflagb, "yval", yval,"eta",eta
	!stressIntegrand1D =- eta**2*(-R**2+3.0*(yval**2+eta**2))/((yval**2+eta**2)**(1.5));
	call trapz1(I2IntegrandX, 0.0, eta, numtrapz, stressIntegrand1DI2)
	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = yvalI2;
	cinterpx(cinternalcount) = eta;
end

function I2INtegrandY(yg)
	use globalinfo
	implicit none
	real (kind=8), intent(in) ::  yi
	real (kind=8) :: I2IntegrandY
	
	I2IntegrandY = (-3./(4.*R**2*pi*mu))*(-3.*R*yi**2/(2.*Sqrt(xvalI2**2+yi**2)**3)+R**3*yi**2/(2.*Sqrt(xvalI2**2+yi**2)**5))*xvalI2
	return
end function I2IntegrandY

function I2IntegrandX(xg)
	use globalinfo
	implicit none
	real (kind=8), intent(in) ::  xi
	real (kind=8) :: I2IntegrandX
	
	I2IntegrandX =  (-3./(4.*R**2*pi*mu))*(-3.*R*yvalI2**2/(2.*Sqrt(xi**2+yvalI2**2)**3)+R**3*yvalI2**2/(2.*Sqrt(xi**2+yvalI2**2)**5))*xi
	return
end function I2IntegrandX
