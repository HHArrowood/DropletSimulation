

function stresstail1D ( xval )  !has been analytically integrated in y already. I'll be doing this numerically.
	use globalinfo
	implicit none
	
	real (kind=8) eta, xval, stresstail1D
	integer (kind=4) startingi
	
	startingi = max(flagb-2, 1)
	call interpbridge(XN+2-startingi, sx(startingi:XN+1), sy(startingi:&
	&XN+1), xval, eta)
	
 	stresstail1D =-xval*(eta*(R**2-3.0*(eta**2+xval**2))/sqrt(eta**2+xv&
 	&al**2)**3)+xval*(yend*(R**2-3.0*(xval**2+yend**2)))/sqrt(xval**2+y&
 	&end**2)**3-xval*6.0*log(eta+sqrt(xval**2+eta**2))+xval*6.0*log(yen&
 	&d+sqrt(xval**2+yend**2))

cinternalcount = cinternalcount +1;
cinterpy(cinternalcount) = eta;
cinterpx(cinternalcount) = xval;

end
 
function stressIntegrandFlat1D ( xval )
	use globalinfo
	implicit none

	real (kind=8) xval, eta, stressIntegrandFlat1D
	integer (kind=4) endingi

	endingi = min(flagb+2, XN+1)
	call interpbridge( endingi, sx(1:endingi), sy(1:endingi), xval, eta)
 	stressIntegrandFlat1D = -xval*(yend*(R**2-3.0*(yend**2+xval**2))/sq&
 	&rt(yend**2+xval**2)**3)+xval*(eta*(R**2-3.0*(xval**2+eta**2)))/sqr&
 	&t(xval**2+eta**2)**3-6.0*log((yend+sqrt(xval**2+yend**2))**xval)+6&
 	&.0*log((eta+sqrt(xval**2+eta**2))**xval)

cinternalcount = cinternalcount +1;
cinterpy(cinternalcount) = eta;
cinterpx(cinternalcount) = xval;
end

function stressIntegrandsphere1D(yval)
	use globalinfo
	implicit none
	
	real (kind=8) yval, eta, stressIntegrandsphere1D	
	integer (kind=4) tempflagb

	tempflagb = min(flagb+2, XN+1)

	call interpbridge( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yva&
	&l, eta )
	stressIntegrandsphere1D = 2.0/R*(R**2-yval**2)-eta**2*(-R**2+3.0*(y&
	&val**2+eta**2))/(yval**2+eta**2)**(1.5)

cinternalcount = cinternalcount +1;
cinterpy(cinternalcount) = yval;
cinterpx(cinternalcount) = eta;
end


function stressIntegrand1D(yval)
	use globalinfo
	implicit none	
	
	integer (kind=4) tempflagb
	real (kind=8) yval, eta, stressIntegrand1D
	
	tempflagb = min(flagb+2, XN+1)
	call interpbridge ( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yv&
	&al, eta )
!print *, "tempflagb",tempflagb, "yval", yval,"eta",eta
	stressIntegrand1D =- eta**2*(-R**2+3.0*(yval**2+eta**2))/((yval**2+&
	&eta**2)**(1.5));
cinternalcount = cinternalcount +1;
cinterpy(cinternalcount) = yval;
cinterpx(cinternalcount) = eta;
end
