

function stresstail1DG2 ( xvalG2temp)
	use globalinfo
	implicit none
	
	real (kind=8) :: eta, stresstail1DG2, xvalG2temp
	real (kind=8), external :: G2IntegrandY
	integer (kind=4) :: startingi
	
	xvalG2=xvalG2temp
	startingi = max(flagb-2, 1)
	call interpbridge(XN+2-startingi, sx(startingi:XN+1), sy(startingi:&
	&XN+1), xvalG2, eta)
	!print*, "eta in stresstail1DG2", eta
 	!stresstail1DG2 =-xval*(eta*(R**2-3.0*(eta**2+xval**2))/sqrt(eta**2+xval**2)**3)&
 	!	+xval*(yend*(R**2-3.0*(xval**2+yend**2)))/sqrt(xval**2+yend**2)**3&
 	!	-xval*6.0*log(eta+sqrt(xval**2+eta**2))+xval*6.0*log(yend+sqrt(xval**2+yend**2))
	call trapz1(G2IntegrandY, yend, eta, numtrapz, stresstail1DG2)

	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = eta;
	cinterpx(cinternalcount) = xvalG2;

end
 
function stressIntegrandFlat1DG2 ( xvalG2temp ) !already integrated in y, altered to be done numerically.
	use globalinfo
	implicit none

	real (kind=8) ::  eta, xvalG2temp, stressIntegrandFlat1DG2
	real (kind=8), external :: G2IntegrandY
	integer (kind=4) :: endingi

	xvalG2=xvalG2temp
	endingi = min(flagb+2, XN+1)
	call interpbridge(endingi,sx(1:endingi),sy(1:endingi),xvalG2,eta)
	!	print*, "eta in stressIntegradFlat1DG2", eta
 !	stressIntegrandFlat1D = -xval*(yend*(R**2-3.0*(yend**2+xval**2))/sqrt(yend**2+xval**2)**3)&
 !		+xval*(eta*(R**2-3.0*(xval**2+eta**2)))/sqrt(xval**2+eta**2)**3&
 	!	-6.0*log((yend+sqrt(xval**2+yend**2))**xval)+6.0*log((eta+sqrt(xval**2+eta**2))**xval)
	call trapz1(G2IntegrandY, yend,eta,numtrapz,stressIntegrandFlat1DG2)
	
	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = eta;
	cinterpx(cinternalcount) = xvalG2;
end

function stressIntegrandsphere1DG2(yvalG2temp) !this time already integrated once in x.
	use globalinfo
	implicit none
	
	real (kind=8) eta, yvalG2temp, stressIntegrandsphere1DG2, xmin
	real (kind=8), external :: G2IntegrandX	
	integer (kind=4) tempflagb
	
	yvalG2=yvalG2temp
	tempflagb = min(flagb+2, XN+1)

	call interpbridge(tempflagb,sy(1:tempflagb),sx(1:tempflagb),yvalG2,&
	&eta)
	!print*, "eta in StressIntegrandsphere1DG2", eta
	xmin = Sqrt(R**2-yvalG2**2)
	if(isnan(xmin))then
		xmin = 0.
	endif

	!stressIntegrandsphere1D = 2.0/R*(R**2-yval**2)-eta**2*(-R**2+3.0*(yval**2+eta**2))/(yval**2+eta**2)**(1.5)
!print*, "xmin in stressIntegrandsphere1DG2", xmin
	call trapz1(G2IntegrandX,eta,xmin,numtrapz,stressIntegrandsphere1DG&
	&2)
!print*, "stressIntegrandsphere1DG2", stressIntegrandSphere1DG2
	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = yvalG2;
	cinterpx(cinternalcount) = eta;
end


function stressIntegrand1DG2(yvalG2temp) !Again, first integration is in X
	use globalinfo
	implicit none	
	
	integer (kind=4) tempflagb
	real (kind=8) eta, yvalG2temp, stressIntegrand1DG2
	real (kind=8), external :: G2IntegrandX
	
	yvalG2 = yvalG2temp
	tempflagb = min(flagb+2, XN+1)
	call interpbridge ( tempflagb, sy(1:tempflagb), sx(1:tempflagb), yv&
	&alG2, eta )
	!print*, " eta in stressIntegrand1DG2", eta
!print *, "tempflagb",tempflagb, "yval", yval,"eta",eta
	!stressIntegrand1D =- eta**2*(-R**2+3.0*(yval**2+eta**2))/((yval**2+eta**2)**(1.5));
	call trapz1(G2IntegrandX,  eta, 0.0, numtrapz, stressIntegrand1DG2)
	cinternalcount = cinternalcount +1;
	cinterpy(cinternalcount) = yvalG2;
	cinterpx(cinternalcount) = eta;
end

function G2INtegrandY(yg)
	use globalinfo
	implicit none
	
	real (kind=8), intent(in) ::  yg
	real (kind=8) :: thet, rt, G2IntegrandY
	real (kind=8) :: T2, T4, T6, T8, T10, T12, T14, T16, T18, T20
	real (kind=8) :: LegendreP, GegenbauerC
	real (kind=8) :: taur, taur2, taur4, taur6, taur8, taur10,taur12, &
	&taur14, taur16, taur18,taur20, taut, taut2, taut4, taut6, taut8, &
	&taut10,taut12, taut14, taut16, taut18, taut20

	rt = sqrt(xvalG2**2+yg**2)
	thet = acos(yg/rt)

	T2 = 3.*pi/4.
	T4 = 21.*pi/32.
	T6 = 165.*pi/256.
	T8 = 2625.*pi/4096.
	T10 = 41895.*pi/65536.
	T12 = (334719.*pi)/524288.
	T14 = (2625673.*pi)/(4194304.)
	T16 =(85579065.*pi)/134217728.
	T18 =(2737609875.*pi)/4294967296.
	T20 =(21895664505.*pi)/34359738368.

	taur2 = -LegendreP(1,thet)*T2/2.*(R**(2+1)/rt**(2+1)-R**(2-1)/rt**&
	&(2-1))
	taur4 = -LegendreP(3,thet)*T4/2.*(R**(4+1)/rt**(4+1)-R**(4-1)/rt**&
	&(4-1))
	taur6 = -LegendreP(5,thet)*T6/2.*(R**(6+1)/rt**(6+1)-R**(6-1)/rt**&
	&(6-1))
	taur8 = -LegendreP(7,thet)*T8/2.*(R**(8+1)/rt**(8+1)-R**(8-1)/rt**&
	&(8-1))
	taur10 = -LegendreP(9,thet)*T10/2.*(R**(10+1)/rt**(10+1)-R**(10-1)&
	&/rt**(10-1))
	taur12 = -LegendreP(11,thet)*T12/2.*(R**(12+1)/rt**(12+1)-R**(12-1)&
	&/rt**(12-1))
	taur14 = -LegendreP(13,thet)*T14/2.*(R**(14+1)/rt**(14+1)-R**(14-1)&
	&/rt**(14-1))
	taur16 = -LegendreP(15,thet)*T16/2.*(R**(16+1)/rt**(16+1)-R**(16-1)&
	&/rt**(16-1))
	taur18 = -LegendreP(17,thet)*T18/2.*(R**(18+1)/rt**(18+1)-R**(18-1)&
	&/rt**(18-1))
	taur20 = -LegendreP(19,thet)*T20/2.*(R**(20+1)/rt**(20+1)-R**(20-1)&
	&/rt**(20-1))
	taur=taur2+taur4+taur6+taur8+taur10+taur12+taur14+taur16+taur18+tau&
	&r20

	if(thet==0.0)then
		taut2 = 0.0
		taut4 = 0.0
		taut6 = 0.0
		taut8 = 0.0
		taut10 = 0.0
		taut12 = 0.0
		taut14 = 0.0
		taut16 = 0.0
		taut18 = 0.0
		taut20 = 0.0
	else
		taut2 = GegenbauerC(2,thet)/sin(thet)*(T2/2.)*(-(2-1)*R**(2+1)/&
		&rt**(2+1)+(2-3)*R**(2-1)/rt**(2-1))
		taut4 = GegenbauerC(4,thet)/sin(thet)*(T4/2.)*(-(4-1)*R**(4+1)/&
		&rt**(4+1)+(4-3)*R**(4-1)/rt**(4-1))
		taut6 = GegenbauerC(6,thet)/sin(thet)*(T6/2.)*(-(6-1)*R**(6+1)/&
		&rt**(6+1)+(6-3)*R**(6-1)/rt**(6-1))
		taut8 = GegenbauerC(8,thet)/sin(thet)*(T8/2.)*(-(8-1)*R**(8+1)/&
		&rt**(8+1)+(8-3)*R**(8-1)/rt**(8-1))
		taut10 = GegenbauerC(10,thet)/sin(thet)*(T10/2.)*(-(10-1)*R**(&
		&10+1)/rt**(10+1)+(10-3)*R**(10-1)/rt**(10-1))
		taut12 = GegenbauerC(12,thet)/sin(thet)*(T12/2.)*(-(12-1)*R**(&
		&12+1)/rt**(12+1)+(12-3)*R**(12-1)/rt**(12-1))
		taut14 = GegenbauerC(14,thet)/sin(thet)*(T14/2.)*(-(14-1)*R**(1&
		&4+1)/rt**(14+1)+(14-3)*R**(14-1)/rt**(14-1))
		taut16 = GegenbauerC(16,thet)/sin(thet)*(T16/2.)*(-(16-1)*R**(1&
		&6+1)/rt**(16+1)+(16-3)*R**(16-1)/rt**(16-1))
		taut18 = GegenbauerC(18,thet)/sin(thet)*(T18/2.)*(-(18-1)*R**(1&
		&8+1)/rt**(18+1)+(18-3)*R**(18-1)/rt**(18-1))
		taut20 = GegenbauerC(20,thet)/sin(thet)*(T20/2.)*(-(20-1)*R**(2&
		&0+1)/rt**(20+1)+(20-3)*R**(20-1)/rt**(20-1))
	endif
	taut = taut2+taut4+taut6+taut8+taut10+taut12+taut14+taut16+taut18+&
	&taut20



	G2IntegrandY = -xvalG2*(cos(thet)*taur+sin(thet)*taut)*g!xvalG2**2*g/R!(-3.*R*xvalG2**2/(4.*Sqrt(xvalG2**2+yg**2)**3)-R**3*xvalG2**2/(4.*Sqrt(xvalG2**2+yg**2)**5))*xvalG2
	!print*, "G2Y"
	return
end function G2IntegrandY

function G2IntegrandX(xg)
	use globalinfo
	implicit none
	real (kind=8), intent(in) ::  xg
	real (kind=8) :: G2IntegrandX
	real (kind=8) ::rt,thet, T2, T4, T6, T8, T10, T12, T14, T16, T18,T20
	real (kind=8) :: LegendreP, GegenbauerC
	real (kind=8) :: taur, taur2, taur4, taur6, taur8, taur10,taur12, &
	&taur14, taur16, taur18,taur20, taut, taut2, taut4, taut6, taut8, &
	&taut10,taut12, taut14, taut16, taut18, taut20

	rt = sqrt(xg**2+yvalG2**2)
	thet = acos(yvalG2/rt)
 !print*, "thet", thet
	!print*, "rt", rt
	T2 = 3.*pi/4.
	T4 = 21.*pi/32.
	T6 = 165.*pi/256.
	T8 = 2625.*pi/4096.
	T10 = 41895.*pi/65536.
	T12 = (334719.*pi)/524288.
	T14 = (2625673.*pi)/(4194304.)
	T16 =(85579065.*pi)/134217728.
	T18 =(2737609875.*pi)/4294967296.
	T20 =(21895664505.*pi)/34359738368.

	taur2 = -LegendreP(1,thet)*T2/2.*(R**(2+1)/rt**(2+1)-R**(2-1)/rt**&
	&(2-1))
	taur4 = -LegendreP(3,thet)*T4/2.*(R**(4+1)/rt**(4+1)-R**(4-1)/rt**&
	&(4-1))
	taur6 = -LegendreP(5,thet)*T6/2.*(R**(6+1)/rt**(6+1)-R**(6-1)/rt**&
	&(6-1))
	taur8 = -LegendreP(7,thet)*T8/2.*(R**(8+1)/rt**(8+1)-R**(8-1)/rt**&
	&(8-1))
	taur10 = -LegendreP(9,thet)*T10/2.*(R**(10+1)/rt**(10+1)-R**(10-1)/&
	&rt**(10-1))
	taur12 = -LegendreP(11,thet)*T12/2.*(R**(12+1)/rt**(12+1)-R**(12-1)&
	&/rt**(12-1))
	taur14 = -LegendreP(13,thet)*T14/2.*(R**(14+1)/rt**(14+1)-R**(14-1)&
	&/rt**(14-1))
	taur16 = -LegendreP(15,thet)*T16/2.*(R**(16+1)/rt**(16+1)-R**(16-1)&
	&/rt**(16-1))
	taur18 = -LegendreP(17,thet)*T18/2.*(R**(18+1)/rt**(18+1)-R**(18-1)&
	&/rt**(18-1))
	taur20 = -LegendreP(19,thet)*T20/2.*(R**(20+1)/rt**(20+1)-R**(20-1)&
	&/rt**(20-1))
	taur=taur2+taur4+taur6+taur8+taur10+taur12+taur14+taur16+taur18+&
	&taur20
	!print*, "taurvals", taur2, taur4, taur6, taur8, taur10, taur12, taur14, taur16, taur18, taur20

	if(thet==0.0)then
		taut2 = 0.0
		taut4 = 0.0
		taut6 = 0.0
		taut8 = 0.0
		taut10 = 0.0
		taut12 = 0.0
		taut14 = 0.0
		taut16 = 0.0
		taut18 = 0.0
		taut20 = 0.0
	else
		taut2 = GegenbauerC(2,thet)/sin(thet)*(T2/2.)*(-(2-1)*R**(2+1)/&
		&rt**(2+1)+(2-3)*R**(2-1)/rt**(2-1))
		taut4 = GegenbauerC(4,thet)/sin(thet)*(T4/2.)*(-(4-1)*R**(4+1)/&
		&rt**(4+1)+(4-3)*R**(4-1)/rt**(4-1))
		taut6 = GegenbauerC(6,thet)/sin(thet)*(T6/2.)*(-(6-1)*R**(6+1)/&
		&rt**(6+1)+(6-3)*R**(6-1)/rt**(6-1))
		taut8 = GegenbauerC(8,thet)/sin(thet)*(T8/2.)*(-(8-1)*R**(8+1)/&
		&rt**(8+1)+(8-3)*R**(8-1)/rt**(8-1))
		taut10 = GegenbauerC(10,thet)/sin(thet)*(T10/2.)*(-(10-1)*R**&
		&(10+1)/rt**(10+1)+(10-3)*R**(10-1)/rt**(10-1))
		taut12 = GegenbauerC(12,thet)/sin(thet)*(T12/2.)*(-(12-1)*R**&
		&(12+1)/rt**(12+1)+(12-3)*R**(12-1)/rt**(12-1))
		taut14 = GegenbauerC(14,thet)/sin(thet)*(T14/2.)*(-(14-1)*R**&
		&(14+1)/rt**(14+1)+(14-3)*R**(14-1)/rt**(14-1))
		taut16 = GegenbauerC(16,thet)/sin(thet)*(T16/2.)*(-(16-1)*R**&
		&(16+1)/rt**(16+1)+(16-3)*R**(16-1)/rt**(16-1))
		taut18 = GegenbauerC(18,thet)/sin(thet)*(T18/2.)*(-(18-1)*R**&
		&(18+1)/rt**(18+1)+(18-3)*R**(18-1)/rt**(18-1))
		taut20 = GegenbauerC(20,thet)/sin(thet)*(T20/2.)*(-(20-1)*R**&
		&(20+1)/rt**(20+1)+(20-3)*R**(20-1)/rt**(20-1))
	endif
	taut = taut2+taut4+taut6+taut8+taut10+taut12+taut14+taut16+taut18+&
	&taut20
	!print*, "tautvals", taut2, taut4, taut6, taut8, taut10, taut12, taut14, taut16, taut18, taut20

	G2IntegrandX = -xg*(cos(thet)*taur+sin(thet)*taut)*g !xg**2*g/R!(-3.*R*xg**2/(4.*Sqrt(xg**2+yvalG2**2)**3)-R**3*xg**2/(4.*Sqrt(xg**2+yvalG2**2)**5))*xg
	!print*, "G2X"
	return
end function G2IntegrandX

