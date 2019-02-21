!----------------------------------------------------------------------- 
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Computes the first few Gegenbauer polynomials of degree n. 
!> @brief
!> Compute \f$ C^{-1/2}_n(\dcos(\theta) \f$. 
!
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation, standardized - H. Arrowood
!
!> @param[in] n  
!> @param[in] theta             
!> @return GegenbauerC
!-----------------------------------------------------------------------  
function GegenbauerC(n,theta)
    use globalinfo
	implicit none
	integer, intent(in)		  :: n !>degree of Gegenbauer polynomial
	real(kind=8), intent(in)  :: theta    !>angle
	real(kind=8)			  :: GegenbauerC  !>function value
 
	if(n==2) then
		GegenbauerC=.5*(1.0-dcos(theta)**2)
	 	return
	end if
	if(n==4) then
		GegenbauerC=.125*(1.0-dcos(theta)**2)*(5.0*dcos(theta)**2-1.0)
		return
	end if

	if(n==6) then
		GegenbauerC=.0625*(1.0-dcos(theta)**2)*(21.0*dcos(theta)**4-14.&
			&*dcos(theta)**2+1.0)
		return
	end if

	if(n==8) then
		GegenbauerC=1.0/128.0*(1.0-dcos(theta)**2)*(429.*dcos(theta)**6&
			&-495.0*dcos(theta)**4+135.0*dcos(theta)**2-5.0)
		return
	end if
	if(n==10) then
		GegenbauerC=0.02734375 - (315*Cos(theta)**2)/256. + (1155*Cos(t&
		&heta)**4)/128. - (3003*Cos(theta)**6)/128. + (6435*Cos(theta)*&
		&*8)/256. - (2431*Cos(theta)**10)/256.
		return
	end if
	if(n==12) then
		GegenbauerC=-0.0205078125 + (693*Cos(theta)**2)/512. - (15015*C&
		&os(theta)**4)/1024. + (15015*Cos(theta)**6)/256. - (109395*Cos&
		&(theta)**8)/1024. + (46189*Cos(theta)**10)/512. - (29393*Cos(t&
		&heta)**12)/1024.
		return
	end if
	if(n==14) then
		GegenbauerC=0.01611328125 - (3003*Cos(theta)**2)/2048. +(45045.&
		&*Cos(theta)**4)/2048. - (255255*Cos(theta)**6)/2048. +(692835.&
		&*Cos(theta)**8)/2048. - (969969*Cos(theta)**10)/2048. +& 
        &(676039*Cos(theta)**12)/2048. - (185725*Cos(theta)**14)/2048.
		return
	end if
	if(n==16) then
		GegenbauerC=        -0.013092041015625 + (6435*Cos(theta)**2)/4&
		&096. - (255255*Cos(theta)**4)/8192. + (969969*Cos(theta)**6)/4&
		&096. - (14549535*Cos(theta)**8)/16384. + (7436429*Cos(theta)**&
		&10)/4096. - (16900975*Cos(theta)**12)/8192. + (5014575*Cos(the&
		&ta)**14)/4096. - (9694845*Cos(theta)**16)/32768.
		return
	end if
	if(n==18) then
		GegenbauerC=0.0109100341796875 - (109395*Cos(theta)**2)/65536. &
		&+ (692835*Cos(theta)**4)/16384. - (6789783*Cos(theta)**6)/1638&
		&4. + (66927861*Cos(theta)**8)/32768. - (185910725*Cos(theta)**&
		&10)/32768. + (152108775*Cos(theta)**12)/16384. - (145422675*Co&
		&s(theta)**14)/16384. + (300540195*Cos(theta)**16)/65536. - (64&
		&822395*Cos(theta)**18)/65536.
		return
	end if
	if(n==20) then
		GegenbauerC=  -0.009273529052734375 + (230945*Cos(theta)**2)/13&
		&1072. - (14549535*Cos(theta)**4)/262144. + (22309287*Cos(theta&
		)**6)/32768. - (557732175*Cos(theta)**8)/131072. + (1003917915*&
		&Cos(theta)**10)/65536. - (4411154475*Cos(theta)**12)/131072. +&
		& (1502700975*Cos(theta)**14)/32768. - (9917826435*Cos(theta)**&
		&16)/262144. + (2268783825*Cos(theta)**18)/131072. - &
        &(883631595*Cos(theta)**20)/262144.
		return
	end if


end function GegenbauerC

!----------------------------------------------------------------------- 
!> @author 
!> H. Arrowood, UNC-CH
!
! DESCRIPTION: 
!> Computes the first few Legendre polynomials of degree n, evaluated
!!at dcos(theta). 
!> @brief
!> Compute \f$ P_n(\dcos(\theta) \f$. 
!
! REVISION HISTORY:
! 07_Nov_2017 - Added documentation, standardized - H. Arrowood
!
!> @param[in] n  
!> @param[in] theta             
!> @return LegendreP
!-----------------------------------------------------------------------

function LegendreP(n,theta)
    implicit none

	integer, intent(in)		  :: n         !>degree of polynomial
	real(kind=8), intent(in)  :: theta     !>angle
	real(kind=8)			  :: LegendreP !>value of polynomial

	if(n==1) then
		LegendreP=dcos(theta)
	 	return
	end if
	if(n==3) then
		LegendreP=.5*(-3.0*dcos(theta)+5.0*dcos(theta)**3)
		return
	end if

	if(n==5) then
		LegendreP=.125*(15.0*dcos(theta)-70.0*dcos(theta)**3+60.0* &
			&dcos(theta)**5)
		return
	end if

	if(n==7) then
		LegendreP=.0625*(-35.0*dcos(theta)+315.0*dcos(theta)**3-693.0* &
			&dcos(theta)**5+429.0*dcos(theta)**7)
		return
	end if
	if(n==9) then
		LegendreP=(315*Cos(theta) - 4620*Cos(theta)**3 + 18018*Cos(thet&
		&a)**5 - 25740*Cos(theta)**7 + 12155*Cos(theta)**9)/128.
		return
	end if
	if(n==11) then
		LegendreP=(-693*Cos(theta) + 15015*Cos(theta)**3 - 90090*Cos(th&
		&eta)**5 + 218790*Cos(theta)**7 - 230945*Cos(theta)**9 + 88179*&
		&Cos(theta)**11)/256.
		return
	end if
	if(n==13) then
		LegendreP=(3003*Cos(theta) - 90090*Cos(theta)**3 + 765765*Cos(t&
		&heta)**5 - 2771340*Cos(theta)**7 + 4849845*Cos(theta)**9 - 405&
		&6234*Cos(theta)**11 + 1300075*Cos(theta)**13)/1024.
		return
	end if
	if(n==15) then
		LegendreP=        (-6435*Cos(theta) + 255255*Cos(theta)**3 - 29&
		&09907*Cos(theta)**5 + 14549535*Cos(theta)**7 - 37182145*Cos(th&
		&eta)**9 + 50702925*Cos(theta)**11 - 35102025*Cos(theta)**13 + &
        &9694845*Cos(theta)**15)/2048.
		return
	end if
	if(n==17) then
		LegendreP=(109395*Cos(theta) - 5542680*Cos(theta)**3 + 81477396&
		&*Cos(theta)**5 - 535422888*Cos(theta)**7 + 1859107250*Cos(thet&
		&a)**9 - 3650610600*Cos(theta)**11+4071834900*Cos(theta)**13 -& 
        &2404321560*Cos(theta)**15 + 583401555*Cos(theta)**17)/32768.
		return
	end if
	if(n==19) then
		LegendreP=  (-230945*Cos(theta) + 14549535*Cos(theta)**3 - 2677&
		&11444*Cos(theta)**5 + 2230928700*Cos(theta)**7 - 10039179150*C&
		&os(theta)**9 + 26466926850*Cos(theta)**11 - 42075627300*Cos(th&
		&eta)**13 + 39671305740*Cos(theta)**15 - 20419054425*Cos(theta)&
		&**17 + 4418157975*Cos(theta)**19)/65536.
		return
	end if

end function LegendreP


