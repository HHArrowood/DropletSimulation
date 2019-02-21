!***********************************************************
!* This program tests the subroutine BESSK to calculate the*
!* modified Bessel function of the third kind of order N   *
!* for any real positive argument X.                       *
!* ------------------------------------------------------- *
!* SAMPLE RUN:                                             *
!*                                                         *
!*  X =  1.20970000000000                                  *
!*                                                         *
!*  N =          0                                         *
!*  Y = 0.314324491956902                                  *
!*                                                         *
!*  N =          1                                         *
!*  Y = 0.428050751380123                                  *
!*                                                         *
!*  N =          2                                         *
!*  Y =  1.02202185722122                                  *
!*                                                         *
!*  N =          3                                         *
!*  Y =  3.80747327670449                                  *
!*                                                         *
!*  N =          4                                         *
!*  Y =  19.9067367949966                                  *
!*                                                         *
!* ------------------------------------------------------- *
!* Reference:   From Numath Library By Tuan Dang Trong     *
!* in Fortran 77 [BIBLI 18].                               *
!*                                                         *
!*                      F90 Version By J-P Moreau, Paris.  *
!***********************************************************     
   !PROGRAM TBESSK
!
!   REAL*8 X, Y
!
!   N=5           !number of integer orders
!   X=1.2097d0    !real argument
!
!   print *,' '
!   print *,' X =', X
!   print *,' '
!
!   do i=0, N-1
!     Y=BESSK(i,X)
!     print *,' N =', i
!     print *,' Y =', Y
!     print *,' '
!   end do
!
!   stop
!   END

      FUNCTION BESSK(N,X)
      REAL *8 X,BESSK,BESSK0,BESSK1,TOX,BK,BKM,BKP
! ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
! ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
! ------------------------------------------------------------------------
      IF (N.EQ.0) THEN
      BESSK = BESSK0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSK = BESSK1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.D0) THEN
      BESSK = 1.D30
      RETURN
      ENDIF
      TOX = 2.D0/X
      BK  = BESSK1(X)
      BKM = BESSK0(X)
      DO 11 J=1,N-1
      BKP = BKM+DFLOAT(J)*TOX*BK
      BKM = BK
      BK  = BKP
   11 CONTINUE
      BESSK = BK
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK0(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      REAL*8 X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7, &
      BESSI0
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0&
      &,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1&
      &,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF(X.EQ.0.D0) THEN
      BESSK0=1.D30
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=-LOG(X/2.D0)*BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/SQRT(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSK1(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
! ----------------------------------------------------------------------
      REAL*8 X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,BE&
      &SSI1
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,  &
      -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1&
      &,0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF(X.EQ.0.D0) THEN
      BESSK0=1.D32
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=LOG(X/2.D0)*BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=EXP(-X)/SQRT(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END

! End of file Tbessk.f90

