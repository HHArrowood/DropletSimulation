  DOUBLE PRECISION FUNCTION DRF (X, Y, Z, IER)
!
!! DRF computes the incomplete or complete elliptic integral of 1st kind.
!
!***PURPOSE  Compute the incomplete or complete elliptic integral of the
!            1st kind.  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them  is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling sequence
!          DRF( X, Y, Z, IER )
!
!          Parameters On entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!          Z      - Double precision, nonnegative variable
!
!
!
!          On Return    (values assigned by the DRF routine)
!
!          DRF     - Double precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!
!         Value of IER assigned by the DRF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z)  <  0.0D0
!                      = 2                MIN(X+Y,X+Z,Y+Z)  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable values for:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0D-78     1.0D+75
!                     CDC 6000/7000 SERIES :   1.0D-292    1.0D+321
!                     UNIVAC 1100 SERIES   :   1.0D-307    1.0D+307
!                     CRAY                 :   2.3D-2466   1.09D+2465
!                     VAX 11 SERIES        :   1.5D-38     3.0D+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!
!          Sample choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0D-3    3.0D-19
!                           3.0D-3    2.0D-16
!                           1.0D-2    3.0D-13
!                           3.0D-2    2.0D-10
!                           1.0D-1    3.0D-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRF Special Comments
!
!
!
!          Check by addition theorem: DRF(X,X+Z,X+W) + DRF(Y,Y+Z,Y+W)
!          = DRF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special double precision functions via DRF
!
!
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!
!                                             2         2   2
!                  F(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2
!                  K(K) = DRF(0,1-K ,1)
!
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!
!                  -----------------------------------------
!
!
!                                          22    2
!                  EL1(X,KC) = X DRF(1,1+KC X ,1+X )
!
!
!                  Lemniscate constant A
!
!                  -----------------------------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = DRF(0,1,2) = DRF(0,2,1)
!                       0
!
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRF
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  DOUBLE PRECISION LOLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION C1, C2, C3, E2, E3, LAMDA
  DOUBLE PRECISION MU, S, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRF
!
  if (FIRST) THEN
     ERRTOL = (4.0D0*D1MACH(3))**(1.0D0/6.0D0)
     LOLIM  = 5.0D0 * D1MACH(1)
     UPLIM  = D1MACH(2)/5.0D0
!
     C1 = 1.0D0/24.0D0
     C2 = 3.0D0/44.0D0
     C3 = 1.0D0/14.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRF = 0.0D0
  if (MIN(X,Y,Z) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MAX(X,Y,Z) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'DRF', &
        'MIN(X+Y,X+Z,Y+Z) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // &
        XERN6, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0D0
  XNDEV = 2.0D0 - (MU+XN)/MU
  YNDEV = 2.0D0 - (MU+YN)/MU
  ZNDEV = 2.0D0 - (MU+ZN)/MU
  EPSLON = MAX(ABS(XNDEV),ABS(YNDEV),ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT = SQRT(XN)
  YNROOT = SQRT(YN)
  ZNROOT = SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  go to 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
  E3 = XNDEV*YNDEV*ZNDEV
  S  = 1.0D0 + (C1*E2-0.10D0-C2*E3)*E2 + C3*E3
  DRF = S/SQRT(MU)
!
  return
end

FUNCTION RF (X, Y, Z, IER)
!
!! RF computes the incomplete or complete elliptic integral of the 1st kind.  
!
!  For X, Y, and Z non-negative and at most one of
!            them zero, RF(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -1/2
!                      (1/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X, Y or Z is zero, the integral is complete.
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      SINGLE PRECISION (RF-S, DRF-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE FIRST KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     RF
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the first kind
!          Standard FORTRAN function routine
!          Single precision version
!          The routine calculates an approximation result to
!          RF(X,Y,Z) = Integral from zero to infinity of
!
!                               -1/2     -1/2     -1/2
!                     (1/2)(t+X)    (t+Y)    (t+Z)    dt,
!
!          where X, Y, and Z are nonnegative and at most one of them
!          is zero.  If one of them is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!          RF( X, Y, Z, IER )
!
!          Parameters on Entry
!          Values assigned by the calling routine
!
!          X      - Single precision, nonnegative variable
!
!          Y      - Single precision, nonnegative variable
!
!          Z      - Single precision, nonnegative variable
!
!
!
!          On Return     (values assigned by the RF routine)
!
!          RF     - Single precision approximation to the integral
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine.  It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!          X, Y, Z are unaltered.
!
!
!   3.    Error Messages
!
!         Value of IER assigned by the RF routine
!
!                  Value assigned         Error Message Printed
!                  IER = 1                MIN(X,Y,Z)  <  0.0E0
!                      = 2                MIN(X+Y,X+Z,Y+Z)  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                   Not less than 5 * (machine minimum).
!
!          UPLIM  - Upper limit of valid arguments
!
!                   Not greater than (machine maximum) / 5.
!
!
!                     Acceptable Values For:   LOLIM      UPLIM
!                     IBM 360/370 SERIES   :   3.0E-78     1.0E+75
!                     CDC 6000/7000 SERIES :   1.0E-292    1.0E+321
!                     UNIVAC 1100 SERIES   :   1.0E-37     1.0E+37
!                     CRAY                 :   2.3E-2466   1.09E+2465
!                     VAX 11 SERIES        :   1.5E-38     3.0E+37
!
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!
!
!          ERRTOL - Relative error due to truncation is less than
!                   ERRTOL ** 6 / (4 * (1-ERRTOL)  .
!
!
!
!              The accuracy of the computed approximation to the inte-
!              gral can be controlled by choosing the value of ERRTOL.
!              Truncation of a Taylor series after terms of fifth order
!              introduces an error less than the amount shown in the
!              second column of the following table for each value of
!              ERRTOL in the first column.  In addition to the trunca-
!              tion error there will be round-off error, but in prac-
!              tice the total error from both sources is usually less
!              than the amount given in the table.
!
!
!
!
!
!          Sample Choices:  ERRTOL   Relative Truncation
!                                    error less than
!                           1.0E-3    3.0E-19
!                           3.0E-3    2.0E-16
!                           1.0E-2    3.0E-13
!                           3.0E-2    2.0E-10
!                           1.0E-1    3.0E-7
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   RF Special Comments
!
!
!
!          Check by addition theorem: RF(X,X+Z,X+W) + RF(Y,Y+Z,Y+W)
!          = RF(0,Z,W), where X,Y,Z,W are positive and X * Y = Z * W.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral RF(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, and Z are unaltered.
!
!
!
!          ********************************************************
!
!          Warning: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!   Special Functions via RF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 1st kind
!                  ----------------------------------------------
!
!
!                                            2         2   2
!                  F(PHI,K) = SIN(PHI) RF(COS (PHI),1-K SIN (PHI),1)
!
!
!                                 2
!                  K(K) = RF(0,1-K ,1)
!
!                         PI/2     2   2      -1/2
!                       = INT  (1-K SIN (PHI) )   D PHI
!                          0
!
!
!
!
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 1st kind
!                  ----------------------------------------------
!
!
!                                         22    2
!                  EL1(X,KC) = X RF(1,1+KC X ,1+X )
!
!
!
!
!                  Lemniscate constant A
!                  ---------------------
!
!
!                       1      4 -1/2
!                  A = INT (1-S )    DS = RF(0,1,2) = RF(0,2,1)
!                       0
!
!
!    -------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Changed calls to XERMSG to standard form, and some
!           editorial changes.  (RWC))
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RF
  real RF
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  REAL LOLIM, UPLIM, EPSLON, ERRTOL
  REAL C1, C2, C3, E2, E3, LAMDA
  REAL MU, S, X, XN, XNDEV
  REAL XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL,LOLIM,UPLIM,C1,C2,C3,FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  RF
!
  if (FIRST) THEN
     ERRTOL = (4.0E0*R1MACH(3))**(1.0E0/6.0E0)
     LOLIM  = 5.0E0 * R1MACH(1)
     UPLIM  = R1MACH(2)/5.0E0
!
     C1 = 1.0E0/24.0E0
     C2 = 3.0E0/44.0E0
     C3 = 1.0E0/14.0E0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  RF = 0.0E0
  if (MIN(X,Y,Z) < 0.0E0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     call XERMSG ('SLATEC', 'RF', &
        'MIN(X,Y,Z) < 0 WHERE X = ' // XERN3 // ' Y = ' // XERN4 // &
        ' AND Z = ' // XERN5, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'RF', &
        'MAX(X,Y,Z) > UPLIM WHERE X = '  // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, 3, 1)
     return
  end if
!
  if (MIN(X+Y,X+Z,Y+Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'RF', &
        'MIN(X+Y,X+Z,Y+Z) < LOLIM WHERE X = ' // XERN3 // &
        ' Y = ' // XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // &
        XERN6, 2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
!
   30 MU = (XN+YN+ZN)/3.0E0
  XNDEV = 2.0E0 - (MU+XN)/MU
  YNDEV = 2.0E0 - (MU+YN)/MU
  ZNDEV = 2.0E0 - (MU+ZN)/MU
  EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT =  SQRT(XN)
  YNROOT =  SQRT(YN)
  ZNROOT =  SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  XN = (XN+LAMDA)*0.250E0
  YN = (YN+LAMDA)*0.250E0
  ZN = (ZN+LAMDA)*0.250E0
  go to 30
!
   40 E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
  E3 = XNDEV*YNDEV*ZNDEV
  S  = 1.0E0 + (C1*E2-0.10E0-C2*E3)*E2 + C3*E3
  RF = S/SQRT(MU)
!
  return
end

DOUBLE PRECISION FUNCTION DRD (X, Y, Z, IER)
!
!! DRD computes the incomplete or complete elliptic integral of 2nd kind.
!
!***PURPOSE  Compute the incomplete or complete elliptic integral of
!            the 2nd kind. For X and Y nonnegative, X+Y and Z positive,
!            DRD(X,Y,Z) = Integral from zero to infinity of
!                                -1/2     -1/2     -3/2
!                      (3/2)(t+X)    (t+Y)    (t+Z)    dt.
!            If X or Y is zero, the integral is complete.
!
!***LIBRARY   SLATEC
!***CATEGORY  C14
!***TYPE      DOUBLE PRECISION (RD-S, DRD-D)
!***KEYWORDS  COMPLETE ELLIPTIC INTEGRAL, DUPLICATION THEOREM,
!             INCOMPLETE ELLIPTIC INTEGRAL, INTEGRAL OF THE SECOND KIND,
!             TAYLOR SERIES
!***AUTHOR  Carlson, B. C.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Notis, E. M.
!             Ames Laboratory-DOE
!             Iowa State University
!             Ames, IA  50011
!           Pexton, R. L.
!             Lawrence Livermore National Laboratory
!             Livermore, CA  94550
!***DESCRIPTION
!
!   1.     DRD
!          Evaluate an INCOMPLETE (or COMPLETE) ELLIPTIC INTEGRAL
!          of the second kind
!          Standard FORTRAN function routine
!          Double precision version
!          The routine calculates an approximation result to
!          DRD(X,Y,Z) = Integral from zero to infinity of
!                              -1/2     -1/2     -3/2
!                    (3/2)(t+X)    (t+Y)    (t+Z)    dt,
!          where X and Y are nonnegative, X + Y is positive, and Z is
!          positive.  If X or Y is zero, the integral is COMPLETE.
!          The duplication theorem is iterated until the variables are
!          nearly equal, and the function is then expanded in Taylor
!          series to fifth order.
!
!   2.     Calling Sequence
!
!          DRD( X, Y, Z, IER )
!
!          Parameters On Entry
!          Values assigned by the calling routine
!
!          X      - Double precision, nonnegative variable
!
!          Y      - Double precision, nonnegative variable
!
!                   X + Y is positive
!
!          Z      - Double precision, positive variable
!
!
!
!          On Return    (values assigned by the DRD routine)
!
!          DRD     - Double precision approximation to the integral
!
!
!          IER    - Integer
!
!                   IER = 0 Normal and reliable termination of the
!                           routine. It is assumed that the requested
!                           accuracy has been achieved.
!
!                   IER >  0 Abnormal termination of the routine
!
!
!          X, Y, Z are unaltered.
!
!   3.    Error Messages
!
!         Value of IER assigned by the DRD routine
!
!                  Value assigned         Error message printed
!                  IER = 1                MIN(X,Y)  <  0.0D0
!                      = 2                MIN(X + Y, Z )  <  LOLIM
!                      = 3                MAX(X,Y,Z)  >  UPLIM
!
!
!   4.     Control Parameters
!
!                  Values of LOLIM, UPLIM, and ERRTOL are set by the
!                  routine.
!
!          LOLIM and UPLIM determine the valid range of X, Y, and Z
!
!          LOLIM  - Lower limit of valid arguments
!
!                    Not less  than 2 / (machine maximum) ** (2/3).
!
!          UPLIM  - Upper limit of valid arguments
!
!                 Not greater than (0.1D0 * ERRTOL / machine
!                 minimum) ** (2/3), where ERRTOL is described below.
!                 In the following table it is assumed that ERRTOL will
!                 never be chosen smaller than 1.0D-5.
!
!
!                    Acceptable values for:   LOLIM      UPLIM
!                    IBM 360/370 SERIES   :   6.0D-51     1.0D+48
!                    CDC 6000/7000 SERIES :   5.0D-215    2.0D+191
!                    UNIVAC 1100 SERIES   :   1.0D-205    2.0D+201
!                    CRAY                 :   3.0D-1644   1.69D+1640
!                    VAX 11 SERIES        :   1.0D-25     4.5D+21
!
!
!          ERRTOL determines the accuracy of the answer
!
!                 The value assigned by the routine will result
!                 in solution precision within 1-2 decimals of
!                 "machine precision".
!
!          ERRTOL    Relative error due to truncation is less than
!                    3 * ERRTOL ** 6 / (1-ERRTOL) ** 3/2.
!
!
!
!        The accuracy of the computed approximation to the integral
!        can be controlled by choosing the value of ERRTOL.
!        Truncation of a Taylor series after terms of fifth order
!        introduces an error less than the amount shown in the
!        second column of the following table for each value of
!        ERRTOL in the first column.  In addition to the truncation
!        error there will be round-off error, but in practice the
!        total error from both sources is usually less than the
!        amount given in the table.
!
!
!
!
!          Sample choices:  ERRTOL   Relative truncation
!                                    error less than
!                           1.0D-3    4.0D-18
!                           3.0D-3    3.0D-15
!                           1.0D-2    4.0D-12
!                           3.0D-2    3.0D-9
!                           1.0D-1    4.0D-6
!
!
!                    Decreasing ERRTOL by a factor of 10 yields six more
!                    decimal digits of accuracy at the expense of one or
!                    two more iterations of the duplication theorem.
!
! *Long Description:
!
!   DRD Special Comments
!
!
!
!          Check: DRD(X,Y,Z) + DRD(Y,Z,X) + DRD(Z,X,Y)
!          = 3 / SQRT(X * Y * Z), where X, Y, and Z are positive.
!
!
!          On Input:
!
!          X, Y, and Z are the variables in the integral DRD(X,Y,Z).
!
!
!          On Output:
!
!
!          X, Y, Z are unaltered.
!
!
!
!          ********************************************************
!
!          WARNING: Changes in the program may improve speed at the
!                   expense of robustness.
!
!
!
!    -------------------------------------------------------------------
!
!
!   Special double precision functions via DRD and DRF
!
!
!                  Legendre form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!
!                                             2         2   2
!                  E(PHI,K) = SIN(PHI) DRF(COS (PHI),1-K SIN (PHI),1) -
!
!                     2      3             2         2   2
!                  -(K/3) SIN (PHI) DRD(COS (PHI),1-K SIN (PHI),1)
!
!
!                                  2        2            2
!                  E(K) = DRF(0,1-K ,1) - (K/3) DRD(0,1-K ,1)
!
!                         PI/2     2   2      1/2
!                       = INT  (1-K SIN (PHI) )  D PHI
!                          0
!
!                  Bulirsch form of ELLIPTIC INTEGRAL of 2nd kind
!
!                  -----------------------------------------
!
!                                               22    2
!                  EL2(X,KC,A,B) = AX DRF(1,1+KC X ,1+X ) +
!
!                                              3          22    2
!                                 +(1/3)(B-A) X DRD(1,1+KC X ,1+X )
!
!
!
!
!                  Legendre form of alternative ELLIPTIC INTEGRAL
!                  of 2nd kind
!
!                  -----------------------------------------
!
!
!
!                            Q     2       2   2  -1/2
!                  D(Q,K) = INT SIN P  (1-K SIN P)     DP
!                            0
!
!
!
!                                     3          2     2   2
!                  D(Q,K) = (1/3) (SIN Q) DRD(COS Q,1-K SIN Q,1)
!
!
!
!
!                  Lemniscate constant  B
!
!                  -----------------------------------------
!
!
!
!
!                       1    2    4 -1/2
!                  B = INT  S (1-S )    DS
!                       0
!
!
!                  B = (1/3) DRD (0,2,1)
!
!
!                  Heuman's LAMBDA function
!
!                  -----------------------------------------
!
!
!
!                  (PI/2) LAMBDA0(A,B) =
!
!                                    2                2
!                 = SIN(B) (DRF(0,COS (A),1)-(1/3) SIN (A) *
!
!                            2               2         2       2
!                  *DRD(0,COS (A),1)) DRF(COS (B),1-COS (A) SIN (B),1)
!
!                            2       3             2
!                  -(1/3) COS (A) SIN (B) DRF(0,COS (A),1) *
!
!                           2         2       2
!                   *DRD(COS (B),1-COS (A) SIN (B),1)
!
!
!
!                  Jacobi ZETA function
!
!                  -----------------------------------------
!
!                             2                 2       2   2
!                  Z(B,K) = (K/3) SIN(B) DRF(COS (B),1-K SIN (B),1)
!
!
!                                       2             2
!                             *DRD(0,1-K ,1)/DRF(0,1-K ,1)
!
!                               2       3           2       2   2
!                            -(K /3) SIN (B) DRD(COS (B),1-K SIN (B),1)
!
!
! ---------------------------------------------------------------------
!
!***REFERENCES  B. C. Carlson and E. M. Notis, Algorithms for incomplete
!                 elliptic integrals, ACM Transactions on Mathematical
!                 Software 7, 3 (September 1981), pp. 398-403.
!               B. C. Carlson, Computing elliptic integrals by
!                 duplication, Numerische Mathematik 33, (1979),
!                 pp. 1-16.
!               B. C. Carlson, Elliptic integrals of the first kind,
!                 SIAM Journal of Mathematical Analysis 8, (1977),
!                 pp. 231-242.
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900510  Modify calls to XERMSG to put in standard form.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DRD
  CHARACTER*16 XERN3, XERN4, XERN5, XERN6
  INTEGER IER
  DOUBLE PRECISION LOLIM, TUPLIM, UPLIM, EPSLON, ERRTOL, D1MACH
  DOUBLE PRECISION C1, C2, C3, C4, EA, EB, EC, ED, EF, LAMDA
  DOUBLE PRECISION MU, POWER4, SIGMA, S1, S2, X, XN, XNDEV
  DOUBLE PRECISION XNROOT, Y, YN, YNDEV, YNROOT, Z, ZN, ZNDEV, &
   ZNROOT
  LOGICAL FIRST
  SAVE ERRTOL, LOLIM, UPLIM, C1, C2, C3, C4, FIRST
  DATA FIRST /.TRUE./
!
!***FIRST EXECUTABLE STATEMENT  DRD
  if (FIRST) THEN
     ERRTOL = (D1MACH(3)/3.0D0)**(1.0D0/6.0D0)
     LOLIM  = 2.0D0/(D1MACH(2))**(2.0D0/3.0D0)
     TUPLIM = D1MACH(1)**(1.0E0/3.0E0)
     TUPLIM = (0.10D0*ERRTOL)**(1.0E0/3.0E0)/TUPLIM
     UPLIM  = TUPLIM**2.0D0
!
     C1 = 3.0D0/14.0D0
     C2 = 1.0D0/6.0D0
     C3 = 9.0D0/22.0D0
     C4 = 3.0D0/26.0D0
  end if
  FIRST = .FALSE.
!
!         call ERROR HANDLER if NECESSARY.
!
  DRD = 0.0D0
  if (  MIN(X,Y) < 0.0D0) THEN
     IER = 1
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     call XERMSG ('SLATEC', 'DRD', &
        'MIN(X,Y) < 0 WHERE X = ' // XERN3 // ' AND Y = ' // &
        XERN4, 1, 1)
     return
  end if
!
  if (MAX(X,Y,Z) > UPLIM) THEN
     IER = 3
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') UPLIM
     call XERMSG ('SLATEC', 'DRD', &
        'MAX(X,Y,Z) > UPLIM WHERE X = ' // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND UPLIM = ' // XERN6, &
        3, 1)
     return
  end if
!
  if (MIN(X+Y,Z) < LOLIM) THEN
     IER = 2
     WRITE (XERN3, '(1PE15.6)') X
     WRITE (XERN4, '(1PE15.6)') Y
     WRITE (XERN5, '(1PE15.6)') Z
     WRITE (XERN6, '(1PE15.6)') LOLIM
     call XERMSG ('SLATEC', 'DRD', &
        'MIN(X+Y,Z) < LOLIM WHERE X = ' // XERN3 // ' Y = ' // &
        XERN4 // ' Z = ' // XERN5 // ' AND LOLIM = ' // XERN6, &
        2, 1)
     return
  end if
!
  IER = 0
  XN = X
  YN = Y
  ZN = Z
  SIGMA = 0.0D0
  POWER4 = 1.0D0
!
   30 MU = (XN+YN+3.0D0*ZN)*0.20D0
  XNDEV = (MU-XN)/MU
  YNDEV = (MU-YN)/MU
  ZNDEV = (MU-ZN)/MU
  EPSLON = MAX(ABS(XNDEV), ABS(YNDEV), ABS(ZNDEV))
  if (EPSLON < ERRTOL) go to 40
  XNROOT = SQRT(XN)
  YNROOT = SQRT(YN)
  ZNROOT = SQRT(ZN)
  LAMDA = XNROOT*(YNROOT+ZNROOT) + YNROOT*ZNROOT
  SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
  POWER4 = POWER4*0.250D0
  XN = (XN+LAMDA)*0.250D0
  YN = (YN+LAMDA)*0.250D0
  ZN = (ZN+LAMDA)*0.250D0
  go to 30
!
   40 EA = XNDEV*YNDEV
  EB = ZNDEV*ZNDEV
  EC = EA - EB
  ED = EA - 6.0D0*EB
  EF = ED + EC + EC
  S1 = ED*(-C1+0.250D0*C3*ED-1.50D0*C4*ZNDEV*EF)
  S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
  DRD = 3.0D0*SIGMA + POWER4*(1.0D0+S1+S2)/(MU*SQRT(MU))
!
  return
end

subroutine XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)

!! XERMSG processes XERROR messages.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERMSG-A)
!***KEYWORDS  ERROR MESSAGE, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
!   XERMSG processes a diagnostic message in a manner determined by the
!   value of LEVEL and the current value of the library error control
!   flag, KONTRL.  See subroutine XSETF for details.
!
!    LIBRAR   A character constant (or character variable) with the name
!             of the library.  This will be 'SLATEC' for the SLATEC
!             Common Math Library.  The error handling package is
!             general enough to be used by many libraries
!             simultaneously, so it is desirable for the routine that
!             detects and reports an error to identify the library name
!             as well as the routine name.
!
!    SUBROU   A character constant (or character variable) with the name
!             of the routine that detected the error.  Usually it is the
!             name of the routine that is calling XERMSG.  There are
!             some instances where a user callable library routine calls
!             lower level subsidiary routines where the error is
!             detected.  In such cases it may be more informative to
!             supply the name of the routine the user called rather than
!             the name of the subsidiary routine that detected the
!             error.
!
!    MESSG    A character constant (or character variable) with the text
!             of the error or warning message.  In the example below,
!             the message is a character constant that contains a
!             generic message.
!
!                   call XERMSG ('SLATEC', 'MMPY',
!                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
!                  *3, 1)
!
!             It is possible (and is sometimes desirable) to generate a
!             specific message--e.g., one that contains actual numeric
!             values.  Specific numeric values can be converted into
!             character strings using formatted WRITE statements into
!             character variables.  This is called standard Fortran
!             internal file I/O and is exemplified in the first three
!             lines of the following example.  You can also catenate
!             substrings of characters to construct the error message.
!             Here is an example showing the use of both writing to
!             an internal file and catenating character strings.
!
!                   CHARACTER*5 CHARN, CHARL
!                   WRITE (CHARN,10) N
!                   WRITE (CHARL,10) LDA
!                10 FORMAT(I5)
!                   call XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
!                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
!                  *   CHARL, 3, 1)
!
!             There are two subtleties worth mentioning.  One is that
!             the // for character catenation is used to construct the
!             error message so that no single character constant is
!             continued to the next line.  This avoids confusion as to
!             whether there are trailing blanks at the end of the line.
!             The second is that by catenating the parts of the message
!             as an actual argument rather than encoding the entire
!             message into one large character variable, we avoid
!             having to know how long the message will be in order to
!             declare an adequate length for that large character
!             variable.  XERMSG calls XERPRN to print the message using
!             multiple lines if necessary.  If the message is very long,
!             XERPRN will break it into pieces of 72 characters (as
!             requested by XERMSG) for printing on multiple lines.
!             Also, XERMSG asks XERPRN to prefix each line with ' *  '
!             so that the total line length could be 76 characters.
!             Note also that XERPRN scans the error message backwards
!             to ignore trailing blanks.  Another feature is that
!             the substring '$$' is treated as a new line sentinel
!             by XERPRN.  If you want to construct a multiline
!             message without having to count out multiples of 72
!             characters, just use '$$' as a separator.  '$$'
!             obviously must occur within 72 characters of the
!             start of each line to have its intended effect since
!             XERPRN is asked to wrap around at 72 characters in
!             addition to looking for '$$'.
!
!    NERR     An integer value that is chosen by the library routine's
!             author.  It must be in the range -99 to 999 (three
!             printable digits).  Each distinct error should have its
!             own error number.  These error numbers should be described
!             in the machine readable documentation for the routine.
!             The error numbers need be unique only within each routine,
!             so it is reasonable for each routine to start enumerating
!             errors from 1 and proceeding to the next integer.
!
!    LEVEL    An integer value in the range 0 to 2 that indicates the
!             level (severity) of the error.  Their meanings are
!
!            -1  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.  An attempt is made to only print this
!                message once.
!
!             0  A warning message.  This is used if it is not clear
!                that there really is an error, but the user's attention
!                may be needed.
!
!             1  A recoverable error.  This is used even if the error is
!                so serious that the routine cannot return any useful
!                answer.  If the user has told the error package to
!                return after recoverable errors, then XERMSG will
!                return to the Library routine which can then return to
!                the user's routine.  The user may also permit the error
!                package to terminate the program upon encountering a
!                recoverable error.
!
!             2  A fatal error.  XERMSG will not return to its caller
!                after it receives a fatal error.  This level should
!                hardly ever be used; it is much better to allow the
!                user a chance to recover.  An example of one of the few
!                cases in which it is permissible to declare a level 2
!                error is a reverse communication Library routine that
!                is likely to be called repeatedly until it integrates
!                across some interval.  If there is a serious error in
!                the input such that another step cannot be taken and
!                the Library routine is called again without the input
!                error having been corrected by the caller, the Library
!                routine will probably be called forever with improper
!                input.  In this case, it is reasonable to declare the
!                error to be fatal.
!
!    Each of the arguments to XERMSG is input; none will be modified by
!    XERMSG.  A routine may make multiple calls to XERMSG with warning
!    level messages; however, after a call to XERMSG with a recoverable
!    error, the routine should return to the user.  Do not try to call
!    XERMSG with a second recoverable error after the first recoverable
!    error because the error package saves the error number.  The user
!    can retrieve this error number by calling another entry point in
!    the error handling package and then clear the error number when
!    recovering from the error.  Calling XERMSG in succession causes the
!    old error number to be overwritten by the latest error number.
!    This is considered harmless for error numbers associated with
!    warning messages but must not be done for error numbers of serious
!    errors.  After a call to XERMSG with a recoverable error, the user
!    must be given a chance to call NUMXER or XERCLR to retrieve or
!    clear the error number.
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
!***REVISION HISTORY  (YYMMDD)
!   880101  DATE WRITTEN
!   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
!           THERE ARE TWO BASIC CHANGES.
!           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
!               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
!               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
!               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
!               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
!               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
!               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
!               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
!           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
!               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
!               OF LOWER CASE.
!   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
!           THE PRINCIPAL CHANGES ARE
!           1.  CLARIFY COMMENTS IN THE PROLOGUES
!           2.  RENAME XRPRNT TO XERPRN
!           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
!               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
!               CHARACTER FOR NEW RECORDS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           CLEAN UP THE CODING.
!   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
!           PREFIX.
!   891013  REVISED TO CORRECT COMMENTS.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
!           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
!           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
!           XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERMSG
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8 XLIBR, XSUBR
  CHARACTER*72  TEMP
  CHARACTER*20  LFIRST
!***FIRST EXECUTABLE STATEMENT  XERMSG
  LKNTRL = J4SAVE (2, 0, .FALSE.)
  MAXMES = J4SAVE (4, 0, .FALSE.)
!
!       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
!       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
!          SHOULD BE PRINTED.
!
!       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
!          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
!          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
!
  if (NERR < -9999999 .OR. NERR > 99999999 .OR. NERR == 0 .OR. &
     LEVEL < -1 .OR. LEVEL > 2) THEN
     call XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' // &
        'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '// &
        'JOB ABORT DUE TO FATAL ERROR.', 72)
     call XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
     call XERHLT (' ***XERMSG -- INVALID INPUT')
     return
  end if
!
!       RECORD THE MESSAGE.
!
  I = J4SAVE (1, NERR, .TRUE.)
  call XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
!
!       HANDLE PRINT-ONCE WARNING MESSAGES.
!
  if (LEVEL == -1 .AND. KOUNT > 1) RETURN
!
!       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
!
  XLIBR  = LIBRAR
  XSUBR  = SUBROU
  LFIRST = MESSG
  LERR   = NERR
  LLEVEL = LEVEL
  call XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
!
  LKNTRL = MAX(-2, MIN(2,LKNTRL))
  MKNTRL = ABS(LKNTRL)
!
!       SKIP PRINTING if THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
!       ZERO AND THE ERROR IS NOT FATAL.
!
  if (LEVEL < 2 .AND. LKNTRL == 0) go to 30
  if (LEVEL == 0 .AND. KOUNT > MAXMES) go to 30
  if (LEVEL == 1 .AND. KOUNT > MAXMES .AND. MKNTRL == 1) go to 30
  if (LEVEL == 2 .AND. KOUNT > MAX(1,MAXMES)) go to 30
!
!       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
!       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
!       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY if CONTROL FLAG
!       IS NOT ZERO.
!
  if (LKNTRL  /=  0) THEN
     TEMP(1:21) = 'MESSAGE FROM ROUTINE '
     I = MIN(LEN(SUBROU), 16)
     TEMP(22:21+I) = SUBROU(1:I)
     TEMP(22+I:33+I) = ' IN LIBRARY '
     LTEMP = 33 + I
     I = MIN(LEN(LIBRAR), 16)
     TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
     TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
     LTEMP = LTEMP + I + 1
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       if LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
!       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
!       FROM EACH OF THE FOLLOWING THREE OPTIONS.
!       1.  LEVEL OF THE MESSAGE
!              'INFORMATIVE MESSAGE'
!              'POTENTIALLY RECOVERABLE ERROR'
!              'FATAL ERROR'
!       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
!              'PROG CONTINUES'
!              'PROG ABORTED'
!       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
!           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
!           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
!              'TRACEBACK REQUESTED'
!              'TRACEBACK NOT REQUESTED'
!       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
!       EXCEED 74 CHARACTERS.
!       WE SKIP THE NEXT BLOCK if THE INTRODUCTORY LINE IS NOT NEEDED.
!
  if (LKNTRL  >  0) THEN
!
!       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
!
     if (LEVEL  <=  0) THEN
        TEMP(1:20) = 'INFORMATIVE MESSAGE,'
        LTEMP = 20
     ELSEIF (LEVEL  ==  1) THEN
        TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
        LTEMP = 30
     ELSE
        TEMP(1:12) = 'FATAL ERROR,'
        LTEMP = 12
     ENDIF
!
!       THEN WHETHER THE PROGRAM WILL CONTINUE.
!
     if ((MKNTRL == 2 .AND. LEVEL >= 1) .OR. &
         (MKNTRL == 1 .AND. LEVEL == 2)) THEN
        TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
        LTEMP = LTEMP + 14
     ELSE
        TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
        LTEMP = LTEMP + 16
     ENDIF
!
!       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
!
     if (LKNTRL  >  0) THEN
        TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
        LTEMP = LTEMP + 20
     ELSE
        TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
        LTEMP = LTEMP + 24
     ENDIF
     call XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
  end if
!
!       NOW SEND OUT THE MESSAGE.
!
  call XERPRN (' *  ', -1, MESSG, 72)
!
!       if LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
!          TRACEBACK.
!
  if (LKNTRL  >  0) THEN
     WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
     DO 10 I=16,22
        if (TEMP(I:I)  /=  ' ') go to 20
   10    CONTINUE
!
   20    call XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
     call FDUMP
  end if
!
!       if LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
!
  if (LKNTRL  /=  0) THEN
     call XERPRN (' *  ', -1, ' ', 72)
     call XERPRN (' ***', -1, 'END OF MESSAGE', 72)
     call XERPRN ('    ',  0, ' ', 72)
  end if
!
!       if THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
!       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
!
   30 if (LEVEL <= 0 .OR. (LEVEL == 1 .AND. MKNTRL <= 1)) RETURN
!
!       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
!       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
!       SUMMARY if THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
!
  if (LKNTRL > 0 .AND. KOUNT < MAX(1,MAXMES)) THEN
     if (LEVEL  ==  1) THEN
        call XERPRN &
           (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
     ELSE
        call XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
     ENDIF
     call XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
     call XERHLT (' ')
  ELSE
     call XERHLT (MESSG)
  end if
  return
end

subroutine XERPRN (PREFIX, NPREF, MESSG, NWRAP)
!
!! XERPRN prints XERROR error messages processed by XERMSG.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERPRN-A)
!***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
!***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
!***DESCRIPTION
!
! This routine sends one or more lines to each of the (up to five)
! logical units to which error messages are to be sent.  This routine
! is called several times by XERMSG, sometimes with a single line to
! print and sometimes with a (potentially very long) message that may
! wrap around into multiple lines.
!
! PREFIX  Input argument of type CHARACTER.  This argument contains
!         characters to be put at the beginning of each line before
!         the body of the message.  No more than 16 characters of
!         PREFIX will be used.
!
! NPREF   Input argument of type INTEGER.  This argument is the number
!         of characters to use from PREFIX.  If it is negative, the
!         intrinsic function LEN is used to determine its length.  If
!         it is zero, PREFIX is not used.  If it exceeds 16 or if
!         LEN(PREFIX) exceeds 16, only the first 16 characters will be
!         used.  If NPREF is positive and the length of PREFIX is less
!         than NPREF, a copy of PREFIX extended with blanks to length
!         NPREF will be used.
!
! MESSG   Input argument of type CHARACTER.  This is the text of a
!         message to be printed.  If it is a long message, it will be
!         broken into pieces for printing on multiple lines.  Each line
!         will start with the appropriate prefix and be followed by a
!         piece of the message.  NWRAP is the number of characters per
!         piece; that is, after each NWRAP characters, we break and
!         start a new line.  In addition the characters '$$' embedded
!         in MESSG are a sentinel for a new line.  The counting of
!         characters up to NWRAP starts over for each new line.  The
!         value of NWRAP typically used by XERMSG is 72 since many
!         older error messages in the SLATEC Library are laid out to
!         rely on wrap-around every 72 characters.
!
! NWRAP   Input argument of type INTEGER.  This gives the maximum size
!         piece into which to break MESSG for printing on multiple
!         lines.  An embedded '$$' ends a line, and the count restarts
!         at the following character.  If a line break does not occur
!         on a blank (it would split a word) that word is moved to the
!         next line.  Values of NWRAP less than 16 will be treated as
!         16.  Values of NWRAP greater than 132 will be treated as 132.
!         The actual line length will be NPREF + NWRAP after NPREF has
!         been adjusted to fall between 0 and 16 and NWRAP has been
!         adjusted to fall between 16 and 132.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   880621  DATE WRITTEN
!   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
!           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
!           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
!           SLASH CHARACTER IN FORMAT STATEMENTS.
!   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
!           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
!           LINES TO BE PRINTED.
!   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
!           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
!   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
!   891214  Prologue converted to Version 4.0 format.  (WRB)
!   900510  Added code to break messages between words.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERPRN
  CHARACTER*(*) PREFIX, MESSG
  INTEGER NPREF, NWRAP
  CHARACTER*148 CBUFF
  INTEGER IU(5), NUNIT
  CHARACTER*2 NEWLIN
  PARAMETER (NEWLIN = '$$')
!***FIRST EXECUTABLE STATEMENT  XERPRN
  call XGETUA(IU,NUNIT)
!
!       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
!       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
!       ERROR MESSAGE UNIT.
!
  N = I1MACH(4)
  DO 10 I=1,NUNIT
     if (IU(I)  ==  0) IU(I) = N
   10 CONTINUE
!
!       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
!       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
!       THE REST OF THIS ROUTINE.
!
  if ( NPREF  <  0 ) THEN
     LPREF = LEN(PREFIX)
  ELSE
     LPREF = NPREF
  end if
  LPREF = MIN(16, LPREF)
  if (LPREF  /=  0) CBUFF(1:LPREF) = PREFIX
!
!       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
!       TIME FROM MESSG TO PRINT ON ONE LINE.
!
  LWRAP = MAX(16, MIN(132, NWRAP))
!
!       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
!
  LENMSG = LEN(MESSG)
  N = LENMSG
  DO 20 I=1,N
     if (MESSG(LENMSG:LENMSG)  /=  ' ') go to 30
     LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
!
!       if THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
!
  if (LENMSG  ==  0) THEN
     CBUFF(LPREF+1:LPREF+1) = ' '
     DO 40 I=1,NUNIT
        WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
     return
  end if
!
!       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
!       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
!       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
!       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
!
!       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
!       INDEX INTRINSIC FUNCTION RETURNS ZERO if THERE IS NO OCCURRENCE
!       OR if THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
!       OF THE SECOND ARGUMENT.
!
!       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
!       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
!       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
!       POSITION NEXTC.
!
!       LPIECE  ==  0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
!                       REMAINDER OF THE CHARACTER STRING.  LPIECE
!                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
!                       WHICHEVER IS LESS.
!
!       LPIECE  ==  1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
!                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
!                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
!                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
!                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
!                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
!                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
!                       SHOULD BE INCREMENTED BY 2.
!
!       LPIECE  >  LWRAP+1  REDUCE LPIECE TO LWRAP.
!
!       ELSE            THIS LAST CASE MEANS 2  <=  LPIECE  <=  LWRAP+1
!                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
!                       PROPERLY HANDLES THE END CASE WHERE LPIECE  ==
!                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
!                       AT THE END OF A LINE.
!
  NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
  if (LPIECE  ==  0) THEN
!
!       THERE WAS NO NEW LINE SENTINEL FOUND.
!
     IDELTA = 0
     LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
     if (LPIECE  <  LENMSG+1-NEXTC) THEN
        DO 52 I=LPIECE+1,2,-1
           if (MESSG(NEXTC+I-1:NEXTC+I-1)  ==  ' ') THEN
              LPIECE = I-1
              IDELTA = 1
              GOTO 54
           ENDIF
   52       CONTINUE
     ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSEIF (LPIECE  ==  1) THEN
!
!       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
!       DON'T PRINT A BLANK LINE.
!
     NEXTC = NEXTC + 2
     go to 50
  ELSEIF (LPIECE  >  LWRAP+1) THEN
!
!       LPIECE SHOULD BE SET DOWN TO LWRAP.
!
     IDELTA = 0
     LPIECE = LWRAP
     DO 56 I=LPIECE+1,2,-1
        if (MESSG(NEXTC+I-1:NEXTC+I-1)  ==  ' ') THEN
           LPIECE = I-1
           IDELTA = 1
           GOTO 58
        ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC = NEXTC + LPIECE + IDELTA
  ELSE
!
!       if WE ARRIVE HERE, IT MEANS 2  <=  LPIECE  <=  LWRAP+1.
!       WE SHOULD DECREMENT LPIECE BY ONE.
!
     LPIECE = LPIECE - 1
     CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
     NEXTC  = NEXTC + LPIECE + 2
  end if
!
!       PRINT
!
  DO 60 I=1,NUNIT
     WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
!
  if (NEXTC  <=  LENMSG) go to 50
  return
end
subroutine XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, &
     ICOUNT)
!
!! XERSVE records that an XERROR error has occurred.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (XERSVE-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
!        CHARACTER * (len) LIBRAR, SUBROU, MESSG
!
!        call XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
!
! *Arguments:
!
!        LIBRAR :IN    is the library that the message is from.
!        SUBROU :IN    is the subroutine that the message is from.
!        MESSG  :IN    is the message to be saved.
!        KFLAG  :IN    indicates the action to be performed.
!                      when KFLAG > 0, the message in MESSG is saved.
!                      when KFLAG=0 the tables will be dumped and
!                      cleared.
!                      when KFLAG < 0, the tables will be dumped and
!                      not cleared.
!        NERR   :IN    is the error number.
!        LEVEL  :IN    is the error severity.
!        ICOUNT :OUT   the number of times this message has been seen,
!                      or zero if the table has overflowed and does not
!                      contain this message specifically.  When KFLAG=0,
!                      ICOUNT will not be altered.
!
! *Description:
!
!   Record that this error occurred and possibly dump and clear the
!   tables.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  I1MACH, XGETUA
!***REVISION HISTORY  (YYMMDD)
!   800319  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900413  Routine modified to remove reference to KFLAG.  (WRB)
!   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
!           sequence, use IF-THEN-ELSE, make number of saved entries
!           easily changeable, changed routine name from XERSAV to
!           XERSVE.  (RWC)
!   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERSVE
  PARAMETER (LENTAB=10)
  INTEGER LUN(5)
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
  CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
  CHARACTER*20 MESTAB(LENTAB), MES
  DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
  SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
  DATA KOUNTX/0/, NMSG/0/
!***FIRST EXECUTABLE STATEMENT  XERSVE
!
  if (KFLAG <= 0) THEN
!
!        Dump the table.
!
     if (NMSG == 0) RETURN
!
!        Print to each unit.
!
     call XGETUA (LUN, NUNIT)
     DO 20 KUNIT = 1,NUNIT
        IUNIT = LUN(KUNIT)
        if (IUNIT == 0) IUNIT = I1MACH(4)
!
!           Print the table header.
!
        WRITE (IUNIT,9000)
!
!           Print body of table.
!
        DO 10 I = 1,NMSG
           WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I), &
              NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
!
!           Print number of other errors.
!
        if (KOUNTX /= 0) WRITE (IUNIT,9020) KOUNTX
        WRITE (IUNIT,9030)
   20    CONTINUE
!
!        Clear the error tables.
!
     if (KFLAG == 0) THEN
        NMSG = 0
        KOUNTX = 0
     ENDIF
  ELSE
!
!        PROCESS A MESSAGE...
!        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
!        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
!
     LIB = LIBRAR
     SUB = SUBROU
     MES = MESSG
     DO 30 I = 1,NMSG
        if (LIB == LIBTAB(I) .AND. SUB == SUBTAB(I) .AND. &
           MES == MESTAB(I) .AND. NERR == NERTAB(I) .AND. &
           LEVEL == LEVTAB(I)) THEN
              KOUNT(I) = KOUNT(I) + 1
              ICOUNT = KOUNT(I)
              return
        ENDIF
   30    CONTINUE
!
     if (NMSG < LENTAB) THEN
!
!           Empty slot found for new message.
!
        NMSG = NMSG + 1
        LIBTAB(I) = LIB
        SUBTAB(I) = SUB
        MESTAB(I) = MES
        NERTAB(I) = NERR
        LEVTAB(I) = LEVEL
        KOUNT (I) = 1
        ICOUNT    = 1
     ELSE
!
!           Table is full.
!
        KOUNTX = KOUNTX+1
        ICOUNT = 0
     ENDIF
  end if
  return
!
!     Formats.
!
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' / &
     ' LIBRARY    SUBROUTINE MESSAGE START             NERR', &
     '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
end

subroutine XERHLT (MESSG)
!
!! XERHLT aborts program execution after printing XERROR error message.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERHLT-A)
!***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        ***Note*** machine dependent routine
!        XERHLT aborts the execution of the program.
!        The error message causing the abort is given in the calling
!        sequence, in case one needs it for printing on a dayfile,
!        for example.
!
!     Description of Parameters
!        MESSG is as in XERMSG.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to delete length of character
!           and changed routine name from XERABT to XERHLT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERHLT
  CHARACTER*(*) MESSG
!***FIRST EXECUTABLE STATEMENT  XERHLT
  STOP
end

subroutine XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
!
!! XERCNT allows user control over handling of XERROR errors.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XERCNT-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        Allows user control over handling of individual errors.
!        Just after each message is recorded, but before it is
!        processed any further (i.e., before it is printed or
!        a decision to abort is made), a call is made to XERCNT.
!        If the user has provided his own version of XERCNT, he
!        can then override the value of KONTROL used in processing
!        this message by redefining its value.
!        KONTRL may be set to any value from -2 to 2.
!        The meanings for KONTRL are the same as in XSETF, except
!        that the value of KONTRL changes only for this message.
!        If KONTRL is set to a value outside the range from -2 to 2,
!        it will be moved back into that range.
!
!     Description of Parameters
!
!      --Input--
!        LIBRAR - the library that the routine is in.
!        SUBROU - the subroutine that XERMSG is being called from
!        MESSG  - the first 20 characters of the error message.
!        NERR   - same as in the call to XERMSG.
!        LEVEL  - same as in the call to XERMSG.
!        KONTRL - the current value of the control flag as set
!                 by a call to XSETF.
!
!      --Output--
!        KONTRL - the new value of KONTRL.  If KONTRL is not
!                 defined, it will remain at its original value.
!                 This changed value of control affects only
!                 the current occurrence of the current message.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900206  Routine changed from user-callable to subsidiary.  (WRB)
!   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
!           names, changed routine name from XERCTL to XERCNT.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XERCNT
  CHARACTER*(*) LIBRAR, SUBROU, MESSG
!***FIRST EXECUTABLE STATEMENT  XERCNT
  return
end

subroutine FDUMP
!
!! FDUMP makes a symbolic dump (should be locally written).
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3
!***TYPE      ALL (FDUMP-A)
!***KEYWORDS  ERROR, XERMSG
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!        ***Note*** Machine Dependent Routine
!        FDUMP is intended to be replaced by a locally written
!        version which produces a symbolic dump.  Failing this,
!        it should be replaced by a version which prints the
!        subprogram nesting list.  Note that this dump must be
!        printed on each of up to five files, as indicated by the
!        XGETUA routine.  See XSETUA and XGETUA for details.
!
!     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  FDUMP
!***FIRST EXECUTABLE STATEMENT  FDUMP
  return
end

subroutine XGETUA (IUNITA, N)
!
!! XGETUA returns unit numbers to which XERROR messages are sent.
!
!***LIBRARY   SLATEC (XERROR)
!***CATEGORY  R3C
!***TYPE      ALL (XGETUA-A)
!***KEYWORDS  ERROR, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        XGETUA may be called to determine the unit number or numbers
!        to which error messages are being sent.
!        These unit numbers may have been set by a call to XSETUN,
!        or a call to XSETUA, or may be a default value.
!
!     Description of Parameters
!      --Output--
!        IUNIT - an array of one to five unit numbers, depending
!                on the value of N.  A value of zero refers to the
!                default unit, as defined by the I1MACH machine
!                constant routine.  Only IUNIT(1),...,IUNIT(N) are
!                defined by XGETUA.  The values of IUNIT(N+1),...,
!                IUNIT(5) are not defined (for N  <  5) or altered
!                in any way by XGETUA.
!        N     - the number of units to which copies of the
!                error messages are being sent.  N will be in the
!                range from 1 to 5.
!
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  J4SAVE
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  XGETUA
  DIMENSION IUNITA(5)
!***FIRST EXECUTABLE STATEMENT  XGETUA
  N = J4SAVE(5,0,.FALSE.)
  DO 30 I=1,N
     INDEX = I+4
     if (I == 1) INDEX = 3
     IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
  return
end

function J4SAVE (IWHICH, IVALUE, ISET)
!
!! J4SAVE saves or recalls global variables needed by error handling routines.
!
!***LIBRARY   SLATEC (XERROR)
!***TYPE      INTEGER (J4SAVE-I)
!***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        J4SAVE saves and recalls several global variables needed
!        by the library error handling routines.
!
!     Description of Parameters
!      --Input--
!        IWHICH - Index of item desired.
!                = 1 Refers to current error number.
!                = 2 Refers to current error control flag.
!                = 3 Refers to current unit number to which error
!                    messages are to be sent.  (0 means use standard.)
!                = 4 Refers to the maximum number of times any
!                     message is to be printed (as set by XERMAX).
!                = 5 Refers to the total number of units to which
!                     each error message is to be written.
!                = 6 Refers to the 2nd unit for error messages
!                = 7 Refers to the 3rd unit for error messages
!                = 8 Refers to the 4th unit for error messages
!                = 9 Refers to the 5th unit for error messages
!        IVALUE - The value to be set for the IWHICH-th parameter,
!                 if ISET is .TRUE. .
!        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
!                 given the value, IVALUE.  If ISET=.FALSE., the
!                 IWHICH-th parameter will be unchanged, and IVALUE
!                 is a dummy parameter.
!      --Output--
!        The (old) value of the IWHICH-th parameter will be returned
!        in the function value, J4SAVE.
!
!***SEE ALSO  XERMSG
!***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
!                 Error-handling Package, SAND82-0800, Sandia
!                 Laboratories, 1982.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900205  Minor modifications to prologue.  (WRB)
!   900402  Added TYPE section.  (WRB)
!   910411  Added KEYWORDS section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  J4SAVE
  LOGICAL ISET
  INTEGER IPARAM(9)
  SAVE IPARAM
  DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
  DATA IPARAM(5)/1/
  DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
!***FIRST EXECUTABLE STATEMENT  J4SAVE
  J4SAVE = IPARAM(IWHICH)
  if (ISET) IPARAM(IWHICH) = IVALUE
  return
end

FUNCTION D1MACH (I)
!
!! D1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        D = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of D above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH( 3) = B**(-T), the smallest relative spacing.
!   D1MACH( 4) = B**(1-T), the largest relative spacing.
!   D1MACH( 5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(14) = T, the number of base-B digits.
!   I1MACH(15) = EMIN, the smallest exponent E.
!   I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   900911  Added SUN 386i constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  D1MACH
!
  double precision d1mach
  INTEGER SMALL(4)
  INTEGER LARGE(4)
  INTEGER RIGHT(4)
  INTEGER DIVER(4)
  INTEGER LOG10(4)
!
  DOUBLE PRECISION DMACH(5)
  SAVE DMACH
!
  EQUIVALENCE (DMACH(1),SMALL(1))
  EQUIVALENCE (DMACH(2),LARGE(1))
  EQUIVALENCE (DMACH(3),RIGHT(1))
  EQUIVALENCE (DMACH(4),DIVER(1))
  EQUIVALENCE (DMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
!     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
!     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
!     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA SMALL(1) / ZC00800000 /
!     DATA SMALL(2) / Z000000000 /
!     DATA LARGE(1) / ZDFFFFFFFF /
!     DATA LARGE(2) / ZFFFFFFFFF /
!     DATA RIGHT(1) / ZCC5800000 /
!     DATA RIGHT(2) / Z000000000 /
!     DATA DIVER(1) / ZCC6800000 /
!     DATA DIVER(2) / Z000000000 /
!     DATA LOG10(1) / ZD00E730E7 /
!     DATA LOG10(2) / ZC77800DC0 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O0000000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O0007777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA SMALL(1) / O1771000000000000 /
!     DATA SMALL(2) / O7770000000000000 /
!     DATA LARGE(1) / O0777777777777777 /
!     DATA LARGE(2) / O7777777777777777 /
!     DATA RIGHT(1) / O1461000000000000 /
!     DATA RIGHT(2) / O0000000000000000 /
!     DATA DIVER(1) / O1451000000000000 /
!     DATA DIVER(2) / O0000000000000000 /
!     DATA LOG10(1) / O1157163034761674 /
!     DATA LOG10(2) / O0006677466732724 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA SMALL(1) / Z"3001800000000000" /
!     DATA SMALL(2) / Z"3001000000000000" /
!     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
!     DATA LARGE(2) / Z"4FFE000000000000" /
!     DATA RIGHT(1) / Z"3FD2800000000000" /
!     DATA RIGHT(2) / Z"3FD2000000000000" /
!     DATA DIVER(1) / Z"3FD3800000000000" /
!     DATA DIVER(2) / Z"3FD3000000000000" /
!     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
!     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA SMALL(1) / 00564000000000000000B /
!     DATA SMALL(2) / 00000000000000000000B /
!     DATA LARGE(1) / 37757777777777777777B /
!     DATA LARGE(2) / 37157777777777777777B /
!     DATA RIGHT(1) / 15624000000000000000B /
!     DATA RIGHT(2) / 00000000000000000000B /
!     DATA DIVER(1) / 15634000000000000000B /
!     DATA DIVER(2) / 00000000000000000000B /
!     DATA LOG10(1) / 17164642023241175717B /
!     DATA LOG10(2) / 16367571421742254654B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn OR -pd8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CC0000000000000' /
!     DATA DMACH(4) / Z'3CD0000000000000' /
!     DATA DMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA SMALL(1) / 201354000000000000000B /
!     DATA SMALL(2) / 000000000000000000000B /
!     DATA LARGE(1) / 577767777777777777777B /
!     DATA LARGE(2) / 000007777777777777774B /
!     DATA RIGHT(1) / 376434000000000000000B /
!     DATA RIGHT(2) / 000000000000000000000B /
!     DATA DIVER(1) / 376444000000000000000B /
!     DATA DIVER(2) / 000000000000000000000B /
!     DATA LOG10(1) / 377774642023241175717B /
!     DATA LOG10(2) / 000007571421742254654B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC DMACH(5)
!
!     DATA SMALL /    20K, 3*0 /
!     DATA LARGE / 77777K, 3*177777K /
!     DATA RIGHT / 31420K, 3*0 /
!     DATA DIVER / 32020K, 3*0 /
!     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA DMACH(1) / '0000000000000010'X /
!     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
!     DATA DMACH(3) / '0000000000003CC0'X /
!     DATA DMACH(4) / '0000000000003CD0'X /
!     DATA DMACH(5) / '79FF509F44133FF3'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FORMAT
!
      DATA DMACH(1) / Z'0010000000000000' /
      DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
      DATA DMACH(3) / Z'3CA0000000000000' /
      DATA DMACH(4) / Z'3CB0000000000000' /
      DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
!     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
!     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
!     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /        128,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
!     DATA DIVER(1), DIVER(2) /       9472,           0 /
!     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
!
!     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1), SMALL(2) /         16,           0 /
!     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
!     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
!     DATA DIVER(1), DIVER(2) /      15568,           0 /
!     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
!
!     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
!     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
!
!     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
!     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
!     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
!     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
!     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
!     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
!     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
!     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
!     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) /  40000B,       0 /
!     DATA SMALL(3), SMALL(4) /       0,       1 /
!     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
!     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
!     DATA RIGHT(3), RIGHT(4) /       0,    225B /
!     DATA DIVER(1), DIVER(2) /  40000B,       0 /
!     DATA DIVER(3), DIVER(4) /       0,    227B /
!     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
!     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
!     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
!     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
!     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
!     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
!     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
!     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
!     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
!     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
!     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
!
!     DATA SMALL(1) / 2.23D-308  /
!     DATA LARGE(1) / 1.79D+308  /
!     DATA RIGHT(1) / 1.11D-16   /
!     DATA DIVER(1) / 2.22D-16   /
!     DATA LOG10(1) / 0.301029995663981195D0 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
!     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
!     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
!     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
!     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
!     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    8388608,           0 /
!     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
!     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
!     DATA DIVER(1), DIVER(2) /  620756992,           0 /
!     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
!
!     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
!     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
!     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
!     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
!     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /    128,      0 /
!     DATA SMALL(3), SMALL(4) /      0,      0 /
!     DATA LARGE(1), LARGE(2) /  32767,     -1 /
!     DATA LARGE(3), LARGE(4) /     -1,     -1 /
!     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
!     DATA RIGHT(3), RIGHT(4) /      0,      0 /
!     DATA DIVER(1), DIVER(2) /   9472,      0 /
!     DATA DIVER(3), DIVER(4) /      0,      0 /
!     DATA LOG10(1), LOG10(2) /  16282,   8346 /
!     DATA LOG10(3), LOG10(4) / -31493, -12296 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA SMALL(3), SMALL(4) / O000000, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA LARGE(3), LARGE(4) / O177777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
!     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
!     DATA DIVER(1), DIVER(2) / O022400, O000000 /
!     DATA DIVER(3), DIVER(4) / O000000, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020232 /
!     DATA LOG10(3), LOG10(4) / O102373, O147770 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
!     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
!     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
!     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA DMACH(1) / Z'0010000000000000' /
!     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3CA0000000000000' /
!     DATA DMACH(4) / Z'3CB0000000000000' /
!     DATA DMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA DMACH(1) / Z'00010000000000000000000000000000' /
!     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
!     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
!     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
!     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
!
!     MACHINE CONSTANTS FOR THE SUN 386i
!
!     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
!     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
!     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
!     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
!     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
!     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
!     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
!     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
!     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  D1MACH = DMACH(I)

  return
end


FUNCTION I1MACH (I)
!
!! I1MACH returns integer machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0  <=  X(I)  <  A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0  <=  X(I)  <  B for I=1,...,T,
!                0  <  X(1), and EMIN  <=  E  <=  EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   891012  Added VAX G-floating constants.  (WRB)
!   891012  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   901009  Correct I1MACH(7) for IBM Mainframes. Should be 2 not 16.
!           (RWC)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added Convex -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!   930618  Corrected I1MACH(5) for Convex -p8 and -pd8 compiler
!           options.  (DWL, RWC and WRB).
!***END PROLOGUE  I1MACH
!
  integer i1mach
  INTEGER IMACH(16),OUTPUT
  SAVE IMACH
  EQUIVALENCE (IMACH(4),OUTPUT)
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        129 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1025 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA IMACH( 1) /          7 /
!     DATA IMACH( 2) /          2 /
!     DATA IMACH( 3) /          2 /
!     DATA IMACH( 4) /          2 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         33 /
!     DATA IMACH( 9) / Z1FFFFFFFF /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -256 /
!     DATA IMACH(13) /        255 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /       -256 /
!     DATA IMACH(16) /        255 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /        -50 /
!     DATA IMACH(16) /         76 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         48 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         39 /
!     DATA IMACH( 9) / O0007777777777777 /
!     DATA IMACH(10) /          8 /
!     DATA IMACH(11) /         13 /
!     DATA IMACH(12) /        -50 /
!     DATA IMACH(13) /         76 /
!     DATA IMACH(14) /         26 /
!     DATA IMACH(15) /     -32754 /
!     DATA IMACH(16) /      32780 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -4095 /
!     DATA IMACH(13) /       4094 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -4095 /
!     DATA IMACH(16) /       4094 /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /    6LOUTPUT/
!     DATA IMACH( 5) /         60 /
!     DATA IMACH( 6) /         10 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         48 /
!     DATA IMACH( 9) / 00007777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /       -929 /
!     DATA IMACH(13) /       1070 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /       -929 /
!     DATA IMACH(16) /       1069 /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / Z'7FFFFFFF' /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16383 /
!     DATA IMACH(16) /      16383 /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -pd8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 9223372036854775807 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1023 /
!     DATA IMACH(13) /       1023 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 46 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         46 /
!     DATA IMACH( 9) / 1777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE CRAY
!     USING THE 64 BIT INTEGER COMPILER OPTION
!
!     DATA IMACH( 1) /        100 /
!     DATA IMACH( 2) /        101 /
!     DATA IMACH( 3) /        102 /
!     DATA IMACH( 4) /        101 /
!     DATA IMACH( 5) /         64 /
!     DATA IMACH( 6) /          8 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         63 /
!     DATA IMACH( 9) / 777777777777777777777B /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         47 /
!     DATA IMACH(12) /      -8189 /
!     DATA IMACH(13) /       8190 /
!     DATA IMACH(14) /         94 /
!     DATA IMACH(15) /      -8099 /
!     DATA IMACH(16) /       8190 /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!
!     DATA IMACH( 1) /         11 /
!     DATA IMACH( 2) /         12 /
!     DATA IMACH( 3) /          8 /
!     DATA IMACH( 4) /         10 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA IMACH( 1) /          5 /
      DATA IMACH( 2) /          6 /
      DATA IMACH( 3) /          6 /
      DATA IMACH( 4) /          6 /
      DATA IMACH( 5) /         32 /
      DATA IMACH( 6) /          4 /
      DATA IMACH( 7) /          2 /
      DATA IMACH( 8) /         31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /          2 /
      DATA IMACH(11) /         24 /
      DATA IMACH(12) /       -125 /
      DATA IMACH(13) /        128 /
      DATA IMACH(14) /         53 /
      DATA IMACH(15) /      -1021 /
      DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING D_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     USING G_FLOATING
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1023 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1022 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         24 /
!     DATA IMACH( 6) /          3 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         23 /
!     DATA IMACH( 9) /    8388607 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         38 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /         43 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          6 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         63 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         39 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION OPTION WITH FTN4
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          4 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         23 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         55 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          7 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         32 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -126 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1015 /
!     DATA IMACH(16) /       1017 /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          7 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) /  Z7FFFFFFF /
!     DATA IMACH(10) /         16 /
!     DATA IMACH(11) /          6 /
!     DATA IMACH(12) /        -64 /
!     DATA IMACH(13) /         63 /
!     DATA IMACH(14) /         14 /
!     DATA IMACH(15) /        -64 /
!     DATA IMACH(16) /         63 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          0 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         54 /
!     DATA IMACH(15) /       -101 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          5 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / "377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         62 /
!     DATA IMACH(15) /       -128 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGER ARITHMETIC.
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          5 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -125 /
!     DATA IMACH(13) /        128 /
!     DATA IMACH(14) /         53 /
!     DATA IMACH(15) /      -1021 /
!     DATA IMACH(16) /       1024 /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          6 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         32 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         31 /
!     DATA IMACH( 9) / 2147483647 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         53 /
!     DATA IMACH(12) /      -1021 /
!     DATA IMACH(13) /       1024 /
!     DATA IMACH(14) /        113 /
!     DATA IMACH(15) /     -16381 /
!     DATA IMACH(16) /      16384 /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
!
!     DATA IMACH( 1) /          5 /
!     DATA IMACH( 2) /          6 /
!     DATA IMACH( 3) /          1 /
!     DATA IMACH( 4) /          6 /
!     DATA IMACH( 5) /         36 /
!     DATA IMACH( 6) /          4 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         35 /
!     DATA IMACH( 9) / O377777777777 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         27 /
!     DATA IMACH(12) /       -128 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         60 /
!     DATA IMACH(15) /      -1024 /
!     DATA IMACH(16) /       1023 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA IMACH( 1) /          1 /
!     DATA IMACH( 2) /          1 /
!     DATA IMACH( 3) /          0 /
!     DATA IMACH( 4) /          1 /
!     DATA IMACH( 5) /         16 /
!     DATA IMACH( 6) /          2 /
!     DATA IMACH( 7) /          2 /
!     DATA IMACH( 8) /         15 /
!     DATA IMACH( 9) /      32767 /
!     DATA IMACH(10) /          2 /
!     DATA IMACH(11) /         24 /
!     DATA IMACH(12) /       -127 /
!     DATA IMACH(13) /        127 /
!     DATA IMACH(14) /         56 /
!     DATA IMACH(15) /       -127 /
!     DATA IMACH(16) /        127 /
!
!***FIRST EXECUTABLE STATEMENT  I1MACH
!
  if ( I < 1 .OR. I > 16 ) then
    WRITE (UNIT = OUTPUT, FMT = 9000)
 9000 FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
    STOP
  end if

  I1MACH = IMACH(I)

  return
end

FUNCTION R1MACH (I)
!
!! R1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   R1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = R1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   R1MACH(3) = B**(-T), the smallest relative spacing.
!   R1MACH(4) = B**(1-T), the largest relative spacing.
!   R1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!***END PROLOGUE  R1MACH
!
  real r1mach
  INTEGER SMALL(2)
  INTEGER LARGE(2)
  INTEGER RIGHT(2)
  INTEGER DIVER(2)
  INTEGER LOG10(2)
!
  REAL RMACH(5)
  SAVE RMACH
!
  EQUIVALENCE (RMACH(1),SMALL(1))
  EQUIVALENCE (RMACH(2),LARGE(1))
  EQUIVALENCE (RMACH(3),RIGHT(1))
  EQUIVALENCE (RMACH(4),DIVER(1))
  EQUIVALENCE (RMACH(5),LOG10(1))
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE AMIGA
!     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7EFFFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE APOLLO
!
!     DATA SMALL(1) / 16#00800000 /
!     DATA LARGE(1) / 16#7FFFFFFF /
!     DATA RIGHT(1) / 16#33800000 /
!     DATA DIVER(1) / 16#34000000 /
!     DATA LOG10(1) / 16#3E9A209B /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
!
!     DATA RMACH(1) / Z400800000 /
!     DATA RMACH(2) / Z5FFFFFFFF /
!     DATA RMACH(3) / Z4E9800000 /
!     DATA RMACH(4) / Z4EA800000 /
!     DATA RMACH(5) / Z500E730E8 /
!
!     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
!
!     DATA RMACH(1) / O1771000000000000 /
!     DATA RMACH(2) / O0777777777777777 /
!     DATA RMACH(3) / O1311000000000000 /
!     DATA RMACH(4) / O1301000000000000 /
!     DATA RMACH(5) / O1157163034761675 /
!
!     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
!
!     DATA RMACH(1) / Z"3001800000000000" /
!     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
!     DATA RMACH(3) / Z"3FD2800000000000" /
!     DATA RMACH(4) / Z"3FD3800000000000" /
!     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
!
!     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
!
!     DATA RMACH(1) / 00564000000000000000B /
!     DATA RMACH(2) / 37767777777777777776B /
!     DATA RMACH(3) / 16414000000000000000B /
!     DATA RMACH(4) / 16424000000000000000B /
!     DATA RMACH(5) / 17164642023241175720B /
!
!     MACHINE CONSTANTS FOR THE CELERITY C1260
!
!     DATA SMALL(1) / Z'00800000' /
!     DATA LARGE(1) / Z'7F7FFFFF' /
!     DATA RIGHT(1) / Z'33800000' /
!     DATA DIVER(1) / Z'34000000' /
!     DATA LOG10(1) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fn COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7FFFFFFF' /
!     DATA RMACH(3) / Z'34800000' /
!     DATA RMACH(4) / Z'35000000' /
!     DATA RMACH(5) / Z'3F9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -fi COMPILER OPTION
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE CONVEX
!     USING THE -p8 OR -pd8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CC0000000000000' /
!     DATA RMACH(4) / Z'3CD0000000000000' /
!     DATA RMACH(5) / Z'3FF34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE CRAY
!
!     DATA RMACH(1) / 200034000000000000000B /
!     DATA RMACH(2) / 577767777777777777776B /
!     DATA RMACH(3) / 377224000000000000000B /
!     DATA RMACH(4) / 377234000000000000000B /
!     DATA RMACH(5) / 377774642023241175720B /
!
!     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
!     STATIC RMACH(5)
!
!     DATA SMALL /    20K,       0 /
!     DATA LARGE / 77777K, 177777K /
!     DATA RIGHT / 35420K,       0 /
!     DATA DIVER / 36020K,       0 /
!     DATA LOG10 / 40423K,  42023K /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING G_FLOAT
!
!     DATA RMACH(1) / '00000080'X /
!     DATA RMACH(2) / 'FFFF7FFF'X /
!     DATA RMACH(3) / '00003480'X /
!     DATA RMACH(4) / '00003500'X /
!     DATA RMACH(5) / '209B3F9A'X /
!
!     MACHINE CONSTANTS FOR THE DEC ALPHA
!     USING IEEE_FLOAT
!
      DATA RMACH(1) / Z'00800000' /
      DATA RMACH(2) / Z'7F7FFFFF' /
      DATA RMACH(3) / Z'33800000' /
      DATA RMACH(4) / Z'34000000' /
      DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC RISC
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE DEC VAX
!     (EXPRESSED IN INTEGER AND HEXADECIMAL)
!     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
!     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
!
!     DATA SMALL(1) /       128 /
!     DATA LARGE(1) /    -32769 /
!     DATA RIGHT(1) /     13440 /
!     DATA DIVER(1) /     13568 /
!     DATA LOG10(1) / 547045274 /
!
!     DATA SMALL(1) / Z00000080 /
!     DATA LARGE(1) / ZFFFF7FFF /
!     DATA RIGHT(1) / Z00003480 /
!     DATA DIVER(1) / Z00003500 /
!     DATA LOG10(1) / Z209B3F9A /
!
!     MACHINE CONSTANTS FOR THE ELXSI 6400
!     (ASSUMING REAL*4 IS THE DEFAULT REAL)
!
!     DATA SMALL(1) / '00800000'X /
!     DATA LARGE(1) / '7F7FFFFF'X /
!     DATA RIGHT(1) / '33800000'X /
!     DATA DIVER(1) / '34000000'X /
!     DATA LOG10(1) / '3E9A209B'X /
!
!     MACHINE CONSTANTS FOR THE HARRIS 220
!
!     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
!     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
!     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
!     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
!     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
!
!     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
!
!     DATA RMACH(1) / O402400000000 /
!     DATA RMACH(2) / O376777777777 /
!     DATA RMACH(3) / O714400000000 /
!     DATA RMACH(4) / O716400000000 /
!     DATA RMACH(5) / O776464202324 /
!
!     MACHINE CONSTANTS FOR THE HP 730
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     3 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 2100
!     4 WORD DOUBLE PRECISION WITH FTN4
!
!     DATA SMALL(1), SMALL(2) / 40000B,       1 /
!     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
!     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
!     DATA DIVER(1), DIVER(2) / 40000B,    327B /
!     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
!
!     MACHINE CONSTANTS FOR THE HP 9000
!
!     DATA SMALL(1) / 00004000000B /
!     DATA LARGE(1) / 17677777777B /
!     DATA RIGHT(1) / 06340000000B /
!     DATA DIVER(1) / 06400000000B /
!     DATA LOG10(1) / 07646420233B /
!
!     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
!     THE PERKIN ELMER (INTERDATA) 7/32.
!
!     DATA RMACH(1) / Z00100000 /
!     DATA RMACH(2) / Z7FFFFFFF /
!     DATA RMACH(3) / Z3B100000 /
!     DATA RMACH(4) / Z3C100000 /
!     DATA RMACH(5) / Z41134413 /
!
!     MACHINE CONSTANTS FOR THE IBM PC
!
!     DATA SMALL(1) / 1.18E-38      /
!     DATA LARGE(1) / 3.40E+38      /
!     DATA RIGHT(1) / 0.595E-07     /
!     DATA DIVER(1) / 1.19E-07      /
!     DATA LOG10(1) / 0.30102999566 /
!
!     MACHINE CONSTANTS FOR THE IBM RS 6000
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE INTEL i860
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
!
!     DATA RMACH(1) / "000400000000 /
!     DATA RMACH(2) / "377777777777 /
!     DATA RMACH(3) / "146400000000 /
!     DATA RMACH(4) / "147400000000 /
!     DATA RMACH(5) / "177464202324 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1) /    8388608 /
!     DATA LARGE(1) / 2147483647 /
!     DATA RIGHT(1) /  880803840 /
!     DATA DIVER(1) /  889192448 /
!     DATA LOG10(1) / 1067065499 /
!
!     DATA RMACH(1) / O00040000000 /
!     DATA RMACH(2) / O17777777777 /
!     DATA RMACH(3) / O06440000000 /
!     DATA RMACH(4) / O06500000000 /
!     DATA RMACH(5) / O07746420233 /
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
!     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
!
!     DATA SMALL(1), SMALL(2) /   128,     0 /
!     DATA LARGE(1), LARGE(2) / 32767,    -1 /
!     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
!     DATA DIVER(1), DIVER(2) / 13568,     0 /
!     DATA LOG10(1), LOG10(2) / 16282,  8347 /
!
!     DATA SMALL(1), SMALL(2) / O000200, O000000 /
!     DATA LARGE(1), LARGE(2) / O077777, O177777 /
!     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
!     DATA DIVER(1), DIVER(2) / O032400, O000000 /
!     DATA LOG10(1), LOG10(2) / O037632, O020233 /
!
!     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!
!     DATA RMACH(1) / Z'00800000' /
!     DATA RMACH(2) / Z'7F7FFFFF' /
!     DATA RMACH(3) / Z'33800000' /
!     DATA RMACH(4) / Z'34000000' /
!     DATA RMACH(5) / Z'3E9A209B' /
!
!     MACHINE CONSTANTS FOR THE SUN
!     USING THE -r8 COMPILER OPTION
!
!     DATA RMACH(1) / Z'0010000000000000' /
!     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
!     DATA RMACH(3) / Z'3CA0000000000000' /
!     DATA RMACH(4) / Z'3CB0000000000000' /
!     DATA RMACH(5) / Z'3FD34413509F79FF' /
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
!
!     DATA RMACH(1) / O000400000000 /
!     DATA RMACH(2) / O377777777777 /
!     DATA RMACH(3) / O146400000000 /
!     DATA RMACH(4) / O147400000000 /
!     DATA RMACH(5) / O177464202324 /
!
!     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
!
!     DATA SMALL(1), SMALL(2) /     0,    256/
!     DATA LARGE(1), LARGE(2) /    -1,   -129/
!     DATA RIGHT(1), RIGHT(2) /     0,  26880/
!     DATA DIVER(1), DIVER(2) /     0,  27136/
!     DATA LOG10(1), LOG10(2) /  8347,  32538/
!
!***FIRST EXECUTABLE STATEMENT  R1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'R1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

  R1MACH = RMACH(I)

  return
end