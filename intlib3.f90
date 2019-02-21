
subroutine simp ( func, a, b, eps, result )
!
!***********************************************************************
!
!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    J N Lyness,
!    Algorithm 379,
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the Association for Computing Machinery,
!    Volume 13 (1970), pages 260-263.
!
!    W M McKeeman and L Tesler,
!    Algorithm 182,
!    Nonrecursive adaptive integration,
!    Communications of the Association for Computing Machinery,
!    Volume 6 (1963), page 315.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!    FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real A, the lower limit of integration.
!
!    Input, real B, the upper limit of integration.
!
!    Input, real EPS, the requested error tolerance.
!
!    Output, real RESULT, the approximation to the integral.
!
  implicit none
!
  integer(kind=4), parameter :: maxlev = 30
!
  real (kind=8) a
  real (kind=8) a1
  real (kind=8) absar
  real (kind=8) b
  real (kind=8) da
  real (kind=8) dx(maxlev)
  real (kind=8) ep
  real (kind=8) eps
  real (kind=8) epsp(maxlev)
  real (kind=8) est
  real (kind=8) est1
  real (kind=8) est2(maxlev)
  real (kind=8) est3(maxlev)
  real (kind=8) f1
  real (kind=8) f2(maxlev)
  real (kind=8) f3(maxlev)
  real (kind=8) f4(maxlev)
  real (kind=8) fa
  real (kind=8) fb
  real (kind=8) fbp(maxlev)
  real (kind=8) fm
  real (kind=8) fmp(maxlev)
  real (kind=8), external :: func
  integer (kind=4) i
  integer (kind=4) j
  integer (kind=4) l
  integer (kind=4) lvl
  integer (kind=4)nrtr(maxlev)
  real (kind=8) pval(maxlev,3)
  real (kind=8) result
  real (kind=8) sum1
  real (kind=8) sx
  real (kind=8) x2(maxlev)
  real (kind=8) x3(maxlev)
!
  result = 0.0E+00
  if ( a == b ) then
    return
  end if
 
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0E+00
 
  lvl = 0
  absar = 0.0E+00
  est = 0.0E+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0E+00 * func ( (a1+b) * 0.5E+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0E+00
  sx = dx(lvl)/6.0E+00
  f1 = 4.0E+00 * func(0.5*dx(lvl)+a1)
  x2(lvl) = a1+dx(lvl)
  f2(lvl) = func(x2(lvl))
  x3(lvl) = x2(lvl)+dx(lvl)
  f3(lvl) = func(x3(lvl))
  epsp(lvl) = ep
  f4(lvl) = 4.0E+00 * func(dx(lvl)*0.5E+00+x3(lvl))
  fmp(lvl) = fm
  est1 = sx*(fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
  sum1 = est1+est2(lvl)+est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )
  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
  if ( lvl >= maxlev ) go to 50
!
!  2 = UP
!
40 continue
 
  if ( lvl > 1 ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) go to 60
  if ( l == 2 ) go to 70
  if ( l == 3 ) go to 80
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
80 continue

  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)

  if ( lvl > 1 ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end

subroutine simp2 ( func, a, b, eps, result )
!
!***********************************************************************
!
!! SIMP approximates the integral of a function using an adaptive Simpson's rule.
!
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    J N Lyness,
!    Algorithm 379,
!    SQUANK - Simpson quadrature used adaptively, noise killed,
!    Communications of the Association for Computing Machinery,
!    Volume 13 (1970), pages 260-263.
!
!    W M McKeeman and L Tesler,
!    Algorithm 182,
!    Nonrecursive adaptive integration,
!    Communications of the Association for Computing Machinery,
!    Volume 6 (1963), page 315.
!
!  Modified:
!
!    30 October 2000
!
!  Parameters:
!
!    Input, real, external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!
!    FUNCTION FUNC(X)
!
!    which evaluates the function at the point X.
!
!    Input, real A, the lower limit of integration.
!
!    Input, real B, the upper limit of integration.
!
!    Input, real EPS, the requested error tolerance.
!
!    Output, real RESULT, the approximation to the integral.
!
implicit none
!
  integer(kind=4), parameter :: maxlev = 30
!
  real (kind=8) a
  real (kind=8) a1
  real (kind=8) absar
  real (kind=8) b
  real (kind=8) da
  real (kind=8) dx(maxlev)
  real (kind=8) ep
  real (kind=8) eps
  real (kind=8) epsp(maxlev)
  real (kind=8) est
  real (kind=8) est1
  real (kind=8) est2(maxlev)
  real (kind=8) est3(maxlev)
  real (kind=8) f1
  real (kind=8) f2(maxlev)
  real (kind=8) f3(maxlev)
  real (kind=8) f4(maxlev)
  real (kind=8) fa
  real (kind=8) fb
  real (kind=8) fbp(maxlev)
  real (kind=8) fm
  real (kind=8) fmp(maxlev)
  real (kind=8), external :: func
  integer (kind=4) i
  integer (kind=4) j
  integer (kind=4) l
  integer (kind=4) lvl
  integer (kind=4)nrtr(maxlev)
  real (kind=8) pval(maxlev,3)
  real (kind=8) result
  real (kind=8) sum1
  real (kind=8) sx
  real (kind=8) x2(maxlev)
  real (kind=8) x3(maxlev)
!
  result = 0.0E+00
  if ( a == b ) then
    return
  end if
 !print *, "limits inside intlib", a, b
  ep = eps
  a1 = a
  nrtr(1:maxlev) = 0
  pval(1:maxlev,1:3) = 0.0E+00
 
  lvl = 0
  absar = 0.0E+00
  est = 0.0E+00
  da = b - a1

  fa = func ( a1 )
  fm = 4.0E+00 * func ( (a1+b) * 0.5E+00 )
  fb = func ( b )
!
!  1 = RECUR
!
   30 continue
 
  lvl = lvl + 1
  dx(lvl) = da / 3.0E+00
  sx = dx(lvl)/6.0E+00
  f1 = 4.0E+00 * func(0.5*dx(lvl)+a1)
  x2(lvl) = a1+dx(lvl)
  f2(lvl) = func(x2(lvl))
  x3(lvl) = x2(lvl)+dx(lvl)
  f3(lvl) = func(x3(lvl))
  epsp(lvl) = ep
  f4(lvl) = 4.0E+00 * func(dx(lvl)*0.5E+00+x3(lvl))
  fmp(lvl) = fm
  est1 = sx*(fa+f1+f2(lvl))
  fbp(lvl) = fb
  est2(lvl) = sx * (f2(lvl)+f3(lvl)+fm)
  est3(lvl) = sx * (f3(lvl)+f4(lvl)+fb)
  sum1 = est1+est2(lvl)+est3(lvl)
  absar = absar - abs ( est ) + abs ( est1 ) + abs ( est2(lvl) ) &
    + abs ( est3(lvl) )
  if ( abs ( est - sum1 ) <= epsp(lvl) * absar ) go to 40
  if ( lvl >= maxlev ) go to 50
!
!  2 = UP
!
40 continue
 
  if ( lvl > 1 ) then
    lvl = lvl-1
  end if

  l = nrtr(lvl)

  if ( l == 0 ) then
    go to 50
  end if

  pval(lvl,l) = sum1

  if ( l == 1 ) go to 60
  if ( l == 2 ) go to 70
  if ( l == 3 ) go to 80
 
50 continue
 
  nrtr(lvl) = 1
  est = est1
  fm = f1
  fb = f2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
60 continue
 
  nrtr(lvl) = 2
  fa = f2(lvl)
  fm = fmp(lvl)
  fb = f3(lvl)
  est = est2(lvl)
  a1 = x2(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
70 continue
 
  nrtr(lvl) = 3
  fa = f3(lvl)
  fm = f4(lvl)
  fb = fbp(lvl)
  est = est3(lvl)
  a1 = x3(lvl)
  ep = epsp(lvl) / 1.7E+00
  da = dx(lvl)
  go to 30
 
80 continue

  sum1 = pval(lvl,1)+pval(lvl,2)+pval(lvl,3)

  if ( lvl > 1 ) then
    go to 40
  end if
 
90 continue
 
  result = sum1
 
  return
end


