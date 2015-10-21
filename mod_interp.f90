    module mod_interp
    implicit none
! This module defines standalone functions and subroutines useful for interpolation of functions, 
! both in one and two dimensions.

    contains
! ------------------------------------------------------------------------------------------------------------------------
    real(8) function neville(points,N,z)
    implicit none
! The Neville algorithm is used for interpolation: the array points holds N sets of (x,f(x)) values,
! and the value of the function at z is required.
     integer N,i,j,near
     real(8) points(2,N),z
     real(8) diff1,diff2
     real(8) C(N),D(N)
     real(8) ho,hp,den,diff_est,w

     near = 1

     diff1 = dabs(z - points(1,1))
! The index, near, of the nearest table entry is found:
     do i=1,N
      diff2 = dabs(z - points(1,i))
      if (diff2 .lt. diff1) then
       near = i
       diff1 = diff2
      endif
      C(i) = points(2,i)
      D(i) = points(2,i)
     enddo

! The initial approximation is found:
     neville = points(2,near)
     near = near - 1
! Loop over current Cs and Ds, and update them:
     do i=1,N-1
      do j=1,N-i
       ho = points(1,j) - z
       hp = points(1,j+i) - z
       w = C(j+1) - D(j)
       den = ho - hp
       if (den .eq. 0.d0) then
        write(*,*) 'Two input values of z are identical!'
        stop
       endif
       den = w / den
       D(j) = hp * den
       C(j) = ho * den
      enddo
      if (2*near .lt. N-i) then
       diff_est = C(near + 1)
      else
       diff_est = D(near)
       near = near - 1
      endif
      neville = neville + diff_est
     enddo

    end function
! ------------------------------------------------------------------------------------------------------------------------
    subroutine bin_interp(surf_points,dfdx,dfdy,d2fdxdy,gspace,point,ans,anszy,anszx)
    implicit none

! This subroutine carries out bicubic interpolation.

    real(8) gspace(2),surf_points(4),dfdx(4),dfdy(4),d2fdxdy(4),coeffs(4,4),point(2)
    real(8) tval,uval,ans,anszx,anszy
    integer i

    call bin_cubic_coeffs(surf_points,dfdx,dfdy,d2fdxdy,gspace,coeffs)

    tval = point(2)/gspace(2)
    uval = point(1)/gspace(1)

!	if (tval .gt. 1.d0) then
!		write(*,*) tval
!		write(*,*) point
!		write(*,*) gspace
!		pause
!	endif

    ans = 0.d0
    anszx = 0.d0
    anszy = 0.d0
    do i=4,1,-1
     ans = tval*ans + ((coeffs(i,4)*uval+coeffs(i,3))*uval+coeffs(i,2))*uval+coeffs(i,1)
     anszx = tval*anszx + (3.d0*coeffs(i,4)*uval + 2.d0*coeffs(i,3))*uval + coeffs(i,2)
     anszy = uval*anszy + (3.d0*coeffs(4,i)*tval + 2.d0*coeffs(3,i))*tval + coeffs(2,i)
    enddo

! Find the unit vector for the surface normal at this point:
    anszy = anszy / gspace(2)
    anszx = anszx / gspace(1)

    end subroutine
! ------------------------------------------------------------------------------------------------------------------------
    subroutine bin_cubic_coeffs(surf_points,dfdx,dfdy,d2fdxdy,gspace,coeffs)
    implicit none

! This subroutine finds the coefficients used for bicubic interpolation.

    real(8) gspace(2),surf_points(4),dfdx(4),dfdy(4),d2fdxdy(4),coeffs(4,4)
    real(8) weights(16,16),tempvec(16)
    integer i

    data weights/1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4 &
     &     ,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4 &
     &     ,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,2 &
     &     ,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2 &
     &     ,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2 &
     &     ,10*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2 &
     &     ,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,-1,0,1,-2,1 &
     &     ,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/

    do i=1,4
     tempvec(i) = surf_points(i)
     tempvec(i+4) = dfdy(i) * gspace(2)
     tempvec(i+8) = dfdx(i) * gspace(1)
     tempvec(i+12) = d2fdxdy(i) * gspace(1) * gspace(2)
    enddo
    tempvec = matmul(weights,tempvec)
    coeffs = reshape(tempvec,(/4,4/),order=(/2,1/))

    end subroutine
! ------------------------------------------------------------------------------------------------------------------------
    real(8) function brent(z1,z2,tol,x,y,fxyz)
    implicit none

    real(8) z1,z2,tol,x,y,fxyz
    integer i,maxit
    real(8) za,zb,zc,zd,ze,funca,funcb,funcc,tol1,xm,p,q,r,s
    parameter(maxit = 100)

    za = z1
    zb = z2
    funca = fxyz(x,y,za)
    funcb = fxyz(x,y,zb)
    if (funca * funcb .gt. 0.d0) then
     write(*,*) 'Root is not bracketed!'
     stop
    endif
    zc = zb
    funcc = funcb
    do i=1,maxit
     if (funcb * funcc .gt. 0.d0) then
      zc = za
      funcc = funca
      zd = zb - za
      ze = zd
     endif
     if (dabs(funcc) .lt. dabs(funcb)) then
      za = zb
      zb = zc
      zc = za
      funca = funcb
      funcb = funcc
      funcc = funca
     endif
     tol1 = 2.d0 * 3.d-8 * abs(zb) + 0.5d0 * tol
     xm = 0.5d0 * (zc - zb)
     if ((dabs(xm) .le. tol1) .or. (dabs(funcb) .eq. 0.d0)) then
      brent = zb
      return
     endif
     if ((dabs(ze) .ge. tol1) .and. (dabs(funca) .gt. dabs(funcb))) then
      s = funcb / funca
      if (za .eq. zc) then
       p = 2.d0 * xm * s
       q = 1.d0 - s
      else
       q = funca / funcc
       r = funcb / funcc
       p = s*(2.d0 * xm * q * (q-r)-(zb-za)*(r-1.d0))
       q = (q-1.d0) * (r-1.d0) * (s-1.d0)
      endif
      if (p .gt. 0.d0) q = -q
      p = dabs(p)
      if (2.d0 * p .lt. min(3.d0*xm*q-dabs(tol1*q),abs(ze*q))) then
       ze = zd
       zd = p / q
      else
       zd = xm
       ze = zd
      endif
     else
      zd = xm
      ze = zd
     endif
     za = zb
     funca = funcb
     if (dabs(zd) .gt. tol1) then
      zb = zb + zd
     else
      zb = zb + sign(tol1,xm)
     endif
    funcb = fxyz(x,y,zb)
    enddo

    write(*,*) 'Brent algorithm has exceeded maximum number of iterations!'
    read(*,*)
    brent = zb
    return

    end function
! ------------------------------------------------------------------------------------------------------------------------
    subroutine bracket(x,y,z1,z2,fxyz,factor)
    implicit none

    integer i
    real(8) x,y,z1,z2,factor,func1,func2,fxyz
    integer numtry
    parameter(numtry = 20)

! Expand the range (z1,z2) until a root of density_field is bracketed by them.
    if (z1 .eq. z2) then
     write(*,*) 'Bracketing algorithm requires a nonzero initial range'
     stop
    endif
    func1 = fxyz(x,y,z1)
    func2 = fxyz(x,y,z2)
    do i=1,numtry
     if (func1 * func2 .lt. 0.d0) return
     if (dabs(func1) .lt. dabs(func2)) then
      z1 = z1 + factor*(z1-z2)
      func1 = fxyz(x,y,z1)
     else
      z2 = z2 + factor*(z2-z1)
      func2 = fxyz(x,y,z2)
     endif
    enddo
    write(*,*) 'Root cannot be bracketed - consider changing numtry parameter.'
    read(*,*)
    end subroutine
! ------------------------------------------------------------------------------------------------------------------------
    end module
