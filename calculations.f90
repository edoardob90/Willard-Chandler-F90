! .........................................................................................................DMW.
! ...................................... SUBROUTINES ..........................................................
! .............................................................................................................
! --------------------------------------------------------------------------------------------------------------

    subroutine find_surface()
    use mod_surf
    use mod_interp
    implicit none

    integer u,v,w,x,y
    logical allfound,found(coord(1),coord(2),2)
    real(DP) lowbound,upbound
    integer xmin,xmax,ymin,ymax

! Firstly, take the centre of the box (box_length(1)/2,box_length(2)/2,0), and use Brent's algorithm
! to find the root of density_field(X,Y,z) - const = 0.d0 for this X and Y - thus find both surfaces.

      found = .false.
      surf = 0.d0
      do x=1,coord(1)
       do y=1,coord(2)
        surf(x,y,:,1) = (x-1)*gspacing(1)
        surf(x,y,:,2) = (y-1)*gspacing(2)
       enddo
      enddo

      x = nint(box_length(1)/2)
      y = nint(box_length(2)/2)
      
      !print *, "Line 31, calculations.f90"
      
      surf(x,y,1,3) = brent( 0.d0, box_length(3)*0.5d0, 1.d-6, surf(x,y,1,1), surf(x,y,1,2), density_field_const )
      surf(x,y,2,3) = brent(box_length(3)*0.5d0, box_length(3), 1.d-9, surf(x,y,2,1), surf(x,y,2,2), density_field_const)
      found(x,y,:) = .true.
      
      !print *, "Line 37, calculations.f90"

! Then, spreading out from this point, find the surface at different (X,Y) values.

      xmin = x
      xmax = x
      ymin = y
      ymax = y
      allfound = .false.

      do while (.not. allfound)
       do w=1,2
        do u = xmin-1,xmax+1
         do v=ymin-1,ymax+1
          if ((u .ge. 1) .and. (u .le. coord(1)) .and. (v .ge. 1) .and. (v .le. coord(2))) then
           if (.not. found(u,v,w)) then

! Find the lower and upper bounds for the z-value:
            lowbound = box_length(3)
            upbound = 0.d0
            do x=u-1,u+1
             do y=v-1,v+1
              if ((x .ge. 1) .and. (x .le. coord(1)) .and. (y .ge. 1) .and. (y .le. coord(2))) then
               if (found(x,y,w))  then
                lowbound = min(lowbound,surf(x,y,w,3))
                upbound = max(upbound,surf(x,y,w,3))
               endif
              endif
             enddo
            enddo
            lowbound = lowbound - gspacing(3)
            upbound = upbound + gspacing(3)
            call bracket(surf(u,v,w,1),surf(u,v,w,2),lowbound,upbound,density_field_const,1.2d0)

! Use Brent's algorithm to find the z-value:
            surf(u,v,w,3) = brent(lowbound,upbound,1.d-9,surf(u,v,w,1),surf(u,v,w,2),density_field_const)
            found(u,v,w) = .true.

           endif
          endif
         enddo
        enddo
       enddo

       xmin = xmin - 1
       xmax = xmax + 1
       ymin = ymin - 1
       ymax = ymax + 1

! Check to see if all of the surface has been found:
       allfound = .true.
       do x=1,coord(1)
        do y=1,coord(2)
         allfound = allfound .and. found(x,y,1) .and. found(x,y,2)
        enddo
       enddo
      enddo

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine int_surface(N,x,y,surf2point,xypoint,find_grad)
    use mod_surf
    use mod_interp
    implicit none

    integer x,y,N
    real(DP) surf_points(4)
    real(DP) xypoint(2),find_grad(3)
    real(DP) dfdx(4),dfdy(4),d2fdxdy(4),anszy,anszx
    real(DP) surf2point

! This subroutine finds the z-value of the surface between 4 (x,y) grid points.
! With all of the relevant data known at the four surrounding grid points, the value of z at which the surface occurs
! is found by a bicubic interpolation:
    surf_points(1) = surf(x,y,N,3)
    surf_points(2) = surf(x,y+1,N,3)
    surf_points(3) = surf(x+1,y+1,N,3)
    surf_points(4) = surf(x+1,y,N,3)
    dfdx(1) = -1.d0 * gradient(x,y,N,1) / gradient(x,y,N,3)
    dfdx(2) = -1.d0 * gradient(x,y+1,N,1) / gradient(x,y+1,N,3)
    dfdx(3) = -1.d0 * gradient(x+1,y+1,N,1) / gradient(x+1,y+1,N,3)
    dfdx(4) = -1.d0 * gradient(x+1,y,N,1) / gradient(x+1,y,N,3)
    dfdy(1) = -1.d0 * gradient(x,y,N,2) / gradient(x,y,N,3)
    dfdy(2) = -1.d0 * gradient(x,y+1,N,2) / gradient(x,y+1,N,3)
    dfdy(3) = -1.d0 * gradient(x+1,y+1,N,2) / gradient(x+1,y+1,N,3)
    dfdy(4) = -1.d0 * gradient(x+1,y,N,2) / gradient(x+1,y,N,3)
    d2fdxdy(1)   = mixed(x,y,N)
    d2fdxdy(2)   = mixed(x,y+1,N)
    d2fdxdy(3)   = mixed(x+1,y+1,N)
    d2fdxdy(4)   = mixed(x+1,y,N)
    call bin_interp(surf_points,dfdx,dfdy,d2fdxdy,gspacing(1:2),xypoint,surf2point,anszy,anszx)

    find_grad(1) = -1.d0 * anszx
    find_grad(2) = -1.d0 * anszy
    find_grad(3) = 1.d0
    find_grad = find_grad / dsqrt(find_grad(1)**2 + find_grad(2)**2 + find_grad(3)**2)

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine calculate_grid()
    use mod_surf
    implicit none

    real(DP) :: grid_spacing
    integer :: i

! Calculate grid

     !read(input,*) (box_length(x),x=1,3,1)

     coord = 0
     grid_spacing = 0.5d0
     do i=1,3
      coord(i) = 1 + nint(box_length(i) / grid_spacing)
      gspacing(i) = box_length(i)/(coord(i) - 1)
     enddo

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine find_terms_int()
    use mod_surf
    implicit none

! Find the gradient and the mixed terms, for interpolation:
    integer v,x,y
    real(DP) dens_field

      do x=1,coord(1)
       do y=1,coord(2)
        do v=1,2
         dens_field = cal_density_field(surf(x,y,v,1),surf(x,y,v,2),surf(x,y,v,3))
         call cal_gradient(gradient(x,y,v,:),surf(x,y,v,1),surf(x,y,v,2),surf(x,y,v,3))
         mixed(x,y,v) = cal_mixed_term(surf(x,y,v,1),surf(x,y,v,2),surf(x,y,v,3), &
     &          dens_field,gradient(x,y,v,:))
        enddo
       enddo
      enddo

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine triangles()
    use mod_surf
    implicit none
! Interpolate between grid points, so that the surface can be drawn with triangles as a smooth function.

    real(DP) :: xypoint(2), grad_int(3)
    integer :: x,y,j

! Draw the surface at the grid points:
10 format('Cl',3x,f20.10,3x,f20.10,3x,f20.10)

     do x=1,coord(1)
      do y=1,coord(2)
        if (.not. ((surf(x,y,1,1) .eq. 0.d0) .and. (surf(x,y,1,2) .eq. 0.d0) &
     &      .and. (surf(x,y,1,3) .eq. 0.d0))) then
              write(f_sur,10) ( surf(x,y,1,j), j=1,3 )
        endif
        if (.not. ((surf(x,y,2,1) .eq. 0.d0) .and. (surf(x,y,2,2) .eq. 0.d0) &
     &      .and. (surf(x,y,2,3) .eq. 0.d0))) then
              write(f_sur,10) (surf(x,y,2,j), j=1,2), surf(x,y,2,3) - box_length(3)
        endif
      enddo
     enddo

    write(*,*) gspacing

! Find and draw the surface 'in-between' the grid points:
     xypoint(1) = 0.5d0 * gspacing(1)
     xypoint(2) = 0.5d0 * gspacing(2)
     do x=1,coord(1)-1
      do y=1,coord(2)-1
        surf2(x,y,1,1) = (x-1)*gspacing(1) + 0.5d0*gspacing(1)
        surf2(x,y,1,2) = (y-1)*gspacing(2) + 0.5d0*gspacing(2)
        call int_surface(1,x,y,surf2(x,y,1,3),xypoint,grad_int)
        surf2(x,y,2,1) = (x-1)*gspacing(1) + 0.5d0*gspacing(1)
        surf2(x,y,2,2) = (y-1)*gspacing(2) + 0.5d0*gspacing(2)
        call int_surface(2,x,y,surf2(x,y,2,3),xypoint,grad_int)

        write(f_sur,10) ( surf2(x,y,1,j), j=1,3 )
        write(f_sur,10) ( surf2(x,y,2,j), j=1,2 ), surf2(x,y,2,3) - box_length(3)

      enddo
     enddo

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
!!!    subroutine find_atom_height_grad(w)
!!!    use mod_surf
!!!    implicit none
!!!
!!!    integer w,gridx,gridy
!!!    real(DP) xyzpoint(2),grad_interp1(3),grad_interp2(3)
!!!
!!!! What is the atom's height compared to the two surfaces? Folding the atom back into the box allows
!!!! this comparison to be carried out:
!!!
!!!    xscratch = atom_position(w,1)
!!!    do while (xscratch .lt. 0.d0 .or. xscratch .gt. box_length(1))
!!!     xscratch = xscratch - sign(box_length(1),xscratch)
!!!    enddo
!!!    yscratch = atom_position(w,2)
!!!    do while (yscratch .lt. 0.d0 .or. yscratch .gt. box_length(2))
!!!     yscratch = yscratch - sign(box_length(2),yscratch)
!!!    enddo
!!!
!!!! Find the (x,y) grid point such that (x,y),(x,y+1),(x+1,y),(x+1,y+1) bound this atom's (x,y) position:
!!!
!!!    gridx = floor(xscratch / gspacing(1))
!!!    gridy = floor(yscratch / gspacing(2))
!!!    xyzpoint(1) = (xscratch/gspacing(1)) - gridx
!!!    xyzpoint(2) = (yscratch/gspacing(2)) - gridy
!!!    xyzpoint(1) = xyzpoint(1) * gspacing(1)
!!!    xyzpoint(2) = xyzpoint(2) * gspacing(2)
!!!
!!!! Find the z positions of the two surfaces for this pair of (x,y) values, and thus the height
!!!! of the atom below both the upper and lower surface:
!!!
!!!    call int_surface(1,gridx+1,gridy+1,zheight(w,1),xyzpoint,grad_interp1)
!!!    call int_surface(2,gridx+1,gridy+1,zheight(w,2),xyzpoint,grad_interp2)
!!!    zheight(w,1) = zheight(w,1) - atom_position(w,3)
!!!    zheight(w,2) = atom_position(w,3) - zheight(w,2)
!!!!    if (dabs(zheight(w,1)) .gt. (0.5d0*box_length(3))) zheight(w,1) = zheight(w,1) - dsign(box_length(3),zheight(w,1))
!!!!    if (dabs(zheight(w,2)) .gt. (0.5d0*box_length(3))) zheight(w,2) = zheight(w,2) - dsign(box_length(3),zheight(w,2))
!!!    do while (dabs(zheight(w,1)) .gt. (0.5d0*box_length(3)))
!!!     zheight(w,1) = zheight(w,1) - dsign(box_length(3),zheight(w,1))
!!!    enddo
!!!    do while (dabs(zheight(w,2)) .gt. (0.5d0*box_length(3)))
!!!     zheight(w,2) = zheight(w,2) - dsign(box_length(3),zheight(w,2))
!!!    enddo
!!!
!!!! Find the unit vector normal to the surface at the (x,y) position of each atom.
!!!
!!!    if (zheight(w,1) .lt. zheight(w,2)) then
!!!     grad_interp = -1.d0 * grad_interp1
!!!    else if (zheight(w,2) .lt. zheight(w,1)) then
!!!     grad_interp = grad_interp2
!!!    else
!!!     write(*,*) 'Molecule is equidistant from the two surfaces -- more programming is needed!'
!!!     stop
!!!    endif
!!!
!!!    end subroutine
