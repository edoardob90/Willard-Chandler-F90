    module mod_surf
    implicit none

! Define global variables:
! ========================
! Constants and input parameters:
    real(8) xi,ximin2,ximin4,ninxisq,phi3xi,const,prefac
    real(8) :: pi=3.1415926535897931d0
    real(8) box_length(3),gspacing(3)
    real(8) xscratch,yscratch,zscratch
    integer coord(3)

! Surface variables:
    real(8), allocatable :: surf2(:,:,:,:),surf(:,:,:,:)
    real(8), allocatable :: gradient(:,:,:,:), mixed(:,:,:)
    real(8), allocatable :: zheight(:,:)
    real(8) grad_interp(3)
    real(8) vec1(3),vec2(3),vec3(3)

! Atoms, their names and positions:
    integer num_frames,num_atom
    character(len=3), allocatable :: atom_name(:)
    real(8), allocatable :: atom_position(:,:)
    integer, allocatable :: oxygen_list(:)
    real(8) unit_cell(3,3),unit_cell_inv(3,3)

! Input and output:
    character(len=32) dname,file_water,file_surface,file_dist
    integer dat,f_wat,f_sur,f_dist
    character(len=100) input,input1

! Functions and subroutines:
    contains
! ------------------------------------------------------------------------------------------------------------------------
    subroutine init_consts()
    implicit none

     xi = 2.4d0
     ximin2 = xi**(-2.d0)
     ximin4 = xi**(-4.d0)
     ninxisq = 9.d0 * xi**2
     const = 0.016d0
     prefac = ((2 * pi * xi**2)**(-1.5d0))
     phi3xi = ( prefac * dexp(-4.5d0))

     end subroutine
! ------------------------------------------------------------------------------------------------------------------------
     real(8) function phi(rsq)
     implicit none

     real(8) rsq

     if (rsq .ge. ninxisq) then
      phi = 0.d0
      return
     endif

     phi = ( prefac * dexp(-0.5d0 * rsq / xi**2 )) - phi3xi

     end function
! ------------------------------------------------------------------------------------------------------------------------
    real(8) function density_field_const(x,y,z)
    implicit none

    real(8) x,y,z

    density_field_const = cal_density_field(x,y,z) - const

    end function
! ------------------------------------------------------------------------------------------------------------------------
    real(8) function cal_density_field(x,y,z)
    implicit none
! Calculate the density field at the position (x,y,z)

    integer w
    real(8) x,y,z,r(3),rsq

    cal_density_field = 0.d0

    do w=1,num_atom/3
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
         r(1) = x - atom_position(oxygen_list(w),1)
         r(2) = y - atom_position(oxygen_list(w),2)
         r(3) = z - atom_position(oxygen_list(w),3)
         do while (dabs(r(1)) .gt. 0.5d0*box_length(1))
          r(1) = r(1) - dsign(box_length(1),r(1))
         enddo
         do while (dabs(r(2)) .gt. 0.5d0*box_length(2))
          r(2) = r(2) - dsign(box_length(2),r(2))
         enddo
         do while (dabs(r(3)) .gt. 0.5d0*box_length(3))
          r(3) = r(3) - dsign(box_length(3),r(3))
         enddo

! Calculate the density field:
         rsq = r(1)**2 + r(2)**2 + r(3)**2
!         rsq = mic(r)
         cal_density_field = cal_density_field + phi(rsq)
    enddo

    end function
! ------------------------------------------------------------------------------------------------------------------------
    subroutine cal_gradient(grad,x,y,z)
    implicit none

    integer w
    real(8) grad(3),x,y,z,r(3),rsq,expterm

    grad(:) = 0.d0

    do w=1,num_atom/3
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
         r(1) = x - atom_position(oxygen_list(w),1)
         r(2) = y - atom_position(oxygen_list(w),2)
         r(3) = z - atom_position(oxygen_list(w),3)
         do while (dabs(r(1)) .gt. 0.5d0*box_length(1))
          r(1) = r(1) - dsign(box_length(1),r(1))
         enddo
         do while (dabs(r(2)) .gt. 0.5d0*box_length(2))
          r(2) = r(2) - dsign(box_length(2),r(2))
         enddo
         do while (dabs(r(3)) .gt. 0.5d0*box_length(3))
          r(3) = r(3) - dsign(box_length(3),r(3))
         enddo

! Calculate the gradient:
         rsq = r(1)**2 + r(2)**2 + r(3)**2
!         rsq = mic(r)
         expterm = phi(rsq)
         grad(1) = grad(1) - ( ximin2 * r(1) * expterm)
         grad(2) = grad(2) - ( ximin2 * r(2) * expterm)
         grad(3) = grad(3) - ( ximin2 * r(3) * expterm)
    enddo

    end subroutine
! ------------------------------------------------------------------------------------------------------------------------
    real(8) function cal_mixed_term(x,y,z,density_field,grad)
    implicit none

    integer w
    real(8) x,y,z,r(3),rsq
    real(8) mixed_xy,mixed_xz,mixed_yz,mixed_zz,density_field,grad(3),expterm

    cal_mixed_term = 0.d0
    mixed_xy = 0.d0
    mixed_xz = 0.d0
    mixed_yz = 0.d0
    mixed_zz = 0.d0

    do w=1,num_atom/3
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
         r(1) = x - atom_position(oxygen_list(w),1)
         r(2) = y - atom_position(oxygen_list(w),2)
         r(3) = z - atom_position(oxygen_list(w),3)
         do while (dabs(r(1)) .gt. 0.5d0*box_length(1))
          r(1) = r(1) - dsign(box_length(1),r(1))
         enddo
         do while (dabs(r(2)) .gt. 0.5d0*box_length(2))
          r(2) = r(2) - dsign(box_length(2),r(2))
         enddo
         do while (dabs(r(3)) .gt. 0.5d0*box_length(3))
          r(3) = r(3) - dsign(box_length(3),r(3))
         enddo

! Calculate the mixed term:
         rsq = r(1)**2 + r(2)**2 + r(3)**2
!         rsq = mic(r)
         expterm = phi(rsq)
         mixed_xy = mixed_xy - ( ximin4 * r(1) * r(2) * expterm)
         mixed_xz = mixed_xz + ( ximin4 * r(1) * r(3) * expterm)
         mixed_yz = mixed_yz + ( ximin4 * r(2) * r(3) * expterm)
         mixed_zz = mixed_zz - ( ximin4 * r(3) * r(3) * expterm)
    enddo

    cal_mixed_term = (mixed_xy / grad(3)) + (mixed_yz * grad(1) / (grad(3)**2)) + &
     &     (mixed_xz * grad(2) / (grad(3)**2)) + &
     &     (ximin2 * grad(1) * grad(2) * density_field / (grad(3)**3)) + &
     &     (mixed_zz * grad(1) * grad(2) / (grad(3)**3))

    end function
! ------------------------------------------------------------------------------------------------------------------------
!    real(8) function mic(r)
    subroutine mic(r)
    implicit none

    real(8) r(3)

    do while (dabs(r(1)) .gt. 0.5d0*box_length(1))
     r(1) = r(1) - dsign(box_length(1),r(1))
    enddo
    do while (dabs(r(2)) .gt. 0.5d0*box_length(2))
     r(2) = r(2) - dsign(box_length(2),r(2))
    enddo
    do while (dabs(r(3)) .gt. 0.5d0*box_length(3))
     r(3) = r(3) - dsign(box_length(3),r(3))
    enddo

    end subroutine

!    mic = r(1)**2 + r(2)**2 + r(3)**2

!    end function
! ------------------------------------------------------------------------------------------------------------------------
    end module
