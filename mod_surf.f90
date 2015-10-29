module mod_surf
    use kinds, only: DP
    use constants, only: pi
    implicit none
    save

! =========================
! Define global variables
! =========================
! Constants and input parameters:
    real(DP) :: xi,ximin2,ximin4,ninxisq,phi3xi,const,prefac
    real(DP), dimension(3) :: box_length, gspacing
    real(DP) :: xscratch,yscratch,zscratch
    integer :: coord(3)

! Surface variables:
    real(DP), allocatable :: surf2(:,:,:,:),surf(:,:,:,:)
    real(DP), allocatable :: gradient(:,:,:,:), mixed(:,:,:)
    real(DP), allocatable :: zheight(:,:)
    real(DP), dimension(3) :: grad_interp(3), vec1(3),vec2(3),vec3(3)

! Atoms, their names and positions:
    integer :: num_frames,num_atom

    ! ALREADY DEFINED in the new type below
    !character(len=3), allocatable :: atom_name(:)
    !real(DP), allocatable :: atom_position(:,:)
    !real(DP), allocatable :: op_values(:) ! order parameter; fourth column of trajectory file
    
    real(DP), dimension(3,3) :: unit_cell, unit_cell_inv

! Input and output:
    character(len=32) :: dname,file_water,file_surface,file_dist
    integer :: dat,f_wat,f_sur,f_dist,stride
    character(len=100) :: input,input1
    real(DP) :: opref(2) ! the two reference values for the order parameter

! NEW: the Atom type
    type Atom
        real(DP), dimension(3) :: xyz
        real(DP) :: op_value
        character(len=2) :: symbol
    end type Atom

! =============================
! Functions and subroutines
! =============================
    contains

    subroutine init_consts()

    implicit none

     !xi = 2.4d0
     ximin2 = 1.d0/xi**2.d0
     ximin4 = 1.d0/xi**4.d0
     ninxisq = 9.d0 * xi**2
     const = 0.016d0
     prefac = ((2 * pi * xi**2)**(-1.5d0))
     phi3xi = ( prefac * dexp(-4.5d0))

     end subroutine
     
! ------------------------------------------------------------------------------------------------------------------------
     
     real(DP) function phi(rsq)
     implicit none

     real(DP) :: rsq

     if (rsq .ge. ninxisq) then
      phi = 0.d0
      return
     endif

     phi = ( prefac * dexp(-0.5d0 * rsq / xi**2 )) - phi3xi

     return

     end function
    
     
! ------------------------------------------------------------------------------------------------------------------------
     
     
    real(DP) function density_field_const(x,y,z)
    implicit none

    real(DP) :: x,y,z

    density_field_const = cal_density_field(x,y,z) - const

    return

    end function

    
! ------------------------------------------------------------------------------------------------------------------------
    

    real(DP) function cal_density_field(x,y,z) result(density_field)
    implicit none
! Calculate the density field at the position (x,y,z)

    integer :: w
    real(DP) :: x,y,z,r(3),rsq

    density_field = 0.d0

    do w=1,num_atom
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
         r(1) = x - atoms(w)%xyz(1)
         r(2) = y - atoms(w)%xyz(2)
         r(3) = z - atoms(w)%xyz(3)
         call mic(r)

! Calculate the density field
        rsq = SUM(r**2)
    
        density_field = density_field + phi(rsq)        
    enddo
        
    return

    end function
! ------------------------------------------------------------------------------------------------------------------------
    real(DP) function cal_op_field(x, y, z) result(op_field)
        implicit none

        real(DP), intent(in) :: x, y, z

        integer :: i

        op_field = 0.d0

        atom_loop: do i=1,num_atom
            r(1) = x - atoms(i)%xyz(1)
            r(2) = y - atoms(i)%xyz(2)
            r(3) = z - atoms(i)%xyz(3)
            call mic(r)

            rsq = SUM(r**2)

            op_field = op_field + ( atoms(i)%op_value * phi(rsq) )

        end do atom_loop

! The real value of op_field at (x,y,z) is: op sum / density_field(x,y,z)
        op_field = op_field / cal_density_field(x,y,z)

        return
        
    end function cal_op_field

!--------------------------------------------------------------------------------------------------------------------------

    subroutine calc_density_gradient(grad,x,y,z)
    implicit none

    real(DP), intent(in) :: x,y,z
    real(DP), intent(inout) :: r(3), grad(3)
    real(DP) :: expterm, rsq
    integer :: w

    grad = 0.d0
    expsum = 0.0d0
    opsum = 0.0d0

    atom_loop: do w=1,num_atom
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
         r(1) = x - atoms(i)%xyz(1)
         r(2) = y - atoms(i)%xyz(2)
         r(3) = z - atoms(i)%xyz(3)
         call mic(r)

! Calculate the gradient:
         rsq = SUM(r**2)
         expterm = phi(rsq)

         grad(1) = grad(1) - ( ximin2 * r(1) * expterm)
         grad(2) = grad(2) - ( ximin2 * r(2) * expterm)
         grad(3) = grad(3) - ( ximin2 * r(3) * expterm)

    end do atom_loop

    return

    end subroutine calc_density_gradient
! ------------------------------------------------------------------------------------------------------------------------
    subroutine cal_gradient(grad, x, y, z)
        implicit none

        real(DP), intent(in) :: x,y,z
        real(DP), intent(out), dimension(3) :: grad

        real(DP), dimension(3) :: op_grad, density_grad
        real(DP) :: density_field, op_field
        integer :: i

        grad = 0.d0
        op_grad = 0.d0

        atom_loop: do i=1,num_atom
            r(1) = x - atoms(i)%xyz(1)
            r(2) = y - atoms(i)%xyz(2)
            r(3) = z - atoms(i)%xyz(3)
            call mic(r)

            rsq = SUM(r**2)

            op_grad(1) = op_grad(1) - (atoms(i)%op_value * ximin2 * r(1) * phi(rsq) )
            op_grad(2) = op_grad(2) - (atoms(i)%op_value * ximin2 * r(2) * phi(rsq) )
            op_grad(3) = op_grad(3) - (atoms(i)%op_value * ximin2 * r(3) * phi(rsq) )

        end do atom_loop

! The gradient of our order-parameter field is a bit more complicated:
!   we need the density field AND the op_field at the point of interest and then we combine all
        density_field = calc_density_field(x,y,z)
        
        op_field = cal_op_field(x,y,z)
        
        call calc_density_grad(density_grad,x,y,z)

! The gradient of our field, finally:
        grad(:) = ( density_field * op_grad(:) - density_grad(:) * op_field ) / ( density_field**2.d0 )
            
        
        return

    end subroutine cal_gradient
! ------------------------------------------------------------------------------------------------------------------------
    real(DP) function cal_mixed_term(x,y,z,density_field,grad)
    implicit none

    integer :: w
    real(DP) :: x,y,z,r(3),rsq
    real(DP) :: mixed_xy,mixed_xz,mixed_yz,mixed_zz,density_field,grad(3),expterm

    cal_mixed_term = 0.d0
    mixed_xy = 0.d0
    mixed_xz = 0.d0
    mixed_yz = 0.d0
    mixed_zz = 0.d0

    do w=1,num_atom
! Find the distance between this position and the position of the w-th oxygen atom by MIC:
            r(1) = x - atoms(i)%xyz(1)
            r(2) = y - atoms(i)%xyz(2)
            r(3) = z - atoms(i)%xyz(3)
            call mic(r)

! Calculate the mixed term:
         rsq = r(1)**2 + r(2)**2 + r(3)**2
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
    subroutine mic(r)
    implicit none

    real(DP), intent(inout) :: r(3)

    do while (dabs(r(1)) .gt. 0.5d0*box_length(1))
     r(1) = r(1) - dsign(box_length(1),r(1))
    enddo
    do while (dabs(r(2)) .gt. 0.5d0*box_length(2))
     r(2) = r(2) - dsign(box_length(2),r(2))
    enddo
    do while (dabs(r(3)) .gt. 0.5d0*box_length(3))
     r(3) = r(3) - dsign(box_length(3),r(3))
    enddo

    return

    end subroutine
! ------------------------------------------------------------------------------------------------------------------------
end module
