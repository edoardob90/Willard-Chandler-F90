!---------------!
! MAIN PROGRAM
!---------------!
    program surface
    use mod_surf
    implicit none
! This program calculates the Chandler surface for liquid water from an inputted XYZ file.

     integer :: w, frame

! For now, unit cell is orthorhombic, so that the matrix is unity:

     unit_cell(:,:) = 0.d0
     unit_cell(1,1) = 1.d0
     unit_cell(2,2) = 1.d0
     unit_cell(3,3) = 1.d0
     unit_cell_inv = unit_cell

! Parse command line for input

    call parse_arguments()

     
! Open input and output files and begin the I/O:
     call start_io()
     
     
! Initialize constants:

     call init_consts()

! Calculate grid:

     call calculate_grid()

     write(*,*) box_length

! Initialize arrays:

     !allocate( atom_name(num_atom) )
     !allocate( atom_position(num_atom,3) )

     allocate( atoms(num_atom) )
     
     allocate( zheight(num_atom,2) )

     allocate( gradient(coord(1),coord(2),2,3) )
     allocate( mixed(coord(1),coord(2),2) )

     allocate( surf2(coord(1)-1,coord(2)-1,2,3) )
     allocate( surf(coord(1),coord(2),2,3) )

! Loop through each frame:

     frame_loop: do frame=1,num_frames

! Write to standard output, write the frame header to the output files, read-in frame
! from XYZ file, and list the number of HOH, DOD and HOD molecules:

      call startio_frame(frame)
!      call list_molecules()

! Find the constant, const, and find the surface in terms of (x,y) grid points:

! const should be: opref(liquid)+opref(solid)/2
      const = SUM(opref)/2
      call find_surface()

! Find the gradient and the mixed terms, and interpolate between grid points to draw the
! surface as a smooth function, in terms of triangles:

      call find_terms_int()
      call triangles()

!!!! Now that processing has been carried out on this frame, the height of an atom compared to each surface is found, and the number
!!!! of H and D atoms above the surface are counted.
!!!
!!!      do w=1,num_atom
!!!
!!!! Find the height of this atom below both the upper and lower surface and the unit vector normal
!!!! to the corresponding point on the surface:
!!!
!!!       call find_atom_height_grad(w)
!!!
!!!! Calculate the three vectors vec1, vec2 and vec3, and print out the atom's height above the surface:
!!!
!!!       call write_distances(mod(w,3),w,min(zheight(w,1),zheight(w,2)))
!!!
!!!      enddo

     end do frame_loop

     write(*,*) 'Processing complete.'

    close(dat)
    close(f_wat)
    close(f_sur)
    deallocate( atoms,surf,surf2,zheight,gradient,mixed )

    contains

        subroutine print_help()
          print '(a)', 'Usage: Surface.x [OPTIONS]'
          print '(a)', ''
          print '(a)', 'Surface.x modes:'
          print '(a)', ''
          print '(a)', '  -i, --interactive        Ask for input interactively'
          print '(a)', ''
          print '(a)', '  -b, --batch        In-line input. The order MUST be the following:'
          print '(a)', '    <input_traj> <wrapped_traj> <surface_file> &
                               & <stride> <box: lx ly lz> <op ref values: 1 2> <xi>'

          return
        end subroutine print_help

        subroutine set_default(filename, def_filename)
            implicit none
            character(len=*) :: filename, def_filename

            if (filename == ' ') filename=def_filename

            return

        end subroutine set_default

        subroutine parse_arguments()

            implicit none
                
                character(len=100) :: arg
                integer :: narg

            ! Parse command line arguments
                narg = command_argument_count()
                call get_command_argument(1, arg)
                arg = trim(adjustl(arg))
            
                ! No arguments is not good
                if (narg == 0) then
                    print *, "No arguments supplied. Are you kidding?!"
                    call print_help()
                    stop
                end if
            
                ! Interactive mode if '-i' is passed as only option
                if (narg == 1 .and. (arg =='-i' .or. arg == '--interactive')) then
                    write(*,'(a,/)') "Surface.x is in interactive mode ..."
                    write(*,'(3x,a)', advance='no') "Trajectory file [traj.xyz] >> "
                    read(*,'(a)') dname
                        call set_default(dname, 'traj.xyz')
                    write(*,'(3x,a)', advance='no') "Output files [wrap_traj.xyz | surface.out] >> "
                    read(*,'(2a)') file_water, file_surface
                        call set_default(file_water, 'wrap_traj.xyz')
                        call set_default(file_surface, 'surface.out')
                    write(*,'(3x,a)', advance='no') "Stride [0 for default] >> "
                    read(*,*) stride
                        if (stride == 0) stride = 1
                    write(*,'(3x,a)', advance='no') "Box size (lx ly lz) >> "
                    read(*,*) box_length(:)
                    write(*,'(3x,a)', advance='no') "Order parameter reference values (op1 op2) >> "
                    read(*,*) opref(:)
                    write(*,'(3x,a)', advance='no') "Gaussian variance [xi] >> "
                    read(*,*) xi
                else if (narg == 1 .and. (arg == '-b' .or. arg == '--batch')) then
                    print *, "Surface.x is in batch mode... but no arguments supplied!"
                    call print_help()
                    stop
                else if (narg > 1) then
                    if (narg < 11) then
                        print *, "Surface.x is in batch mode... but more input is needed!"
                        call print_help()
                        stop
                    else
                        call get_command_argument(1,arg)
                        batch: if (arg == '-b' .or. arg == '--batch') then
                            ! Batch mode
                            !write(6,'(a,/,a)') "Batch mode: order must be supplied in THIS exact order", &
                            !    & "<input_traj> <wrapped_traj> <surface_file> &
                            !    & <stride> <box: lx ly lz> <op ref values: 1 2> <xi>"
                            arguments: do w=2,narg
                                call get_command_argument(w, arg)
                                select case (w)
                                case (2)
                                    dname = arg
                                case (3)
                                    file_water = arg
                                case (4)
                                    file_surface = arg
                                case (5)
                                    !arg = trim(adjustl(arg))
                                    read(arg,'(i)') stride
                                case (6,7,8) ! box size
                                    !arg = trim(adjustl(arg))
                                    read(arg,*) box_length(w-5)
                                case (9,10) ! op ref values
                                    !arg = trim(adjustl(arg))
                                    read(arg,*) opref(w-8)
                                case (11)
                                    !arg = trim(adjustl(arg))
                                    read(arg,*) xi
                                end select
                            end do arguments
                        else
                            write(*,*) "Problem in parsing arguments. Did you supplied them correctly?"
                            call print_help()
                            stop
                        end if batch
                    end if
                end if


            return
        end subroutine parse_arguments
    
end program surface
