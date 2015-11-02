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

    !write(6,*) trim(dname),' ', normal_is, normal_along_z, trim(file_water), trim(file_surface), stride, box_length(:), opref(:), xi
    !stop
     
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
     
     !allocate( zheight(num_atom,2) )

     allocate( gradient(coord(1),coord(2),2,3) )
     allocate( mixed(coord(1),coord(2),2) )

     allocate( surf2(coord(1)-1,coord(2)-1,2,3) )
     allocate( surf(coord(1),coord(2),2,3) )

! Loop through each frame:

     frame_loop: do frame=1,num_frames
     
         ! Check if we need to process every frame or not
         frame_skip: if ( mod(frame,stride) /= 0) then
             
             call skip_frame()
         
         else
             ! Write to standard output, write the frame header to the output files, read-in frame
             ! from XYZ file, and list the number of HOH, DOD and HOD molecules:
             
                   call startio_frame(frame)
                   
             
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
         end if frame_skip

     end do frame_loop

     write(*,*) 'Processing complete.'

    close(dat)
    close(f_wat)
    close(f_sur)
    
    deallocate( atoms,surf,surf2,gradient,mixed )

    contains

        subroutine print_help()
          
          print '(a)', 'Usage: Surface.x [-b] [-i] <ARGUMENTS>'
          print '(a)', ''
          print '(a)', 'Surface.x modes:'
          print '(a)', ''
          print '(a)', '  -i, --interactive        Ask for input interactively'
          print '(a)', ''
          print '(a)', '  -b, --batch        Read from standard input. The order MUST be the following:'
          print '(a)', '    <input_traj> <interface normal (x,y,z)> <wrapped_traj> <surface_file> &
                               & <stride> <box: lx ly lz> <op ref values: 1 2> <Gaussian xi>'

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
                real(DP) :: temp

            ! Parse command line arguments
                narg = command_argument_count()
                call get_command_argument(1, arg)
                arg = trim(adjustl(arg))
            
        parser: if (narg == 0) then
                ! If no arguments are passed, then read from standard input
                
                    print *, "No arguments supplied!"
                    call print_help()
                    stop
            
                ! Interactive mode if '-i' is passed as only option
                else if (narg == 1 .and. (arg =='-i' .or. arg == '--interactive')) then
                
                    write(*,'(a,/)') "Surface.x is in interactive mode ..."
                    
                    write(*,'(3x,a)', advance='no') "Trajectory file [traj.xyz] >> "
                    read(*,'(a)') dname
                        call set_default(dname, 'traj.xyz')
                    ! After reading traj file, we need to know if the expected interface is perpendicular to
                    ! Z-axis (default for this code). If not, we must record this to swap coordinates later
                    write(*,'(3x,a)', advance='no') "Interface normal along which axis? [x, y or z; default is z] >> "
                    read(*,'(a)') normal_is
                        if (normal_is == 'z' .or. normal_is == '') then
                            normal_along_z = .true.
                            normal_is = 'z'
                        else
                            normal_along_z = .false.
                        end if
                    
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
                
                ! If '-b' is passed, batch mode: read from standard input
                else if (narg == 1 .and. (arg == '-b' .or. arg == '--batch')) then
                
                    print *, "Surface.x is in batch mode... reading from standard input"
                    read(5,*) dname, normal_is, file_water, file_surface, stride, box_length(:), opref(:), xi
                    if (normal_is == 'T' .or. normal_is == '') then
                            normal_along_z = .true.
                            normal_is = 'z'
                    else
                            normal_along_z = .false.
                    end if
                
                else
                    print *, "Wrong number of arguments or option not recognized. Stop."
                    call print_help()
                    stop
                                        
                end if parser
                
                ! If normal is not along Z, swap the box lengths
                    if (.not. normal_along_z) then
                        !print *, "Enter in the swap check"
                        if (normal_is == 'x') then
                            temp = box_length(1)
                            box_length(1) = box_length(3)
                            box_length(3) = temp
                        else if (normal_is == 'y') then
                            temp = box_length(2)
                            box_length(2) = box_length(3)
                            box_length(3) = temp
                        end if
                    end if


            return
            
        end subroutine parse_arguments

    
end program surface
