!---------------!
! MAIN PROGRAM
!---------------!
    program surface
    use mod_surf
    implicit none
! This program calculates the Chandler surface for liquid water from an inputted XYZ file.

     integer :: w, frame, i, j, ntimes, L, M

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


! Initialize arrays:
     ! Shortcuts of arrays extents
     L = coord(1)
     M = coord(2)

     allocate( atoms(num_atom) )
     
     allocate( gradient(L,M,2,3) )
     allocate( mixed(L,M,2) )

     allocate( surf2(L-1,M-1,2,3) )
     allocate( surf(L,M,2,3) )

     allocate( h_xy(L-1,M-1))
     allocate( ck_xy(((L-1)/2+1),M-1), ck_xy_r(((L-1)/2+1),M-1), ck_averaged(((L-1)/2+1),M-1) )
     allocate( k_x(L-1), k_y(M-1) )
    


! Loop through each frame:

     ntimes = 0 ! Number of times the FT is computed (for the time average)
     ck_averaged = 0.0d0 ! Averages of the coefficients, set initially to zero

     frame_loop: do frame=1,num_frames
     
         ! Check if we need to process every frame or not
         frame_skip: if ( mod(frame,stride) /= 0) then
             
                          call skip_frame()
                          cycle frame_loop
         
                     else

                         !!! PART ONE: find the Willard-Chandler surface !!!
                         
                         ! Write to standard output, write the frame header to the output files, read-in frame
                         ! from XYZ file, and list the number of HOH, DOD and HOD molecules:
                         
                               call startio_frame(frame)
                               
                         
                         ! Find the constant, const, and find the surface in terms of (x,y) grid points:
                         
                         !!! Value of the CONST already set by CONTOUR
                         ! const should be: opref(liquid)+opref(solid)/2
                         !      const = SUM(opref)/2
                               
                               call find_surface()
                         
                         ! Find the gradient and the mixed terms, and interpolate between grid points to draw the
                         ! surface as a smooth function, in terms of triangles:
                         
                               call find_terms_int()
                               
                               call triangles()

                               
                         !!! PART TWO: Fourier transform of the height profile of the interface !!!
                            
                         ! Reshape the matrix surface to be a 2-dimensional matrix
                         ! Each element z=(i,j) is:
                         !  i = point of the grid along X
                         !  j = point of the grid along Y
                         !  z = height of the surface at the corresponding grid point

                         fft_compute: if (compute_fft) then

                               plan = fftw_plan_dft_r2c_2d(M-1, L-1, h_xy, ck_xy, FFTW_MEASURE)

                               ! Wave vectors of the Fourier sum
                               ! In general: k_(x,y) = (i,j) * 2*pi/L_(x,y)
                               do i=1,L-1
                                  k_x(i) = (i-1) * 2*pi/box_length(1)
                               end do
                               do j=1,M-1
                                  k_y(j) = (j-1) * 2*pi/box_length(2)
                               end do

                               ! Profile height of the surface at current frame
                               forall (i=1:L-1, j=1:M-1)
                                   h_xy(i,j) = surf2(i,j,1,3) 
                               end forall
                        
                               !write(999,*) "At frame ", frame, ":"
                               !write(999,*) h_xy(1,1), surf(1,1,1,3), "ref: ", h_xy_t0(1,1)
                               
                               write(6,'(a,i8)') "Computing Fourier transform at frame ", frame
 
                               call fftw_execute_dft_r2c(plan, h_xy, ck_xy)

                               ! Compute the square modulus of each Fourier coefficients
                               ck_xy_r(:,:) = ABS(ck_xy(:,:))**2
     
                               !!!if (frame /= stride) then
                               !!!   write(999,*) "*************************"
                               !!!   write(999,*) "Frame: ", frame
                               !!!   write(999,*) "*************************"
                               !!!   do i=1,L-1
                               !!!      do j=1,M-1
                               !!!         write(999,*) k_x(i), k_y(j), ck_xy_r(i,j)
                               !!!      end do
                               !!!   end do
                               !!!   write(999,*) "*************************"
                               !!!end if

                               ck_averaged(:,:) = ck_averaged(:,:) + ck_xy_r(:,:)
                               

                               ntimes = ntimes + 1
                         
                           end if fft_compute
                         
                     end if frame_skip

     end do frame_loop

     ! Compute the average of the Fourier coefficients and write output
     if (compute_fft) then
         ck_averaged = ck_averaged / ntimes

         write(output_unit,*) "Writing averaged Fourier coefficients ..."

         do i=1,(L-1)/2+1
            do j=1,(M-1)/2+1
               write(f_fft,'(2f10.5,x,e15.7)') k_x(i), k_y(j), ck_averaged(i,j)
            end do
         end do
     end if

     write(output_unit,*) 'Processing complete.'

    close(dat)
    close(f_wat)
    close(f_sur)
    close(f_fft)
    
    deallocate( atoms, surf, surf2, gradient, mixed, ck_averaged, h_xy, ck_xy, ck_xy_r, k_x, k_y)
    
    call fftw_destroy_plan(plan)

    contains

        subroutine print_help()
          
          write(error_unit,'(a)') 'Usage: Surface.x [-b < INPUT] [-i] [-h]'
          write(error_unit,*)
          write(error_unit,'(a)') 'Surface.x switches:'
          write(error_unit,*)
          write(error_unit,'(a)') '  -i, --interactive        Ask for input interactively'
          write(error_unit,*)
          write(error_unit,'(a)') '  -b, --batch              Read from standard input. Type -h to know what arguments are accepted'
          write(error_unit,*)
          write(error_unit,'(a)') '  -h, --help               Detailed description of arguments and their order'

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
                character(len=3) :: ext
                integer :: narg, llen
                real(DP) :: temp

            ! Parse command line arguments
                narg = command_argument_count()
                call get_command_argument(1, arg)
                arg = trim(adjustl(arg))
            
        parser: if (narg == 0) then
                ! If no arguments are passed, then read from standard input
                
                    write(error_unit,*) "No arguments supplied!"
                    call print_help()
                    stop
            
                ! Interactive mode if '-i' is passed as only option
                else if (narg == 1 .and. (arg =='-i' .or. arg == '--interactive')) then
                
                    write(output_unit,'(a,/)') "Surface.x is in interactive mode ..."
                    
                    write(output_unit,'(3x,a)', advance='no') "Trajectory file [traj.xyz] >> "
                    read(input_unit,'(a)') dname
                        call set_default(dname, 'traj.xyz')
                    ! After reading traj file, we need to know if the expected interface is perpendicular to
                    ! Z-axis (default for this code). If not, we must record this to swap coordinates later
                    write(output_unit,'(3x,a)', advance='no') "Interface normal along which axis? [x, y or z; default is z] >> "
                    read(input_unit,'(a)') normal_is
                        if (normal_is == 'z' .or. normal_is == '') then
                            normal_along_z = .true.
                            normal_is = 'z'
                        else
                            normal_along_z = .false.
                        end if
                    
                    write(output_unit,'(3x,a)', advance='no') "Output files [wrap_traj.xyz | surface.out] >> "
                    read(input_unit,'(2a)') file_water, file_surface
                        call set_default(file_water, 'wrap_traj.xyz')
                        call set_default(file_surface, 'surface.out')
                    
                    write(output_unit,'(3x,a)', advance='no') "Stride [0 for default] >> "
                    read(input_unit,*) stride
                        if (stride == 0) stride = 1
                    
                    write(output_unit,'(3x,a)', advance='no') "Box size (lx ly lz) >> "
                    read(input_unit,*) box_length(:)
                    
                    write(output_unit,'(3x,a)', advance='no') "Order parameter reference values (op1 op2) >> "
                    read(input_unit,*) contour
                    
                    write(output_unit,'(3x,a)', advance='no') "Bandwidth of Gaussian kernel [xi] >> "
                    read(input_unit,*) xi
                    
                    write(output_unit,'(3x,a)', advance='no') "Grid spacing >> "
                    read(input_unit,*) grid_spacing
                
                    write(output_unit,'(3x,a)', advance='no') "Compute FT of height profile? [Y/n] >> "
                    read(input_unit,'(a)') fft_answer
                        if (fft_answer(1:1) == "y" .or. fft_answer(1:1) == "Y") then
                            compute_fft = .true.
                            fft_answer = "yes"
                        else
                            compute_fft = .false.
                            fft_answer = "no"
                        end if
                    
                    write(output_unit,'(3x,a)', advance='no') "File of REAL Fourier coefficients [fourier.dat] >> "
                    read(input_unit,'(a)') file_fft
                        call set_default(file_fft, 'fourier.dat')

                ! If '-b' is passed, batch mode: read from standard input
                else if (narg == 1 .and. (arg == '-b' .or. arg == '--batch')) then
                
                    write(output_unit,'(a,/)') "Surface.x is in batch mode... reading from standard input"
                    
                    read(input_unit,*) dname, normal_is, file_water, file_surface, stride, box_length(:), contour, xi, grid_spacing, fft_answer, file_fft

20 format(a,/,3x,2a,/,3x,3a,/,3x,4a,/,3x,a,i10,a,/,3x,a,3f10.5,/,3x,a,f10.5,/,3x,a,f5.3,/,3x,a,f5.3,a,/,3x,2a,/,3x,2a,/)

                    write(output_unit,20) & 
                        & "Willard-Chandler surface will be computed according to the following input:", &
                        & "Input trajectory: ", to_upper(trim(dname)), &
                        & "Interface normal along ", to_upper(normal_is), " axis", &
                        & "Output files: ", to_upper(trim(file_water)), " and ", to_upper(trim(file_surface)), &
                        & "Processing every ", stride, " frames", &
                        & "Simulation box size: ", box_length(:), &
                        & "Value of the CONTOUR where to draw the surface: ", contour, &
                        & "Bandwidth of the Gaussian kernel: ", xi, & 
                        & "Grid spacing: ", grid_spacing, " (in appropriate distance units)", &
                        & "Computing FT of the interface height profile: ", to_upper(fft_answer), &
                        & "REAL Fourier coefficients will be written to: ", to_upper(trim(file_fft))

                    ! Check surface normal
                    if (normal_is == 'z' .or. normal_is == '') then
                            normal_along_z = .true.
                            normal_is = 'z'
                    else
                            normal_along_z = .false.
                    end if

                    ! Check if FFT is to be computed
                    if (fft_answer(1:1) == 'y' .or. fft_answer(1:1) == 'Y') then
                        compute_fft = .true.
                    else
                        compute_fft = .false.
                    end if

                
                else if (narg == 1 .and. (arg == '-h' .or. arg == '--help')) then
                    call print_arguments()
                    stop
                
                else
                    write(error_unit,*) "Wrong number of arguments or option not recognized. Stop."
                    call print_help()
                    stop
                                        
                end if parser
                
                ! If normal is not along Z, swap the box lengths
                    if (.not. normal_along_z) then
                        !print *, "Enter in the swap check"
                        if (normal_is == 'x') then
                            call swap_coord(box_length, 1, 3)
                        else if (normal_is == 'y') then
                            call swap_coord(box_length, 2, 3)
                        end if
                    end if
            return
            
        end subroutine parse_arguments


        function to_upper(strIn) result(strOut)
        
             implicit none
        
             character(len=*), intent(in) :: strIn
             character(len=len(strIn)) :: strOut
             integer :: i,j
        
             do i = 1, len(strIn)
                  j = iachar(strIn(i:i))
                  if (j>= iachar("a") .and. j<=iachar("z") ) then
                       strOut(i:i) = achar(iachar(strIn(i:i))-32)
                  else
                       strOut(i:i) = strIn(i:i)
                  end if
             end do
        
        end function to_upper


        subroutine print_arguments()

            implicit none

            write(error_unit, '(a)') "Surface.x detailed list of arguments for batch mode ..."
            write(error_unit, '(a)') "ARGUMENTS. Their order MUST be the following: ", &
                        & "   -> Input trajectory filename (XYZ file)", &
                        & "   -> Interface normal along [x/y/z] axis", &
                        & "   -> Output files: wrapped trajectory, surface plot", &
                        & "   -> Stride of processing input trajectory", &
                        & "   -> Simulation box size (3 numbers, SPACE separated)", &
                        & "   -> Value of the CONTOUR where to draw the surface", &
                        & "   -> Bandwidth of the Gaussian kernel", & 
                        & "   -> Spacing of the regular grid", &
                        & "   -> Compute the FT of the surface profile? [Y/N]", &
                        & "   -> Output file where to write time-averaged, real Fourier coefficients (unnormalized)"


            return

        end subroutine print_arguments

    
end program surface
