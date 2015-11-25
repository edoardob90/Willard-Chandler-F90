! .........................................................................................................DMW.
! ...................................... SUBROUTINES ..........................................................
! .............................................................................................................
    subroutine start_io()
        use mod_surf
        
        implicit none
        
        integer :: ios

        dat = 7
        f_wat = 8
        f_sur = 9
        f_dist = 10
        f_fft = 11
    ! Find number of lines in file, store as num_frames:
        open(f_wat, file=trim(adjustl(file_water)), status='unknown', action='readwrite')
        open(f_sur, file=trim(adjustl(file_surface)), status='unknown', action='readwrite')     
        open(f_fft, file=trim(adjustl(file_fft)), status='unknown', action='readwrite')
        open(dat, file=trim(adjustl(dname)), status='old', action='read')
        
        num_frames = 0
        do
            read(dat,*, iostat=ios)
            if (ios /= 0) exit
            num_frames = num_frames + 1
        end do
        
        rewind(dat)

    ! Read in number of atoms and use it to calculate number of frames:
        read(dat,*) num_atom
        num_frames = num_frames / (num_atom + 2)
        read(dat,*) ! Throw away the second line of frame 1 containing box size (the code already knows it from input)

        return
        
    end subroutine start_io
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine swap_coord(vector, i, j)
        
        implicit none

        real(DP), dimension(3), intent(inout) :: vector
        integer, intent(in) :: i,j
        real(DP) :: temp

        temp = vector(i)
        vector(i) = vector(j)
        vector(j) = temp

        return
    end subroutine swap_coord
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine skip_frame(f)

        use mod_surf

        implicit none

        integer :: i, f

        !write(6,*) "Skipping frame ", f
        !flush(6)

        ! Skip header lines (first 2)
        read(dat,*)
        read(dat,*)
        
        ! Skip coordinates
        do i=1,num_atom
           read(dat,*)
        end do

        return

    end subroutine skip_frame
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine startio_frame(frame)
    
    use mod_surf
    
    implicit none

    integer :: frame, x, i

! Write to standard output at the beginning of the frame:

      write(6,'(a,i6,a,i6)') 'Processing frame ',frame,' of ',num_frames

! Write frame headers to output files:

10 format(3f10.5)
11 format(a3,x,f12.5,x,f12.5,x,f12.5)
      
      write(f_wat,*) num_atom
      write(f_wat,10) box_length(:)
      
      write(f_sur,'(i6)') coord(1) * coord(2) + (coord(1)-1) * (coord(2)-1)
      write(f_sur,10) box_length(:)
      
      !write(f_dist,*) num_atom
      !write(f_dist,10) box_length(:)

! Read-in frame from input file:

      if (frame > 1) then
          read(dat,*)
          ! Update the box dimensions if we are in NPT ensemble
          read(dat,*) box_length(:)

          ! Swap box dimensions if the interface is NOT perp. to Z
          if (.not. normal_along_z) then
             if (normal_is == 'x') then
                 call swap_coord(box_length, 1, 3)
             else if (normal_is == 'y') then
                 call swap_coord(box_length, 2, 3)
             end if
          end if
      end if

      atom_loop: do x=1,num_atom
! Read from trajectory
    ! First: read current record
        read(dat,*,iostat=i) atoms(x)%symbol, atoms(x)%xyz(:), atoms(x)%op_value
        
    ! If the interface is perpendicular to Z, nothing to be done
    ! Otherwise, we must swap coordinates to get a correct interface orientation along Z
        if ( .not. normal_along_z ) then
            if ( normal_is == 'x' ) then
                ! Swap X and Z
                call swap_coord(atoms(x)%xyz, 1, 3)
            
            else if ( normal_is == 'y') then
                ! Swap Y and Z
                call swap_coord(atoms(x)%xyz, 2, 3)
            end if
        end if

! Apply periodic boundary conditions:
        xscratch = atoms(x)%xyz(1)
        do while (xscratch < 0.d0 .or. xscratch > box_length(1))
            xscratch = xscratch - sign(box_length(1),xscratch)
        enddo
        yscratch = atoms(x)%xyz(2)
        do while (yscratch < 0.d0 .or. yscratch > box_length(2))
            yscratch = yscratch - sign(box_length(2),yscratch)
        enddo
        zscratch = atoms(x)%xyz(3)
        do while (zscratch < 0.d0 .or. zscratch > box_length(3))
            zscratch = zscratch - sign(box_length(3),zscratch)
        enddo

! Write the (possibly shifted) positions of the atoms to the output file:
       write(f_wat,11) atoms(x)%symbol, xscratch, yscratch, zscratch

      end do atom_loop

      return

    end subroutine
