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
    ! Find number of lines in file, store as num_frames:
        open(f_wat, file=trim(adjustl(file_water)), status='unknown', action='readwrite')
        open(f_sur, file=trim(adjustl(file_surface)), status='unknown', action='readwrite')     
        !open(f_dist,file=file_dist)
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
    subroutine skip_frame()

        use mod_surf

        implicit none

        integer :: i

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
    !real(8) :: com(3),mass(3)

! Write to standard output at the beginning of the frame:

      write(*,*) 'Processing frame ',frame,' of ',num_frames

! Write frame headers to output files:

10 format(3f10.5)
11 format(a3,x,e12.5,x,e12.5,x,e12.5)
      
      write(f_wat,*) num_atom
      write(f_wat,10) box_length(:)
      
      write(f_sur,*) 2 * coord(1) * coord(2) + 2 * (coord(1)-1) * (coord(2)-1)
      write(f_sur,10) box_length(:)
      
      !write(f_dist,*) num_atom
      !write(f_dist,10) box_length(:)

! Read-in frame from input file:

      if (frame .gt. 1) then
       read(dat,*)
       read(dat,*)
      endif

      atom_loop: do x=1,num_atom
! Read from trajectory
    ! If the interface is perpendicular to Z, nothing to be done
    ! Otherwise, we must swap the correct coordinate to get a correct interface orientation
        if ( normal_along_z ) then
            read(dat,*) atoms(x)%symbol, atoms(x)%xyz(:), atoms(x)%op_value
        
        else
            
            if ( normal_is == 'x' ) then
                ! Swap X and Z
                read(dat,*) atoms(x)%symbol, (atoms(x)%xyz(i), i=3,1,-1), atoms(x)%op_value
            
            else if ( normal_is == 'y') then
                ! Swap Y and Z
                read(dat,*) atoms(x)%symbol, atoms(x)%xyz(1), atoms(x)%xyz(3), atoms(x)%xyz(2), atoms(x)%op_value
            
            end if
            
        end if

! Apply periodic boundary conditions:
        xscratch = atoms(x)%xyz(1)
        do while (xscratch .lt. 0.d0 .or. xscratch .gt. box_length(1))
            xscratch = xscratch - sign(box_length(1),xscratch)
        enddo
        yscratch = atoms(x)%xyz(2)
        do while (yscratch .lt. 0.d0 .or. yscratch .gt. box_length(2))
            yscratch = yscratch - sign(box_length(2),yscratch)
        enddo

! Write the (possibly shifted) positions of the atoms to the output file:
       write(f_wat,11) atoms(x)%symbol, atoms(x)%xyz(1)+xscratch, atoms(x)%xyz(2)+yscratch, atoms(x)%xyz(3)

      end do atom_loop

      return

    end subroutine

    
! -----------------------------------------------------------------------------------------------------------------------------
!!!    subroutine list_molecules()
!!!    use mod_surf
!!!    implicit none
!!!
!!!    integer x
!!!    integer num_HOH,num_DOD,num_HOD
!!!
!!!    num_HOH = 0
!!!    num_DOD = 0
!!!    num_HOD = 0
!!!
!!!! Count number of HOH,HOD and DOD molecules:
!!!
!!!       do x=1,num_atom/3
!!!        if (atom_name(3*x - 2) .ne. 'O') then
!!!         write(*,*) 'Not a catastrophic error, but more programming needs to be done!'
!!!         stop
!!!        endif
!!!        if (atom_name(3*x - 1) .eq. 'H' .and. atom_name(3*x) .eq. 'H') then
!!!         num_HOH = num_HOH + 1
!!!        else if (atom_name(3*x - 1) .eq. 'D' .and. atom_name(3*x) .eq. 'D') then
!!!         num_DOD = num_DOD + 1
!!!        else if ((atom_name(3*x - 1) .eq. 'H' .and. atom_name(3*x) .eq. 'D') .or. &
!!!     &          (atom_name(3*x - 1) .eq. 'D' .and. atom_name(3*x) .eq. 'H')) then
!!!         num_HOD = num_HOD + 1
!!!        else
!!!         write(*,*) 'Should not get here! -- molecule',x,' is neither HOH, HOD or DOD.'
!!!         write(*,*) atom_name(3*x - 2),atom_name(3*x - 1),atom_name(3*x),3*x - 2
!!!         stop
!!!        endif
!!!       enddo
!!!
!!!       write(*,*) 'Number of molecules of each type => fraction:'
!!!       write(*,*) 'HOH: ',num_HOH,'=>',num_HOH/float(num_atom/3)
!!!       write(*,*) 'DOD: ',num_DOD,'=>',num_DOD/float(num_atom/3)
!!!       write(*,*) 'HOD: ',num_HOD,'=>',num_HOD/float(num_atom/3)
!!!
!!!    end subroutine
!!!! -----------------------------------------------------------------------------------------------------------------------------
!!!    subroutine write_distances(atom_n,w,height)
!!!    use mod_surf
!!!    implicit none
!!!
!!!    integer w,atom_n
!!!    real(8) height
!!!
!!!    if (atom_n .eq. 0) atom_n = 3
!!!
!!!! Calculate the three vectors vec1 (from the oxygen atom to the midpoint
!!!! between the two hydrogens), vec2 and vec3 (the two O-H/D bond vectors):
!!!
!!!    if (atom_name(w) .eq. 'O') then
!!!     vec1 = 0.5d0*(atom_position(w+1,:)+atom_position(w+2,:))-atom_position(w,:)
!!!     vec2 = atom_position(w+1,:) - atom_position(w,:)
!!!     vec3 = atom_position(w+2,:) - atom_position(w,:)
!!!     vec1 = vec1 / dsqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
!!!     vec2 = vec2 / dsqrt(vec2(1)**2 + vec2(2)**2 + vec2(3)**2)
!!!     vec3 = vec3 / dsqrt(vec3(1)**2 + vec3(2)**2 + vec3(3)**2)
!!!    endif
!!!
!!!! Write to the distances output file the distance from each atom to the Chandler surface in the z direction, as well
!!!! as the projections of the three important vectors on the unit surface normal, and the projections in the z direction:
!!!
!!!    select case (atom_n)
!!!     case (1)
!!!      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
!!!     &         -1.d0*((grad_interp(1)*vec1(1)) + (grad_interp(2)*vec1(2)) + (grad_interp(3)*vec1(3))), &
!!!     &         ' ', vec1(3)
!!!     case (2)
!!!      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
!!!     &          -1.d0*((grad_interp(1)*vec2(1)) + (grad_interp(2)*vec2(2)) + (grad_interp(3)*vec2(3))), &
!!!     &         ' ', vec2(3)
!!!     case (3)
!!!      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
!!!     &          -1.d0*((grad_interp(1)*vec3(1)) + (grad_interp(2)*vec3(2)) + (grad_interp(3)*vec3(3))), &
!!!     &         ' ', vec2(3)
!!!    end select
!!!
!!!    end subroutine
