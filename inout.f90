! .........................................................................................................DMW.
! ...................................... SUBROUTINES ..........................................................
! .............................................................................................................
    subroutine start_io()
    use mod_surf
    implicit none

     dat = 7
     f_wat = 8
     f_sur = 9
     f_dist = 10
! Find number of lines in file, store as num_frames:
     open(f_wat,file=file_water)
     open(f_sur,file=file_surface)
     open(f_dist,file=file_dist)
     open(dat,file=dname,status='old')
     num_frames = 0
25   read(dat,'(a)',err=50,end=50)
     num_frames = num_frames + 1
     goto 25
50   rewind(dat)

! Read in number of atoms and use it to calculate number of frames:
     read(dat,*) num_atom
     num_frames = num_frames / (num_atom + 2)
     read(dat,'(A)') input
     input1 = input
     input = adjustl(input(2:len(input)))

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine startio_frame(frame)
    use mod_surf
    implicit none

    integer frame,w,x,y
    real(8) com(3),mass(3)

! Write to standard output at the beginning of the frame:

      write(*,*) 'Processing frame ',frame,' of ',num_frames

! Write frame headers to output files:

      write(f_wat,*) num_atom
      write(f_wat,*) input1
      write(f_sur,*) 2 * coord(1) * coord(2) + 2 * (coord(1)-1) * (coord(2)-1)
!      write(f_sur,*) 2 * coord(1) * coord(2)
      write(f_sur,*) input1
      write(f_dist,*) num_atom
      write(f_dist,*) input1

! Read-in frame from input file:

      if (frame .gt. 1) then
       read(dat,'(A)')
       read(dat,'(A)')
      endif

      w = 1
      do x=1,num_atom
       read(dat,*) atom_name(x),(atom_position(x,y),y=1,3,1)
       if (atom_name(x) .eq. 'O') then
        oxygen_list(w) = x
        w = w + 1
       endif
      enddo

! Write atom position to output file; water molecules are folded back into the simulation box,
! for aesthetic reasons:

!       xscratch = atom_position(x,1)
!       do while (xscratch .lt. 0.d0 .or. xscratch .gt. box_length(1))
!        xscratch = xscratch - sign(box_length(1),xscratch)
!       enddo
!       yscratch = atom_position(x,2)
!       do while (yscratch .lt. 0.d0 .or. yscratch .gt. box_length(2))
!        yscratch = yscratch - sign(box_length(2),yscratch)
!       enddo

!       write(f_wat,*) atom_name(x),'    ',xscratch,'    ',yscratch,'    ',atom_position(x,3)
!      enddo

! Write atom positions to output file -- the centre of mass of each water molecule is folded back into
! the simulation box, for aesthetic reasons:

      do x=1,num_atom,3
       mass(1) = 16.d0
       mass(2) = 1.d0
       mass(3) = 1.d0
       com(:) = 0.d0

       if (atom_name(x+1) .eq.'D') mass(2) = 2.d0
       if (atom_name(x+2) .eq. 'D') mass(3) = 2.d0
       do y=0,2
        com(:) = com(:) + mass(y+1)*atom_position(x+y,:)
       enddo
       com(:) = com(:) / sum(mass)

! Apply periodic boundary conditions:
       xscratch = com(1)
       do while (xscratch .lt. 0.d0 .or. xscratch .gt. box_length(1))
        xscratch = xscratch - sign(box_length(1),xscratch)
       enddo
       yscratch = com(2)
       do while (yscratch .lt. 0.d0 .or. yscratch .gt. box_length(2))
        yscratch = yscratch - sign(box_length(2),yscratch)
       enddo

! Write the (possibly shifted) positions of the atoms to the output file:

       write(f_wat,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") &
     &      atom_name(x),' ',atom_position(x,1)+(xscratch-com(1)),' ', &
     &     atom_position(x,2)+(yscratch-com(2)),' ',atom_position(x,3)
       write(f_wat,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") &
     &     atom_name(x+1),' ',atom_position(x+1,1)+(xscratch-com(1)),' ', &
     &     atom_position(x+1,2)+(yscratch-com(2)),' ',atom_position(x+1,3)
       write(f_wat,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") &
     &     atom_name(x+2),' ',atom_position(x+2,1)+(xscratch-com(1)),' ', &
     &     atom_position(x+2,2)+(yscratch-com(2)),' ',atom_position(x+2,3)

      enddo

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine list_molecules()
    use mod_surf
    implicit none

    integer x
    integer num_HOH,num_DOD,num_HOD

    num_HOH = 0
    num_DOD = 0
    num_HOD = 0

! Count number of HOH,HOD and DOD molecules:

       do x=1,num_atom/3
        if (atom_name(3*x - 2) .ne. 'O') then
         write(*,*) 'Not a catastrophic error, but more programming needs to be done!'
         stop
        endif
        if (atom_name(3*x - 1) .eq. 'H' .and. atom_name(3*x) .eq. 'H') then
         num_HOH = num_HOH + 1
        else if (atom_name(3*x - 1) .eq. 'D' .and. atom_name(3*x) .eq. 'D') then
         num_DOD = num_DOD + 1
        else if ((atom_name(3*x - 1) .eq. 'H' .and. atom_name(3*x) .eq. 'D') .or. &
     &          (atom_name(3*x - 1) .eq. 'D' .and. atom_name(3*x) .eq. 'H')) then
         num_HOD = num_HOD + 1
        else
         write(*,*) 'Should not get here! -- molecule',x,' is neither HOH, HOD or DOD.'
         write(*,*) atom_name(3*x - 2),atom_name(3*x - 1),atom_name(3*x),3*x - 2
         stop
        endif
       enddo

       write(*,*) 'Number of molecules of each type => fraction:'
       write(*,*) 'HOH: ',num_HOH,'=>',num_HOH/float(num_atom/3)
       write(*,*) 'DOD: ',num_DOD,'=>',num_DOD/float(num_atom/3)
       write(*,*) 'HOD: ',num_HOD,'=>',num_HOD/float(num_atom/3)

    end subroutine
! -----------------------------------------------------------------------------------------------------------------------------
    subroutine write_distances(atom_n,w,height)
    use mod_surf
    implicit none

    integer w,atom_n
    real(8) height

    if (atom_n .eq. 0) atom_n = 3

! Calculate the three vectors vec1 (from the oxygen atom to the midpoint
! between the two hydrogens), vec2 and vec3 (the two O-H/D bond vectors):

    if (atom_name(w) .eq. 'O') then
     vec1 = 0.5d0*(atom_position(w+1,:)+atom_position(w+2,:))-atom_position(w,:)
     vec2 = atom_position(w+1,:) - atom_position(w,:)
     vec3 = atom_position(w+2,:) - atom_position(w,:)
     vec1 = vec1 / dsqrt(vec1(1)**2 + vec1(2)**2 + vec1(3)**2)
     vec2 = vec2 / dsqrt(vec2(1)**2 + vec2(2)**2 + vec2(3)**2)
     vec3 = vec3 / dsqrt(vec3(1)**2 + vec3(2)**2 + vec3(3)**2)
    endif

! Write to the distances output file the distance from each atom to the Chandler surface in the z direction, as well
! as the projections of the three important vectors on the unit surface normal, and the projections in the z direction:

    select case (atom_n)
     case (1)
      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
     &         -1.d0*((grad_interp(1)*vec1(1)) + (grad_interp(2)*vec1(2)) + (grad_interp(3)*vec1(3))), &
     &         ' ', vec1(3)
     case (2)
      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
     &          -1.d0*((grad_interp(1)*vec2(1)) + (grad_interp(2)*vec2(2)) + (grad_interp(3)*vec2(3))), &
     &         ' ', vec2(3)
     case (3)
      write(f_dist,"(A3,A1,E12.5,A1,E12.5,A1,E12.5)") atom_name(w),' ',-1.d0*height,' ', &
     &          -1.d0*((grad_interp(1)*vec3(1)) + (grad_interp(2)*vec3(2)) + (grad_interp(3)*vec3(3))), &
     &         ' ', vec2(3)
    end select

    end subroutine
