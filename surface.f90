! .........................................................................................................DMW.
! .................................... MAIN PROGRAM ...........................................................
! .............................................................................................................
    program surface
    use mod_surf
    implicit none
! This program calculates the Chandler surface for liquid water from an inputted XYZ file.

     integer :: w,frame
     
! For now, unit cell is orthorhombic, so that the matrix is unity:

     unit_cell(:,:) = 0.d0
     unit_cell(1,1) = 1.d0
     unit_cell(2,2) = 1.d0
     unit_cell(3,3) = 1.d0
     unit_cell_inv = unit_cell

     if (iargc() .ne. 4) then
      write(*,*) 'USAGE: surface [input file] [water output file] [surface output file] [distances output file]'
      stop
     endif 



! Open input and output files and begin the I/O:

     !call getarg(1,dname)
     !call getarg(2,file_water)
     !call getarg(3,file_surface)
     !call getarg(4,file_dist)
     
     open(unit=100, file='SURFACE_INPUT', status='old', action='read')
     ! Read in arguments: dname traj file_surface, ...
     read(100,'(3a,i)') dname, file_water, file_surface, stride
     read(100,*) box_length(:)
     read(100,*) xi, opref1, opref2
     
     
     call start_io()
     
     
! Initialize constants:

     call init_consts()

! Calculate grid:

     call calculate_grid()

     write(*,*) box_length

! Initialize arrays:

     allocate( atom_name(num_atom) )
     allocate( atom_position(num_atom,3) )
     !allocate( oxygen_list(num_atom / 3) )
     allocate( zheight(num_atom,2) )

     allocate( gradient(coord(1),coord(2),2,3) )
     allocate( mixed(coord(1),coord(2),2) )

     allocate( surf2(coord(1)-1,coord(2)-1,2,3) )
     allocate( surf(coord(1),coord(2),2,3) )

! Loop through each frame:

     do frame=1,num_frames

! Write to standard output, write the frame header to the output files, read-in frame
! from XYZ file, and list the number of HOH, DOD and HOD molecules:

      call startio_frame(frame)
      call list_molecules()

! Find the constant, const, and find the surface in terms of (x,y) grid points:

      const = cal_density_field(box_length(1)*0.5d0,box_length(2)*0.5d0,0.d0) * 0.5d0
      call find_surface()

! Find the gradient and the mixed terms, and interpolate between grid points to draw the
! surface as a smooth function, in terms of triangles:

      call find_terms_int()
      call triangles()

! Now that processing has been carried out on this frame, the height of an atom compared to each surface is found, and the number
! of H and D atoms above the surface are counted.

      do w=1,num_atom

! Find the height of this atom below both the upper and lower surface and the unit vector normal
! to the corresponding point on the surface:

       call find_atom_height_grad(w)

! Calculate the three vectors vec1, vec2 and vec3, and print out the atom's height above the surface:

       call write_distances(mod(w,3),w,min(zheight(w,1),zheight(w,2)))

      enddo

     enddo

     write(*,*) 'Processing complete.'

    close(dat)
    close(f_wat)
    close(f_sur)
    deallocate( atom_name,atom_position,surf,surf2,zheight,gradient,mixed )
    end program
! ---------------------------------------------------------------------------------------------------------------------
