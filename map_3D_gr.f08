! USAGE: this_file.x -crd [COORDINATE FILE] -fit [FIT PARAMETER FILE] -cfg [CONFIGURATION FILE] -out [OUTPUT FILE]
!
!
! CRDfile = [solute coordinate file name]
! FitFile = [fitting parameter file name]
! CfgFile = [python config file]
! OutFile = [output file name]
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!   Modules   !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize the xyz-coordinates of the solute molecule, the array of atom types corresponding to those coordinates as well as the
! number or atoms and number of atom types.
module crdData
	double precision, allocatable :: x(:,:)
	integer, allocatable :: aType(:)
	integer nAtoms, nTypes

endmodule crdData

! initialize the array of fit parameter types (x1, x2, x0), an array containing the value of each parameter, and the total number
! of parameters given the number of atom types.
module fitData
	integer, allocatable :: pType(:)
	double precision, allocatable :: pValue(:)
	integer nVars

endmodule fitData

! initialize the values to be read from the gr.2D.py script configuration file.
module cfgData
	double precision :: hist_dist_min, hist_dist_max, bin_dist_size, hist_ang_min, hist_ang_max, bin_ang_size

endmodule cfgData

! initialize the arrays and values to be used in calculating the g(r) values from each atom
module grData
	double precision, allocatable :: gr(:,:,:), xAxis(:), yAxis(:), zAxis(:)
	double precision :: R(3), rSolv(3)
	double precision mag_x, x1, x2, x0, gx_out

endmodule grData


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program map_3D_gr
	implicit none
	character*64 crdFile, fitFile, cfgFile, outFile
	double precision omp_get_wtime, ti, tf

	ti = omp_get_wtime()

	! read the crd and fit file names from the command line.
	call parse_command_line(crdFile,fitFile,cfgFile,outFile)

	! read the solute coordinates and atom types from the crd file.
	call read_crd(crdFile)

	! read the fitting parameters for the IS-SPA fit from fit file.
	call read_fit(fitFile)

	! read the python config parameters.
	call read_cfg(cfgFile)

	! construct IS-SPA 3D g(r).
	call map_gr(outFile)
	
	! Print time taken to finish calculation.
	tf = omp_get_wtime()
	write(*,*) "Total time elapsed: ", tf-ti, "seconds"

endprogram map_3D_gr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Subroutines !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! read coordinate file, fit parameter file, config file, and output file from command line
subroutine parse_command_line(crdFile,fitFile,cfgFile,outFile)
	implicit none
	character*16 arg
	character*64 crdFile, fitFile, cfgFile, outFile
	integer i
	logical crdFileFlag, fitFileFlag, cfgFileFlag, outFileFlag

	crdFileFlag = .false.
	fitFileFlag = .false.
	cfgFileFlag = .false.
	outFileFlag = .false.

	i=1
	do
		call get_command_argument(i,arg)
		select case (arg)

		case ('-crd')
			i = i+1
			call get_command_argument(i,crdFile)
			crdFileFlag=.true.
			print*, "CRD File: ", crdFile
		case ('-fit')
			i = i+1
			call get_command_argument(i,fitFile)
			fitFileFlag=.true.
			print*, "Fit File: ", fitFile
		case ('-cfg')
			i = i+1
			call get_command_argument(i,cfgFile)
			cfgFileFlag=.true.
			print*, "cfg File: ", cfgFile
		case ('-out')
			i = i+1
			call get_command_argument(i,outFile)
			outFileFlag=.true.
			print*, "out File: ", outFile
		case default
			print*, 'Unrecognized command-line option: ', arg
			print*, 'Usage: map_3D_gr.x -crd [crd file]'
			print*, 'Usage: map_3D_gr.x -fit [fit file]'
			print*, 'Usage: map_3D_gr.x -cfg [cfg file]'
			print*, 'Usage: map_3D_gr.x -out [out file]'
			stop

		end select
		i = i+1
		if (i.ge.command_argument_count()) exit
	enddo

	if (crdFileFlag.eqv..false.) then
		write(*,*) "Must provide a crd file using command line argument -crd [crd file name]"
		stop
	endif

	if (fitFileFlag.eqv..false.) then
		write(*,*) "Must provide a fit file using command line argument -fit [fit file name]"
		stop
	endif

	if (cfgFileFlag.eqv..false.) then
		write(*,*) "Must provide a config file using command line argument -cfg [cfg file name]"
		stop
	endif

	if (outFileFlag.eqv..false.) then
		write(*,*) "Must provide a config file using command line argument -out [out file name]"
		stop
	endif

endsubroutine parse_command_line


! read coordinate information from file
subroutine read_crd(crdFile)
	use crdData
	implicit none
	character*64 crdFile
	integer i, j

	open(20,file=crdFile)

	! read first line and assign values. Then increment to the next line.
	read(20,999) nAtoms, nTypes
	write(*,*) "Number of Atoms:	", nAtoms
	write(*,*) "Number of Atom Types:	", nTypes

	allocate( aType(nAtoms), x(3,nAtoms) )

	! read the rest of the lines and make arrays of types and corresponding xyz-coordinates.
	do i = 1, nAtoms
		read(20,998) aType(i), (x(j,i),j=1,3)
!		write(*,998) aType(i), (x(j,i),j=1,3) ! for testing purposes
	enddo
	
	close(20)

	! Input Formats
999		format (i4, 1x, i4) ! first line of CRD file. 
998		format (i3,3(1x, f12.6)) ! atom lines of CRD file.

endsubroutine read_crd


! read fitting parameters from file
subroutine read_fit(fitFile)
	use crdData
	use fitData
	implicit none
	character*64 fitFile
	character*32 junk
	integer i

	nVars = 3 * nTypes

	open(20,file=fitFile)

	! allocate parameter-type and parameter-value arrays.
	allocate( pType(nVars), pValue(nVars) )

	do i = 1, nVars
		read(20,*) pType(i), junk, junk, pValue(i), junk
!		write(*,*) pType(i), pValue(i) ! for testing purposes
	enddo

	close(20)

endsubroutine read_fit


! read python cfg file for histogram values.
subroutine read_cfg(cfgFile)
	use cfgData
	implicit none
	character*64 cfgFile
	character*128 line
	character*32 firstWord, sep
	integer ios
	logical distMinFlag, distMaxFlag, distSizeFlag, angMinFlag, angMaxFlag, angSizeFlag

	distMinFlag = .false.
	distMaxFlag = .false.
	distSizeFlag = .false.
	angMinFlag = .false.
	angMaxFlag = .false.
	angSizeFlag = .false.
	
	ios = 0

	open(20,file=cfgFile)
	do while(ios>=0)
		read(20,'(a)',IOSTAT=ios) line
		call split(line,'=',firstWord, sep)
		if (line .ne. "") then
			if (firstWord .eq. "hist_dist_min") then
				read(line,*) hist_dist_min							! read the character 'line' and cast it as whatever format
				write(*,*) "Distance Minimum:	", hist_dist_min	! 'hist_dist_min' was initialized as which in this case is a
				distMinFlag = .true.								! dble.
			else if (firstWord .eq. "hist_dist_max") then
				read(line,*) hist_dist_max
				write(*,*) "Distance Maximum:	", hist_dist_max
				distMaxFlag = .true.
			else if (firstWord .eq. "bin_dist_size") then
				read(line,*) bin_dist_size
				write(*,*) "Distance Step Size:	", bin_dist_size
				distSizeFlag = .true.
			else if (firstWord .eq. "hist_ang_min") then
				read(line,*) hist_ang_min
				write(*,*) "Angle Minimum:		", hist_ang_min
				angMinFlag = .true.
			else if (firstWord .eq. "hist_ang_max") then
				read(line,*) hist_ang_max
				write(*,*) "Angle Maximum:		", hist_ang_max
				angMaxFlag = .true.
			else if (firstWord .eq. "bin_ang_size") then
				read(line,*) bin_ang_size
				write(*,*) "Angle Step Size:	", bin_ang_size
				angSizeFlag = .true.
			endif
		endif
	enddo
	close(20)

	if (distMinFlag.eqv..false.) then
		write(*,*) "Config file must have a 'hist_dist_min' value"
		stop
	endif
	if (distMaxFlag.eqv..false.) then
		write(*,*) "Config file must have a 'hist_dist_max' value"
		stop
	endif
	if (distSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'bin_dist_size' value"
		stop
	endif
	if (angMinFlag.eqv..false.) then
		write(*,*) "Config file must have a 'hist_ang_min' value"
		stop
	endif
	if (angMaxFlag.eqv..false.) then
		write(*,*) "Config file must have a 'hist_ang_max' value"
		stop
	endif
	if (angSizeFlag.eqv..false.) then
		write(*,*) "Config file must have a 'bin_ang_size' value"
		stop
	endif

endsubroutine read_cfg


! construct 3D grid of superimposed g(r) values.
subroutine map_gr(outFile)
	use crdData
	use fitData
	use cfgData
	use grData
	implicit none
	character*64 outFile
	integer num_dist_bins, a, i, j, k, l

	num_dist_bins = int( 2 * ( (hist_dist_max - hist_dist_min)/bin_dist_size ) + 1 )
	write(*,*) "Number of Distance Bins:", num_dist_bins

	allocate( gr(num_dist_bins,num_dist_bins,num_dist_bins), xAxis(num_dist_bins), yAxis(num_dist_bins), zAxis(num_dist_bins) )

	! Set equal to 1 because I'm multiplying the original/old values by the new ones. 0 would make everything 0.
	gr = 1

	write(6,*) "Generating 3D g(r) Array..."
	do a = 1, nAtoms
		! Pick the parameter values from pValue array corresponding to the atom type and save them as scalars.
		x1 = pValue( 3*(aType(a)-1) + 1 )
		x2 = pValue( 3*(aType(a)-1) + 2 )
		x0 = pValue( 3*(aType(a)-1) + 3 )

		! Set solute position vector
		do l = 1, 3
			R(l) = x(l,a)
		enddo

		! Assign values to 3D g(r) array.
		do k = 1, num_dist_bins
			rSolv(3) = ((k + 0.5) * bin_dist_size) + hist_dist_min

			do j = 1, num_dist_bins
				rSolv(2) = ((j + 0.5) * bin_dist_size) + hist_dist_min

				do i = 1, num_dist_bins
					rSolv(1) = ((i + 0.5) * bin_dist_size) + hist_dist_min

					! Calculate the g(r) value at 'rSolv' for atom a.
					call g_of_x(R, rSolv, x1, x2, x0, gx_out)

					! Multiply current value at this 'rSolv' by new value.
					gr(i,j,k) = gr(i,j,k) * gx_out

				enddo
			enddo
		enddo
	enddo

	write(6,*) "Generating Plot Axes..."
	do i = 1, num_dist_bins
		xAxis(i) = ((i + 0.5) * bin_dist_size) + hist_dist_min
		yAxis(i) = ((i + 0.5) * bin_dist_size) + hist_dist_min
		zAxis(i) = ((i + 0.5) * bin_dist_size) + hist_dist_min
	enddo

	! Write output file
	write(6,*) "Writing Output File:	", outFile
	open(25,file=outFile)

	! write first line with labels
	write(25,899) "#", "X", "Y", "Z", "g(r)"

	! write xyz-coordinates and cooresponding g(r) value.
	do k = 1, num_dist_bins
		do j = 1, num_dist_bins
			do i = 1, num_dist_bins
				!write(6,*) xAxis(i),yAxis(j),zAxis(k)
				write(25,898) xAxis(i), yAxis(j), zAxis(k), gr(i,j,k)
			enddo
		enddo
	enddo

	close(25)

	! Output Formats
899		format (a1,4(a12)) ! First line of output file, labeling the columns.
898		format (4(1x,f12.6)) ! XYZ-coordinates and g(r) value.

endsubroutine map_gr


! Functional form of parabolic potential around a given atom.
! R		= position vector of solute atom
! rSolv	= position vector of solvent atom
! x1	= distance at which parabola == 0 
! x2	= distance at which parabola == 1 and changes to flat value of 1
! x0	= distance at which parabola peaks
subroutine g_of_x(R, rSolv, x1, x2, x0, gx_out)
	implicit none
	double precision :: x(3), R(3), rSolv(3)	! Using the grData module would overwrite these variables and give a compilation
	double precision mag_x, x1, x2, x0, gx_out	! error. Initializing them here, outside of the module, makes it so that they are
	integer i									! separate variables and can be used concurrently in this nested subroutine.

	! Create vector 'x' that points from solute atom to solvent volume bin.
	do i = 1, 3
		x(i) = rSolv(i) - R(i)
	enddo

	mag_x = 0
	do i = 1, 3
		mag_x = mag_x + x(i)**2
	enddo
	mag_x = sqrt(mag_x)

	if (mag_x .le. x1) then
		gx_out = 0
	elseif ((mag_x .gt. x1) .and. (mag_x .lt. x2)) then
		gx_out = ( (x1 - x0)**2 - (mag_x - x0)**2 ) / ( (x1 - x0)**2 - (x2 - x0)**2 )
	elseif (mag_x .ge. x2) then
		gx_out = 1
	endif

endsubroutine g_of_x
