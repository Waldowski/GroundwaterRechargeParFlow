! I/O subroutines for reading and writing files for the Parflow analysis
module IO
    
    use declaration
    
    contains 
    
    ! reads in grid data from .pfb files
    subroutine pfb_read_grid(fname)
        
        implicit none
        
        ! declaration
        real*8 :: x1, y1, z1        ! lower boundary coordinates      
        character*200 :: fname      ! file name
        
        ! open file
        open(100,file=trim(adjustl(fname)),form='unformatted',access='stream',convert='BIG_ENDIAN',status='old')
        
        ! Start: reading of domain spatial information
        read(100) x1 !X
        read(100) y1 !Y  
        read(100) z1 !Z

        read(100) nix !NX 
        read(100) niy !NY 
        read(100) niz !NZ 

        read(100) meandx !DX 
        read(100) meandy !DY 
        read(100) meandz !DZ 
        
        print *,trim(adjustl(fname)),': nx,ny,nz ',nix,niy,niz
        ! End: reading of domain spatial information
        
        ! close file
        close(100)
        
    end subroutine
    
	subroutine const_value(matrix,const,nx,ny,nz)
	
	implicit none
	
	! declaration
	integer*4 :: nx, ny, nz                                ! number of grid cells
	real*8 :: const                                        ! const value to save
	real*8 :: matrix(nx,ny,nz)                             ! matrix where constant value is stored
	
	! save const value into 3d matrix
		do  k= 1, nz 
			do  j= 1 , ny 
				do  i= 1 , nx 
					matrix(i,j,k) = const
				end do
			end do
		end do
		
	end subroutine
	
    ! reads in data from .pfb files
    subroutine pfb_read(value,fname) 
        
        implicit none
        
        ! declaration
        real*8 :: value(nix,niy,niz)                           ! variable to be read in
        integer*4 :: nx, ny, nz                                ! number of grid cells
        real*8 :: dx, dy, dz                                   ! size of grid cells
        real*8 :: x1, y1, z1                                   ! lower boundary coordinates                      
        integer*4 :: i, j, k, is                               ! loop indices
        integer*4 :: ns                                        ! number subgrids
        integer*4 :: ix, iy, iz, rx, ry, rz, nnx, nny, nnz     ! subgrid information    
        character*200 :: fname                                 ! file name
        
        ! print file name
        print *,trim(adjustl(fname))    
        
        ! open file
        open(100,file=trim(adjustl(fname)),form='unformatted',access='stream',convert='BIG_ENDIAN',status='old')

        ! Start: reading of domain spatial information
        read(100) x1 !X
        read(100) y1 !Y  
        read(100) z1 !Z

        read(100) nx !NX 
        read(100) ny !NY 
        read(100) nz !NZ 

        read(100) dx !DX 
        read(100) dy !DY 
        read(100) dz !DZ 
        
        read(100) ns !num_subgrids
        ! End: reading of domain spatial information

        ! Start: loop over number of sub grids
        do is = 0, (ns-1)

            ! Start: reading of sub-grid spatial information
            read(100) ix
            read(100) iy
            read(100) iz

            read(100) nnx 
            read(100) nny 
            read(100) nnz 
            
            read(100) rx
            read(100) ry
            read(100) rz
            ! End: reading of sub-grid spatial information

            ! read in data from each individual subgrid
            do  k=iz +1 , iz + nnz 
                do  j=iy +1 , iy + nny 
                    do  i=ix +1 , ix + nnx 
                        read(100) value(i,j,k)
                    end do
                end do
            end do
            
        end do
        ! End: loop over number of sub grids
  
        ! close file
        close(100)
        
    end subroutine
    
    ! reads in all needed data that is constant over simulation time
    subroutine read_stationary_data
        
        implicit none

			call pfb_read(mask,namemask_curr)
     	
        if (doGrossRecharge .or. doNetRecharge) then
		
			if (constsx) then
				call const_value(slopex,sx_const,nix,niy,1)
			else
				call pfb_read(slopex,nameslopex_curr)
			endif
			
			if (constsy) then
				call const_value(slopey,sy_const,nix,niy,1)
			else
				call pfb_read(slopey,nameslopey_curr)
			endif	
			if (constporo) then
				call const_value(porosity,poro_const,nix,niy,niz)
			else
				call pfb_read(porosity,nameporo_curr)
			endif
			if (constspecstor) then
				call const_value(specstor,specstor_const,nix,niy,niz)
			else
				call pfb_read(specstor,namespecstor_curr)
			endif
			
        endif
		if (doGrossRecharge .or. doNetRecharge .or. doFluxCrossing .or. doVirtualLysimeter) then
		
			if (constksat) then
				call const_value(ksat,ksat_const,nix,niy,niz)
			else
				call pfb_read(ksat,nameksat_curr)			
			endif
			if (constalpha) then
				call const_value(alpha,alpha_const,nix,niy,niz)
			else
				call pfb_read(alpha,namealpha_curr) 				
			endif
			if (constn) then
				call const_value(n,n_const,nix,niy,niz)
			else
				call pfb_read(n,namen_curr)				
			endif
			if (constmaskgwr) then
				call const_value(mask_gwr, maskgwr_const,nix,niy,1)
			else
				call pfb_read(mask_gwr,namemaskgwr_curr)
			endif
			
		endif
        
    end subroutine
    
    ! reads in all needed data that changes with time
    subroutine read_instationary_data
        
        implicit none
        
        call pfb_read(pressure,namepress_curr) 
        if (doGrossRecharge .or. doNetRecharge) then
            call pfb_read(saturation,namesatur_curr) 
        endif  
        
    end subroutine
  
    ! writes data to binary files
    subroutine write_binary(value,fname,nx,ny,nz,irec)
      
        implicit none
      
        ! declaration
        real*8 :: value(nx,ny,nz)       ! data to be saved
        integer*4 :: nx, ny, nz         ! number of cells 
        integer*4 :: reclen, irec       ! size and position of data to be saved    
        character*200 :: fname          ! file name
      
        ! check data size
        reclen = nx*ny*nz
        inquire(iolength=reclen) value
        
        ! open file
        open(200,file=trim(adjustl(fname)),access='direct',recl=reclen,form='unformatted')
        
        ! write data to file
        write(200,rec=irec) value
        
        ! close file
        close(200)
        
    end subroutine
    
    ! writes scalar data to binary files
    subroutine write_binary_scalar(value,fname,irec)
      
        implicit none
      
        ! declaration
        real*8 :: value                 ! data to be saved
        integer*4 :: reclen, irec       ! size and position of data to be saved    
        character*200 :: fname          ! file name
      
        ! check data size
        reclen = 1
        inquire(iolength=reclen) value
        
        ! open file
        open(200,file=trim(adjustl(fname)),access='direct',recl=reclen,form='unformatted')
        
        ! write data to file
        write(200,rec=irec) value
        
        ! close file
        close(200)
        
    end subroutine
    
    ! writes all calculated data into output files
    subroutine writeOutput
        
        implicit none 
		
		if (doGWTable) then
            call write_binary(watertab,trim(nameoutputpath_curr)//'watertab_'//&
            nameofcase_curr,nix,niy,1,irec)
        endif
		        
        if (doGrossRecharge .and. irec > 1) then
            call write_binary(recharge_gross,trim(nameoutputpath_curr)//'recharge_gross_'//&
			nameofcase_curr,nix,niy,1,irec-1)
        endif
		
        if (doNetRecharge .and. irec > 1) then
            call write_binary(recharge_net,trim(nameoutputpath_curr)//'recharge_net_'//&
			nameofcase_curr,nix,niy,1,irec-1)
        endif
		
		if (doFluxCrossing .and. irec > 1) then
			call write_binary(recharge_crossing,trim(nameoutputpath_curr)//'flux_crossing_'//&
			nameofcase_curr,nix,niy,1,irec-1)
		endif
		
		if (doRechargeSources .and. doNetRecharge .and. irec > 1) then
			call write_binary(recharge_net_sat,trim(nameoutputpath_curr)//'recharge_loc_pond_net_'//& ! recharge at local ponding spots
			nameofcase_curr,nix,niy,1,irec-1)
			call write_binary(riv_recharge_net,trim(nameoutputpath_curr)//'recharge_river_net_'//& ! recharge at rivers
			nameofcase_curr,nix,niy,1,irec-1)
			call write_binary(gwcomp,trim(nameoutputpath_curr)//'recharge_ss_component_'//& ! specific storage component of recharge
			nameofcase_curr,nix,niy,1,irec-1)
		endif
		
		if (doRechargeSources .and. doGrossRecharge .and. irec > 1) then
			call write_binary(recharge_gross_sat,trim(nameoutputpath_curr)//'recharge_loc_pond_gross_'//&
			nameofcase_curr,nix,niy,1,irec-1)
			call write_binary(surf_wat_exch,trim(nameoutputpath_curr)//'surf_wat_exch_flux_'//&
			nameofcase_curr,nix,niy,1,irec-1)
			call write_binary(riv_recharge_gross,trim(nameoutputpath_curr)//'recharge_river_gross_'//&
			nameofcase_curr,nix,niy,1,irec-1)
		endif
		
		if (doVirtualLysimeter .and. irec > 1) then
            call write_binary(vl,trim(nameoutputpath_curr)//'vl_'//&
			nameofcase_curr,nix,niy,1,irec-1)
        endif
		
    end subroutine
       
    ! reads in user input from .txt file
    subroutine read_userinput(fname)
        
        implicit none
        
        ! declaration
        character*200 :: fname          ! file name
        logical :: doRead = .true.      ! continue reading?
        real*8 :: val                   ! current variable value (numeric))
        character*200 :: valChar        ! current variable value (string) 
        character*20 :: var             ! current variable name
        integer*4 :: i                  ! index for vardz loop
        integer*4 :: n                  ! status of reading
        
        ! open file
        open(300,file=trim(adjustl(fname)),status='old')
        
        ! read until end of file
        do while(doRead)
            
            ! read line
            read(300,*,iostat=n) var, val
            ! no numeric value?
            if (index(var,'name') > 0 .and. index(var,'#') == 0) then
                backspace(300)
                read(300,'(A)',iostat=n) valChar
                var = valChar(1:index(valChar,' ')-1)
                valChar = adjustl(valChar(index(valChar,' '):200))
                valChar = valChar(1:index(valChar,' ')-1)
            endif
            
            ! assign current variable
            if (var == 'istep') then
                istep = int(val)
            else if (var == 'nstep') then
                nstep = int(val)
            else if (var == 'dt') then
                dt = val
            else if (var == 'nvardz') then
                allocate(vardz(int(val)))
                ! loop over all vardz entries
                do i=1,int(val)
                    read(300,*) var, val
                    vardz(i) = val
                end do
			else if (var == 'rootend') then
                rootend = int(val)
            else if (var == 'ensstart') then
                ensstart = int(val)
            else if (var == 'ensend') then
                ensend = int(val)
            else if (var == 'doGrossRecharge') then
                doGrossRecharge = int(val) /= 0
            else if (var == 'doNetRecharge') then
                doNetRecharge = int(val) /= 0
			else if (var == 'doFluxCrossing') then
                doFluxCrossing = int(val) /= 0
			else if (var == 'doRechargeSources') then
                doRechargeSources = int(val) /= 0
			else if (var == 'doVirtualLysimeter') then
                doVirtualLysimeter = int(val) /= 0
            else if (var == 'doGWTable') then
                doGWTable = int(val) /= 0
			else if (var == 'constsx') then
				constsx = int(val) /= 0
			else if (var == 'constsy') then
				constsy = int(val) /= 0
			else if (var == 'constporo') then
				constporo = int(val) /= 0
			else if (var == 'constksat') then
				constksat = int(val) /= 0
			else if (var == 'constalpha') then
				constalpha = int(val) /= 0
			else if (var == 'constn') then
				constn = int(val) /= 0
			else if (var == 'constmaskgwr') then
				constmaskgwr = int(val) /= 0
			else if (var == 'constspecstor') then
				constspecstor = int(val) /= 0			
            else if (var == 'namesatur') then
                namesatur = valChar
            else if (var == 'namepress') then
                namepress = valChar
            else if (var == 'nameslopex') then
				nameslopex = valChar
            else if (var == 'nameslopey') then
                nameslopey = valChar
            else if (var == 'nameporo') then
                nameporo = valChar
            else if (var == 'nameksat') then
                nameksat = valChar
            else if (var == 'namealpha') then
                namealpha = valChar
            else if (var == 'namen') then
                namen = valChar
            else if (var == 'namemask') then
                namemask = valChar
            else if (var == 'namemaskgwr') then
                namemaskgwr = valChar
            else if (var == 'namespecstor') then
                namespecstor = valChar
            else if (var == 'nameoutputpath') then
                nameoutputpath = valChar
			else if (var == 'nameofcase') then
                nameofcase = valChar
			else if (var == 'sx_const') then
				sx_const = val
			else if (var == 'sy_const') then
				sy_const = val
			else if (var == 'poro_const') then
				poro_const = val
			else if (var == 'ksat_const') then
				ksat_const = val
			else if (var == 'alpha_const') then
				alpha_const = val
			else if (var == 'n_const') then
				n_const = val
			else if (var == 'maskgwr_const') then
				maskgwr_const = val
			else if (var == 'specstor_const') then
				specstor_const = val	
            else if (scan(var,'#') .gt. 0) then
                continue
            else if (var == 'end') then
                doRead = .false.
            else
                print *,'Variable ', var,' does not exist!'
            endif
            
        end do
        
    end subroutine
    
end module IO

