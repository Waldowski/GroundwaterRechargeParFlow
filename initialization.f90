! Init and update subroutines for Parflow analysis
module initialization
    
    use declaration
    use IO
    
    contains
    
    ! allocates all needed arrays
    subroutine allocateArrays
        
        implicit none
        
        allocate(mask(nix,niy,niz))
        allocate(pressure(nix,niy,niz))
        
		if (doGrossRecharge) then
			allocate(riv_recharge_gross(nix,niy))
			allocate(recharge_gross_sat(nix,niy))
			allocate(recharge_gross(nix,niy))
         endif
        if (doNetRecharge) then
			allocate(riv_recharge_net(nix,niy))
			allocate(recharge_net_sat(nix,niy))
			allocate(gwcomp(nix,niy))
			allocate(recharge_net(nix,niy))	
        endif
		if (doVirtualLysimeter) then
			allocate(vl(nix,niy))
		endif	
        if (doGrossRecharge .or. doNetRecharge) then
			allocate(cumrech_stor(nix,niy))
			allocate(cumrech_stor_old(nix,niy))
			allocate(watertab_height_init(nix,niy))
			allocate(watertab_height(nix,niy))
			allocate(theta_init(nix,niy,niz))
			allocate(theta_watertab_init(nix,niy))
			allocate(theta_watertab(nix,niy))
			allocate(recharge_wt(nix,niy))
			allocate(watertab_cell_init(nix,niy))
			allocate(saturation(nix,niy,niz))
            allocate(porosity(nix,niy,niz))
        endif 
		if (doGWTable .or. doGrossRecharge .or. doNetRecharge) then
            allocate(watertab(nix,niy))
			allocate(watertab_switch(nix,niy))
		endif
		if (doGrossRecharge .or. doFluxCrossing) then
			allocate(recharge_crossing(nix,niy))
			allocate(surf_wat_exch(nix,niy))
			allocate(topcomp(nix,niy))
		endif	
        if (doGrossRecharge .or. doNetRecharge .or. doFluxCrossing) then
            allocate(slopex(nix,niy,1))
            allocate(slopey(nix,niy,1))
        endif
		if (doGrossRecharge .or. doNetRecharge .or. doFluxCrossing .or. doVirtualLysimeter) then
			allocate(mask_gwr(nix,niy))
			allocate(specstor(nix,niy,niz))
			allocate(p_old(nix,niy,niz))
			allocate(ksat(nix,niy,niz))
            allocate(alpha(nix,niy,niz))
            allocate(n(nix,niy,niz))
		endif
		
    end subroutine
    
    ! replaces part of the filename with another string
    subroutine updateFileNames(str_old,str_new)
        
        implicit none
        
        integer*4 :: pos                     ! position of old string
        character(len=*) :: str_old, str_new ! strings to be switched
        
        pos = index(namesatur_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namesatur_curr = namesatur_curr(1:pos-1)//trim(adjustl(str_new))// &
                namesatur_curr(pos+len(trim(adjustl(str_old))):len(namesatur_curr))
			pos = index(namesatur_curr,trim(adjustl(str_old)))
        end do
        pos = index(namepress_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namepress_curr = namepress_curr(1:pos-1)//trim(adjustl(str_new))// &
                namepress_curr(pos+len(trim(adjustl(str_old))):len(namepress_curr))
			pos = index(namepress_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameslopex_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameslopex_curr = nameslopex_curr(1:pos-1)//trim(adjustl(str_new))// & 
                nameslopex_curr(pos+len(trim(adjustl(str_old))):len(nameslopex_curr))
			pos = index(nameslopex_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameslopey_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameslopey_curr = nameslopey_curr(1:pos-1)//trim(adjustl(str_new))// &
                nameslopey_curr(pos+len(trim(adjustl(str_old))):len(nameslopey_curr))
			pos = index(nameslopey_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameporo_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameporo_curr = nameporo_curr(1:pos-1)//trim(adjustl(str_new))// &
                nameporo_curr(pos+len(trim(adjustl(str_old))):len(nameporo_curr))
			pos = index(nameporo_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameksat_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameksat_curr = nameksat_curr(1:pos-1)//trim(adjustl(str_new))// &
                nameksat_curr(pos+len(trim(adjustl(str_old))):len(nameksat_curr))
			pos = index(nameksat_curr,trim(adjustl(str_old)))
        end do
        pos = index(namealpha_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namealpha_curr = namealpha_curr(1:pos-1)//trim(adjustl(str_new))// &
                namealpha_curr(pos+len(trim(adjustl(str_old))):len(namealpha_curr))
				pos = index(namealpha_curr,trim(adjustl(str_old)))
        end do
        pos = index(namen_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namen_curr = namen_curr(1:pos-1)//trim(adjustl(str_new))// &
                namen_curr(pos+len(trim(adjustl(str_old))):len(namen_curr))
				pos = index(namen_curr,trim(adjustl(str_old)))
        end do
        pos = index(namemask_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namemask_curr = namemask_curr(1:pos-1)//trim(adjustl(str_new))// &
                namemask_curr(pos+len(trim(adjustl(str_old))):len(namemask_curr))
				pos = index(namemask_curr,trim(adjustl(str_old)))
        end do
        pos = index(namemaskgwr_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namemaskgwr_curr = namemaskgwr_curr(1:pos-1)//trim(adjustl(str_new))// &
                namemaskgwr_curr(pos+len(trim(adjustl(str_old))):len(namemaskgwr_curr))
				pos = index(namemaskgwr_curr,trim(adjustl(str_old)))
        end do
        pos = index(namespecstor_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            namespecstor_curr = namespecstor_curr(1:pos-1)//trim(adjustl(str_new))// &
                namespecstor_curr(pos+len(trim(adjustl(str_old))):len(namespecstor_curr))
				pos = index(namespecstor_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameoutputpath_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameoutputpath_curr = nameoutputpath_curr(1:pos-1)//trim(adjustl(str_new))// &
                nameoutputpath_curr(pos+len(trim(adjustl(str_old))):len(nameoutputpath_curr))
				pos = index(nameoutputpath_curr,trim(adjustl(str_old)))
        end do
        pos = index(nameofcase_curr,trim(adjustl(str_old)))
        do while (pos > 0)
            nameofcase_curr = nameofcase_curr(1:pos-1)//trim(adjustl(str_new))// &
                nameofcase_curr(pos+len(trim(adjustl(str_old))):len(nameofcase_curr))
				pos = index(nameofcase_curr,trim(adjustl(str_old)))
        end do		
        
    end subroutine
    
    ! initialization of current file names
    subroutine initCurrentFileNames
        
        implicit none
        
        namesatur_curr = namesatur
        namepress_curr = namepress
        nameslopex_curr = nameslopex
        nameslopey_curr = nameslopey
        nameporo_curr = nameporo
        nameksat_curr = nameksat
        namealpha_curr = namealpha
        namen_curr = namen
        namemask_curr = namemask
        namemaskgwr_curr = namemaskgwr
        namespecstor_curr = namespecstor
        nameoutputpath_curr = nameoutputpath
		nameofcase_curr = nameofcase
        
    end subroutine
    
    ! initialization of anaylsis setup
    subroutine init
        
        implicit none
        
        ! input file
        narg = command_argument_count()
        if (narg > 0) then
            do arg=1,narg
                call get_command_argument(arg,argName)
                nameinput = argName
            end do
        else
            nameinput = ''
        endif

        ! read in user input
        call read_userinput(nameinput)

        ! init file names
        call initCurrentFileNames

        ! read in grid data
        write(ensname,"(I5.5)") 0
        call updateFileNames('ENSE',ensname)
        call pfb_read_grid(namemask_curr) 
        call initCurrentFileNames

        ! update adaptive grid sizes
        if (size(vardz) /= niz) then
            print *,'Entries for vardz do not match number of grid cells!'
        endif
        vardz = vardz * meandz

        ! allocation
        call allocateArrays
        
    end subroutine
       
    ! initialization of instationary output data
    subroutine initInstationaryData
        
        implicit none
    
        ! initialize temporally varying variables
        if (doGWTable .or. doNetRecharge .or. doGrossRecharge) then
            watertab(:,:)=-9999999.0
		! reset watertab_switch
			watertab_switch(:,:) = int(1) /= 0
        endif
		! initialize lateral fluxes
		if (doNetRecharge .or. doGrossRecharge) then ! for recharge
			f_l_sum = 0.0
			f_r_sum = 0.0
			f_b_sum = 0.0
			f_f_sum = 0.0	
			f_l = 0.0
            f_r = 0.0
            f_b = 0.0
            f_f = 0.0
			f_l_out = 0.0
            f_r_out = 0.0
            f_b_out = 0.0
            f_f_out = 0.0
		endif

        ! get pressure from last time step
        if (i > istep .and. doNetRecharge .or. i > istep .and. doGrossRecharge) then
            p_old = pressure
        endif
            
    end subroutine
    
end module initialization

