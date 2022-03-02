program analysis
    
    ! load modules
    use declaration
    use IO
    use initialization
    use calculations
    
    implicit none
	
	logical :: file_exists
	    
    ! init setup for analysis
    call init
	  
    ! start calculations
    do ensrun=ensstart,ensend ! ENS loop
        
        ! current ensemble numbers
        write(ensname,"(I5.5)")ensrun-1
        write(ensnumber,"(I3.1)")ensrun-1
        
        ! update file names with current ensemble number
        call updateFileNames('ENSE',ensname)
        call updateFileNames('NUM',ensnumber)
        
        ! read in needed data
        call read_stationary_data
        
        ! reset file names
        call initCurrentFileNames
        
        ! data from first time step 
        irec = 1
        do i=istep,istep+nstep ! time loop
			            
            ! current timestep number
            write(nameruntime,"(I5.5)") i
            
            ! update file names with current timestep number
            call updateFileNames('ENSE',ensname)
            call updateFileNames('NUM',ensnumber)
            call updateFileNames('TIME',nameruntime)
			
				INQUIRE(FILE=namepress_curr, EXIST=file_exists)
				
				if (file_exists .and. doBoreholesDA .and. irec > 1) then
					if (dt==1) then
						INQUIRE(FILE=nameupdate_curr, EXIST=file_exists)
						else
						print *,'Please note that doBoreholesDA currently only works with timestepsize of 1 hour.'
						file_exists = 1==0
					endif
				endif
			
			if (file_exists) then
			
				! initialize temporally varying variables
				call initInstationaryData
				! read in needed data
				call read_instationary_data
				
				do j=1,nix ! x loop
				do k=1,niy ! y loop
					
					! if location part of subcatchment of interest
					if (mask(j,k,1) > 0) then
					
					! calculations that need to be done prior to XYZ calculations
					call PreXYZ
						
						do l=1,niz ! z loop
													
							! calculate everything that does also depend on Z
							call calculateXYZDependentData

						end do ! z loop
											
						! calculate everything that does not depend on Z
						call calculateXYDependentData					
						
						! calculations needed after z loop is finished
						call finalizeZLoop						

					endif		
					
				end do ! y loop
				end do ! x loop

				! save data
				call writeOutput	

				endif				
				
				! update data position
				irec = irec + 1

				! reset file names
				call initCurrentFileNames	
		
        end do ! time loop

        ! update file names
        call updateFileNames('ENSE',ensname)
        call updateFileNames('NUM',ensnumber)

        ! reset file names
        call initCurrentFileNames
        
    end do ! ENS loop
    
end program
