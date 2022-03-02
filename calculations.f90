module calculations
    
    use declaration
	use lat_flux
    use gw_recharge
    use groundwatertable
    use runoff
	use surfacewaterlevel
    
    contains
    
    ! calculates everything that does not depend on Z
    subroutine calculateXYDependentData
        
        implicit none
		
        if (doGrossRecharge .or. doNetRecharge) then
			! calculate watertable height
			watertab_height(j,k) = sum(vardz(1:niz))-watertab(j,k)
			if (i == istep) then
				watertab_height_init(j,k) = watertab_height(j,k)
			endif
        endif
		
		! calculate recharge as vertical flux
        if (doFluxCrossing .and. i > istep .or. doGrossRecharge .and. i > istep) then
			call verticalFlux
        endif
		
		! calculate flux below roots
		if (doVirtualLysimeter) then
			call VirtualLysimeter
		endif

    end subroutine
    
    ! calculates everything that does also depend on Z
    subroutine calculateXYZDependentData
        
        implicit none
		        
        ! calculate lateral fluxes in groundwater and sum them up
		if (i > istep .and. doNetRecharge) then
            call lateralFluxesSum
        endif
        
        ! determine depth to groundwater table
        if (doNetRecharge .or. doGrossRecharge) then
            call depthToGroundwater
        endif
		
		if (i == istep .and. doNetRecharge .or. i == istep .and. doGrossRecharge) then
			theta_init(j,k,l) = saturation(j,k,l) * porosity(j,k,l)
		endif
		
    end subroutine
	
	subroutine PreXYZ
	
		if (doGrossRecharge .or. doNetRecharge .or. doFluxCrossing) then
			call WatertabCell ! get positions of last groundwater cell
		endif

	end subroutine

    ! all calculations needed after z loop is finished
    subroutine finalizeZLoop
        
        implicit none
		
		! Groundwater recharge
        				
		if (doNetRecharge .or. doGrossRecharge) then
			call localgwstorchange ! saves storage change as recharge_wt
			call latgwdis ! lateral component of groundwater discharge
		endif
		if (i > istep) then
			if (doGrossRecharge .or. doFluxCrossing) then
			surf_wat_exch(j,k) = 0.0 ! reset value
				if (pressure(j,k,niz) >=0) then
					call lateralfluxes(niz)
					call TopverticalFlux ! saves top vertical flux as surf_wat_exch
					call Toplayerspecificstorage
					surf_wat_exch(j,k) =  surf_wat_exch(j,k) + (f_r - f_l + f_f - f_b)/(meandx*meandy) + topcomp(j,k)
					if (i_sat_new == niz) then !one consistent waterbody
						recharge_crossing(j,k) = surf_wat_exch(j,k)
					endif
				endif
			endif
			if (doGrossRecharge) then
				recharge_gross_sat(j,k) = 0.0 ! reset value
				recharge_gross(j,k) = recharge_wt(j,k) + recharge_crossing(j,k) + &  ! add flux which crosses the water table to groundwater storage changes and also account for lateral gw-discharge
				(-f_r_out + f_l_out - f_f_out + f_b_out)/(meandx*meandy)
				recharge_gross_unsat(j,k) = recharge_gross(j,k)
				if (i_sat_new == niz) then ! one consistent waterbody
					recharge_gross_unsat(j,k) = 0.0
					recharge_gross_sat(j,k) = surf_wat_exch(j,k) + &
					(-f_r_out + f_l_out - f_f_out + f_b_out)/(meandx*meandy)
				endif
				! mask out rivers for gross recharge
				if (mask_gwr(j,k) == 0) then ! river
					if (i_sat_new < niz) then
						riv_recharge_gross_unsat(j,k) = recharge_gross(j,k)
					else
						riv_recharge_gross_unsat(j,k) = 0.0
					endif
					riv_recharge_gross(j,k) = recharge_gross(j,k)
					riv_exch_flux(j,k) = surf_wat_exch(j,k)	
					recharge_gross_sat(j,k) = 0.0
					recharge_gross_unsat(j,k) = 0.0
					recharge_crossing(j,k) = 0.0
					recharge_gross(j,k) = 0.0
				endif
			endif
			! mask out rivers for flux crossing
			if (doFluxCrossing) then
				if(mask_gwr(j,k) == 0) then ! river	
					recharge_crossing(j,k) = 0.0
				endif
			endif					
			if (doNetRecharge) then
				call Groundwaterspecificstorage
				recharge_net(j,k) = recharge_wt(j,k) + (f_r_sum- f_r_out - f_l_sum + f_l_out + f_f_sum - f_f_out - f_b_sum + f_b_out) / & ! add lateral gw outflow + Ss term to recharge estimate based on watertable fluctuations - also account for lateral gw-discharge
				(meandx*meandy) + gwcomp(j,k)	
				recharge_net_unsat(j,k) = recharge_net(j,k)
				recharge_net_sat(j,k) = 0.0 ! reset value
				if (pressure(j,k,niz) >=0) then
					recharge_net_sat(j,k) = recharge_net(j,k) - recharge_wt(j,k) ! recharge_wt is not a part of local ponding
					if (i_sat_new == niz) then !one consistent waterbody
						recharge_net_unsat(j,k) = 0.0
					endif
				endif
				! mask out rivers
				if (mask_gwr(j,k) == 0) then ! river
					if (i_sat_new < niz) then
						riv_recharge_net_unsat(j,k) = recharge_net(j,k)
					else
						riv_recharge_net_unsat(j,k) = 0.0
					endif
					riv_recharge_net(j,k) = recharge_net(j,k)
					recharge_net_sat(j,k) = 0.0
					recharge_net_unsat(j,k) = 0.0
					recharge_net(j,k) = 0.0
				endif
			endif					
		endif

		! reset summed up lateral fluxes
		f_l_sum = 0.0
		f_r_sum = 0.0
		f_b_sum = 0.0
		f_f_sum = 0.0
		f_l_out = 0.0
		f_r_out = 0.0
		f_b_out = 0.0
		f_f_out = 0.0	
		
	end subroutine
	
    
end module calculations
