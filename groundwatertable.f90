! Subroutines for calculating depth to groundwater table
module groundwatertable
    
    use declaration
    
    contains
    
    ! calculates depth to groundwater table
    ! determine groundwatertable as location at which the pressure == 0 via linear interpolation
    subroutine depthToGroundwater
        
        implicit none
        if (l /= niz) then
            if (pressure(j,k,niz-l)/pressure(j,k,niz+1-l) < 0) then
                watertab(j,k) = sum(vardz(1:l-1)) + vardz(l)/2 - &
                pressure(j,k,niz-l+1)*(vardz(l)/2 + vardz(l+1)/2)/(pressure(j,k,niz-l)-pressure(j,k,niz-l+1))
				watertab_switch(j,k) = int(0) /= 0
            endif
		else if (l == niz .and. watertab_switch(j,k)) then
			watertab(j,k) = 0
        endif

        
    end subroutine
    
end module groundwatertable

