! Subroutines for calculating depth to groundwater table
module surfacewaterlevel
    
    use declaration
    
    contains
    
    ! calculates surface water above land surface
	subroutine calculateWaterlevel
	
        implicit none
		
		if (pressure(j,k,niz) > 0) then
			waterlevel(j,k) = pressure(j,k,niz)
		else
			waterlevel(j,k) = 0
		endif
        
    end subroutine
    
end module surfacewaterlevel

