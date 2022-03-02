! Subroutines for calculating overland flow
module runoff
    
    use declaration
    
    contains
    
    ! calculates runoff
    subroutine calculateRunoff

        implicit none

        if (slopey(j,k,1) > 0 .and. pressure(j,k,niz) > 0) then
            runoffy(j,k) =  (((sqrt(abs(slopey(j,k,1)))/mannings(j,k,1))*(pressure(j,k,niz)**(5.0/3.0)))*meandx)
			runoff2d(j,k) = runoff2d(j,k) + runoffy(j,k)
			! how much surface water leaves the start of the domain in y direction?
			if (k==1) then
				runoff_leaving(j,k) = runoff_leaving(j,k) + runoffy(j,k)
			endif			
        endif
		if (slopey(j,k,1) < 0 .and. pressure(j,k,niz) > 0) then
            runoffy(j,k) =  (((sqrt(abs(slopey(j,k,1)))/mannings(j,k,1))*(pressure(j,k,niz)**(5.0/3.0)))*meandx)
			runoff2d(j,k) = runoff2d(j,k) + runoffy(j,k)
			! how much surface water leaves the end of the domain in y direction?
			if (k==niy) then
				runoff_leaving(j,k) = runoff_leaving(j,k) + runoffy(j,k)
			endif	
        endif
		
        
        if (slopex(j,k,1) > 0 .and. pressure(j,k,niz) > 0) then
            runoffx(j,k) = (((sqrt(abs(slopex(j,k,1)))/mannings(j,k,1))*(pressure(j,k,niz)**(5.0/3.0)))*meandy)
			runoff2d(j,k) = runoff2d(j,k) + runoffx(j,k)
			! how much surface water leaves the start of the domain in x direction?
			if (j==1) then
				runoff_leaving(j,k) = runoff_leaving(j,k) + runoffx(j,k)
			endif
        endif
		
        if (slopex(j,k,1) < 0 .and. pressure(j,k,niz) > 0) then
            runoffx(j,k) = (((sqrt(abs(slopex(j,k,1)))/mannings(j,k,1))*(pressure(j,k,niz)**(5.0/3.0)))*meandy)
			runoff2d(j,k) = runoff2d(j,k) + runoffx(j,k)
			! how much surface water leaves the end of the domain in x direction?
			if (j==nix) then
				runoff_leaving(j,k) = runoff_leaving(j,k) + runoffx(j,k)
			endif
        endif
        
			
		
    end subroutine
    
end module runoff

