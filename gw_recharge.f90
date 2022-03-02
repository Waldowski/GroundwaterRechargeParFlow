! Subroutines and functions needed for groundwater recharge calculation
module gw_recharge
    
    use declaration
	use groundwatertable
	use lat_flux
    
    contains		
	
	! calculates lateral groundwater discharge
	subroutine latgwdis
	
		implicit none
		
		integer*4 :: m                                                  ! loop index
		integer*4 :: i_sat_right,i_sat_left,i_sat_front,i_sat_behind    ! water table cells of neighbours
		real*8 :: gc_n, g_n                                             ! averaged slopes
        real*8 :: kr_i, kr_n                                            ! relative hydraulic conductivies
        real*8 :: up, down                                              ! up and downstream gradient
	
	! Only calculate water table cell of neighbour, if there is a neighbour	
	if (j < nix) then
		i_sat_right = findWaterTable(pressure(j+1,k,:))
	else
		i_sat_right = i_sat_new
	endif
	if (j > 1) then
		i_sat_left = findWaterTable(pressure(j-1,k,:))
	else
		i_sat_left = i_sat_new
	endif
	if (k < niy) then
		i_sat_front = findWaterTable(pressure(j,k+1,:))
	else
		i_sat_front = i_sat_new
	endif
	if (k > 1) then
		i_sat_behind = findWaterTable(pressure(j,k-1,:))
	else
		i_sat_behind = i_sat_new
	endif
	
	! right
	if (i_sat_new .ne. i_sat_right) then
		do m = min(i_sat_new,i_sat_right)+1,max(i_sat_new,i_sat_right) ! loop for every cell in z direction, which is unsaturated in neighbour, but saturated in current cell
			gc_n = g*0.5*(cos(atan(slopex(j,k,1))) + cos(atan(slopex(j+1,k,1))))
			g_n = g*0.5*(sin(atan(slopex(j,k,1))) + sin(atan(slopex(j+1,k,1))))
			kr_i = relativeConductivity(alpha(j,k,m),pressure(j,k,m),n(j,k,m),(1-1/n(j,k,m)))
			kr_n = relativeConductivity(alpha(j+1,k,m),pressure(j+1,k,m),n(j+1,k,m),(1-1/n(j+1,k,m)))
			up = (pressure(j,k,m)-pressure(j+1,k,m))/meandx * gc_n - g_n 
			down = 0.
			f_r_out = f_r_out + meandy*vardz(niz-m+1)*harmmean(ksat(j,k,m),ksat(j+1,k,m)) * up * &
				upstreammean(up,down,kr_i*rho,kr_n*rho)/mu
		end do
	endif
	
	! left 
	if (i_sat_new .ne. i_sat_left) then
		do m = min(i_sat_new,i_sat_left)+1,max(i_sat_new,i_sat_left) ! loop for every cell in z direction, which is unsaturated in neighbour, but saturated in current cell
			gc_n = g*0.5*(cos(atan(slopex(j,k,1))) + cos(atan(slopex(j-1,k,1))))
			g_n = g*0.5*(sin(atan(slopex(j,k,1))) + sin(atan(slopex(j-1,k,1))))
			kr_i = relativeConductivity(alpha(j,k,m),pressure(j,k,m),n(j,k,m),(1-1/n(j,k,m)))
			kr_n = relativeConductivity(alpha(j-1,k,m),pressure(j-1,k,m),n(j-1,k,m),(1-1/n(j-1,k,m)))
			up = (pressure(j-1,k,m)-pressure(j,k,m))/meandx * gc_n - g_n 
			down = 0.
			f_l_out = f_l_out + meandy*vardz(niz-m+1)*harmmean(ksat(j,k,m),ksat(j-1,k,m)) * up * &
				upstreammean(up,down,kr_n*rho,kr_i*rho)/mu
		end do	
	endif
	
	! front
	if (i_sat_new .ne. i_sat_front) then
		do m = min(i_sat_new,i_sat_front)+1,max(i_sat_new,i_sat_front) ! loop for every cell in z direction, which is unsaturated in neighbour, but saturated in current cell
			gc_n = g*0.5*(cos(atan(slopey(j,k,1))) + cos(atan(slopey(j,k+1,1))))
			g_n = g*0.5*(sin(atan(slopey(j,k,1))) + sin(atan(slopey(j,k+1,1))))
			kr_i = relativeConductivity(alpha(j,k,m),pressure(j,k,m),n(j,k,m),(1-1/n(j,k,m)))
			kr_n = relativeConductivity(alpha(j,k+1,m),pressure(j,k+1,m),n(j,k+1,m),(1-1/n(j,k+1,m)))
			up = (pressure(j,k,m)-pressure(j,k+1,m))/meandy * gc_n - g_n 
			down = 0.
			f_f_out = f_f_out + meandx*vardz(niz-m+1)*harmmean(ksat(j,k,m),ksat(j,k+1,m)) * up * &
				upstreammean(up,down,kr_i*rho,kr_n*rho)/mu
		end do
	endif
	
	! back
	if (i_sat_new .ne. i_sat_behind) then
		do m = min(i_sat_new,i_sat_behind)+1,max(i_sat_new,i_sat_behind) ! loop for every cell in z direction, which is unsaturated in neighbour, but saturated in current cell
			gc_n = g*0.5*(cos(atan(slopey(j,k,1))) + cos(atan(slopey(j,k-1,1))))
			g_n = g*0.5*(sin(atan(slopey(j,k,1))) + sin(atan(slopey(j,k-1,1))))
			kr_i = relativeConductivity(alpha(j,k,m),pressure(j,k,m),n(j,k,m),(1-1/n(j,k,m)))
			kr_n = relativeConductivity(alpha(j,k-1,m),pressure(j,k-1,m),n(j,k-1,m),(1-1/n(j,k-1,m)))
			up = (pressure(j,k-1,m)-pressure(j,k,m))/meandy * gc_n - g_n 
			down = 0.
			f_b_out = f_b_out + meandx*vardz(niz-m+1)*harmmean(ksat(j,k,m),ksat(j,k-1,m)) * up * &
				upstreammean(up,down,kr_n*rho,kr_i*rho)/mu
		end do
	endif
	
	end subroutine
    
    ! calculates the vertical flux from the first unsaturated into the last saturated cell
    subroutine verticalFlux
        
        implicit none
        
        real*8 :: kr, kr_u      ! relative hydraulic conductivity
        real*8 :: lower, upper  ! lower and upper relative hydraulic head
        
		if (i_sat_new < niz) then
	   
            kr = relativeConductivity(alpha(j,k,i_sat_new),pressure(j,k,i_sat_new),n(j,k,i_sat_new), &
                (1-1/n(j,k,i_sat_new)))
            kr_u = relativeConductivity(alpha(j,k,i_sat_new+1),pressure(j,k,i_sat_new+1),n(j,k,i_sat_new+1), &
                (1-1/n(j,k,i_sat_new+1)))
               
            lower = pressure(j,k,i_sat_new)/(rho*g)
            upper = pressure(j,k,i_sat_new+1)/(rho*g)+(vardz(niz-i_sat_new+1)+vardz(niz-i_sat_new))*0.5
            recharge_crossing(j,k) = &
			DZmean(ksat(j,k,i_sat_new),ksat(j,k,i_sat_new+1),vardz(niz-i_sat_new+1),vardz(niz-i_sat_new)) * &
			((pressure(j,k,i_sat_new+1)-pressure(j,k,i_sat_new))/(0.5*(vardz(niz-i_sat_new+1)+vardz(niz-i_sat_new)))+rho*g) * &
			upstreamMean(lower,upper,kr*rho,kr_u*rho)/mu                

        else  ! almost fully or fully saturated
		
            recharge_crossing(j,k) = 0 ! it is calculated later!
			
        endif

		
	end subroutine
	
	subroutine localgwstorchange
		
		implicit none
		
		real*8 :: rest_dist_bot, rest_dist_top ! distance between water table and end/start of the cell
		real*8 :: dh ! water table change (rise is positive)
		integer*4 :: watertab_cell_top    ! upper water table cell
		integer*4 :: watertab_cell_bot    ! lower water table cell
		real*8 :: watertab_height_top  ! upper water table height
		real*8 :: watertab_height_bot  ! lower water table height
		
		! get exact cell in which watertable is	
		if (watertab(j,k) > sum(vardz(1:niz-i_sat_new)) .or. i_sat_new == niz) then
			watertab_cell = i_sat_new
		else
			watertab_cell = i_sat_new + 1
		endif
		
		! save initial watertab_cell
		if (i == istep) then
			watertab_cell_init(j,k) = watertab_cell
		endif
		
		! calculate water table change from t = 0 to current timestep
		dh = watertab_height(j,k) - watertab_height_init(j,k)
		
		if (dh .NE. 0) then
		
			! calculate water contents at parts where water table fluctuates - take different scenarios into account
			
			! water table is in same cell as initial		
			if (watertab_cell == watertab_cell_init(j,k)) then 
			rest_dist_bot = 0
			rest_dist_top = 0			
			
			theta_watertab(j,k) = saturation(j,k,watertab_cell)*porosity(j,k,watertab_cell)
			theta_watertab_init(j,k) = theta_init(j,k,watertab_cell_init(j,k))
			
			! water table is one cell above initial
			else if (watertab_cell == watertab_cell_init(j,k)+1) then
			
			watertab_cell_top = watertab_cell
			watertab_height_top = watertab_height(j,k)
			watertab_cell_bot = watertab_cell_init(j,k)
			watertab_height_bot = watertab_height_init(j,k)
			
			! calculate distances from water tables to cell edges
			rest_dist_bot = sum(vardz(niz-watertab_cell_bot+1:niz))-watertab_height_bot
			rest_dist_top = watertab_height_top - sum(vardz(niz-watertab_cell_top+2:niz))
			
			theta_watertab(j,k) = (saturation(j,k,watertab_cell_bot) * porosity(j,k,watertab_cell_bot) * rest_dist_bot + &
			saturation(j,k,watertab_cell_top)*porosity(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			theta_watertab_init(j,k) = (theta_init(j,k,watertab_cell_bot) * rest_dist_bot + &
			theta_init(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			if(abs((rest_dist_bot+rest_dist_top) /dh)-1 > 1.0E-6) then
				print *,'fatal error - water table distance for calculations does not match dh'
				stop
			endif
			
			! water table is more than one cell above initial
			else if (watertab_cell > watertab_cell_init(j,k)+1) then
			
			watertab_cell_top = watertab_cell
			watertab_height_top = watertab_height(j,k)
			watertab_cell_bot = watertab_cell_init(j,k)
			watertab_height_bot = watertab_height_init(j,k)
			
			! calculate distances from water tables to cell edges
			rest_dist_bot = sum(vardz(niz-watertab_cell_bot+1:niz))-watertab_height_bot
			rest_dist_top = watertab_height_top - sum(vardz(niz-watertab_cell_top+2:niz))
		
			theta_watertab(j,k) = (saturation(j,k,watertab_cell_bot) * porosity(j,k,watertab_cell_bot) * rest_dist_bot + &
			sum(saturation(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			porosity(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1)) + &
			saturation(j,k,watertab_cell_top)*porosity(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			theta_watertab_init(j,k) = (theta_init(j,k,watertab_cell_bot) * rest_dist_bot + &
			sum(theta_init(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1)) + &
			theta_init(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			if(abs((rest_dist_bot+rest_dist_top+sum(vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1))) /dh)-1 > 1.0E-6) then
				print *,'fatal error - water table distance for calculations does not match dh'
				stop
			endif
			
			! water table is one cell below initial
			else if (watertab_cell_init(j,k) == watertab_cell+1) then 
			
			watertab_cell_top = watertab_cell_init(j,k)
			watertab_height_top = watertab_height_init(j,k)
			watertab_cell_bot = watertab_cell
			watertab_height_bot = watertab_height(j,k)
			
			! calculate distances from water tables to cell edges
			rest_dist_bot = sum(vardz(niz-watertab_cell_bot+1:niz))-watertab_height_bot
			rest_dist_top = watertab_height_top - sum(vardz(niz-watertab_cell_top+2:niz))
			
			theta_watertab(j,k) = (saturation(j,k,watertab_cell_bot) * porosity(j,k,watertab_cell_bot) * rest_dist_bot + &
			saturation(j,k,watertab_cell_top)*porosity(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			theta_watertab_init(j,k) = (theta_init(j,k,watertab_cell_bot) * rest_dist_bot + &
			theta_init(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			if(abs((rest_dist_bot+rest_dist_top) /dh)-1 > 1.0E-6) then
				print *,'fatal error - water table distance for calculations does not match dh'
				stop
			endif
			
			! water table is more than one cell below initial
			else if (watertab_cell_init(j,k) > watertab_cell+1) then 
			
			watertab_cell_top = watertab_cell_init(j,k)
			watertab_height_top = watertab_height_init(j,k)
			watertab_cell_bot = watertab_cell
			watertab_height_bot = watertab_height(j,k)
						
			! calculate distances from water tables to cell edges
			rest_dist_bot = sum(vardz(niz-watertab_cell_bot+1:niz))-watertab_height_bot
			rest_dist_top = watertab_height_top - sum(vardz(niz-watertab_cell_top+2:niz))	

			theta_watertab(j,k) = (saturation(j,k,watertab_cell_bot) * porosity(j,k,watertab_cell_bot) * rest_dist_bot + &
			sum(saturation(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			porosity(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1)) + &
			saturation(j,k,watertab_cell_top)*porosity(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)
			
			theta_watertab_init(j,k) = (theta_init(j,k,watertab_cell_bot) * rest_dist_bot + &
			sum(theta_init(j,k,watertab_cell_bot+1:watertab_cell_top-1)* &
			vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1)) + &
			theta_init(j,k,watertab_cell_top)*rest_dist_top)/abs(dh)			
			
			if(abs((rest_dist_bot+rest_dist_top+sum(vardz(niz-watertab_cell_bot:niz-watertab_cell_top+2:-1))) /dh)-1 > 1.0E-6) then
				print *,'fatal error - water table distance for calculations does not match dh'
				stop
			endif
			
			endif		
					
		endif
		
		! calculate cumulative storage part of recharge based on watertab fluctuations - take absolute values when subtracting water contents, so that groundwater discharge has a negative sign
		cumrech_stor(j,k) = abs(theta_watertab(j,k)-theta_watertab_init(j,k))*(dh)
		
		! save watertab fluct. part in recharge flux
		recharge_wt(j,k) = (cumrech_stor(j,k)-cumrech_stor_old(j,k))/dt
		
		cumrech_stor_old(j,k) = cumrech_stor(j,k)	

	end subroutine
		
	subroutine TopverticalFlux
        
        implicit none
		       
        real*8 :: kr, kr_u      ! relative hydraulic conductivity
        real*8 :: lower, upper  ! lower and upper relative hydraulic head
        
            kr = relativeConductivity(alpha(j,k,niz-1),pressure(j,k,niz-1),n(j,k,niz-1), &
                (1-1/n(j,k,niz-1)))
            kr_u = relativeConductivity(alpha(j,k,niz),pressure(j,k,niz),n(j,k,niz), &
                (1-1/n(j,k,niz)))
            lower = pressure(j,k,niz-1)/(rho*g)
            upper = pressure(j,k,niz)/(rho*g)+(vardz(1)+vardz(2))*0.5
            surf_wat_exch(j,k) = DZmean(ksat(j,k,niz-1),ksat(j,k,niz),vardz(2),vardz(1)) * &
			((pressure(j,k,niz)-pressure(j,k,niz-1))/(0.5*(vardz(2)+vardz(1)))+rho*g) * &
			upstreamMean(lower,upper,kr*rho,kr_u*rho)/mu
				
    end subroutine
	
	subroutine Toplayerspecificstorage
		
		implicit none
		
		topcomp(j,k) = specstor(j,k,niz) * vardz(1) * (pressure(j,k,niz)-p_old(j,k,niz))/dt ! saturation = 1, unit: m/h!

	end subroutine
	
	subroutine Groundwaterspecificstorage
		
		implicit none
		
		gwcomp(j,k) = sum(specstor(j,k,i_sat_new:1:-1) * vardz(niz-i_sat_new+1:niz) * &
		(pressure(j,k,i_sat_new:1:-1)-p_old(j,k,i_sat_new:1:-1)))/dt
		
	end subroutine
	
	! Calculate flux below roots
	subroutine VirtualLysimeter
        
        implicit none
		       
        real*8 :: kr, kr_u      ! relative hydraulic conductivity
        real*8 :: lower, upper  ! lower and upper relative hydraulic head
        
            kr = relativeConductivity(alpha(j,k,rootend-1),pressure(j,k,rootend-1),n(j,k,rootend-1), &
                (1-1/n(j,k,rootend-1)))
            kr_u = relativeConductivity(alpha(j,k,rootend),pressure(j,k,rootend),n(j,k,rootend), &
                (1-1/n(j,k,rootend)))
            lower = pressure(j,k,rootend-1)/(rho*g)
            upper = pressure(j,k,rootend)/(rho*g)+(vardz(niz-rootend+1)+vardz(niz-rootend+2))*0.5
            vl(j,k) = DZmean(ksat(j,k,rootend-1),ksat(j,k,rootend),vardz(niz-rootend+2),vardz(niz-rootend+1)) * &
			((pressure(j,k,rootend)-pressure(j,k,rootend-1))/(0.5*(vardz(niz-rootend+2)+vardz(niz-rootend+1)))+rho*g) * &
			upstreamMean(lower,upper,kr*rho,kr_u*rho)/mu
				
    end subroutine
    
end module gw_recharge

