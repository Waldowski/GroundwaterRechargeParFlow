! Subroutines and functions needed to calculate lateral fluxes
module lat_flux
    
    use declaration
	use groundwatertable
    
    contains
    
	! General comments on all lateral flux calculations:
	! 1. All fluxes point in positive x- and y-direction.
	! 2. The negative sign from the darcy flow equation is included in the pressure difference by calculating start_pressure - end_pressure. 
	! Therefore we need to calculate - sin(theta) instad of + sin(theta) and no additional negative sign is needed in the equation.
	
    ! find last saturated cell
    subroutine WatertabCell
        
        implicit none
        
        ! get positions of last groundwater cell
        i_sat_new = findWaterTable(pressure(j,k,:))
		
    end subroutine 
	
    ! calculates the lateral fluxes between a grid cell and its four neighbours in groundwater domain and sum them up
    subroutine lateralFluxesSum
        
		if (l <= i_sat_new) then
			call lateralfluxes(l)
			f_r_sum = f_r_sum + f_r
			f_l_sum = f_l_sum + f_l
			f_f_sum = f_f_sum + f_f
			f_b_sum = f_b_sum + f_b
		endif

    end subroutine
	
	! Calculate the volume weighted mean lateral fluxes in whole subsurface and vadose zone
	subroutine lateralFluxesMeanBH
	
	    implicit none
		
		integer*4 :: b    ! loop index for boreholes

		call lateralfluxes(l)
		
		! print*,'l: ',l
		! print*,'j: ',j
		! print*,'k: ',k
		! print*,'i_sat_new: ',i_sat_new
		! print*,'f_r',f_r
		! print*,'f_l',f_l
		! print*,'f_f',f_f
		! print*,'f_b',f_b
		
		! calculate resulting flux
		f_res = sqrt((f_r-f_l)**2+(f_f-f_b)**2)
		
		! print*,'f_res: ',f_res
		
		latfluxmean(j,k) = latfluxmean(j,k) + f_res * vardz(niz-l+1)/sum(vardz(1:niz))
		if (l > i_sat_new) then
			latfluxvz(j,k) = latfluxvz(j,k) + f_res * vardz(niz-l+1)/sum(vardz(i_sat_new+1:niz))
		endif
		! print*,'latfluxmean: ',latfluxmean(j,k)
		
		! Calculate resulting lateral fluxes at every depth at artificial boreholes
		do b=1,size(latflux_bh,1) ! borehole loop  
			if (cfact == 1) then
				! print*,'checkpoint 1'
				if (j == (positions_bh(b,1)-1)*cfact+1 .and. k == (positions_bh(b,2)-1)*cfact+1) then
					! print*,'checkpoint 2'
					! print*,'f_r',f_r
					! print*,'f_l',f_l
					! print*,'f_f',f_f
					! print*,'f_b',f_b
					
					latflux_x = (f_r-f_l);
					latflux_y = (f_f-f_b);
					latflux_bh(b,l) = sqrt(latflux_x**2+latflux_y**2);
					! print*,'latflux_x: ',latflux_x
					! print*,'latflux_y',latflux_y
					! print*,'latflux_bh(1,l) inside',latflux_bh(1,l)
				endif
			else
				if (j .GE. (positions_bh(b,1)-1)*cfact+1 .and. j .LE. positions_bh(b,1)*cfact .and. &
				k .GE. (positions_bh(b,2)-1)*cfact+1 .and. k .LE. positions_bh(b,2)*cfact)  then
					print*,'this should not have happened'
					latflux_x = (f_r-f_l);
					latflux_y = (f_f-f_b);

					latflux_x_mean = latflux_x_mean + latflux_x/(cfact**2);
					latflux_y_mean = latflux_y_mean + latflux_y/(cfact**2);
					latflux_bh_resmean(b,:) =  latflux_bh_resmean(b,:) + sqrt(latflux_x**2+latflux_y**2)/(cfact**2); ! all parts of the mean are summed up
					if (j == positions_bh(b,1)*cfact .and. k == positions_bh(b,2)*cfact) then
						latflux_bh(b,l) = sqrt(latflux_x_mean**2+latflux_y_mean**2); ! calculate resulting flux from summed up means
						latflux_x_mean = 0.0; ! reset values
						latflux_y_mean = 0.0;
					endif
				endif
			endif
		end do
	
	end subroutine
	
	! general subroutine to calculate lateral fluxes
	subroutine lateralfluxes(z)
	
	    real*8 :: gc_n, g_n     ! averaged slopes
        real*8 :: kr_i, kr_n    ! relative hydraulic conductivies
        real*8 :: up, down      ! up and downstream gradient
		integer*4 :: z          ! height coordinate
		
		! print*,'slopex',slopex(j,k,1)
		! print*,'slopey',slopey(j,k,1)
		! print*,'ks',ksat(j,k,z)
		
		! print*,'alpha',slopey(j,k,z)
		! print*,'pressure',pressure(j,k,z)
		
		! at right boundary?
        if (j < nix) then        
			gc_n = g*0.5*(cos(atan(slopex(j,k,1))) + cos(atan(slopex(j+1,k,1))))
			g_n = g*0.5*(sin(atan(slopex(j,k,1))) + sin(atan(slopex(j+1,k,1))))
			kr_i = relativeConductivity(alpha(j,k,z),pressure(j,k,z),n(j,k,z),(1-1/n(j,k,z)))
			kr_n = relativeConductivity(alpha(j+1,k,z),pressure(j+1,k,z),n(j+1,k,z),(1-1/n(j+1,k,z)))
			up = (pressure(j,k,z)-pressure(j+1,k,z))/meandx * gc_n - g_n 
			down = 0.
			f_r = meandy*vardz(niz-z+1)*harmmean(ksat(j,k,z),ksat(j+1,k,z)) * up * &
				upstreammean(up,down,kr_i*rho,kr_n*rho)/mu
		else
			f_r = 0.0
        endif
        
        ! at left boundary?
        if (j > 1) then
			gc_n = g*0.5*(cos(atan(slopex(j,k,1))) + cos(atan(slopex(j-1,k,1))))
			g_n = g*0.5*(sin(atan(slopex(j,k,1))) + sin(atan(slopex(j-1,k,1))))
			kr_i = relativeConductivity(alpha(j,k,z),pressure(j,k,z),n(j,k,z),(1-1/n(j,k,z)))
			kr_n = relativeConductivity(alpha(j-1,k,z),pressure(j-1,k,z),n(j-1,k,z),(1-1/n(j-1,k,z)))
			up = (pressure(j-1,k,z)-pressure(j,k,z))/meandx * gc_n - g_n 
			down = 0.
			f_l = meandy*vardz(niz-z+1)*harmmean(ksat(j,k,z),ksat(j-1,k,z)) * up * &
				upstreammean(up,down,kr_n*rho,kr_i*rho)/mu
		else
			f_l = 0.0
        endif
        
        ! at front boundary?
        if (k < niy) then
			gc_n = g*0.5*(cos(atan(slopey(j,k,1))) + cos(atan(slopey(j,k+1,1))))
			g_n = g*0.5*(sin(atan(slopey(j,k,1))) + sin(atan(slopey(j,k+1,1))))
			kr_i = relativeConductivity(alpha(j,k,z),pressure(j,k,z),n(j,k,z),(1-1/n(j,k,z)))
			kr_n = relativeConductivity(alpha(j,k+1,z),pressure(j,k+1,z),n(j,k+1,z),(1-1/n(j,k+1,z)))
			up = (pressure(j,k,z)-pressure(j,k+1,z))/meandy * gc_n - g_n 
			down = 0.
			f_f = meandx*vardz(niz-z+1)*harmmean(ksat(j,k,z),ksat(j,k+1,z)) * up * &
				upstreammean(up,down,kr_i*rho,kr_n*rho)/mu
		else
			f_f = 0.0
        endif
        
        ! at back boundary?
        if (k > 1) then     
			gc_n = g*0.5*(cos(atan(slopey(j,k,1))) + cos(atan(slopey(j,k-1,1))))
			g_n = g*0.5*(sin(atan(slopey(j,k,1))) + sin(atan(slopey(j,k-1,1))))
			kr_i = relativeConductivity(alpha(j,k,z),pressure(j,k,z),n(j,k,z),(1-1/n(j,k,z)))
			kr_n = relativeConductivity(alpha(j,k-1,z),pressure(j,k-1,z),n(j,k-1,z),(1-1/n(j,k-1,z)))
			up = (pressure(j,k-1,z)-pressure(j,k,z))/meandy * gc_n - g_n 
			down = 0.
			f_b = meandx*vardz(niz-z+1)*harmmean(ksat(j,k,z),ksat(j,k-1,z)) * up * &
				upstreammean(up,down,kr_n*rho,kr_i*rho)/mu
			else
			f_b = 0.0
        endif
		
		end subroutine
	
    ! returns the index of the last groundwater cell in the column
    function findWaterTable(p) result (i_sat)

        implicit none

        integer*4 :: i_sat, k   ! index of positions and loop index
        real*8 :: p(:)          ! pressure column
        
        ! move upwards through column
        do k=1,size(p,1)
            ! first unsaturated cell
            if (p(k)<0) then
                i_sat = k-1
                exit
            endif
            ! fully saturated
            if (k==size(p,1)) then
                i_sat = k
                exit
            endif
        end do

    end function findWaterTable
    
    ! calculates the relative conductivity according to van Genuchten-Mualem
    function relativeConductivity(alpha,p,n,m) result (kr)

        implicit none

        real*8 :: alpha, p, n, m, a, se, kr
        
        if (p .GE. 0) then
            kr = 1
        else
            a = alpha*abs(p)
            se = 1 / ((1 + a**n) **m)
            kr = (1 - a**(n-1) * se)**2 * se**0.5
        endif

        

    end function relativeConductivity
    
    ! calculates the harmonic mean of two scalars
    function harmmean(a, b) result (h)

        implicit none 

        real*8 :: a, b, h

        h = 2 / (1/a + 1/b);

    end function harmmean
    
    ! calculates the upstream weighted mean of two scalars
    function upstreammean(a,b,c,d) result (h)

        implicit none

        real*8 :: a, b, c, d, h

        if ((a-b)>=0) then
            h = c
        else
            h = d
        endif

    end function upstreammean
    
    ! calculates the DZ weighted mean of two scalars
    function DZmean(a,b,c,d) result (h)

        implicit none 

        real*8 :: a, b, c, d, h

        if ((b*c + a*d) /= 0) then 
            h = (c+d)*a*b / (b*c+a*d)
        else
            h = 0
        endif

    end function DZmean
    
end module lat_flux

