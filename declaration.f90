!     
! File:   declaration.f90
! Author: Brandhorst
!
! Created on 24. Juli 2018, 15:29
!

! Declaration of all variables possibly needed for the Parflow analysis
module declaration
    
    ! constant parameters
    parameter(g = 1.0)                                                        ! gravity
    parameter(rho = 1.0)                                                      ! density
    parameter(mu = 1.0)                                                       ! viscosity
	real*8 :: specstor_const                                                  ! specific storage
    real*8 :: poro_const                                                      ! porosity
    real*8 :: ksat_const                                                      ! saturated hydraulic conductivity
    real*8 :: alpha_const, n_const                                            ! van Genuchten model parameters
    real*8 :: sx_const                                                        ! slope in x-direction
    real*8 :: sy_const                                                        ! slope in y-direction
	real*8 :: maskgwr_const                                                   ! groundwater recharge mask
    
    ! timing
    integer*4 :: istep                                                        ! starting timestep
    integer*4 :: nstep                                                        ! number of timesteps
    real*8 :: dt                                                              ! timestep size
    integer*4 :: i                                                            ! time loop index
    
    ! grid
    integer*4 :: nix                                                          ! number of X-Cells
    integer*4 :: niy                                                          ! number of Y-Cells
    integer*4 :: niz                                                          ! number of Z-Cells
    real*8 :: meandx                                                          ! mean size of X-Cells
    real*8 :: meandy                                                          ! mean size of Y-Cells
    real*8 :: meandz                                                          ! mean size of Z-Cells
    real*8,allocatable :: vardz(:)                                            ! actual sizes of Z-Cells
    integer*4 :: j, k, l                                                      ! loop indices for loops over X,Y,Z
    real*8,allocatable :: mask(:,:,:)                                         ! defines which cells are used for analysis
    
    ! ensemble
    integer*4 :: ensend                                                       ! number of last ensemble member
    integer*4 :: ensstart                                                     ! number of first ensemble member
    integer*4 :: ensrun                                                       ! ensemble loop index
    
    ! parameters
    real*8,allocatable :: specstor(:,:,:)                                     ! specific storage
	real*8,allocatable :: topcomp(:,:)                                        ! specific storage effects at most upper layer
	real*8,allocatable :: gwcomp(:,:)                                         ! specific storage effects summed up for groundwater
    real*8,allocatable :: porosity(:,:,:)                                     ! porosity
    real*8,allocatable :: ksat(:,:,:)                                         ! saturated hydraulic conductivity
    real*8,allocatable :: alpha(:,:,:), n(:,:,:)                              ! van Genuchten model parameters
    real*8,allocatable :: mannings(:,:,:)                                     ! Manning's coefficient
    real*8,allocatable :: slopex(:,:,:)                                       ! slope in x-direction
    real*8,allocatable :: slopey(:,:,:)                                       ! slope in y-direction
	
    ! raw output
    real*8,allocatable :: saturation(:,:,:)                                   ! saturation 
    real*8,allocatable :: pressure(:,:,:)                                     ! pressure	
    
    ! groundwater recharge
    integer*4 :: i_sat_old, i_sat_new                                         ! index of last saturated cell at old and new timestep
	integer*4 :: watertab_cell                                                ! cell which includes water table at current time step
	integer*4,allocatable :: watertab_cell_init(:,:)                          ! cell which includes water table at current initial step
    real*8 :: f_l, f_r, f_f, f_b, f_l_sum, f_r_sum, f_f_sum, f_b_sum,f_res    ! lateral fluxes to four neighbouring cells and their sums
    real*8 :: f_r_out, f_l_out, f_b_out, f_f_out                              ! lateral component of discharge (recharge)
	real*8,allocatable :: p_old(:,:,:)                                        ! pressure at old timestep
    real*8,allocatable :: recharge_net(:,:),recharge_gross(:,:)               ! net and gross groundwater recharge
	real*8,allocatable :: recharge_crossing(:,:),recharge_wt(:,:)             ! flux crossing gwt and recharge due to water table fluctuations
	real*8,allocatable :: recharge_net_sat(:,:),recharge_gross_sat(:,:)       ! saturated recharge parts
	real*8,allocatable :: recharge_net_unsat(:,:), recharge_gross_unsat(:,:) ! unsaturated recharge parts
	real*8,allocatable :: riv_recharge_net(:,:), riv_recharge_gross(:,:)                             ! saturated parts of river recharge
	real*8,allocatable :: riv_recharge_net_unsat(:,:), riv_recharge_gross_unsat(:,:)                       ! unsaturated parts of river recharge
	real*8,allocatable :: riv_exch_flux(:,:),surf_wat_exch(:,:)               ! surface water exchange fluxes
    real*8,allocatable :: mask_gwr(:,:)                                       ! mask for gwr to account for rivers
	real*8,allocatable :: theta_init(:,:,:)                                   ! initial water content
	real*8,allocatable :: cumrech_stor(:,:)                                   ! cumulative storage component of recharge
	real*8,allocatable :: cumrech_stor_old(:,:)                               ! cumulative storage component of recharge (old timestep)
	real*8,allocatable :: theta_watertab(:,:)                                 ! water content at area where water table fluctuates
	real*8,allocatable :: theta_watertab_init(:,:)                            ! initial water content at area where water table fluctuates
	integer*4 :: rootend                                                      ! final cell of rootzone (counted bottom to top)
	real*8,allocatable :: vl(:,:)                                             ! flux below rootzone
	
    ! groundwater table
    real*8,allocatable :: watertab(:,:)                                       ! depth to groundwater table
	real*8,allocatable :: watertab_height(:,:)                                ! height to groundwater table
	real*8,allocatable :: watertab_height_init(:,:)                           ! initial groundwater table height
	logical,allocatable :: watertab_switch(:,:)                               ! switch to make sure water table is only calculated once at each x-y location
    
    ! file names
    character*200 :: namesatur, namesatur_curr                                ! name of saturation file (current)
    character*200 :: namepress, namepress_curr                                ! name of pressure file (current)
    character*200 :: nameslopex, nameslopex_curr                      ! name of x-slope file (current)
    character*200 :: nameslopey, nameslopey_curr                      ! name of y-slope file (current)
    character*200 :: nameporo, nameporo_curr                          ! name of porosity file (current)
    character*200 :: nameksat, nameksat_curr                          ! name of conductivity file (current)
    character*200 :: namealpha, namealpha_curr                        ! name of van Genuchten alpha file (current)
    character*200 :: namen, namen_curr                                ! name of van Genuchten n file (current)
    character*200 :: namemask, namemask_curr                          ! name of mask file (current)
    character*200 :: namemaskgwr, namemaskgwr_curr                    ! name of gwr mask file (current)
    character*200 :: namespecstor, namespecstor_curr                  ! name of specific storage file (current)
    character*200 :: nameoutputpath, nameoutputpath_curr          ! path to generated output files
    character*200 :: nameofcase, nameofcase_curr          ! specific name contained in all generated output
    character*200 :: nameinput                                        ! name of user input file
	character*200 :: nameoutput                                           ! name of output file + path
    character*5 :: nameruntime                                        ! string of current time step
    character*5 :: ensname                                            ! string of current ensemble name (number - 1)
    character*3 :: ensnumber                                          ! string of current ensemble number
	
    ! saving
    integer*4 :: reclen, irec                                         ! size and position of data to be saved
    
    ! flags for calculation 
    logical :: doGrossRecharge                                        ! calculate groundwater recharge (Fastest method! For ParFlow standalone, it is the same as net recharge)             
    logical :: doNetRecharge                                          ! calculate groundwater recharge (Slower. For ParFlow-CLM. Groundwater extracted by plants is not counted as recharge here)
	logical :: doFluxCrossing                                         ! calculate the vertical flux crossing the water table
	logical :: doRechargeSources                                      ! split groundwater recharge into its different sources
	logical :: doVirtualLysimeter                                     ! calculate vertical flux below roots
    logical :: doGWTable                                              ! calculate depth to groundwater table
	
	logical :: constsx
	logical :: constsy
	logical :: constporo
	logical :: constksat
	logical :: constalpha
	logical :: constn
	logical :: constmask
	logical :: constmaskgwr
	logical :: constspecstor
    
    ! command line input arguments
    integer*4 :: narg,arg                                             ! total and current number of input arguments
    character*200 :: argName                                          ! string of input argument
    
end module declaration

