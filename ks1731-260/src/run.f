program run_dStar
	use iso_fortran_env, only : output_unit, error_unit
    use NScool_def
    use NScool_lib
    use constants_def, only : boltzmann
    use argparse
    
    character(len=*), parameter :: default_dStar_dir = '/home/justin/dStar'
    character(len=*), parameter :: default_inlist_file = 'inlist'
    character(len=*), parameter :: default_output_directory = 'LOGS'
    real(dp), parameter :: default_core_mass = 1.4 ! Msun
    real(dp), parameter :: default_core_radius = 10.0 ! km
    real(dp), parameter :: default_core_temperature = 9.35e7 ! K
    real(dp), parameter :: default_Qimp = 4.2 ! dimensionless
    real(dp), parameter :: default_Q_heating_shallow = 1.37 ! MeV
    real(dp), parameter :: default_Q_heating_inner = 1.5 ! MeV
    character(len=64) :: my_dStar_dir, inlist, output_directory, core_mass_arg, &
    	& core_radius_arg, core_temperature_arg, Qimp_arg, &
    	& Q_heating_shallow_arg, Q_heating_inner_arg
    	 
    
    real(dp) :: eV_to_MK
    real(dp) :: core_mass, core_radius, core_temperature, Qimp, &
    	& Q_heating_shallow, Q_heating_inner
    type(NScool_info), pointer :: s=>null()
    integer :: ierr, NScool_id, i
    real(dp), dimension(8) :: pred_Teff, obs_Teff, obs_Teff_sig
    real(dp) :: chi2
    
    ierr = 0
    call command_arg_set( &
        & 'dStar_directory',"sets the main dStar root directory",ierr, &
        & flag='D',takes_parameter=.TRUE.)
    call check_okay('set command argument dStar_directory',ierr)
    
    call command_arg_set( &
        & 'inlist_file','sets the namelist parameter file',ierr, &
        & flag='I',takes_parameter=.TRUE.)
    call check_okay('set command argument inlist file',ierr)
    
    call command_arg_set( &
        & 'output_directory','sets the output directory',ierr, &
        & flag='L',takes_parameter=.TRUE.)
    call check_okay('set command argument output directory',ierr)

    call command_arg_set( &
        & 'Mcore','core mass',ierr, &
        & flag='M',takes_parameter=.TRUE.)
    call check_okay('set command argument core mass',ierr)

    call command_arg_set( &
        & 'Rcore','core radius',ierr, &
        & flag='R',takes_parameter=.TRUE.)
    call check_okay('set command argument core radius',ierr)
    
    call command_arg_set( &
        & 'Tcore','core temperature',ierr, &
        & flag='T',takes_parameter=.TRUE.)
    call check_okay('set command argument core temperature',ierr)
    
    call command_arg_set( &
        & 'Qimp','impurity parameter',ierr, &
        & flag='Q',takes_parameter=.TRUE.)
    call check_okay('set command argument impurity parameter',ierr)
    
    call command_arg_set( &
        & 'Qheating','shallow heating per nucleon',ierr, &
        & flag='s',takes_parameter=.TRUE.)
    call check_okay('set command argument shallow heating',ierr)
    
    call command_arg_set( &
        & 'Qinner','inner heating per nucleon',ierr, &
        & flag='i',takes_parameter=.TRUE.)
    call check_okay('set command argument inner heating',ierr)
    
    call parse_arguments(ierr)
    call check_okay('parse_arguments',ierr)

    my_dStar_dir = trim(get_command_arg('dStar_directory'))
    if (len_trim(my_dStar_dir)==0) my_dStar_dir = default_dStar_dir
    inlist = trim(get_command_arg('inlist_file'))
    if (len_trim(inlist)==0) inlist = default_inlist_file
    output_directory = trim(get_command_arg('output_directory'))
    if (len_trim(output_directory)==0) output_directory = default_output_directory

    core_mass_arg = trim(get_command_arg('Mcore'))
    if (len_trim(core_mass_arg) > 0)  then
        read(core_mass_arg,*) core_mass
    else
        core_mass = default_core_mass
    end if
    
    core_radius_arg = trim(get_command_arg('Rcore'))
    if (len_trim(core_radius_arg) > 0)  then
        read(core_radius_arg,*) core_radius
    else
        core_radius = default_core_radius
    end if
    
    core_temperature_arg = trim(get_command_arg('Tcore'))
    if (len_trim(core_temperature_arg) > 0)  then
        read(core_temperature_arg,*) core_temperature
    else
        core_temperature = default_core_temperature
    end if
    
    Qimp_arg = trim(get_command_arg('Qimp'))
    if (len_trim(Qimp_arg) > 0)  then
        read(Qimp_arg,*) Qimp
    else
        Qimp = default_Qimp
    end if
    
    Q_heating_shallow_arg = trim(get_command_arg('Qheating'))
    if (len_trim(Q_heating_shallow_arg) > 0)  then
        read(Q_heating_shallow_arg,*) Q_heating_shallow
    else
        Q_heating_shallow = default_Q_heating_shallow
    end if
    
    Q_heating_inner_arg = trim(get_command_arg('Qinner'))
    if (len_trim(Q_heating_inner_arg) > 0)  then
        read(Q_heating_inner_arg,*) Q_heating_inner
    else
        Q_heating_inner = default_Q_heating_inner
    end if

    call NScool_init(my_dStar_dir, ierr)
    call check_okay('NScool_init',ierr)
    
    NScool_id = alloc_NScool(ierr)
    call check_okay('NScool_id',ierr)
    
    call NScool_setup(NScool_id,inlist,ierr)
    call check_okay('NScool_setup',ierr)
    
    call get_NScool_info_ptr(NScool_id,s,ierr)
    call check_okay('get_NScool_info_ptr',ierr)
    s% core_mass = core_mass
    s% core_radius = core_radius
    s% core_temperature = core_temperature
    s% output_directory = output_directory
    s% base_profile_filename = trim(s% output_directory)//'/profile'
    s% history_filename = trim(s% output_directory)//'/history.data'
    s% profile_manifest_filename = trim(s% output_directory)//'/profiles'
    s% Qimp = Qimp
    s% Q_heating_shallow = Q_heating_shallow
    s% Q_heating_inner = Q_heating_inner
    s% Mcore = core_mass
    s% Rcore = core_radius
    s% Tcore = core_temperature
    
    call NScool_create_model(NScool_id,ierr)
    call check_okay('NScool_create_model',ierr)

    call NScool_evolve_model(NScool_id,ierr)        
    call check_okay('NScool_evolve_model',ierr)
    
    ! we don't want to compare the effective temp. at t = 0, the end of the 
    ! outburst
    pred_Teff = s% Teff_monitor(2:)/1.0e6
    eV_to_MK = 1.602176565e-12_dp/boltzmann/1.0e6
    ! observed effective temperatures (eV) and uncertainties
    obs_Teff = [104.6,89.5,76.4,73.8,71.7,70.3,64.5,64.4] * eV_to_MK
    obs_Teff_sig = [1.3,1.03,1.8,1.9,1.4,1.9,1.8,1.2] * eV_to_MK
    chi2 = sum((pred_Teff-obs_Teff)**2 / obs_Teff_sig**2 )
    
    write(output_unit,*) chi2
    !write(output_unit,'(f8.2)') chi2
    
    call NScool_shutdown
    
contains
	subroutine check_okay(msg,ierr)
		character(len=*), intent(in) :: msg
		integer, intent(inout) :: ierr
		if (ierr /= 0) then
			write (error_unit,*) trim(msg)//': ierr = ',ierr
			if (ierr < 0) stop
		end if
	end subroutine check_okay
end program run_dStar
