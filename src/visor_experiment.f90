!===============================================================================
! visor_experiment.f90
! Main program for visor aerodynamics experiment
! Saves flow field snapshots for animation and comprehensive drag analysis
!
! Runtime configuration via namelist file (visor.nml) or command-line:
!   ./visor_cfd                  # Uses visor.nml or defaults
!   ./visor_cfd 2000             # Override iterations to 2000
!===============================================================================
program visor_experiment
    use grid_module
    use geometry_module
    use solver_module
    implicit none
    
    ! Grid and solver
    type(Grid3D) :: grid
    type(FVMSolver) :: solver
    
    ! Geometry mask
    logical, allocatable :: mask(:,:,:)
    
    ! Experiment parameters (can be set via namelist)
    real(dp) :: runner_pos(3)
    real(dp) :: inlet_velocity
    integer :: n_iterations
    integer :: report_interval
    integer :: snapshot_interval
    
    ! Namelist for runtime configuration
    namelist /experiment/ n_iterations, report_interval, snapshot_interval, inlet_velocity
    
    ! Visor angles to test (including backward angles)
    real(dp), allocatable :: angles(:)
    integer :: n_angles
    
    ! Results
    real(dp), allocatable :: drag_results(:)
    real(dp) :: baseline_drag
    
    ! Loop variables
    integer :: i_angle, iter, snap_count
    real(dp) :: angle, Fx, Fy, Fz, u_max, dt
    real(dp) :: start_time, end_time, elapsed
    character(len=256) :: filename, angle_str, iter_suffix
    character(len=32) :: arg
    integer :: head_k, visor_cells, ios, nargs
    real(dp) :: dz_at_head
    logical :: nml_exists
    
    ! Initialize
    call cpu_time(start_time)
    
    !---------------------------------------------------------------------------
    ! Set defaults, then read namelist, then check command-line override
    !---------------------------------------------------------------------------
    n_iterations = 500
    report_interval = 50
    snapshot_interval = 50
    inlet_velocity = 5.0_dp
    
    ! Try to read namelist file
    inquire(file='visor.nml', exist=nml_exists)
    if (nml_exists) then
        open(unit=99, file='visor.nml', status='old', iostat=ios)
        if (ios == 0) then
            read(99, nml=experiment, iostat=ios)
            close(99)
        end if
    end if
    
    ! Check for command-line override of n_iterations
    nargs = command_argument_count()
    if (nargs >= 1) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=ios) n_iterations
        if (ios /= 0) then
            print '(A)', 'Warning: Invalid command-line argument, using namelist/default value'
        end if
    end if
    
    ! Create iteration suffix for output files (e.g., "_500" or "_2000")
    write(iter_suffix, '(A,I0)') '_', n_iterations
    
    print '(A)', ''
    print '(A)', '========================================================================'
    print '(A)', 'VISOR AERODYNAMICS EXPERIMENT'
    print '(A)', 'Finite Volume Method with Non-Uniform Grid (Fortran Implementation)'
    print '(A)', '========================================================================'
    print '(A)', ''
    
    !---------------------------------------------------------------------------
    ! Setup grid - 2.5m height, 150 vertical levels, non-uniform
    ! Moderate resolution balances visor resolution with solver stability
    ! dz ~ 4mm at head (5mm visor reliably detected at all angles)
    !---------------------------------------------------------------------------
    call grid%init(nx=60, ny=45, nz=150, &
                   length=2.0_dp, width=1.7_dp, height=2.5_dp, &
                   z_refine_min=1.6_dp, z_refine_max=2.1_dp)
    
    call grid%print_info()
    
    ! Report visor resolution
    head_k = 1
    do while (grid%z(head_k) < 1.84_dp .and. head_k < grid%nz)
        head_k = head_k + 1
    end do
    dz_at_head = grid%dz(head_k)
    visor_cells = nint(0.005_dp / dz_at_head)
    print '(A,F6.2,A)', '  dz at head height (1.84m): ', dz_at_head * 1000.0_dp, ' mm'
    print '(A,I3,A)', '  5mm visor spans approximately ', visor_cells, ' cells'
    print '(A)', ''
    
    !---------------------------------------------------------------------------
    ! Experiment configuration
    !---------------------------------------------------------------------------
    runner_pos = [0.8_dp, 0.85_dp, 0.0_dp]  ! Runner base position
    
    print '(A)', 'Experiment Configuration:'
    print '(A,3F8.3)', '  Runner position: ', runner_pos
    print '(A,F6.2,A,F6.1,A)', '  Inlet velocity: ', inlet_velocity, ' m/s (', &
          inlet_velocity * 3.6_dp, ' km/h)'
    print '(A,I5)', '  Iterations per config: ', n_iterations
    print '(A,I5)', '  Snapshot interval: ', snapshot_interval
    print '(A)', ''
    
    !---------------------------------------------------------------------------
    ! Define visor angles to test
    ! Positive = tilted up (normal wear)
    ! Negative = tilted down (including backward-facing)
    ! Increment by 5 degrees from -45 to +45
    !---------------------------------------------------------------------------
    n_angles = 19
    allocate(angles(n_angles))
    ! -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45
    angles = [-45.0_dp, -40.0_dp, -35.0_dp, -30.0_dp, -25.0_dp, -20.0_dp, -15.0_dp, -10.0_dp, -5.0_dp, &
               0.0_dp, &
               5.0_dp, 10.0_dp, 15.0_dp, 20.0_dp, 25.0_dp, 30.0_dp, 35.0_dp, 40.0_dp, 45.0_dp]
    
    allocate(drag_results(n_angles))
    allocate(mask(grid%nx, grid%ny, grid%nz))
    
    print '(A,I3,A)', 'Testing ', n_angles, ' visor angles plus baseline (no visor)'
    print '(A)', 'Angles: -45 to +45 degrees in 5° increments (negative = tilted down/backward)'
    print '(A)', ''
    
    ! Create output directory for snapshots
    call execute_command_line('mkdir -p output/snapshots', wait=.true.)
    
    !---------------------------------------------------------------------------
    ! Run baseline (no visor)
    !---------------------------------------------------------------------------
    print '(A)', '================================================================'
    print '(A)', 'Running: NO VISOR (baseline)'
    print '(A)', '================================================================'
    
    call create_geometry_mask(grid, runner_pos, 0.0_dp, mask, use_visor=.false.)
    
    call solver%init(grid, inlet_velocity)
    call solver%set_solid_mask(mask)
    
    snap_count = 0
    do iter = 1, n_iterations
        call solver%step(dt)
        
        if (mod(iter, report_interval) == 0 .or. iter == n_iterations) then
            call solver%compute_drag(Fx, Fy, Fz)
            u_max = solver%get_max_velocity()
            print '(A,I5,A,ES10.3,A,F8.3,A,ES12.4,A)', &
                  '  Iter ', iter, ': dt=', dt, 's, u_max=', u_max, ' m/s, Fx=', Fx, ' N'
        end if
        
        ! Save snapshot for animation
        if (mod(iter, snapshot_interval) == 0) then
            snap_count = snap_count + 1
            write(filename, '(A,A,A,I4.4,A)') 'output/snapshots/baseline', &
                  trim(iter_suffix), '_', snap_count, '.bin'
            call save_snapshot(solver, grid, filename)
        end if
    end do
    
    call solver%compute_drag(Fx, Fy, Fz)
    baseline_drag = Fx
    print '(A,ES12.4,A)', 'Baseline drag: ', baseline_drag, ' N'
    
    call solver%cleanup()
    
    !---------------------------------------------------------------------------
    ! Run each visor angle
    !---------------------------------------------------------------------------
    do i_angle = 1, n_angles
        angle = angles(i_angle)
        
        print '(A)', ''
        print '(A)', '================================================================'
        print '(A,F7.1,A)', 'Running: Visor angle = ', angle, ' degrees'
        print '(A)', '================================================================'
        
        call create_geometry_mask(grid, runner_pos, angle, mask, use_visor=.true.)
        
        call solver%init(grid, inlet_velocity)
        call solver%set_solid_mask(mask)
        
        snap_count = 0
        do iter = 1, n_iterations
            call solver%step(dt)
            
            if (mod(iter, report_interval) == 0 .or. iter == n_iterations) then
                call solver%compute_drag(Fx, Fy, Fz)
                u_max = solver%get_max_velocity()
                print '(A,I5,A,ES10.3,A,F8.3,A,ES12.4,A)', &
                      '  Iter ', iter, ': dt=', dt, 's, u_max=', u_max, ' m/s, Fx=', Fx, ' N'
            end if
            
            ! Save snapshot for animation
            if (mod(iter, snapshot_interval) == 0) then
                snap_count = snap_count + 1
                write(filename, '(A,I4.4,A,A,I4.4,A)') 'output/snapshots/visor_', &
                      nint(angle) + 100, trim(iter_suffix), '_', snap_count, '.bin'
                call save_snapshot(solver, grid, filename)
            end if
        end do
        
        call solver%compute_drag(Fx, Fy, Fz)
        drag_results(i_angle) = Fx
        
        call solver%cleanup()
    end do
    
    !---------------------------------------------------------------------------
    ! Print summary
    !---------------------------------------------------------------------------
    call cpu_time(end_time)
    elapsed = end_time - start_time
    
    print '(A)', ''
    print '(A)', '========================================================================'
    print '(A)', 'EXPERIMENT SUMMARY'
    print '(A)', '========================================================================'
    print '(A)', ''
    print '(A)', 'Configuration              Drag (N)        vs Baseline      % Change'
    print '(A)', '------------------------------------------------------------------------'
    print '(A,ES14.4,A)', 'No visor (baseline)    ', baseline_drag, '           ---            ---'
    
    do i_angle = 1, n_angles
        angle = angles(i_angle)
        Fx = drag_results(i_angle)
        print '(A,F6.1,A,ES14.4,A,ES12.4,A,F8.2,A)', &
              'Visor ', angle, ' deg       ', Fx, '    ', &
              Fx - baseline_drag, '     ', &
              100.0_dp * (Fx - baseline_drag) / baseline_drag, '%'
    end do
    
    print '(A)', '------------------------------------------------------------------------'
    print '(A)', ''
    
    ! Find optimal angle (minimum absolute drag)
    i_angle = minloc(abs(drag_results), 1)
    print '(A,F7.1,A)', 'Optimal visor angle: ', angles(i_angle), ' degrees'
    print '(A,ES12.4,A)', 'Drag at optimal angle: ', drag_results(i_angle), ' N'
    print '(A,ES12.4,A,F8.2,A)', 'Change from baseline: ', &
          drag_results(i_angle) - baseline_drag, ' N (', &
          100.0_dp * (drag_results(i_angle) - baseline_drag) / baseline_drag, '%)'
    print '(A)', ''
    print '(A,F10.1,A)', 'Total computation time: ', elapsed, ' seconds'
    
    !---------------------------------------------------------------------------
    ! Save results to CSV file (with iteration count in filename)
    !---------------------------------------------------------------------------
    write(filename, '(A,A,A)') 'output/visor_results', trim(iter_suffix), '.csv'
    open(unit=10, file=trim(filename), status='replace')
    write(10, '(A)') 'angle_degrees,drag_N,drag_diff_N,drag_diff_percent'
    write(10, '(A,ES16.8,A)') 'baseline,', baseline_drag, ',0.0,0.0'
    do i_angle = 1, n_angles
        write(10, '(F8.1,A,ES16.8,A,ES16.8,A,ES16.8)') &
              angles(i_angle), ',', drag_results(i_angle), ',', &
              drag_results(i_angle) - baseline_drag, ',', &
              100.0_dp * (drag_results(i_angle) - baseline_drag) / baseline_drag
    end do
    close(10)
    print '(A,A)', 'Results saved to: ', trim(filename)
    
    ! Save grid info for Python plotting (with iteration suffix for consistency)
    write(filename, '(A,A,A)') 'output/grid_info', trim(iter_suffix), '.csv'
    open(unit=11, file=trim(filename), status='replace')
    write(11, '(A)') 'nx,ny,nz,length,width,height'
    write(11, '(I5,A,I5,A,I5,A,F8.4,A,F8.4,A,F8.4)') &
          grid%nx, ',', grid%ny, ',', grid%nz, ',', &
          grid%length, ',', grid%width, ',', grid%height
    close(11)
    
    ! Save z-coordinates for Python (with iteration suffix)
    write(filename, '(A,A,A)') 'output/z_coords', trim(iter_suffix), '.csv'
    open(unit=12, file=trim(filename), status='replace')
    write(12, '(A)') 'k,z,dz'
    do iter = 1, grid%nz
        write(12, '(I5,A,F10.6,A,F10.6)') iter, ',', grid%z(iter), ',', grid%dz(iter)
    end do
    close(12)
    
    print '(A,A,A)', 'Grid info saved to: output/grid_info', trim(iter_suffix), '.csv'
    print '(A,A,A)', 'Z-coords saved to: output/z_coords', trim(iter_suffix), '.csv'
    
    !---------------------------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------------------------
    call grid%cleanup()
    deallocate(angles)
    deallocate(drag_results)
    deallocate(mask)
    
    print '(A)', ''
    print '(A,A,A)', 'Done! Run: python plot_results.py ', trim(iter_suffix(2:)), &
          ' to generate animations and plots.'
    
contains

    !---------------------------------------------------------------------------
    ! Save flow field snapshot to binary file
    !---------------------------------------------------------------------------
    subroutine save_snapshot(slv, grd, fname)
        type(FVMSolver), intent(in) :: slv
        type(Grid3D), intent(in) :: grd
        character(len=*), intent(in) :: fname
        
        integer :: unit_num
        
        unit_num = 20
        open(unit=unit_num, file=fname, status='replace', access='stream', &
             form='unformatted')
        
        ! Write dimensions
        write(unit_num) grd%nx, grd%ny, grd%nz
        
        ! Write velocity fields
        write(unit_num) slv%u
        write(unit_num) slv%v
        write(unit_num) slv%w
        
        ! Write pressure
        write(unit_num) slv%p
        
        ! Write solid mask
        write(unit_num) slv%solid
        
        close(unit_num)
        
    end subroutine save_snapshot
    
end program visor_experiment
