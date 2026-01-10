!===============================================================================
! visor_experiment.f90
! Main program for visor aerodynamics experiment
! Saves flow field snapshots for animation and comprehensive drag analysis
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
    
    ! Experiment parameters
    real(dp) :: runner_pos(3)
    real(dp) :: inlet_velocity
    integer :: n_iterations
    integer :: report_interval
    integer :: snapshot_interval
    
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
    character(len=256) :: filename, angle_str
    integer :: head_k, visor_cells
    real(dp) :: dz_at_head
    
    ! Initialize
    call cpu_time(start_time)
    
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
    inlet_velocity = 5.0_dp                   ! 5 m/s = 18 km/h running pace
    n_iterations = 500                        ! Iterations per configuration (increased for convergence)
    report_interval = 50
    snapshot_interval = 50                    ! Save snapshot every N iterations
    
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
    !---------------------------------------------------------------------------
    n_angles = 13
    allocate(angles(n_angles))
    ! Backward/down: -45, -30, -15, 0, then forward/up: 15, 30, 45, etc.
    angles = [-45.0_dp, -30.0_dp, -15.0_dp, -10.0_dp, -5.0_dp, 0.0_dp, &
              5.0_dp, 10.0_dp, 15.0_dp, 20.0_dp, 30.0_dp, 40.0_dp, 45.0_dp]
    
    allocate(drag_results(n_angles))
    allocate(mask(grid%nx, grid%ny, grid%nz))
    
    print '(A,I3,A)', 'Testing ', n_angles, ' visor angles plus baseline (no visor)'
    print '(A)', 'Angles: -45 to +45 degrees (negative = tilted down/backward)'
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
            write(filename, '(A,I4.4,A)') 'output/snapshots/baseline_', snap_count, '.bin'
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
                write(filename, '(A,I4.4,A,I4.4,A)') 'output/snapshots/visor_', &
                      nint(angle) + 100, '_', snap_count, '.bin'
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
    ! Save results to CSV file
    !---------------------------------------------------------------------------
    open(unit=10, file='output/visor_results.csv', status='replace')
    write(10, '(A)') 'angle_degrees,drag_N,drag_diff_N,drag_diff_percent'
    write(10, '(A,ES16.8,A)') 'baseline,', baseline_drag, ',0.0,0.0'
    do i_angle = 1, n_angles
        write(10, '(F8.1,A,ES16.8,A,ES16.8,A,ES16.8)') &
              angles(i_angle), ',', drag_results(i_angle), ',', &
              drag_results(i_angle) - baseline_drag, ',', &
              100.0_dp * (drag_results(i_angle) - baseline_drag) / baseline_drag
    end do
    close(10)
    print '(A)', 'Results saved to: output/visor_results.csv'
    
    ! Save grid info for Python plotting
    open(unit=11, file='output/grid_info.csv', status='replace')
    write(11, '(A)') 'nx,ny,nz,length,width,height'
    write(11, '(I5,A,I5,A,I5,A,F8.4,A,F8.4,A,F8.4)') &
          grid%nx, ',', grid%ny, ',', grid%nz, ',', &
          grid%length, ',', grid%width, ',', grid%height
    close(11)
    
    ! Save z-coordinates for Python
    open(unit=12, file='output/z_coords.csv', status='replace')
    write(12, '(A)') 'k,z,dz'
    do iter = 1, grid%nz
        write(12, '(I5,A,F10.6,A,F10.6)') iter, ',', grid%z(iter), ',', grid%dz(iter)
    end do
    close(12)
    
    print '(A)', 'Grid info saved to: output/grid_info.csv, output/z_coords.csv'
    
    !---------------------------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------------------------
    call grid%cleanup()
    deallocate(angles)
    deallocate(drag_results)
    deallocate(mask)
    
    print '(A)', ''
    print '(A)', 'Done! Run: python plot_results.py to generate animations and plots.'
    
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
