!===============================================================================
! test_solver.f90
! Diagnostic test: run a few iterations and check for stability
!===============================================================================
program test_solver
    use grid_module
    use geometry_module
    use solver_module
    implicit none
    
    type(Grid3D) :: grid
    type(FVMSolver) :: solver
    logical, allocatable :: mask(:,:,:)
    
    real(dp) :: runner_pos(3)
    real(dp) :: inlet_velocity
    integer :: iter, n_iterations
    real(dp) :: Fx, Fy, Fz, u_max, dt
    real(dp) :: u_min, u_mean, p_max, p_min
    integer :: n_solid
    
    print '(A)', ''
    print '(A)', '========================================'
    print '(A)', 'SOLVER DIAGNOSTIC TEST'
    print '(A)', '========================================'
    print '(A)', ''
    
    !---------------------------------------------------------------------------
    ! Use nz=150 grid - balanced for visor detection AND stability
    ! dz ~ 4mm at head (5mm visor reliably detected)
    !---------------------------------------------------------------------------
    print '(A)', 'Initializing nz=150 grid (balanced for detection + stability)...'
    call grid%init(nx=60, ny=45, nz=150, &
                   length=2.0_dp, width=1.7_dp, height=2.5_dp, &
                   z_refine_min=1.6_dp, z_refine_max=2.1_dp)
    
    call grid%print_info()
    
    !---------------------------------------------------------------------------
    ! Setup geometry (no visor for simplicity)
    !---------------------------------------------------------------------------
    runner_pos = [0.8_dp, 0.85_dp, 0.0_dp]
    inlet_velocity = 5.0_dp
    
    allocate(mask(grid%nx, grid%ny, grid%nz))
    
    print '(A)', ''
    print '(A)', 'Creating geometry mask (runner only, no visor)...'
    call create_geometry_mask(grid, runner_pos, 0.0_dp, mask, use_visor=.false.)
    
    n_solid = count(mask)
    print '(A,I8)', '  Solid cells: ', n_solid
    print '(A,F8.4)', '  Solid fraction: ', real(n_solid) / real(grid%nx * grid%ny * grid%nz)
    
    !---------------------------------------------------------------------------
    ! Initialize solver
    !---------------------------------------------------------------------------
    print '(A)', ''
    print '(A)', 'Initializing solver...'
    print '(A,F6.2,A)', '  Inlet velocity: ', inlet_velocity, ' m/s'
    print '(A,ES10.3,A)', '  Viscosity: ', 1.5e-5_dp, ' m²/s'
    print '(A,F8.4,A)', '  Density: ', 1.225_dp, ' kg/m³'
    
    call solver%init(grid, inlet_velocity)
    call solver%set_solid_mask(mask)
    
    !---------------------------------------------------------------------------
    ! Run diagnostic iterations
    !---------------------------------------------------------------------------
    n_iterations = 50
    
    print '(A)', ''
    print '(A)', '========================================'
    print '(A,I4,A)', 'Running ', n_iterations, ' diagnostic iterations'
    print '(A)', '========================================'
    print '(A)', ''
    print '(A)', 'Iter |    dt (s)    | u_max (m/s) | u_min (m/s) | p_max (Pa)  | p_min (Pa)  |   Fx (N)'
    print '(A)', '-----|--------------|-------------|-------------|-------------|-------------|------------'
    
    do iter = 1, n_iterations
        call solver%step(dt)
        
        ! Compute diagnostics
        u_max = maxval(solver%u, mask=.not.solver%solid)
        u_min = minval(solver%u, mask=.not.solver%solid)
        p_max = maxval(solver%p, mask=.not.solver%solid)
        p_min = minval(solver%p, mask=.not.solver%solid)
        
        call solver%compute_drag(Fx, Fy, Fz)
        
        print '(I4,A,ES12.4,A,F11.4,A,F11.4,A,ES11.3,A,ES11.3,A,ES10.3)', &
              iter, ' |', dt, ' |', u_max, ' |', u_min, ' |', p_max, ' |', p_min, ' |', Fx
        
        ! Check for blow-up
        if (abs(u_max) > 1000.0_dp .or. abs(u_min) > 1000.0_dp) then
            print '(A)', ''
            print '(A)', '*** DIVERGENCE DETECTED: velocity > 1000 m/s ***'
            print '(A)', 'Stopping early.'
            exit
        end if
        
        if (isnan(u_max) .or. isnan(p_max)) then
            print '(A)', ''
            print '(A)', '*** NaN DETECTED ***'
            print '(A)', 'Stopping early.'
            exit
        end if
    end do
    
    print '(A)', ''
    print '(A)', '========================================'
    print '(A)', 'FINAL STATE'
    print '(A)', '========================================'
    print '(A,ES12.4)', 'u_max: ', maxval(abs(solver%u), mask=.not.solver%solid)
    print '(A,ES12.4)', 'v_max: ', maxval(abs(solver%v), mask=.not.solver%solid)
    print '(A,ES12.4)', 'w_max: ', maxval(abs(solver%w), mask=.not.solver%solid)
    print '(A,ES12.4)', 'p_max: ', maxval(abs(solver%p), mask=.not.solver%solid)
    print '(A)', ''
    
    if (abs(maxval(solver%u)) > 100.0_dp) then
        print '(A)', 'STATUS: UNSTABLE - solver is diverging'
        print '(A)', ''
        print '(A)', 'Likely causes:'
        print '(A)', '  1. CFL number too high for explicit scheme'
        print '(A)', '  2. Grid too fine (small dz) causing viscous instability'
        print '(A)', '  3. Pressure correction needs more iterations'
        print '(A)', ''
        print '(A)', 'Suggested fixes:'
        print '(A)', '  - Reduce CFL from 0.3 to 0.1 or lower'
        print '(A)', '  - Increase pressure iterations from 50 to 200+'
        print '(A)', '  - Add velocity clamping to prevent runaway'
    else
        print '(A)', 'STATUS: STABLE - solver appears to be working'
    end if
    
    call solver%cleanup()
    call grid%cleanup()
    deallocate(mask)
    
    print '(A)', ''
    print '(A)', 'Diagnostic test complete.'
    
end program test_solver
