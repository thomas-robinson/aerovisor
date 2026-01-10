!===============================================================================
! solver_module.f90
! Finite Volume CFD Solver with non-uniform grid support
!===============================================================================
module solver_module
    use grid_module, only: dp, Grid3D
    implicit none
    private
    
    !---------------------------------------------------------------------------
    ! FVM Solver type
    !---------------------------------------------------------------------------
    type, public :: FVMSolver
        ! Grid reference
        type(Grid3D), pointer :: grid => null()
        
        ! Flow properties
        real(dp) :: inlet_velocity = 5.0_dp   ! m/s
        real(dp) :: viscosity = 1.5e-5_dp     ! m²/s (air at 20°C)
        real(dp) :: density = 1.225_dp        ! kg/m³
        
        ! Field variables
        real(dp), allocatable :: u(:,:,:)     ! x-velocity
        real(dp), allocatable :: v(:,:,:)     ! y-velocity
        real(dp), allocatable :: w(:,:,:)     ! z-velocity
        real(dp), allocatable :: p(:,:,:)     ! pressure
        
        ! Solid mask
        logical, allocatable :: solid(:,:,:)
        
        ! Solver parameters - tuned for moderate grids
        real(dp) :: alpha_u = 0.5_dp          ! Velocity under-relaxation
        real(dp) :: alpha_p = 0.1_dp          ! Pressure under-relaxation
        integer :: max_pressure_iters = 100   ! Pressure iterations
        real(dp) :: pressure_tol = 1.0e-4_dp
        
        ! Statistics
        integer :: step_count = 0
        real(dp) :: current_dt = 0.0_dp
        
    contains
        procedure :: init => solver_init
        procedure :: set_solid_mask => solver_set_solid_mask
        procedure :: step => solver_step
        procedure :: compute_drag => solver_compute_drag
        procedure :: get_max_velocity => solver_get_max_velocity
        procedure :: cleanup => solver_cleanup
        
        ! Private procedures
        procedure, private :: apply_bc => solver_apply_bc
        procedure, private :: compute_dt => solver_compute_dt
        procedure, private :: pressure_correction => solver_pressure_correction
        
    end type FVMSolver
    
contains

    !---------------------------------------------------------------------------
    ! Initialize solver
    !---------------------------------------------------------------------------
    subroutine solver_init(self, grid, inlet_velocity)
        class(FVMSolver), intent(inout) :: self
        type(Grid3D), target, intent(in) :: grid
        real(dp), intent(in), optional :: inlet_velocity
        
        integer :: nx, ny, nz
        
        self%grid => grid
        if (present(inlet_velocity)) self%inlet_velocity = inlet_velocity
        
        nx = grid%nx
        ny = grid%ny
        nz = grid%nz
        
        ! Allocate field arrays
        allocate(self%u(nx, ny, nz))
        allocate(self%v(nx, ny, nz))
        allocate(self%w(nx, ny, nz))
        allocate(self%p(nx, ny, nz))
        allocate(self%solid(nx, ny, nz))
        
        ! Initialize fields
        self%u = self%inlet_velocity
        self%v = 0.0_dp
        self%w = 0.0_dp
        self%p = 0.0_dp
        self%solid = .false.
        
        self%step_count = 0
        
    end subroutine solver_init
    
    !---------------------------------------------------------------------------
    ! Set solid geometry mask
    !---------------------------------------------------------------------------
    subroutine solver_set_solid_mask(self, mask)
        class(FVMSolver), intent(inout) :: self
        logical, intent(in) :: mask(:,:,:)
        
        self%solid = mask
        
        ! Zero velocity in solid cells
        where (self%solid)
            self%u = 0.0_dp
            self%v = 0.0_dp
            self%w = 0.0_dp
        end where
        
    end subroutine solver_set_solid_mask
    
    !---------------------------------------------------------------------------
    ! Compute stable time step (CFL condition)
    !---------------------------------------------------------------------------
    real(dp) function solver_compute_dt(self) result(dt)
        class(FVMSolver), intent(in) :: self
        
        real(dp) :: u_max, dx_min, dy_min, dz_min, d_min
        real(dp) :: cfl
        
        cfl = 0.2_dp  ! Moderate CFL for stability
        
        u_max = max(maxval(abs(self%u)), maxval(abs(self%v)), maxval(abs(self%w)), 0.1_dp)
        dx_min = minval(self%grid%dx)
        dy_min = minval(self%grid%dy)
        dz_min = minval(self%grid%dz)
        d_min = min(dx_min, dy_min, dz_min)
        
        ! CFL condition
        dt = cfl * d_min / u_max
        
        ! Viscous stability limit
        dt = min(dt, 0.25_dp * d_min**2 / self%viscosity)
        
    end function solver_compute_dt
    
    !---------------------------------------------------------------------------
    ! Apply boundary conditions
    !---------------------------------------------------------------------------
    subroutine solver_apply_bc(self)
        class(FVMSolver), intent(inout) :: self
        
        integer :: nx, ny, nz
        
        nx = self%grid%nx
        ny = self%grid%ny
        nz = self%grid%nz
        
        ! Inlet (x = 1): Fixed velocity
        self%u(1, :, :) = self%inlet_velocity
        self%v(1, :, :) = 0.0_dp
        self%w(1, :, :) = 0.0_dp
        
        ! Outlet (x = nx): Zero gradient
        self%u(nx, :, :) = self%u(nx-1, :, :)
        self%v(nx, :, :) = self%v(nx-1, :, :)
        self%w(nx, :, :) = self%w(nx-1, :, :)
        self%p(nx, :, :) = 0.0_dp  ! Reference pressure
        
        ! Bottom wall (z = 1): No-slip
        self%u(:, :, 1) = 0.0_dp
        self%v(:, :, 1) = 0.0_dp
        self%w(:, :, 1) = 0.0_dp
        
        ! Top (z = nz): Free slip
        self%u(:, :, nz) = self%u(:, :, nz-1)
        self%v(:, :, nz) = self%v(:, :, nz-1)
        self%w(:, :, nz) = 0.0_dp
        
        ! Side walls (y = 1, ny): Symmetric
        self%u(:, 1, :) = self%u(:, 2, :)
        self%u(:, ny, :) = self%u(:, ny-1, :)
        self%v(:, 1, :) = 0.0_dp
        self%v(:, ny, :) = 0.0_dp
        self%w(:, 1, :) = self%w(:, 2, :)
        self%w(:, ny, :) = self%w(:, ny-1, :)
        
        ! Solid cells: No-slip
        where (self%solid)
            self%u = 0.0_dp
            self%v = 0.0_dp
            self%w = 0.0_dp
        end where
        
    end subroutine solver_apply_bc
    
    !---------------------------------------------------------------------------
    ! Pressure correction (simplified SIMPLE-like)
    !---------------------------------------------------------------------------
    subroutine solver_pressure_correction(self)
        class(FVMSolver), intent(inout) :: self
        
        integer :: i, j, k, iter, nx, ny, nz
        real(dp) :: div, coeff
        real(dp) :: dx, dy, dz
        real(dp), allocatable :: p_new(:,:,:), div_field(:,:,:)
        real(dp) :: residual, max_div
        
        nx = self%grid%nx
        ny = self%grid%ny
        nz = self%grid%nz
        
        allocate(p_new(nx, ny, nz))
        allocate(div_field(nx, ny, nz))
        
        ! Compute divergence field
        div_field = 0.0_dp
        
        !$omp parallel do collapse(3) private(i,j,k,dx,dy,dz,div)
        do k = 2, nz-1
            do j = 2, ny-1
                do i = 2, nx-1
                    if (self%solid(i,j,k)) cycle
                    
                    dx = self%grid%dx(i)
                    dy = self%grid%dy(j)
                    dz = self%grid%dz(k)
                    
                    div = (self%u(i+1,j,k) - self%u(i-1,j,k)) / (2.0_dp * dx) + &
                          (self%v(i,j+1,k) - self%v(i,j-1,k)) / (2.0_dp * dy) + &
                          (self%w(i,j,k+1) - self%w(i,j,k-1)) / (2.0_dp * dz)
                    
                    div_field(i,j,k) = div
                end do
            end do
        end do
        !$omp end parallel do
        
        ! Jacobi iteration for pressure
        do iter = 1, self%max_pressure_iters
            p_new = self%p
            max_div = 0.0_dp
            
            !$omp parallel do collapse(3) private(i,j,k,dx,dy,dz,coeff) reduction(max:max_div)
            do k = 2, nz-1
                do j = 2, ny-1
                    do i = 2, nx-1
                        if (self%solid(i,j,k)) cycle
                        
                        dx = self%grid%dx(i)
                        dy = self%grid%dy(j)
                        dz = self%grid%dz(k)
                        
                        coeff = 2.0_dp/dx**2 + 2.0_dp/dy**2 + 2.0_dp/dz**2
                        
                        p_new(i,j,k) = ( &
                            (self%p(i+1,j,k) + self%p(i-1,j,k)) / dx**2 + &
                            (self%p(i,j+1,k) + self%p(i,j-1,k)) / dy**2 + &
                            (self%p(i,j,k+1) + self%p(i,j,k-1)) / dz**2 - &
                            self%density * div_field(i,j,k) / self%current_dt &
                        ) / coeff
                        
                        max_div = max(max_div, abs(div_field(i,j,k)))
                    end do
                end do
            end do
            !$omp end parallel do
            
            ! Under-relaxation
            self%p = self%alpha_p * p_new + (1.0_dp - self%alpha_p) * self%p
            
            if (max_div < self%pressure_tol) exit
        end do
        
        ! Correct velocities with pressure gradient
        !$omp parallel do collapse(3) private(i,j,k,dx,dy,dz)
        do k = 2, nz-1
            do j = 2, ny-1
                do i = 2, nx-1
                    if (self%solid(i,j,k)) cycle
                    
                    dx = self%grid%dx(i)
                    dy = self%grid%dy(j)
                    dz = self%grid%dz(k)
                    
                    self%u(i,j,k) = self%u(i,j,k) - &
                        self%current_dt / self%density * &
                        (self%p(i+1,j,k) - self%p(i-1,j,k)) / (2.0_dp * dx)
                    
                    self%v(i,j,k) = self%v(i,j,k) - &
                        self%current_dt / self%density * &
                        (self%p(i,j+1,k) - self%p(i,j-1,k)) / (2.0_dp * dy)
                    
                    self%w(i,j,k) = self%w(i,j,k) - &
                        self%current_dt / self%density * &
                        (self%p(i,j,k+1) - self%p(i,j,k-1)) / (2.0_dp * dz)
                end do
            end do
        end do
        !$omp end parallel do
        
        deallocate(p_new, div_field)
        
    end subroutine solver_pressure_correction
    
    !---------------------------------------------------------------------------
    ! Perform one time step
    !---------------------------------------------------------------------------
    subroutine solver_step(self, dt_out)
        class(FVMSolver), intent(inout) :: self
        real(dp), intent(out), optional :: dt_out
        
        integer :: i, j, k, nx, ny, nz
        real(dp) :: dt
        real(dp) :: dx, dy, dz, dz_m, dz_p
        real(dp) :: dudx, dudy, dudz
        real(dp) :: dvdx, dvdy, dvdz
        real(dp) :: dwdx, dwdy, dwdz
        real(dp) :: d2u, d2v, d2w
        real(dp) :: conv_u, conv_v, conv_w
        real(dp) :: diff_u, diff_v, diff_w
        real(dp), allocatable :: u_new(:,:,:), v_new(:,:,:), w_new(:,:,:)
        
        nx = self%grid%nx
        ny = self%grid%ny
        nz = self%grid%nz
        
        ! Compute time step
        dt = self%compute_dt()
        self%current_dt = dt
        
        allocate(u_new(nx, ny, nz))
        allocate(v_new(nx, ny, nz))
        allocate(w_new(nx, ny, nz))
        
        u_new = self%u
        v_new = self%v
        w_new = self%w
        
        ! Update interior cells
        !$omp parallel do collapse(3) private(i,j,k,dx,dy,dz,dz_m,dz_p, &
        !$omp    dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz, &
        !$omp    conv_u,conv_v,conv_w,d2u,d2v,d2w,diff_u,diff_v,diff_w)
        do k = 2, nz-1
            do j = 2, ny-1
                do i = 2, nx-1
                    if (self%solid(i,j,k)) cycle
                    
                    dx = self%grid%dx(i)
                    dy = self%grid%dy(j)
                    dz = self%grid%dz(k)
                    dz_m = self%grid%dz(max(1,k-1))
                    dz_p = self%grid%dz(min(nz,k+1))
                    
                    ! Convective terms (upwind)
                    ! u derivatives
                    if (self%u(i,j,k) > 0) then
                        dudx = (self%u(i,j,k) - self%u(i-1,j,k)) / dx
                    else
                        dudx = (self%u(i+1,j,k) - self%u(i,j,k)) / dx
                    end if
                    
                    if (self%v(i,j,k) > 0) then
                        dudy = (self%u(i,j,k) - self%u(i,j-1,k)) / dy
                    else
                        dudy = (self%u(i,j+1,k) - self%u(i,j,k)) / dy
                    end if
                    
                    if (self%w(i,j,k) > 0) then
                        dudz = (self%u(i,j,k) - self%u(i,j,k-1)) / ((dz + dz_m)/2)
                    else
                        dudz = (self%u(i,j,k+1) - self%u(i,j,k)) / ((dz + dz_p)/2)
                    end if
                    
                    conv_u = self%u(i,j,k)*dudx + self%v(i,j,k)*dudy + self%w(i,j,k)*dudz
                    
                    ! v derivatives (similar)
                    if (self%u(i,j,k) > 0) then
                        dvdx = (self%v(i,j,k) - self%v(i-1,j,k)) / dx
                    else
                        dvdx = (self%v(i+1,j,k) - self%v(i,j,k)) / dx
                    end if
                    
                    if (self%v(i,j,k) > 0) then
                        dvdy = (self%v(i,j,k) - self%v(i,j-1,k)) / dy
                    else
                        dvdy = (self%v(i,j+1,k) - self%v(i,j,k)) / dy
                    end if
                    
                    if (self%w(i,j,k) > 0) then
                        dvdz = (self%v(i,j,k) - self%v(i,j,k-1)) / ((dz + dz_m)/2)
                    else
                        dvdz = (self%v(i,j,k+1) - self%v(i,j,k)) / ((dz + dz_p)/2)
                    end if
                    
                    conv_v = self%u(i,j,k)*dvdx + self%v(i,j,k)*dvdy + self%w(i,j,k)*dvdz
                    
                    ! w derivatives
                    if (self%u(i,j,k) > 0) then
                        dwdx = (self%w(i,j,k) - self%w(i-1,j,k)) / dx
                    else
                        dwdx = (self%w(i+1,j,k) - self%w(i,j,k)) / dx
                    end if
                    
                    if (self%v(i,j,k) > 0) then
                        dwdy = (self%w(i,j,k) - self%w(i,j-1,k)) / dy
                    else
                        dwdy = (self%w(i,j+1,k) - self%w(i,j,k)) / dy
                    end if
                    
                    if (self%w(i,j,k) > 0) then
                        dwdz = (self%w(i,j,k) - self%w(i,j,k-1)) / ((dz + dz_m)/2)
                    else
                        dwdz = (self%w(i,j,k+1) - self%w(i,j,k)) / ((dz + dz_p)/2)
                    end if
                    
                    conv_w = self%u(i,j,k)*dwdx + self%v(i,j,k)*dwdy + self%w(i,j,k)*dwdz
                    
                    ! Diffusion terms (central difference)
                    d2u = (self%u(i+1,j,k) - 2*self%u(i,j,k) + self%u(i-1,j,k)) / dx**2 + &
                          (self%u(i,j+1,k) - 2*self%u(i,j,k) + self%u(i,j-1,k)) / dy**2 + &
                          (self%u(i,j,k+1) - 2*self%u(i,j,k) + self%u(i,j,k-1)) / dz**2
                    
                    d2v = (self%v(i+1,j,k) - 2*self%v(i,j,k) + self%v(i-1,j,k)) / dx**2 + &
                          (self%v(i,j+1,k) - 2*self%v(i,j,k) + self%v(i,j-1,k)) / dy**2 + &
                          (self%v(i,j,k+1) - 2*self%v(i,j,k) + self%v(i,j,k-1)) / dz**2
                    
                    d2w = (self%w(i+1,j,k) - 2*self%w(i,j,k) + self%w(i-1,j,k)) / dx**2 + &
                          (self%w(i,j+1,k) - 2*self%w(i,j,k) + self%w(i,j-1,k)) / dy**2 + &
                          (self%w(i,j,k+1) - 2*self%w(i,j,k) + self%w(i,j,k-1)) / dz**2
                    
                    diff_u = self%viscosity * d2u
                    diff_v = self%viscosity * d2v
                    diff_w = self%viscosity * d2w
                    
                    ! Time integration
                    u_new(i,j,k) = self%u(i,j,k) + dt * (-conv_u + diff_u)
                    v_new(i,j,k) = self%v(i,j,k) + dt * (-conv_v + diff_v)
                    w_new(i,j,k) = self%w(i,j,k) + dt * (-conv_w + diff_w)
                end do
            end do
        end do
        !$omp end parallel do
        
        ! Update velocities with under-relaxation
        self%u = self%alpha_u * u_new + (1.0_dp - self%alpha_u) * self%u
        self%v = self%alpha_u * v_new + (1.0_dp - self%alpha_u) * self%v
        self%w = self%alpha_u * w_new + (1.0_dp - self%alpha_u) * self%w
        
        ! Clamp velocities to prevent runaway (safety net only)
        self%u = max(-100.0_dp, min(100.0_dp, self%u))
        self%v = max(-100.0_dp, min(100.0_dp, self%v))
        self%w = max(-100.0_dp, min(100.0_dp, self%w))
        
        deallocate(u_new, v_new, w_new)
        
        ! Pressure correction
        call self%pressure_correction()
        
        ! Clamp pressure (safety net only)
        self%p = max(-1.0e5_dp, min(1.0e5_dp, self%p))
        
        ! Apply boundary conditions
        call self%apply_bc()
        
        self%step_count = self%step_count + 1
        
        if (present(dt_out)) dt_out = dt
        
    end subroutine solver_step
    
    !---------------------------------------------------------------------------
    ! Compute drag force on solid body
    !---------------------------------------------------------------------------
    subroutine solver_compute_drag(self, Fx, Fy, Fz)
        class(FVMSolver), intent(in) :: self
        real(dp), intent(out) :: Fx, Fy, Fz
        
        integer :: i, j, k, nx, ny, nz
        real(dp) :: dx, dy, dz, area
        real(dp) :: Fx_local, Fy_local, Fz_local
        
        nx = self%grid%nx
        ny = self%grid%ny
        nz = self%grid%nz
        
        Fx = 0.0_dp
        Fy = 0.0_dp
        Fz = 0.0_dp
        
        !$omp parallel do collapse(3) private(i,j,k,dx,dy,dz,area) &
        !$omp reduction(+:Fx,Fy,Fz)
        do k = 2, nz-1
            do j = 2, ny-1
                do i = 2, nx-1
                    if (.not. self%solid(i,j,k)) cycle
                    
                    dx = self%grid%dx(i)
                    dy = self%grid%dy(j)
                    dz = self%grid%dz(k)
                    
                    ! X-faces
                    if (i > 1 .and. .not. self%solid(i-1,j,k)) then
                        area = dy * dz
                        Fx = Fx - self%p(i-1,j,k) * area
                        Fx = Fx - self%density * self%viscosity * &
                             self%u(i-1,j,k) / dx * area
                    end if
                    
                    if (i < nx .and. .not. self%solid(i+1,j,k)) then
                        area = dy * dz
                        Fx = Fx + self%p(i+1,j,k) * area
                        Fx = Fx + self%density * self%viscosity * &
                             self%u(i+1,j,k) / dx * area
                    end if
                    
                    ! Y-faces
                    if (j > 1 .and. .not. self%solid(i,j-1,k)) then
                        area = dx * dz
                        Fy = Fy - self%p(i,j-1,k) * area
                    end if
                    
                    if (j < ny .and. .not. self%solid(i,j+1,k)) then
                        area = dx * dz
                        Fy = Fy + self%p(i,j+1,k) * area
                    end if
                    
                    ! Z-faces  
                    if (k > 1 .and. .not. self%solid(i,j,k-1)) then
                        area = dx * dy
                        Fz = Fz - self%p(i,j,k-1) * area
                    end if
                    
                    if (k < nz .and. .not. self%solid(i,j,k+1)) then
                        area = dx * dy
                        Fz = Fz + self%p(i,j,k+1) * area
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        
    end subroutine solver_compute_drag
    
    !---------------------------------------------------------------------------
    ! Get maximum velocity magnitude
    !---------------------------------------------------------------------------
    real(dp) function solver_get_max_velocity(self) result(u_max)
        class(FVMSolver), intent(in) :: self
        
        real(dp) :: speed_sq
        integer :: i, j, k
        
        u_max = 0.0_dp
        
        !$omp parallel do collapse(3) private(i,j,k,speed_sq) reduction(max:u_max)
        do k = 1, self%grid%nz
            do j = 1, self%grid%ny
                do i = 1, self%grid%nx
                    if (self%solid(i,j,k)) cycle
                    speed_sq = self%u(i,j,k)**2 + self%v(i,j,k)**2 + self%w(i,j,k)**2
                    u_max = max(u_max, sqrt(speed_sq))
                end do
            end do
        end do
        !$omp end parallel do
        
    end function solver_get_max_velocity
    
    !---------------------------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------------------------
    subroutine solver_cleanup(self)
        class(FVMSolver), intent(inout) :: self
        
        if (allocated(self%u)) deallocate(self%u)
        if (allocated(self%v)) deallocate(self%v)
        if (allocated(self%w)) deallocate(self%w)
        if (allocated(self%p)) deallocate(self%p)
        if (allocated(self%solid)) deallocate(self%solid)
        
        self%grid => null()
        
    end subroutine solver_cleanup
    
end module solver_module
