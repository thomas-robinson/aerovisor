!===============================================================================
! geometry_module.f90
! Runner and visor geometry for CFD simulation
!===============================================================================
module geometry_module
    use grid_module, only: dp, Grid3D
    implicit none
    private
    
    !---------------------------------------------------------------------------
    ! Runner geometry type
    !---------------------------------------------------------------------------
    type, public :: RunnerGeometry
        real(dp) :: height = 2.0_dp        ! Total height (m)
        real(dp) :: position(3) = [0.8_dp, 0.85_dp, 0.0_dp]  ! Base position
        
        ! Body dimensions (derived from height)
        real(dp) :: head_radius
        real(dp) :: head_center_z
        real(dp) :: neck_radius, neck_height
        real(dp) :: torso_width, torso_depth, torso_height, torso_bottom
        real(dp) :: leg_radius, leg_length
        real(dp) :: arm_radius, arm_length
        real(dp) :: eye_height
        
    contains
        procedure :: init => runner_init
        procedure :: is_inside => runner_is_inside
        
    end type RunnerGeometry
    
    !---------------------------------------------------------------------------
    ! Visor geometry type
    !---------------------------------------------------------------------------
    type, public :: VisorGeometry
        real(dp) :: angle_rad = 0.0_dp     ! Angle from horizontal (radians)
        real(dp) :: visor_length = 0.08_dp  ! 8cm brim length
        real(dp) :: visor_width = 0.16_dp   ! 16cm width
        real(dp) :: visor_thickness = 0.005_dp  ! 5mm thickness (realistic)
        real(dp) :: attachment_height
        real(dp) :: attachment_x
        real(dp) :: runner_y
        
    contains
        procedure :: init => visor_init
        procedure :: is_inside => visor_is_inside
        
    end type VisorGeometry
    
    ! Public procedures
    public :: create_geometry_mask
    
contains

    !---------------------------------------------------------------------------
    ! Initialize runner geometry
    !---------------------------------------------------------------------------
    subroutine runner_init(self, height, position)
        class(RunnerGeometry), intent(inout) :: self
        real(dp), intent(in), optional :: height
        real(dp), intent(in), optional :: position(3)
        
        if (present(height)) self%height = height
        if (present(position)) self%position = position
        
        ! Calculate body dimensions based on height
        self%head_radius = self%height * 0.06_dp
        self%head_center_z = self%height - self%head_radius + self%position(3)
        
        self%neck_radius = self%head_radius * 0.5_dp
        self%neck_height = self%height * 0.03_dp
        
        self%torso_width = self%height * 0.22_dp
        self%torso_depth = self%height * 0.12_dp
        self%torso_height = self%height * 0.30_dp
        self%torso_bottom = self%height * 0.52_dp + self%position(3)
        
        self%leg_radius = self%height * 0.055_dp
        self%leg_length = self%height * 0.50_dp
        
        self%arm_radius = self%height * 0.035_dp
        self%arm_length = self%height * 0.30_dp
        
        self%eye_height = self%head_center_z
        
    end subroutine runner_init
    
    !---------------------------------------------------------------------------
    ! Check if point is inside runner
    !---------------------------------------------------------------------------
    logical function runner_is_inside(self, x, y, z) result(inside)
        class(RunnerGeometry), intent(in) :: self
        real(dp), intent(in) :: x, y, z
        
        real(dp) :: rel_x, rel_y, rel_z, r2
        real(dp) :: torso_top, torso_center_y
        real(dp) :: leg_sep, leg1_y, leg2_y
        real(dp) :: arm_sep, arm1_y, arm2_y, arm_top
        
        inside = .false.
        
        rel_x = x - self%position(1)
        rel_y = y - self%position(2)
        rel_z = z - self%position(3)
        
        ! Head (sphere)
        r2 = rel_x**2 + rel_y**2 + (z - self%head_center_z)**2
        if (r2 <= self%head_radius**2) then
            inside = .true.
            return
        end if
        
        ! Neck (cylinder)
        if (z >= self%head_center_z - self%head_radius - self%neck_height .and. &
            z <= self%head_center_z - self%head_radius * 0.7_dp) then
            r2 = rel_x**2 + rel_y**2
            if (r2 <= self%neck_radius**2) then
                inside = .true.
                return
            end if
        end if
        
        ! Torso (ellipsoid)
        torso_top = self%torso_bottom + self%torso_height
        torso_center_y = self%position(2)
        if (z >= self%torso_bottom .and. z <= torso_top) then
            r2 = (rel_x / self%torso_depth)**2 + (rel_y / self%torso_width)**2
            if (r2 <= 1.0_dp) then
                inside = .true.
                return
            end if
        end if
        
        ! Legs (cylinders)
        if (z >= 0.0_dp .and. z <= self%torso_bottom + 0.05_dp) then
            leg_sep = self%torso_width * 0.35_dp
            leg1_y = self%position(2) - leg_sep
            leg2_y = self%position(2) + leg_sep
            
            r2 = rel_x**2 + (y - leg1_y)**2
            if (r2 <= self%leg_radius**2) then
                inside = .true.
                return
            end if
            
            r2 = rel_x**2 + (y - leg2_y)**2
            if (r2 <= self%leg_radius**2) then
                inside = .true.
                return
            end if
        end if
        
        ! Arms (cylinders, slightly behind torso)
        arm_top = self%torso_bottom + self%torso_height * 0.95_dp
        if (z >= arm_top - self%arm_length .and. z <= arm_top) then
            arm_sep = self%torso_width * 0.6_dp
            arm1_y = self%position(2) - arm_sep
            arm2_y = self%position(2) + arm_sep
            
            r2 = (rel_x + self%torso_depth * 0.3_dp)**2 + (y - arm1_y)**2
            if (r2 <= self%arm_radius**2) then
                inside = .true.
                return
            end if
            
            r2 = (rel_x + self%torso_depth * 0.3_dp)**2 + (y - arm2_y)**2
            if (r2 <= self%arm_radius**2) then
                inside = .true.
                return
            end if
        end if
        
    end function runner_is_inside
    
    !---------------------------------------------------------------------------
    ! Initialize visor geometry
    !---------------------------------------------------------------------------
    subroutine visor_init(self, runner, angle_degrees)
        class(VisorGeometry), intent(inout) :: self
        type(RunnerGeometry), intent(in) :: runner
        real(dp), intent(in) :: angle_degrees
        
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        
        self%angle_rad = angle_degrees * pi / 180.0_dp
        self%attachment_height = runner%eye_height + 0.02_dp
        self%attachment_x = runner%head_radius * 0.95_dp + runner%position(1)
        self%runner_y = runner%position(2)
        
    end subroutine visor_init
    
    !---------------------------------------------------------------------------
    ! Check if point is inside visor
    !---------------------------------------------------------------------------
    logical function visor_is_inside(self, x, y, z) result(inside)
        class(VisorGeometry), intent(in) :: self
        real(dp), intent(in) :: x, y, z
        
        real(dp) :: rel_x, rel_y, rel_z
        real(dp) :: x_rot, z_rot
        real(dp) :: cos_a, sin_a
        
        inside = .false.
        
        ! Transform to visor-local coordinates
        rel_x = x - self%attachment_x
        rel_y = y - self%runner_y
        rel_z = z - self%attachment_height
        
        ! Rotate around y-axis (visor angle)
        cos_a = cos(self%angle_rad)
        sin_a = sin(self%angle_rad)
        x_rot = rel_x * cos_a + rel_z * sin_a
        z_rot = -rel_x * sin_a + rel_z * cos_a
        
        ! Check if inside visor bounding box (in rotated frame)
        if (x_rot >= 0.0_dp .and. x_rot <= self%visor_length .and. &
            abs(rel_y) <= self%visor_width / 2.0_dp .and. &
            abs(z_rot) <= self%visor_thickness / 2.0_dp) then
            inside = .true.
        end if
        
    end function visor_is_inside
    
    !---------------------------------------------------------------------------
    ! Create geometry mask on grid
    !---------------------------------------------------------------------------
    subroutine create_geometry_mask(grid, runner_pos, visor_angle, mask, use_visor)
        type(Grid3D), intent(in) :: grid
        real(dp), intent(in) :: runner_pos(3)
        real(dp), intent(in) :: visor_angle
        logical, intent(out) :: mask(grid%nx, grid%ny, grid%nz)
        logical, intent(in) :: use_visor
        
        type(RunnerGeometry) :: runner
        type(VisorGeometry) :: visor
        integer :: i, j, k
        integer :: runner_voxels, visor_voxels
        
        ! Initialize runner
        call runner%init(height=2.0_dp, position=runner_pos)
        
        ! Initialize visor
        if (use_visor) then
            call visor%init(runner, visor_angle)
        end if
        
        ! Create mask
        mask = .false.
        runner_voxels = 0
        visor_voxels = 0
        
        !$omp parallel do collapse(3) private(i,j,k) reduction(+:runner_voxels,visor_voxels)
        do k = 1, grid%nz
            do j = 1, grid%ny
                do i = 1, grid%nx
                    if (runner%is_inside(grid%x(i), grid%y(j), grid%z(k))) then
                        mask(i, j, k) = .true.
                        runner_voxels = runner_voxels + 1
                    else if (use_visor) then
                        if (visor%is_inside(grid%x(i), grid%y(j), grid%z(k))) then
                            mask(i, j, k) = .true.
                            visor_voxels = visor_voxels + 1
                        end if
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        
        print '(A,I8)', '  Runner voxels: ', runner_voxels
        if (use_visor) then
            print '(A,I8)', '  Visor voxels:  ', visor_voxels
            if (visor_voxels == 0) then
                print '(A,F6.1,A)', '  WARNING: Visor at ', visor_angle, ' degrees has 0 voxels!'
            end if
        end if
        
    end subroutine create_geometry_mask
    
end module geometry_module
