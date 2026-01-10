!===============================================================================
! grid_module.f90
! Non-uniform 3D grid for CFD simulation
!===============================================================================
module grid_module
    implicit none
    private
    
    ! Double precision
    integer, parameter, public :: dp = selected_real_kind(15, 307)
    
    !---------------------------------------------------------------------------
    ! Non-uniform grid type
    !---------------------------------------------------------------------------
    type, public :: Grid3D
        ! Domain dimensions (meters)
        real(dp) :: length = 2.0_dp    ! x direction
        real(dp) :: width = 1.7_dp     ! y direction  
        real(dp) :: height = 2.5_dp    ! z direction (reduced from 3.0m)
        
        ! Number of cells
        integer :: nx = 60
        integer :: ny = 45
        integer :: nz = 350  ! 350 levels for finer resolution near visor
        
        ! Refinement zone (around head at ~1.8m)
        real(dp) :: z_refine_min = 1.6_dp
        real(dp) :: z_refine_max = 2.1_dp
        
        ! Grid coordinates (cell centers)
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: z(:)
        
        ! Cell sizes
        real(dp), allocatable :: dx(:)
        real(dp), allocatable :: dy(:)
        real(dp), allocatable :: dz(:)
        
        ! Grid info
        integer :: n_below, n_fine, n_above
        real(dp) :: dz_min, dz_max
        
    contains
        procedure :: init => grid_init
        procedure :: print_info => grid_print_info
        procedure :: cell_volume => grid_cell_volume
        procedure :: cleanup => grid_cleanup
        
    end type Grid3D
    
contains

    !---------------------------------------------------------------------------
    ! Initialize grid with non-uniform z spacing
    !---------------------------------------------------------------------------
    subroutine grid_init(self, nx, ny, nz, length, width, height, &
                         z_refine_min, z_refine_max)
        class(Grid3D), intent(inout) :: self
        integer, intent(in), optional :: nx, ny, nz
        real(dp), intent(in), optional :: length, width, height
        real(dp), intent(in), optional :: z_refine_min, z_refine_max
        
        integer :: i, j, k
        real(dp) :: refine_height, below_height, above_height
        real(dp) :: target_dz_fine
        integer :: min_n_fine, remaining
        real(dp) :: below_frac
        real(dp), allocatable :: z_faces(:)
        
        ! Set optional parameters
        if (present(nx)) self%nx = nx
        if (present(ny)) self%ny = ny
        if (present(nz)) self%nz = nz
        if (present(length)) self%length = length
        if (present(width)) self%width = width
        if (present(height)) self%height = height
        if (present(z_refine_min)) self%z_refine_min = z_refine_min
        if (present(z_refine_max)) self%z_refine_max = z_refine_max
        
        ! Allocate arrays
        allocate(self%x(self%nx))
        allocate(self%y(self%ny))
        allocate(self%z(self%nz))
        allocate(self%dx(self%nx))
        allocate(self%dy(self%ny))
        allocate(self%dz(self%nz))
        
        ! Generate uniform x grid
        do i = 1, self%nx
            self%dx(i) = self%length / real(self%nx, dp)
            self%x(i) = (real(i, dp) - 0.5_dp) * self%dx(i)
        end do
        
        ! Generate uniform y grid
        do j = 1, self%ny
            self%dy(j) = self%width / real(self%ny, dp)
            self%y(j) = (real(j, dp) - 0.5_dp) * self%dy(j)
        end do
        
        ! Generate non-uniform z grid
        refine_height = self%z_refine_max - self%z_refine_min
        below_height = self%z_refine_min
        above_height = self%height - self%z_refine_max
        
        ! Target 1.25mm in refined zone to resolve 5mm visor with 4 cells
        target_dz_fine = 1.25e-3_dp
        min_n_fine = ceiling(refine_height / target_dz_fine)
        remaining = self%nz - min_n_fine
        
        if (remaining < 10) then
            ! Not enough cells - favor refined zone
            self%n_fine = int(0.80_dp * real(self%nz, dp))
            remaining = self%nz - self%n_fine
            below_frac = below_height / (below_height + above_height)
            self%n_below = max(3, int(real(remaining, dp) * below_frac))
            self%n_above = max(3, remaining - self%n_below)
            self%n_fine = self%nz - self%n_below - self%n_above
        else
            self%n_fine = min_n_fine
            below_frac = below_height / (below_height + above_height)
            self%n_below = max(3, int(real(remaining, dp) * below_frac))
            self%n_above = max(3, remaining - self%n_below)
            self%n_fine = self%nz - self%n_below - self%n_above
        end if
        
        print '(A,I4,A,I4,A,I4,A,I4)', 'Grid distribution: ', self%n_below, &
              ' below + ', self%n_fine, ' refined + ', self%n_above, ' above = ', self%nz
        
        ! Generate z faces
        allocate(z_faces(self%nz + 1))
        
        ! Below refinement zone
        do k = 1, self%n_below
            z_faces(k) = (real(k-1, dp) / real(self%n_below, dp)) * self%z_refine_min
        end do
        
        ! Refined zone
        do k = 1, self%n_fine
            z_faces(self%n_below + k) = self%z_refine_min + &
                (real(k-1, dp) / real(self%n_fine, dp)) * refine_height
        end do
        
        ! Above refinement zone
        do k = 1, self%n_above + 1
            z_faces(self%n_below + self%n_fine + k) = self%z_refine_max + &
                (real(k-1, dp) / real(self%n_above, dp)) * above_height
        end do
        
        ! Cell centers and sizes
        do k = 1, self%nz
            self%z(k) = 0.5_dp * (z_faces(k) + z_faces(k+1))
            self%dz(k) = z_faces(k+1) - z_faces(k)
        end do
        
        self%dz_min = minval(self%dz)
        self%dz_max = maxval(self%dz)
        
        deallocate(z_faces)
        
    end subroutine grid_init
    
    !---------------------------------------------------------------------------
    ! Print grid information
    !---------------------------------------------------------------------------
    subroutine grid_print_info(self)
        class(Grid3D), intent(in) :: self
        
        print '(A)', '=================================================='
        print '(A)', 'Grid Configuration'
        print '(A)', '=================================================='
        print '(A,I4,A,I4,A,I4)', '  Grid size: ', self%nx, ' x ', self%ny, ' x ', self%nz
        print '(A,F5.2,A,F5.2,A,F5.2,A)', '  Domain: ', self%length, 'm x ', &
              self%width, 'm x ', self%height, 'm'
        print '(A,F6.2,A)', '  dx (uniform): ', self%dx(1)*1000.0_dp, ' mm'
        print '(A,F6.2,A)', '  dy (uniform): ', self%dy(1)*1000.0_dp, ' mm'
        print '(A,F6.2,A,F6.2,A)', '  dz range: ', self%dz_min*1000.0_dp, ' - ', &
              self%dz_max*1000.0_dp, ' mm'
        print '(A,I4,A,I4,A,I4)', '  Z distribution: ', self%n_below, ' below + ', &
              self%n_fine, ' refined + ', self%n_above, ' above'
        print '(A,F5.2,A,F5.2,A)', '  Refinement zone: ', self%z_refine_min, 'm to ', &
              self%z_refine_max, 'm'
        print '(A)', ''
        
    end subroutine grid_print_info
    
    !---------------------------------------------------------------------------
    ! Get cell volume
    !---------------------------------------------------------------------------
    real(dp) function grid_cell_volume(self, i, j, k) result(vol)
        class(Grid3D), intent(in) :: self
        integer, intent(in) :: i, j, k
        
        vol = self%dx(i) * self%dy(j) * self%dz(k)
        
    end function grid_cell_volume
    
    !---------------------------------------------------------------------------
    ! Cleanup
    !---------------------------------------------------------------------------
    subroutine grid_cleanup(self)
        class(Grid3D), intent(inout) :: self
        
        if (allocated(self%x)) deallocate(self%x)
        if (allocated(self%y)) deallocate(self%y)
        if (allocated(self%z)) deallocate(self%z)
        if (allocated(self%dx)) deallocate(self%dx)
        if (allocated(self%dy)) deallocate(self%dy)
        if (allocated(self%dz)) deallocate(self%dz)
        
    end subroutine grid_cleanup
    
end module grid_module
