module params
    implicit none
    integer,  parameter :: dp = kind(1.0d0)
    real(dp), parameter :: PI = acos(-1.0_dp)

    ! Helix (centreline)
    real(dp) :: helix_R = 1.0_dp   ! Helix path radius
    real(dp) :: helix_h = 6.0_dp   ! Total axial height
    real(dp) :: n_turns = 4.0_dp   ! Number of full turns

    !  Tube cross-section
    real(dp) :: tube_r   = 0.3_dp  ! Tube cross-section radius
    integer  :: num_sides = 6       ! Base polygon sides (6 pts/ring)

    ! Spiral discretisation
    integer :: n_s = 0              ! placeholder; set in main before use

    !  Hemisphere end-cap
    integer :: n_hemi = 0          ! (unused — cap uses icosahedral subdivision)

    integer :: refine_level = 3

    !  Output file
    character(len=256) :: outfile = 'spiral.obj'

end module params
