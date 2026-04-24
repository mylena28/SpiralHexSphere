! params.f90 — Global parameters for the SpiralHexSphere mesh generator
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
    integer  :: num_sides = 6       ! Base polygon sides (3 → 6 pts/ring)
                                    !   pts_per_ring = num_sides * 2

    ! Spiral discretisation (Zinchenko Level-1 refinement)
    !   total_rings  = (n_s - 1) * 2 + 1   (doubled rings for midpoints)
    !   pts_per_ring = num_sides * 2         (doubled pts for midpoints)
    !integer :: n_s = 36            ! Number of original cross-sections (>= 2)
    ! n_s is computed at runtime in main.f90 (equilateral-triangle condition)
    integer :: n_s = 0              ! placeholder; set in main before use

    !  Hemisphere end-cap
    !   The cap now uses Zinchenko Level-1 icosahedral subdivision (same method
    !   as the spiral body).  n_hemi is kept for reference but not used by geometry.f90.
    integer :: n_hemi = 0          ! (unused — cap uses icosahedral subdivision)

    !  Mesh refinement level (0 = base mesh; each level splits every triangle into 4)
    integer :: refine_level = 0

    !  Output file
    character(len=256) :: outfile = 'spiral.dat'

end module params
