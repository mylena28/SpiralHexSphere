! main.f90 — SpiralHexSphere mesh generator
program main
    use params
    use geometry
    use io
    implicit none

    integer  :: total_nv, total_nf, nv_spiral
    integer  :: total_rings, pts_per_ring
    integer  :: lev, nv_new, nf_new, nv_sp_new
    !real(dp) :: proj_err
    !integer  :: proj_idx
    real(dp), allocatable :: verts(:,:), verts_new(:,:)
    integer,  allocatable :: faces(:,:), faces_new(:,:)
    real(dp), allocatable :: s_frac(:), s_frac_new(:)
    integer,  allocatable :: cap_flag(:), cap_flag_new(:)

    ! Equilateral-triangle condition: arc length per ring ≈ arc length per side
    ! 2*PI*helix_R*n_turns / n_s  =  2*PI*tube_r / num_sides
    ! => n_s = helix_R * n_turns * num_sides / tube_r
    n_s = nint(helix_R * n_turns * real(num_sides, dp) / tube_r)

    total_rings  = n_s
    pts_per_ring = num_sides

    write(*,'(A)') '======================================================'
    write(*,'(A)') ' Helix + 2 Hemisphere Caps with Zinchenko mesh'
    write(*,'(A)') '======================================================'
    write(*,'(A,F6.2)') '  Helix radius       : ', helix_R
    write(*,'(A,F6.2)') '  Helix height       : ', helix_h
    write(*,'(A,F6.2)') '  Number of turns    : ', n_turns
    write(*,'(A,F6.2)') '  Tube radius        : ', tube_r
    write(*,'(A,I4)')   '  Polygon sides      : ', num_sides
    write(*,'(A,I4)')   '  Pts / ring         : ', pts_per_ring
    write(*,'(A,I4)')   '  Original sections  : ', n_s
    write(*,'(A,I4)')   '  Total spiral rings : ', total_rings
    write(*,'(A,A)')    '  Cap method         : ', 'Level-0 hexagonal (6 interior pts)'
    write(*,'(A,I4)')   '  Cap new verts/cap  : ', 3
    write(*,'(A,I4)')   '  Cap faces/cap      : ', 10
    write(*,'(A,I4)')   '  Refinement level   : ', refine_level
    write(*,'(A)') '------------------------------------------------------'

    call generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces, s_frac, cap_flag)
    !call check_projection_error(total_nv, verts, s_frac, cap_flag, proj_err, proj_idx)
    !write(*,'(A,ES10.3,A,I0,A)') '  Projection error L0: ', proj_err, '  (vertex ', proj_idx, ')'

    ! Apply 4-1 midpoint subdivision refine_level times
    do lev = 1, refine_level
        write(*,'(A,I2,A)') '  Subdividing level ', lev, ' ...'
        call subdivide_mesh(total_nv, total_nf, nv_spiral, verts, faces, &
                            s_frac, cap_flag, &
                            nv_new, nf_new, nv_sp_new, verts_new, faces_new, &
                            s_frac_new, cap_flag_new)
        call move_alloc(verts_new, verts)
        call move_alloc(faces_new, faces)
        call move_alloc(s_frac_new, s_frac)
        call move_alloc(cap_flag_new, cap_flag)
        total_nv  = nv_new
        total_nf  = nf_new
        nv_spiral = nv_sp_new
        !call check_projection_error(total_nv, verts, s_frac, cap_flag, proj_err, proj_idx)
        !write(*,'(A,I2,A,ES10.3,A,I0,A)') &
        !    '  Projection error L', lev, ': ', proj_err, '  (vertex ', proj_idx, ')'
    end do
    write(*,'(A,I8)') '------------------------------------------------------'
    write(*,'(A,I8)') '  Total vertices     : ', total_nv
    write(*,'(A,I8)') '  Total faces        : ', total_nf
    write(*,'(A,I8)') '  Spiral vertices    : ', nv_spiral

    call write_mesh(trim(outfile), total_nv, total_nf, nv_spiral, pts_per_ring, verts, faces)

    deallocate(verts, faces, s_frac, cap_flag)
    write(*,'(A)') '======================================================'
    write(*,'(A)') '  Done.'

end program main
