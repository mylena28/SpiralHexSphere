! main.f90 — SpiralHexSphere mesh generator
program main
    use params
    use geometry
    use io
    implicit none

    integer  :: total_nv, total_nf, nv_spiral
    integer  :: total_rings, pts_per_ring
    real(dp), allocatable :: verts(:,:)
    integer,  allocatable :: faces(:,:)

    ! Equilateral-triangle condition: arc length per ring ≈ arc length per side
    ! 2*PI*helix_R*n_turns / n_s  =  2*PI*tube_r / num_sides
    ! => n_s = helix_R * n_turns * num_sides / tube_r
    n_s = nint(helix_R * n_turns * real(num_sides, dp) / tube_r)

    total_rings  = n_s
    pts_per_ring = num_sides

    write(*,'(A)') '======================================================'
    write(*,'(A)') '  SpiralHexSphere — Zinchenko Helix + 2 Hemisphere Caps'
    write(*,'(A)') '======================================================'
    write(*,'(A,F6.2)') '  Helix radius       : ', helix_R
    write(*,'(A,F6.2)') '  Helix height       : ', helix_h
    write(*,'(A,F6.2)') '  Number of turns    : ', n_turns
    write(*,'(A,F6.2)') '  Tube radius        : ', tube_r
    write(*,'(A,I4)')   '  Polygon sides      : ', num_sides
    write(*,'(A,I4)')   '  Pts / ring         : ', pts_per_ring
    write(*,'(A,I4)')   '  Original sections  : ', n_s
    write(*,'(A,I4)')   '  Total spiral rings : ', total_rings
    write(*,'(A,A)')    '  Cap method         : ', 'Level-0 triangulated (3 interior pts)'
    write(*,'(A,I4)')   '  Cap new verts/cap  : ', 3
    write(*,'(A,I4)')   '  Cap faces/cap      : ', 10
    write(*,'(A)') '------------------------------------------------------'

    call generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces)

    write(*,'(A,I8)') '  Total vertices     : ', total_nv
    write(*,'(A,I8)') '  Total faces        : ', total_nf
    write(*,'(A,I8)') '  Spiral vertices    : ', nv_spiral

    call write_mesh(trim(outfile), total_nv, total_nf, nv_spiral, pts_per_ring, verts, faces)

    deallocate(verts, faces)
    write(*,'(A)') '======================================================'
    write(*,'(A)') '  Done.'

end program main
